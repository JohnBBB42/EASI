##############################################
# analysis_dpv_streamlit.py
# A complete refactor of Analysis_DPV for Streamlit
##############################################

import os
import re
import io
import csv
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from collections import defaultdict
from itertools import combinations
from scipy.stats import t, linregress
from scipy.signal import find_peaks, savgol_filter
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.ndimage import gaussian_filter1d

###########################
# 1) Custom Error
###########################
class CustomValueError(Exception):
    """Custom error for handling specific value errors."""
    def __init__(self, message="An invalid value was provided."):
        super().__init__(message)

###########################
# 2) Core DPV Analysis (DataFrame-based)
###########################
def analysis_dpv_streamlit(dpv_df: pd.DataFrame, blank_df: pd.DataFrame = None) -> dict:
    """
    Refactored DPV analysis that works on a DataFrame instead of a file path.
    It expects the uploaded DPV file (CSV/Excel) to include a few metadata rows
    at the top. We try to detect these rows, remove them, and then use the remaining
    numeric data. If a blank_df is provided, it is used for LOD analysis.

    Returns a dictionary with:
      - results: list of dictionaries (peak detection info)
      - mean_peak_currents, std_peak_currents, cov_peak_currents
      - lod_results, t_test_results
      - headers, data_array, parsed_metadata
    """

    ##############################################
    # Step A: Detect and extract metadata rows
    ##############################################
    def detect_metadata_rows(df: pd.DataFrame, max_check=10):
        metadata_lines = []
        drop_indices = []
        for i in range(min(max_check, len(df))):
            row = df.iloc[i, :].astype(str).tolist()
            # Count non-numeric cells in the row
            non_numeric = sum(1 for cell in row if not re.match(r"^-?\d+\.?\d*$", cell))
            if non_numeric >= len(row) // 2:
                metadata_lines.append(','.join(row))
                drop_indices.append(i)
            else:
                break  # first numeric row encountered
        new_df = df.drop(index=drop_indices).reset_index(drop=True)
        return metadata_lines, new_df

    def extract_metadata_entries(metadata_lines):
        if not metadata_lines:
            return []
        # We assume the last metadata line is our metadata row
        last_line = metadata_lines[-1]
        row_str = re.sub(r'^Metadata row:\s*', '', last_line)
        # Split on the delimiter (here we assume either ',,' or ',')
        splitted = row_str.split(',,') if ',,' in row_str else row_str.split(',')
        return [s.strip() for s in splitted if s.strip()]

    def clean_and_parse_metadata(entries):
        def clean(entry):
            entry = re.sub(r'PP_', '', entry)
            entry = entry.replace('-', '_')
            parts = entry.split('_')
            return '_'.join(parts[:4])
        def parse(entry):
            parts = entry.split('_')
            return {
                'Electrode': parts[0] if len(parts) > 0 else "UnknownElectrode",
                'Analytes': parts[1] if len(parts) > 1 else "UnknownAnalytes",
                'Concentration': parts[2] if len(parts) > 2 else "UnknownConcentration",
                'Method': parts[3] if len(parts) > 3 else "UnknownMethod"
            }
        cleaned = [clean(e) for e in entries]
        return [parse(e) for e in cleaned]

    metadata_lines, numeric_df = detect_metadata_rows(dpv_df)
    metadata_entries = extract_metadata_entries(metadata_lines)
    parsed_metadata = clean_and_parse_metadata(metadata_entries) if metadata_entries else []

    ##############################################
    # Step B: Get headers and numeric data array
    ##############################################
    # Assume that the first row of the remaining df is the header if non-numeric
    header_row = numeric_df.iloc[0, :].astype(str).tolist()
    numeric_count = sum(1 for cell in header_row if re.match(r"^-?\d+\.?\d*$", cell))
    if numeric_count < len(header_row) // 2:
        # Use row 0 as headers and drop it from data
        headers = header_row
        data_df = numeric_df.drop(index=0).reset_index(drop=True)
    else:
        headers = [str(c) for c in range(numeric_df.shape[1])]
        data_df = numeric_df.copy()

    # Convert remaining data to a NumPy array (floats; missing values as NaN)
    data_array = data_df.to_numpy(dtype=float, copy=True)

    ##############################################
    # Step C: Find peaks and compute peak data
    ##############################################
    # Guess voltage and current columns from headers (simple case)
    voltage_columns = [i for i, col in enumerate(headers) if 'v' in col.lower()]
    current_columns = [i for i, col in enumerate(headers) if ('µa' in col.lower() or 'ua' in col.lower() or 'current' in col.lower())]
    if not voltage_columns: 
        voltage_columns = [0]
    if not current_columns:
        current_columns = [1]

    def baseline_als_simplified(y, lam=1e5, p=0.01, niter=10):
        L = len(y)
        D = sparse.csc_matrix(np.diff(np.eye(L), 2))
        w = np.ones(L)
        for i in range(niter):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.transpose())
            z = np.linalg.solve(Z.toarray(), (w * y))
            w = p * (y > z) + (1 - p) * (y < z)
        return z

    def find_peaks_in_data(v_data, c_data):
        window_length = 11
        polyorder = 3
        if len(c_data) < window_length:
            window_length = len(c_data) if len(c_data) % 2 == 1 else len(c_data) - 1
        if window_length <= polyorder or window_length < 3:
            smoothed = gaussian_filter1d(c_data, sigma=2)
        else:
            smoothed = savgol_filter(c_data, window_length=window_length, polyorder=polyorder)
        peaks, _ = find_peaks(smoothed, prominence=0.1)
        return peaks, smoothed

    collected_results = []
    # Loop over each metadata entry (or use a fallback if none)
    if parsed_metadata:
        loop_range = range(len(parsed_metadata))
    else:
        loop_range = range(min(len(voltage_columns), len(current_columns)))
    for idx in loop_range:
        # Use column idx if available; if not, fallback to first pair
        if idx >= len(voltage_columns) or idx >= len(current_columns):
            break
        vcol = voltage_columns[idx]
        ccol = current_columns[idx]
        v_data = data_array[:, vcol]
        c_data = data_array[:, ccol]
        valid = (~np.isnan(v_data)) & (~np.isnan(c_data))
        v_data = v_data[valid]
        c_data = c_data[valid]
        peaks, smoothed = find_peaks_in_data(v_data, c_data)
        baseline = baseline_als_simplified(smoothed)
        adjusted = smoothed - baseline
        for pk in peaks:
            if pk < len(v_data):
                collected_results.append({
                    'Entry': idx,
                    'Electrode': parsed_metadata[idx]['Electrode'] if parsed_metadata else "NA",
                    'Analytes': parsed_metadata[idx]['Analytes'] if parsed_metadata else "NA",
                    'Method': parsed_metadata[idx]['Method'] if parsed_metadata else "DPV",
                    'Concentration': parsed_metadata[idx]['Concentration'] if parsed_metadata else "unknown",
                    'Peak Voltage': v_data[pk],
                    'Peak Current': c_data[pk],
                    'Adjusted Peak Current': adjusted[pk],
                    'Y Offset': c_data[pk] - adjusted[pk],
                    'Baseline at Peak': baseline[pk]
                })

    ##############################################
    # Step D: Compute statistics (mean, std, cov)
    ##############################################
    def compute_mean_peak_currents(results):
        grouped = defaultdict(list)
        for r in results:
            key = (r['Concentration'], r['Analytes'])
            grouped[key].append(r['Adjusted Peak Current'])
        return {k: [np.mean(v)] for k, v in grouped.items()}

    def compute_std_peak_currents(results):
        grouped = defaultdict(list)
        for r in results:
            key = (r['Concentration'], r['Analytes'])
            grouped[key].append(r['Adjusted Peak Current'])
        return {k: [np.std(v)] for k, v in grouped.items()}

    def compute_cov_peaks(mean_peaks, std_peaks):
        covs = {}
        for k in mean_peaks:
            m = mean_peaks[k][0]
            s = std_peaks.get(k, [0])[0]
            covs[k] = [(s / m) * 100 if m != 0 else 0]
        return covs

    mean_peak_currents = compute_mean_peak_currents(collected_results)
    std_peak_currents = compute_std_peak_currents(collected_results)
    cov_peak_currents = compute_cov_peaks(mean_peak_currents, std_peak_currents)

    ##############################################
    # Step E: LOD analysis
    ##############################################
    def compute_lod(mean_peaks, blank_df):
        if blank_df is None or blank_df.empty:
            blank_array = np.array([0.1, 0.15, 0.05, 0.1, 0.08])
        else:
            blank_array = blank_df.to_numpy().flatten()
        std_blank = np.std(blank_array)
        results = {}
        # Group by analytes
        group = defaultdict(list)
        for key, m in mean_peaks.items():
            conc_str, analytes = key
            try:
                conc_val = float(conc_str.rstrip("uM"))
            except Exception:
                conc_val = 0
            group[analytes].append((conc_val, m[0]))
        for analyte, items in group.items():
            items.sort(key=lambda x: x[0])
            if len(items) < 2:
                continue
            concs, means = zip(*items)
            log_concs = np.log10(concs)
            slope, intercept, r_val, p_val, std_err = linregress(log_concs, means)
            lod = (3.3 * std_blank) / slope if slope != 0 else np.nan
            results[analyte] = lod
        return results

    lod_results = None
    if blank_df is not None:
        lod_results = compute_lod(mean_peak_currents, blank_df)

    ##############################################
    # Step F: T-test analysis
    ##############################################
    def perform_t_tests(mean_peaks, std_peaks, sample_size=3):
        t_results = []
        group = defaultdict(lambda: defaultdict(dict))
        for key, m in mean_peaks.items():
            conc_str, analytes = key
            try:
                conc_val = float(conc_str.rstrip("uM"))
            except Exception:
                conc_val = 0
            group[analytes][conc_val]['mean'] = m[0]
        for key, s in std_peaks.items():
            conc_str, analytes = key
            try:
                conc_val = float(conc_str.rstrip("uM"))
            except Exception:
                conc_val = 0
            group[analytes][conc_val]['std'] = s[0]
        for analyte, conc_dict in group.items():
            sorted_concs = sorted(conc_dict.keys())
            for (c1, c2) in combinations(sorted_concs, 2):
                mean1 = conc_dict[c1].get('mean', 0)
                mean2 = conc_dict[c2].get('mean', 0)
                std1 = conc_dict[c1].get('std', 0)
                std2 = conc_dict[c2].get('std', 0)
                n1 = n2 = sample_size
                pooled_std = np.sqrt(((n1 - 1) * std1 ** 2 + (n2 - 1) * std2 ** 2) / (n1 + n2 - 2))
                t_stat = (mean1 - mean2) / (pooled_std * np.sqrt(1 / n1 + 1 / n2)) if pooled_std != 0 else 0
                df = n1 + n2 - 2
                crit_val = t.ppf(1 - 0.025, df)
                t_results.append({
                    'Analyte': analyte,
                    'Concentration1': c1,
                    'Concentration2': c2,
                    'T-statistic': t_stat,
                    'Critical value': crit_val,
                    'Significant': abs(t_stat) > crit_val
                })
        return t_results

    t_test_results = perform_t_tests(mean_peak_currents, std_peak_currents, sample_size=3)

    ##############################################
    # Return final dictionary
    ##############################################
    return {
        'results': collected_results,
        'mean_peak_currents': mean_peak_currents,
        'std_peak_currents': std_peak_currents,
        'cov_peak_currents': cov_peak_currents,
        'lod_results': lod_results,
        't_test_results': t_test_results,
        'headers': [str(h) for h in headers],
        'data_array': data_array,
        'parsed_metadata': parsed_metadata
    }
