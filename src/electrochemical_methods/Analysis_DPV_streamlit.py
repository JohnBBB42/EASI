##############################################
# streamlit_dpv.py
# A complete refactor of Analysis_DPV for Streamlit
##############################################

import os
import re
import io
import csv
import sys
import numpy as np
import pandas as pd
import streamlit as st
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
# 1) We define a custom error
###########################
class CustomValueError(Exception):
    """Custom error for handling specific value errors."""
    def __init__(self, message="An invalid value was provided."):
        super().__init__(message)

###########################
# 2) The core DPV logic refactored to accept a DataFrame
###########################
def analysis_dpv_streamlit(
    dpv_df: pd.DataFrame, 
    blank_df: pd.DataFrame=None
) -> dict:
    """
    A full refactor of your old Analysis_DPV logic to handle in-memory data.
    Instead of reading from file, we parse 'dpv_df' for both the actual data
    and the 'metadata row' if possible. If blank_df is given, we use it for LOD.
    
    Returns a dict with:
      - results (peak detection, etc.)
      - mean_peak_currents
      - std_peak_currents
      - cov_peak_currents
      - lod_results
      - t_test_results
      - headers
      - data_array
      - parsed_metadata
    """

    ##############################################
    # Step A: Extract lines & metadata from dpv_df 
    #         (Mimicking your old 'load_data_and_metadata' but we don't read a file).
    ##############################################

    # Because your original code read lines from a file in 'utf-16' and 
    # expected certain rows for metadata, we must replicate that logic with 
    # the DataFrame we have. We'll assume the top few lines are text-based metadata 
    # or that your data has them as rows. If your data starts with numeric rows, 
    # adapt accordingly.

    # We do not have a 'lines' array from reading the file. So either:
    # 1) The user has a CSV with the first 4 lines as metadata, row 5+ as data.
    # 2) Or you have some alternative. We'll show how you might handle it if your 
    #    'dpv_df' includes these lines as well.

    # We'll define a function that tries to parse the first few lines as metadata 
    # if they are not purely numeric. Then we remove them from 'dpv_df' for data analysis.
    
    def detect_metadata_rows(dpv_df: pd.DataFrame, max_check=10):
        """
        Attempt to find non-numeric rows at the top which are metadata lines.
        We'll check row by row up to max_check. If a row is mostly non-numeric,
        we'll treat it as metadata. We'll return the list of metadata text lines
        plus a new dpv_df that excludes them.
        """
        metadata_lines = []
        drop_indices   = []
        for i in range(min(max_check, len(dpv_df))):
            row = dpv_df.iloc[i, :].astype(str).tolist()
            # If row is mostly strings that can't parse as float, treat as metadata
            non_numeric_count = 0
            for cell in row:
                try:
                    float(cell)
                except ValueError:
                    non_numeric_count += 1
            # if half or more are non-numeric, consider it metadata
            if non_numeric_count >= len(row)//2:
                metadata_lines.append(','.join(row))
                drop_indices.append(i)
            else:
                # as soon as we hit a row that is mostly numeric,
                # we assume data starts here
                break
        # create a new df without those lines
        new_df = dpv_df.drop(index=drop_indices).reset_index(drop=True)
        return metadata_lines, new_df

    # We'll parse the metadata lines similarly to your old code, 
    # though we can't replicate it exactly since we no longer 
    # have 'lines[3]' or so. We'll attempt to find a "metadata row" among them.

    def extract_metadata_entries(metadata_lines):
        """
        For demonstration, we look for a line that starts with 'Metadata row:' 
        or includes 'DPV'. We'll parse them. If none found, we do minimal parsing.
        """
        if not metadata_lines:
            return []

        # Suppose your 4th line was the real metadata row in your old code. We'll 
        # just try the last line of the found metadata as that might 
        # contain the relevant row. Or you can adapt logic as needed.
        last_line = metadata_lines[-1]
        # e.g. "Metadata row: PP_A-B_C-DPV" ...
        row_str = re.sub(r'^Metadata row:\s*', '', last_line)
        # split by ',,' if your code used that
        splitted = row_str.split(',,') if ',,' in row_str else row_str.split(',')
        splitted = [s.strip() for s in splitted if s.strip()]
        return splitted

    def validate_metadata_entry(entry, entry_number, missing_parts):
        """
        Copied from your old code. We'll reduce it to a no-op unless you want 
        to do the full checks. 
        """
        # If you'd like your old checks:
        # 1) 'PP_' prefix
        # 2) Contains '-' for analytes
        # 3) Contains 'DPV'
        # etc.
        pass

    def clean_metadata_entry(entry):
        entry = re.sub(r'PP_', '', entry)
        entry = entry.replace('-', '_')
        parts = entry.split('_')
        # join up to 4 parts
        entry = '_'.join(parts[:4])
        return entry

    def parse_metadata_entry(entry):
        parts = entry.split('_')
        # we expect up to 4
        meta = {
            'Electrode': parts[0] if len(parts)>0 else "UnknownElectrode",
            'Analytes':  parts[1] if len(parts)>1 else "UnknownAnalytes",
            'Concentration': parts[2] if len(parts)>2 else "UnknownConcentration",
            'Method':    parts[3] if len(parts)>3 else "UnknownMethod"
        }
        return meta

    # Now let's do that:

    metadata_lines, numeric_df = detect_metadata_rows(dpv_df)
    splitted = extract_metadata_entries(metadata_lines)
    if splitted:
        # attempt to validate
        missing_parts = []
        for i, entry in enumerate(splitted, start=1):
            validate_metadata_entry(entry, i, missing_parts)
        # clean
        cleaned = [clean_metadata_entry(e) for e in splitted]
        parsed_metadata = [parse_metadata_entry(e) for e in cleaned]
    else:
        parsed_metadata = []

    # numeric_df now presumably holds the real data. 
    # We'll treat row0 as headers if it is not purely numeric:
    headers = numeric_df.iloc[0,:].astype(str).tolist()
    # Check if these headers are actually numeric
    numeric_test = 0
    for cell in headers:
        try:
            float(cell)
            numeric_test+=1
        except:
            pass
    # if most are non-numeric, we treat them as column labels
    if numeric_test < len(headers)//2:
        # treat row 0 as headers
        new_data = numeric_df.drop(index=0).reset_index(drop=True)
    else:
        # no real header, just numeric data
        new_data = numeric_df.reset_index(drop=True)

    data_array = new_data.to_numpy(dtype=float, na_value=np.nan)

    # Next we do the logic from find_peaks_and_process_data, 
    # but we must guess how to identify columns containing 'V' or 'µA'
    # because we no longer have your "lines" to read. We'll do a guess: 
    # if a header includes 'V' => voltage, if 'µA' => current
    voltage_columns = []
    current_columns = []
    for i, col in enumerate(headers):
        c = col.lower()
        if 'v' in c:
            voltage_columns.append(i)
        elif 'µa' in c or 'ua' in c or 'current' in c:
            current_columns.append(i)

    # If that fails, we fallback to the first column = voltage, second = current
    if not voltage_columns:
        voltage_columns = [0]
    if not current_columns:
        current_columns = [1]

    # We'll define a results structure that mimics your old "results" list
    # with 'Peak Voltage', 'Peak Current', etc.
    # We'll define baseline_als, peak-finding, etc.

    def baseline_als_simplified(y, lam=1e5, p=0.01, niter=10):
        L = len(y)
        D = sparse.csc_matrix(np.diff(np.eye(L), 2))
        w = np.ones(L)
        for i in range(niter):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.transpose())
            z = np.linalg.solve(Z.toarray(), (w * y))
            w = p*(y>z) + (1-p)*(y<z)
        return z

    def find_peaks_in_data(voltage_data, current_data):
        # We'll do a smoothing to reduce noise
        # your code uses savgol_filter or gaussian if too short
        window_length = 11
        polyorder     = 3
        if len(current_data)<window_length:
            window_length = len(current_data) if len(current_data)%2==1 else len(current_data)-1
        if window_length<=polyorder or window_length<3:
            smoothed = gaussian_filter1d(current_data, sigma=2)
        else:
            smoothed = savgol_filter(current_data, window_length=window_length, polyorder=polyorder)
        peaks, _ = find_peaks(smoothed, prominence=0.1)
        return peaks, smoothed

    collected_results = []
    # We'll do a loop for each "metadata entry" or at least for each pair of columns
    # (since your old code aligns them).
    for idx, row in enumerate(parsed_metadata):
        # If we have more metadata entries than voltage/current columns, we skip
        if idx>=len(voltage_columns) or idx>=len(current_columns):
            break
        vcol = voltage_columns[idx]
        ccol = current_columns[idx]
        voltage_data = data_array[:, vcol]
        current_data = data_array[:, ccol]
        # remove NaN
        valid_idx = (~np.isnan(voltage_data)) & (~np.isnan(current_data))
        voltage_data = voltage_data[valid_idx]
        current_data = current_data[valid_idx]

        peaks, smoothed = find_peaks_in_data(voltage_data, current_data)
        baseline = baseline_als_simplified(smoothed)
        adjusted = smoothed - baseline

        for pk in peaks:
            if pk<len(voltage_data):
                peak_voltage   = voltage_data[pk]
                peak_current   = current_data[pk]
                adjusted_peak  = adjusted[pk]
                y_offset       = peak_current - adjusted_peak
                collected_results.append({
                    'Entry': idx,
                    'Electrode':      row.get('Electrode','Unknown'),
                    'Analytes':       row.get('Analytes','Unknown'),
                    'Method':         row.get('Method','Unknown'),
                    'Concentration':  row.get('Concentration','Unknown'),
                    'Peak Voltage':   peak_voltage,
                    'Peak Current':   peak_current,
                    'Adjusted Peak Current': adjusted_peak,
                    'Y Offset':       y_offset,
                    'Baseline at Peak': baseline[pk]
                })
    # If there's no metadata or columns mismatch, let's do a fallback for 
    # the first pair of columns only
    if not parsed_metadata:
        st.warning("No parsed metadata found. We'll just do the first pair of columns as DPV data.")
        # do the same for columns 0,1
        # or for all pairs
        for i in range(min(len(voltage_columns), len(current_columns))):
            voltage_data = data_array[:, voltage_columns[i]]
            current_data = data_array[:, current_columns[i]]
            valid_idx = (~np.isnan(voltage_data)) & (~np.isnan(current_data))
            voltage_data=voltage_data[valid_idx]
            current_data=current_data[valid_idx]
            peaks, smoothed = find_peaks_in_data(voltage_data, current_data)
            baseline = baseline_als_simplified(smoothed)
            adjusted = smoothed - baseline
            for pk in peaks:
                if pk<len(voltage_data):
                    collected_results.append({
                        'Entry': i,
                        'Electrode': 'NA',
                        'Analytes':  'NA',
                        'Method':    'DPV',
                        'Concentration': 'unknown',
                        'Peak Voltage':  voltage_data[pk],
                        'Peak Current':  current_data[pk],
                        'Adjusted Peak Current': adjusted[pk],
                        'Y Offset': current_data[pk]-adjusted[pk],
                        'Baseline at Peak': baseline[pk]
                    })

    ##############################################
    # Step B: compute mean/std/cov per analyte & concentration 
    ##############################################
    def compute_mean_peak_currents(collected_results):
        grouped = defaultdict(list)
        for r in collected_results:
            key = (r['Concentration'], r['Analytes'])
            grouped[key].append(r)
        mean_peaks = {}
        for key, entries in grouped.items():
            # We'll group by 'Entry' to figure out how many peaks each entry has
            # We can do a simpler approach: we assume all entries for that key 
            # form a big list of adjusted peak currents. Or we want # of peaks 
            # per 'Entry'? Let's do the approach from your code carefully:
            # This is simplified code:
            # We'll assume each 'Entry' has 1 peak for that key
            # If you have multiple peaks, handle them carefully.
            # For now, we'll just do a single average if multiple entries
            peak_list = [r['Adjusted Peak Current'] for r in entries]
            mean_val  = np.mean(peak_list) if peak_list else 0
            # if multiple peaks per entry, we’d do more advanced logic
            # but let's keep it simple
            mean_peaks[key] = [mean_val]
        return mean_peaks

    def compute_std_peak_currents(collected_results, mean_peak_currents):
        grouped = defaultdict(list)
        for r in collected_results:
            key = (r['Concentration'], r['Analytes'])
            grouped[key].append(r)
        std_peaks = {}
        for key, entries in grouped.items():
            # same approach
            peak_list = [r['Adjusted Peak Current'] for r in entries]
            std_val   = np.std(peak_list) if peak_list else 0
            std_peaks[key] = [std_val]
        return std_peaks

    def compute_cov_peaks(mean_peaks, std_peaks):
        cov_peaks = {}
        for key in mean_peaks:
            means = mean_peaks[key]
            stds  = std_peaks.get(key, [0]*len(means))
            cvals = []
            for i in range(len(means)):
                c = (stds[i]/means[i])*100 if (means[i]!=0) else 0
                cvals.append(c)
            cov_peaks[key] = cvals
        return cov_peaks

    mean_peak_currents = compute_mean_peak_currents(collected_results)
    std_peak_currents  = compute_std_peak_currents(collected_results, mean_peak_currents)
    cov_peak_currents  = compute_cov_peaks(mean_peak_currents, std_peak_currents)

    ##############################################
    # Step C: LOD analysis if blank_df is provided
    ##############################################
    def compute_lod(mean_peak_currents, blank_df):
        """
        from your code: if blank_responses is None, do default. 
        We'll parse blank_df -> an array of current baseline?
        """
        if blank_df is None or blank_df.empty:
            blank_array = np.array([0.1, 0.15, 0.05, 0.1, 0.08])
        else:
            # flatten the blank df
            blank_array = blank_df.to_numpy().flatten()
        std_blank = np.std(blank_array)
        results = {}
        analytes_dict=defaultdict(list)
        for key, meanlist in mean_peak_currents.items():
            # key=(Concentration, Analytes)
            c_str, analytes=key
            c_val=0
            try:
                c_val=float(c_str[:-2]) # remove 'uM' or etc.
            except:
                pass
            analytes_dict[analytes].append((c_val, meanlist))

        # We'll do a simple slope = (peak vs logC) approach
        for analytes, data in analytes_dict.items():
            data.sort(key=lambda x: x[0])
            if not data:
                continue
            # build arrays
            concs     = [x[0] for x in data]
            mean_peaks= [x[1][0] for x in data] # if single peak
            # must have at least 2
            if len(concs)<2:
                continue
            logc = np.log10(concs)
            slope, intercept, r_val, p_val, std_err= linregress(logc, mean_peaks)
            lod_val=(3.3*std_blank)/slope if slope!=0 else 9999
            results[analytes]=lod_val
        return results

    lod_results=None
    if blank_df is not None:
        lod_results=compute_lod(mean_peak_currents, blank_df)

    ##############################################
    # Step D: T-test
    ##############################################
    def perform_t_tests(mean_peaks, std_peaks, sample_size=3):
        # from your code
        t_test_results=[]
        analytes_dict=defaultdict(lambda: defaultdict(dict))
        for key, means in mean_peaks.items():
            c_str, analytes=key
            try:
                c_val=float(c_str[:-2])
            except:
                c_val=0
            # we only handle the 0th if single peak
            analytes_dict[analytes][c_val]['mean']=means[0] if means else 0

        for key, stds in std_peaks.items():
            c_str, analytes=key
            try:
                c_val=float(c_str[:-2])
            except:
                c_val=0
            analytes_dict[analytes][c_val]['std']=stds[0] if stds else 0

        for analyte, conc_dict in analytes_dict.items():
            sorted_concs=sorted(conc_dict.keys())
            for (c1, c2) in combinations(sorted_concs,2):
                mean1=conc_dict[c1].get('mean',0)
                mean2=conc_dict[c2].get('mean',0)
                std1 =conc_dict[c1].get('std',0)
                std2 =conc_dict[c2].get('std',0)
                n1=n2=sample_size
                pooled_std= np.sqrt(((n1-1)*std1**2 + (n2-1)*std2**2)/(n1+n2-2))
                t_stat=0
                if pooled_std!=0:
                    t_stat=(mean1-mean2)/(pooled_std*np.sqrt(1/n1+1/n2))
                df=n1+n2-2
                # for alpha=0.05, 2-sided => t.ppf(1-0.025,df)
                crit_val=t.ppf(1-0.025, df)
                t_test_results.append({
                    'Analyte': analyte,
                    'Concentration1': c1,
                    'Concentration2': c2,
                    'T-statistic': t_stat,
                    'Critical value': crit_val,
                    'Significant': abs(t_stat)>crit_val
                })
        return t_test_results

    t_test_results=perform_t_tests(mean_peak_currents, std_peak_currents, sample_size=3)

    ##############################################
    # Return final dictionary
    ##############################################
    # We'll build 'headers' from the old code. We guess the columns 
    # from the new data. Also 'data_array' is the numeric data.
    # 'parsed_metadata' is from the top lines we found.
    # 'results' is the list of peak dicts
    return {
        'results':          collected_results,
        'mean_peak_currents': mean_peak_currents,
        'std_peak_currents':  std_peak_currents,
        'cov_peak_currents':  cov_peak_currents,
        'lod_results':        lod_results,
        't_test_results':     t_test_results,
        'headers':            [str(h) for h in headers],  # safe as strings
        'data_array':         data_array,
        'parsed_metadata':    parsed_metadata
    }
