#_________________________________________
### Differential Pulse Voltametry (DPV)
#_________________________________________

import os
import re
import sys
import csv
import pylab as p
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from matplotlib.pyplot import *
from scipy.optimize import curve_fit
from collections import defaultdict
from itertools import combinations
from scipy.stats import t, linregress
from scipy.signal import find_peaks, savgol_filter
from scipy import sparse
from scipy.sparse.linalg import spsolve
import warnings

class CustomValueError(Exception):
    """Custom error for handling specific value errors."""
    def __init__(self, message="An invalid value was provided."):
        super().__init__(message)

def Analysis_DPV(file_path, blank_responses=None, metadata_fields=None):
    #_________________________________________
    # DEFINE FUNCTIONS
    #_________________________________________

    #---------------------------
    # Flexible metadata loading function
    #---------------------------
    def load_data_and_metadata(file_path, metadata_fields=None, encoding='utf-16', metadata_header_index=3, metadata_delimiter=",,"):
        """
        Load the file, extract and validate metadata using a flexible naming scheme,
        and load the accompanying numerical data.

        Parameters:
            file_path (str): Path to the file.
            metadata_fields (list, optional): Expected metadata field names.
                Default is ["Electrode", "Analytes", "Concentration", "Method"].
            encoding (str): File encoding. Default 'utf-16'.
            metadata_header_index (int): Row index for the metadata header (default is 3, i.e. the fourth line).
            metadata_delimiter (str): Delimiter between metadata entries (default is ",,").

        Returns:
            tuple: (headers, data_array, parsed_metadata) on success; (None, None, None) if errors occur.
        """
        if metadata_fields is None:
            metadata_fields = ["Electrode", "Analytes", "Concentration", "Method"]

        # ---------------------------
        # Helper functions (unchanged except for flexibility)
        # ---------------------------
        def read_file_lines(fp, encoding=encoding):
            with open(fp, 'r', encoding=encoding) as file:
                return file.readlines()

        def extract_metadata_row(lines, index=metadata_header_index, delimiter=metadata_delimiter):
            metadata_row = lines[index].strip()
            metadata_row = re.sub(r'^Metadata row: ', '', metadata_row)
            metadata_entries = [entry.strip() for entry in metadata_row.split(delimiter) if entry.strip()]
            cleaned_entries = []
            missing_parts = []
            for entry in metadata_entries:
                entry = entry.strip()
                if 'DPV' in entry:
                    entry = entry.split('DPV')[0] + 'DPV'
                else:
                    missing_parts.append("Method")
                cleaned_entries.append(entry)
            return cleaned_entries, missing_parts

        def validate_metadata_entry(entry, entry_number, missing_parts):
            # (Exactly as in your original code)
            entry_to_validate = re.sub(r'PP_', '', entry)
            parts_to_validate = entry_to_validate.split('_')
            if len(parts_to_validate) > 0:
                if '-' in parts_to_validate[0]:
                    missing_parts.append("Electrode")
            if len(parts_to_validate) > 1:
                if '-' not in parts_to_validate[1]:
                    if any(char.isdigit() for char in parts_to_validate[1]):
                        missing_parts.append("Analytes")
                    else:
                        missing_parts.append("Concentration")
            if len(parts_to_validate) > 2:
                if 'DPV' not in entry:
                    missing_parts.append("Method")
            if missing_parts:
                raise CustomValueError(f"Error in Entry #{entry_number}: Entry is missing required components: {', '.join(missing_parts)}. Expected format: {'_'.join(metadata_fields)}.")

        def clean_metadata_entry(entry):
            # Exactly as in your original code
            entry = re.sub(r'PP_', '', entry)
            entry = entry.replace('-', '_')
            parts = entry.split('_')
            # Here, we use only as many parts as defined by metadata_fields
            entry = '_'.join(parts[:len(metadata_fields)])
            return entry

        def parse_metadata_entry(entry):
            parts = entry.split('_')
            # Assume len(parts) is >= len(metadata_fields)
            return {field: part for field, part in zip(metadata_fields, parts)}

        def load_and_clean_data(fp, encoding=encoding, header_index=4, metadata_identifiers=['Date and time measurement:']):
            with open(fp, 'r', encoding=encoding) as file:
                reader = csv.reader(file)
                data = list(reader)[header_index + 1:]
            headers = data[0]
            data_cleaned = [row for row in data if not any(identifier in row[0] for identifier in metadata_identifiers)]
            valid_data = [row for row in data_cleaned[1:] if len(row) == len(headers)]
            return headers, valid_data

        def safe_float_convert(value):
            try:
                return float(value)
            except ValueError:
                return np.nan

        def convert_data_to_array(data_rows):
            return np.array([[safe_float_convert(cell) for cell in row] for row in data_rows])
        
        # --- Execute the metadata and data loading ---
        lines = read_file_lines(file_path)

        metadata_entries, missing_parts = extract_metadata_row(lines)
        # Make a second call to get a separate copy for validation:
        metadata_entries_val, missing_parts = extract_metadata_row(lines)
        
        error_occurred = False
        try:
            for entry_number, entry in enumerate(metadata_entries_val, start=1):
                validate_metadata_entry(entry, entry_number, missing_parts)
        except CustomValueError as e:
            print(e)
            error_occurred = True

        if not error_occurred:
            # Clean each metadata entry and then parse it
            cleaned_metadata_entries = [clean_metadata_entry(entry) for entry in metadata_entries]
            parsed_metadata = [parse_metadata_entry(entry) for entry in cleaned_metadata_entries]
        else:
            print("Stopping due to validation error. No further processing will be performed, check your dataset.")
            return None, None, None

        headers, valid_data = load_and_clean_data(file_path)
        data_array = convert_data_to_array(valid_data)
        return headers, data_array, parsed_metadata

    def find_peaks_and_process_data(headers, data_array, parsed_metadata):
        voltage_columns = [i for i, col in enumerate(headers) if 'V' in col]
        current_columns = [i for i, col in enumerate(headers) if 'µA' in col]

        def baseline_als(y, lam=1e5, p=0.01, niter=10):
            L = len(y)
            D = sparse.csc_matrix(np.diff(np.eye(L), 2))
            w = np.ones(L)
            for i in range(niter):
                W = sparse.spdiags(w, 0, L, L)
                Z = W + lam * D.dot(D.transpose())
                z = spsolve(Z, w*y)
                w = p * (y > z) + (1-p) * (y < z)
            return z
        from scipy.ndimage import gaussian_filter1d

        def find_peaks_in_data(voltage_data, current_data):
            window_length = 11; polyorder = 3
            if len(current_data) < window_length:
                window_length = len(current_data) if len(current_data) % 2 == 1 else len(current_data) - 1 
            if window_length <= polyorder:
                print("Data is too short for the specified window_length and polyorder. Using Gaussian filter instead.")
                smoothed_current = gaussian_filter1d(current_data, sigma=2)  # Use sigma=2 for smoothing
            else:
                smoothed_current = savgol_filter(current_data, window_length=window_length, polyorder=polyorder)
            peaks, _ = find_peaks(smoothed_current, prominence=0.1)
            return peaks, smoothed_current

        results = []
        for index, row in enumerate(parsed_metadata):
            # If the lengths do not match, print a warning and return
            if len(parsed_metadata) > len(voltage_columns) or len(parsed_metadata) > len(current_columns):
                raise ValueError("Mismatch between parsed metadata entries and available columns. Please check your dataset.")

            voltage_data = data_array[:, voltage_columns[index]]
            current_data = data_array[:, current_columns[index]]

            # Remove rows with NaN values
            valid_idx = ~np.isnan(voltage_data) & ~np.isnan(current_data)
            voltage_data = voltage_data[valid_idx]
            current_data = current_data[valid_idx]

            peaks, smoothed_current = find_peaks_in_data(voltage_data, current_data)

            # Apply ALS baseline correction
            baseline = baseline_als(smoothed_current)
            adjusted_current = smoothed_current - baseline

            for peak_idx in peaks:
                peak_voltage = voltage_data[peak_idx]
                peak_current = current_data[peak_idx]
                adjusted_peak_current = adjusted_current[peak_idx]
                y_offset = peak_current - adjusted_peak_current

                results.append({
                    'Entry': index,
                    'Electrode': row['Electrode'],
                    'Analytes': row['Analytes'],
                    'Method': row['Method'],
                    'Concentration': row['Concentration'],
                    'Peak Voltage': peak_voltage,
                    'Peak Current': peak_current,
                    'Adjusted Peak Current': adjusted_peak_current,
                    'Y Offset': y_offset,
                    'Baseline at Peak': baseline[peak_idx]
                })

        return results

    def calculate_cov_peak_currents(results):
        def calculate_mean_peak_currents(results):
            grouped_data = defaultdict(list)
            for result in results:
                key = (result['Concentration'], result['Analytes'])
                grouped_data[key].append(result)
            
            mean_peak_currents = {}
            for key, entries in grouped_data.items():
                max_peaks = max(len([res for res in entries if res['Entry'] == entry['Entry']])
                                for entry in entries)
                peak_sums = [0] * max_peaks
                peak_counts = [0] * max_peaks
                for entry in entries:
                    peaks_in_entry = [res for res in entries if res['Entry'] == entry['Entry']]
                    for i, res in enumerate(peaks_in_entry):
                        if i < max_peaks:
                            peak_sums[i] += res['Adjusted Peak Current']
                            peak_counts[i] += 1
                mean_peak_currents[key] = [peak_sums[i] / peak_counts[i] if peak_counts[i] > 0 else 0
                                           for i in range(max_peaks)]
            return mean_peak_currents
        
        def calculate_std_peak_currents(results):
            # Group results by (Concentration, Analytes)
            grouped_entries = defaultdict(lambda: defaultdict(list))
            for result in results:
                key = (result['Concentration'], result['Analytes'])
                entry = result['Entry']
                grouped_entries[key][entry].append(result['Adjusted Peak Current'])
            
            # For each group, align peaks across entries
            std_peak_currents = {}
            for key, entries in grouped_entries.items():
                # Find the maximum number of peaks across all entries
                max_peaks = max(len(peak_list) for peak_list in entries.values())
                peak_currents = [[] for _ in range(max_peaks)]
                
                # Align peaks by their position
                for peak_list in entries.values():
                    for i, current in enumerate(peak_list):
                        peak_currents[i].append(current)
                
                # Calculate standard deviation for each peak position
                std_list = [np.std(currents) if len(currents) > 1 else 0 for currents in peak_currents]
                std_peak_currents[key] = std_list
            return std_peak_currents
        

        mean_peak_currents = calculate_mean_peak_currents(results)
        std_peak_currents = calculate_std_peak_currents(results)

        cov_peak_currents = {}
        for key in mean_peak_currents.keys():
            means = mean_peak_currents[key]
            stds = std_peak_currents[key]
            covs = [(stds[i] / means[i]) * 100 if means[i] != 0 else 0 for i in range(len(means))]
            cov_peak_currents[key] = covs

        return cov_peak_currents, mean_peak_currents, std_peak_currents
    
    def calculate_lod(mean_peak_currents, blank_responses=None):
        if blank_responses is None:
            blank_responses = [0.1, 0.15, 0.05, 0.1, 0.08]
        
        std_blank = np.std(blank_responses)
        lod_results = {}
        analytes_dict = defaultdict(list)

        # Build a dict: analyte -> list of (concentration_float, [peak1, peak2, ...])
        for (conc_str, analytes), peak_list in mean_peak_currents.items():
            try:
                # strip off unit, e.g. "50µM" → 50.0
                c_val = float(conc_str.rstrip('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZµμ'))
            except ValueError:
                continue
            analytes_dict[analytes].append((c_val, peak_list))

        for analytes, data in analytes_dict.items():
            # sort by concentration
            data.sort(key=lambda x: x[0])
            concs, peaks_matrix = zip(*data)
            n_points = len(concs)

            if n_points < 2:
                warnings.warn(
                    f"Cannot calculate LOD for '{analytes}': only {n_points} concentration point(s).",
                    RuntimeWarning
                )
                # still reserve keys, but set to nan
                n_peaks = len(peaks_matrix[0])
                for i in range(n_peaks):
                    lod_results[(analytes, f"Peak {i+1}")] = np.nan
                continue

            log_concs = np.log10(np.array(concs, dtype=float) + 1e-12)
            if np.allclose(log_concs, log_concs[0]):
                warnings.warn(
                    f"Zero variance in log(conc) for '{analytes}'; skipping LOD.",
                    RuntimeWarning
                )
                for i in range(len(peaks_matrix[0])):
                    lod_results[(analytes, f"Peak {i+1}")] = np.nan
                continue

            # now loop over each peak position
            for peak_idx, peak_vals in enumerate(zip(*peaks_matrix), start=1):
                y = np.array(peak_vals, dtype=float)
                # if all y are the same, slope=0 → skip
                slope, intercept, r_val, p_val, stderr = linregress(log_concs, y)
                if slope == 0 or np.isnan(slope):
                    warnings.warn(
                        f"Slope zero or NaN for '{analytes}' peak {peak_idx}; LOD undefined.",
                        RuntimeWarning
                    )
                    lod_results[(analytes, f"Peak {peak_idx}")] = np.nan
                else:
                    lod = (3.3 * std_blank) / slope
                    lod_results[(analytes, f"Peak {peak_idx}")] = lod

        return lod_results

    def perform_t_tests(mean_peak_currents, std_peak_currents, sample_size=None):
        analytes_dict = defaultdict(lambda: defaultdict(dict))
        for key, mean_peak_current in mean_peak_currents.items():
            concentration, analytes = key
            analyte_list = analytes.split(',')
            for i, analyte in enumerate(analyte_list):
                analyte_name = analyte.strip()
                if len(analyte_list) > 1:
                    analyte_name = '3' + analyte_name
                analytes_dict[analyte_name][float(concentration[:-2])]['mean'] = mean_peak_current[i]
        
        for key, std_peak_current in std_peak_currents.items():
            concentration, analytes = key
            analyte_list = analytes.split(',')
            for i, analyte in enumerate(analyte_list):
                analyte_name = analyte.strip()
                if len(analyte_list) > 1:
                    analyte_name = '3' + analyte_name
                analytes_dict[analyte_name][float(concentration[:-2])]['std'] = std_peak_current[i]
        
        t_test_results = []
        for analyte, conc_dict in analytes_dict.items():
            concentrations = sorted(conc_dict.keys())
            for (concentration1, concentration2) in combinations(concentrations, 2):
                mean1 = conc_dict[concentration1]['mean']
                mean2 = conc_dict[concentration2]['mean']
                std1 = conc_dict[concentration1]['std']
                std2 = conc_dict[concentration2]['std']
                n1 = n2 = sample_size if sample_size is not None else 3
                
                pooled_std = np.sqrt(((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2))
                if pooled_std == 0:
                    t_stat = 0
                else:
                    t_stat = (mean1 - mean2) / (pooled_std * np.sqrt(1/n1 + 1/n2))
                
                df = n1 + n2 - 2
                critical_value = t.ppf(1 - 0.025, df)
                
                t_test_result = {
                    'Analyte': analyte,
                    'Concentration1': concentration1,
                    'Concentration2': concentration2,
                    'T-statistic': t_stat,
                    'Critical value': critical_value,
                    'Significant': abs(t_stat) > critical_value
                }
                
                t_test_results.append(t_test_result)
        
        return t_test_results

    # Main execution using the new loader:
    headers, data_array, parsed_metadata = load_data_and_metadata(file_path, metadata_fields=metadata_fields)

    if headers is None or data_array is None or parsed_metadata is None:
        return None
    else:
        results = find_peaks_and_process_data(headers, data_array, parsed_metadata)
        cov_peak_currents, mean_peak_currents, std_peak_currents = calculate_cov_peak_currents(results)
        lod_results = calculate_lod(mean_peak_currents, blank_responses)
        t_test_results = perform_t_tests(mean_peak_currents, std_peak_currents)
        return {
            'results': results,
            'mean_peak_currents': mean_peak_currents,
            'std_peak_currents': std_peak_currents,
            'cov_peak_currents': cov_peak_currents,
            'lod_results': lod_results,
            't_test_results': t_test_results,
            'headers': headers,
            'data_array': data_array,
            'parsed_metadata': parsed_metadata
        }