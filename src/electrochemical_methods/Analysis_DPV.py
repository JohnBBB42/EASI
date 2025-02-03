#_________________________________________
### Differential Pulse Voltametry (DPV)
#_________________________________________

class CustomValueError(Exception):
    """Custom error for handling specific value errors."""
    def __init__(self, message="An invalid value was provided."):
        super().__init__(message)

def Analysis_DPV(file_path, blank_responses=None):
    #_________________________________________
    # IMPORT PACKAGES
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

    #_________________________________________
    # DEFINE FUNCTIONS
    #_________________________________________

    def load_data_and_metadata(file_path):
        # Function to read lines from the file
        def read_file_lines(file_path, encoding='utf-16'):
            with open(file_path, 'r', encoding=encoding) as file:
                return file.readlines()

        def extract_metadata_row(lines, index=3):
            # Get the fourth row, remove leading 'Metadata row: ' if present
            metadata_row = lines[index].strip()
            metadata_row = re.sub(r'^Metadata row: ', '', metadata_row)
            metadata_entries = metadata_row.split(',,')
            
            # Clean each metadata entry by removing everything after 'DPV'
            cleaned_metadata_entries = []
            missing_parts = []
            for entry in metadata_entries:
                entry = entry.strip()  # Remove leading/trailing whitespace

                # If 'DPV' exists, remove everything after it
                if 'DPV' in entry:
                    entry = entry.split('DPV')[0] + 'DPV'
                else:
                    missing_parts.append("Method")
                cleaned_metadata_entries.append(entry)
            # Filter out any empty strings that might exist due to trailing delimiters
            cleaned_metadata_entries = [entry for entry in cleaned_metadata_entries if entry]
            return cleaned_metadata_entries, missing_parts
        
        def validate_metadata_entry(entry, entry_number, missing_parts):
            for entry_number, entry in enumerate(metadata_entries, start=1):
                entry_to_validate = re.sub(r'PP_', '', entry)
                parts_to_validate = entry_to_validate.split('_')

                # Check if parts_to_validate has at least 1 part to access index [0]
                if len(parts_to_validate) > 0:
                    if '-' in parts_to_validate[0]:
                        missing_parts.append("Electrode")

                # Check if parts_to_validate has at least 2 parts to access index [1]
                if len(parts_to_validate) > 1:
                    if '-' not in parts_to_validate[1]:
                        if any(char.isdigit() for char in parts_to_validate[1]):
                            missing_parts.append("Analytes")
                        else:
                            missing_parts.append("Concentration")
                
                if len(parts_to_validate) > 2:
                    if 'DPV' not in entry:
                        missing_parts.append("Method")

                # Raise an error if there are missing parts
                if missing_parts:
                    raise CustomValueError(f"Error in Entry #{entry_number}: Entry is missing required components: {', '.join(missing_parts)}. Expected format: 'Electrode_Analytes-Concentration_Method'.")

        def clean_metadata_entry(entry):
            entry = re.sub(r'PP_', '', entry)
            entry = entry.replace('-', '_')
            parts = entry.split('_')
            entry = '_'.join(parts[:4])
            return entry

        def parse_metadata_entry(entry):
            parts = entry.split('_')
            return {
                'Electrode': parts[0],
                'Analytes': parts[1],
                'Concentration': parts[2],
                'Method': parts[3]
            }

        def load_and_clean_data(file_path, encoding='utf-16', header_index=4, metadata_identifiers=['Date and time measurement:']):
            with open(file_path, 'r', encoding=encoding) as file:
                reader = csv.reader(file)
                data = list(reader)[header_index + 1:]

            headers = data[0]
            data_cleaned = [row for row in data if not any(identifier in row[0] for identifier in metadata_identifiers)]
            valid_data = [row for row in data_cleaned[1:] if len(row) == len(headers)]

            return headers, valid_data

        def convert_data_to_array(valid_data):
            def safe_float_convert(value):
                try:
                    return float(value)
                except ValueError:
                    return np.nan
            return np.array([[safe_float_convert(cell) for cell in row] for row in valid_data])

        # Read lines from the file
        lines = read_file_lines(file_path)

        # Extract and clean metadata
        metadata_entries, missing_parts = extract_metadata_row(lines)
        metadata_entries_val, missing_parts = extract_metadata_row(lines)
        
        # Validate Metadata Entries - Loop through all entries and stop if all are valid
        entry_number = []
        error_occurred = False

        try:
            for entry_number, entry in enumerate(metadata_entries_val, start=1):
                validate_metadata_entry(entry, entry_number, missing_parts)
        except CustomValueError as e:
            print(e)
            error_occurred = True

        if not error_occurred:
            # Clean metadata entries and parse metadata
            cleaned_metadata_entries = [clean_metadata_entry(entry) for entry in metadata_entries]
            parsed_metadata = [parse_metadata_entry(entry) for entry in cleaned_metadata_entries]

            # Load and clean the data
            headers, valid_data = load_and_clean_data(file_path)
            data_array = convert_data_to_array(valid_data)

            return headers, data_array, parsed_metadata
        else:
            print("Stopping due to validation error. No further processing will be performed, check your dataset.")
            return None, None, None

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
        
        lod_results = {}
        std_blank = np.std(blank_responses)
        analytes_dict = defaultdict(list)
        for key, mean_peak_current in mean_peak_currents.items():
            concentration, analytes = key
            analytes_dict[analytes].append((float(concentration[:-2]), mean_peak_current))
        
        for analytes, data in analytes_dict.items():
            data.sort(key=lambda x: x[0])
            concentrations, mean_peak_currents = zip(*data)
            for i, mean_peak_current in enumerate(zip(*mean_peak_currents)):
                log_concentrations = np.log10(concentrations)
                slope, intercept, r_value, p_value, std_err = linregress(log_concentrations, mean_peak_current)
                lod = (3.3 * std_blank) / slope
                lod_results[(analytes, f'Peak {i+1}')] = lod
        
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

    # Main execution
    headers, data_array, parsed_metadata = load_data_and_metadata(file_path)

    # Check if any of the returned values are None
    if headers is None or data_array is None or parsed_metadata is None:
        results = mean_peak_currents = std_peak_currents = cov_peak_currents = lod_results = t_test_results = None
        return results
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
