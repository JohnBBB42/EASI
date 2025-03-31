import os
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.stats import linregress

class DPVPlotting:
    def plot_data_array(self, data_array, headers, parsed_metadata, save_path, selected_indices=None):
        plt.figure(figsize=(10, 6))

        if selected_indices is None:
            selected_indices = range(0, len(headers), 2)

        for i in selected_indices:
            potential = data_array[:, i]
            current = data_array[:, i + 1]
            
            meta = parsed_metadata[i//2]  # if each pair of columns matches one metadata entry
            analyte = meta['Analytes']
            concentration = meta['Concentration']
            label = f"{analyte} - {concentration}"
            plt.plot(potential, current, label=label)

        plt.xlabel('Potential (V)')
        plt.ylabel('Current (µA)')
        plt.title('Overlapping Plots of Potential vs Current')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close()

    def plot_data_array_with_corrected_baseline(self, data_array, headers, parsed_metadata, save_path, selected_indices=None):
        plt.figure(figsize=(10, 6))

        if selected_indices is None:
            selected_indices = range(0, len(headers), 2)

        def baseline_als(y, lam=1e5, p=0.01, niter=10):
            L = len(y)
            D = sparse.csc_matrix(np.diff(np.eye(L), 2))
            w = np.ones(L)
            for i in range(niter):
                W = sparse.spdiags(w, 0, L, L)
                Z = W + lam * D.dot(D.transpose())
                z = spsolve(Z, w * y)
                w = p * (y > z) + (1 - p) * (y < z)
            return z

        for i in selected_indices:
            potential = data_array[:, i]
            current = data_array[:, i + 1]

            # Filter out NaN values
            valid_idx = ~np.isnan(potential) & ~np.isnan(current)
            potential = potential[valid_idx]
            current = current[valid_idx]

            if len(potential) == 0 or len(current) == 0:
                print(f"Skipping plot for index {i//2} due to insufficient data")
                continue

            # Apply baseline correction
            baseline = baseline_als(current)
            corrected_current = current - baseline

            meta = parsed_metadata[i//2]  # if each pair of columns matches one metadata entry
            analyte = meta['Analytes']
            concentration = meta['Concentration']
            label = f"{analyte} - {concentration}"

            plt.plot(potential, corrected_current, label=label)

        plt.xlabel('Potential (V)')
        plt.ylabel('Current (µA)')
        plt.title('Corrected Plots of Potential vs Current')
        plt.legend(
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        fontsize='small',        # or 'small', 'medium', etc.
        ncol=2                      
        ) 
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close()

        return save_path


    def plot_concentrations_vs_mean_peak_currents(self, mean_peak_currents, base_save_path, selected_indices=None):
        temp_dir = base_save_path
        plot_filenames = []

        analytes_dict = defaultdict(list)  # defaultdict is now imported
        for key, mean_peak_current in mean_peak_currents.items():
            concentration, analytes = key
            analytes_dict[analytes].append((float(concentration[:-2]), mean_peak_current))

        for analytes, data in analytes_dict.items():
            data.sort(key=lambda x: x[0])
            concentrations, mean_peak_currents = zip(*data)
            plt.figure(figsize=(10, 6))
            for i, mean_peak_current in enumerate(zip(*mean_peak_currents)):
                if selected_indices is None or i in selected_indices:
                    label_str = f"{analytes} - Peak {i+1}"
                    plt.scatter(concentrations, mean_peak_current_col, label=label_str)
            plt.xscale('log')
            plt.xlabel('Concentration (uM)')
            plt.ylabel('Mean Peak Current')
            plt.title(f'Mean Peak Currents for Analytes: {analytes}')
            plt.legend()
            plt.tight_layout()

            save_path = os.path.join(temp_dir, f"{analytes}.png")
            plt.savefig(save_path)
            plt.close()
            plot_filenames.append(save_path)

        return plot_filenames

    def analyze_peak_currents(self, mean_peak_currents, std_peak_currents, base_save_path, selected_indices=None):
        # Ensure std_peak_currents is a dictionary
        if not isinstance(std_peak_currents, dict):
            raise TypeError(f"Expected std_peak_currents to be a dictionary, but got {type(std_peak_currents)} instead.")

        temp_dir = base_save_path
        plot_filenames = []

        analytes_dict = defaultdict(list)
        for key, mean_peak_current in mean_peak_currents.items():
            concentration, analytes = key  # Key is a tuple (Concentration, Analytes)

            # Ensure std_peak_currents is accessed with the correct key
            if key not in std_peak_currents:
                std_error = [0] * len(mean_peak_current)  # Default to zeros if key is missing
            else:
                std_error = std_peak_currents[key]

            analytes_dict[analytes].append((float(concentration[:-2]), mean_peak_current, std_error))

        for analytes, data in analytes_dict.items():
            data.sort(key=lambda x: x[0])  # Sort by concentration
            concentrations, mean_peak_currents, std_errors = zip(*data)
            plt.figure(figsize=(10, 6))

            for i, mean_peak_current in enumerate(zip(*mean_peak_currents)):
                if selected_indices is None or i in selected_indices:
                    error = [std[i] for std in std_errors] if std_errors else None
                    # Scatter or errorbar
                    if error:
                        plt.errorbar(
                            concentrations,
                            mean_peak_current,
                            yerr=error,
                            fmt='o',
                            label=f"{analytes} - Peak {i+1}"
                        )
                    else:
                        plt.scatter(concentrations, mean_peak_col, label=f"{analytes} - Peak {i+1}")
                    
                    log_concentrations = np.log10(concentrations)
                    slope, intercept, r_value, p_value, std_err = linregress(log_concentrations, mean_peak_current)
                    fit_line = slope * log_concentrations + intercept
                    plt.plot(concentrations, fit_line, label=f'Fit Peak {i+1} (R²={r_value ** 2:.2f}, Slope={slope:.2f}, Intercept={intercept:.2f})')

            plt.xscale('log')
            plt.xlabel('Concentration (uM)')
            plt.ylabel('Mean Peak Current')
            plt.title(f'Mean Peak Currents and Linear Fit for Analytes: {analytes}')
            plt.legend()
            plt.tight_layout()

            save_path = os.path.join(temp_dir, f"{analytes}.png")
            plt.savefig(save_path)
            plt.close()
            plot_filenames.append(save_path)

        return plot_filenames

    def plot_observed_vs_expected_concentration(self, mean_peak_currents, base_save_path, selected_indices=None):
        """
        Plots observed concentration vs. expected concentration for each analyte.
        Adds R² to the legend labels.
        """
        plot_filenames = []
        analytes_dict = defaultdict(list)

        # Group data by analytes
        for key, mean_peaks in mean_peak_currents.items():
            concentration_str, analytes = key
            try:
                c_val = float(concentration_str[:-2])  # e.g. "50uM" -> 50.0
            except ValueError:
                continue
            # Each entry: (concentration_value, [peak1, peak2, ...])
            analytes_dict[analytes].append((c_val, mean_peaks))

        for analytes, data in analytes_dict.items():
            data.sort(key=lambda x: x[0])  # Sort by ascending concentration
            max_peaks = max(len(peaks) for (_, peaks) in data)

            plt.figure(figsize=(10, 6))

            for peak_idx in range(max_peaks):
                # Gather x_vals (expected) and y_vals (peak current) for this peak index
                x_vals, y_vals = [], []
                for (conc, peaks) in data:
                    if peak_idx < len(peaks):
                        x_vals.append(conc)
                        y_vals.append(peaks[peak_idx])

                # Must have at least 2 data points to do regression
                if len(x_vals) < 2:
                    continue

                # Fit: i = m*C + b
                slope, intercept, r_val, p_val, std_err = linregress(x_vals, y_vals)
                # Observed concentration
                c_obs = [(y - intercept) / slope for y in y_vals]

                # Include R² in label
                label_str = (
                    f"{analytes} - Peak {peak_idx+1}\n"
                    f"R²={r_val**2:.3f}, N={len(x_vals)}"
                )
                plt.scatter(x_vals, c_obs, label=label_str)

            # Identity line
            all_concs = [c for (c, peaks) in data]
            if all_concs:
                min_val, max_val = min(all_concs), max(all_concs)
                plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.4)

            plt.xlabel("Expected Concentration (uM)")
            plt.ylabel("Observed Concentration (uM)")
            plt.title(f"Observed vs Expected: {analytes}")
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()

            filename = os.path.join(base_save_path, f"Obs_vs_Exp_{analytes}.png")
            plt.savefig(filename)
            plt.close()

            plot_filenames.append(filename)

        return plot_filenames
