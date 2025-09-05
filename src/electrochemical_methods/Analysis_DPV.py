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

def Analysis_DPV(file_path, blank_responses=None, metadata_fields=None, peak_regions=None):
    #_________________________________________
    # DEFINE FUNCTIONS
    #_________________________________________

    # ─────────────────────────────────────────────────────────────────────────────
    # FLEXIBLE CSV / .pssession LOADER  (keeps dots *and* commas in concentrations)
    # ─────────────────────────────────────────────────────────────────────────────
    def load_data_and_metadata(
            file_path,
            metadata_fields=None,             # what the user typed in the GUI
            *,                                # keyword‑only from here down
            encoding          ="utf‑16",
            metadata_row_idx  =3,
            metadata_delim    =",,",
    ):
        """
        Return
            headers        – the column header row
            data_array     – numpy array with all numeric data
            parsed_metadata – list[dict] 1‑per‑curve

        • Concentrations such as “29.60µM*”, “29,50uM”, “50mM” are preserved
        verbatim (commas are converted to dots later by numeric_conc()).
        • Works for any ordering of the four tokens.
        """
        import csv, re, numpy as np
        from collections import defaultdict

        # ------------------------------------------------------------------ helpers
        unit_pat  = r"(?:nM|uM|µM|mM|M|(?:g|mg|µg)\/?L)"
        conc_rx   = re.compile(rf"(\d+(?:[.,]\d+)?)\s*{unit_pat}\*?", re.I)
        analy_rx  = re.compile(r"^(HX|UA|XAN|AA|DA)$", re.I)   # extend at will
        method_rx = re.compile(r"DPV", re.I)

        def classify(tok: str):
            """return ('Concentration', cleaned) | ('Analytes', cleaned) | … | None"""
            if conc_rx.fullmatch(tok):
                return "Concentration", tok.replace(",", ".")      # keep dot‑decimal
            if analy_rx.match(tok):
                return "Analytes", tok.upper()
            if method_rx.search(tok):
                return "Method", "DPV"
            return None, None

        def safe_float(x):
            try:
                return float(x.replace(",", "."))  # accept both 0,5 and 0.5
            except ValueError:
                return np.nan

        # ------------------------------------------------- 0) pull the tiny CSV part
        with open(file_path, "r", encoding=encoding, newline="") as fh:
            lines = fh.readlines()

        meta_raw = lines[metadata_row_idx].strip()
        meta_raw = re.sub(r"^Metadata row:\s*", "", meta_raw)
        raw_entries = [e.strip() for e in meta_raw.split(metadata_delim) if e.strip()]

        if not metadata_fields:
            metadata_fields = ["Electrode", "Analytes", "Concentration", "Method"]

        # ------------------------------------------------- 1) build the token table
        parsed = []                                                  # list[dict]
        for entry in raw_entries:
            tokens = re.split(r"[._\-\s]+", entry)
            rec    = defaultdict(str)
            if tokens:
                rec["Electrode"] = tokens[0]                         # always first

            for tok in tokens[1:]:
                k, v = classify(tok)
                if k and not rec[k]:
                    rec[k] = v

            rec.setdefault("Method", "DPV")                          # safety net
            parsed.append({f: rec.get(f, "") for f in metadata_fields})

        # ------------------------------------------------- 2) numeric data section
        with open(file_path, "r", encoding=encoding, newline="") as fh:
            sample = fh.read(2048)
            fh.seek(0)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=[",", ";", "\t"])
            except csv.Error:
                dialect          = csv.excel()
                dialect.delimiter = ","

            rows = list(csv.reader(fh, dialect))

        hdr_row   = metadata_row_idx + 2                # blank + “Date …” rows
        headers   = rows[hdr_row]
        data_rows = rows[hdr_row + 1:]

        data_array = np.array([
            [safe_float(c) for c in row]
            for row in data_rows
            if len(row) == len(headers)
        ])

        return headers, data_array, parsed



    def find_peaks_and_process_data(headers, data_array, parsed_metadata, peak_regions=None):
        #voltage_columns = [i for i, col in enumerate(headers) if 'V' in col]
        #current_columns = [i for i, col in enumerate(headers) if 'µA' in col]

        voltage_columns = [i for i,col in enumerate(headers) if col.upper().startswith("V")]
        current_columns = [i for i,col in enumerate(headers) if "µA" in col]

        # ─── SAFETY NET (do not crash if counts mismatch) ────────────────────
        orig_meta_len = len(parsed_metadata)
        n_pairs = min(len(voltage_columns), len(current_columns), orig_meta_len)
        voltage_columns = voltage_columns[:n_pairs]
        current_columns = current_columns[:n_pairs]
        parsed_metadata = parsed_metadata[:n_pairs]
        if n_pairs < orig_meta_len:
            print(f"⚠️  Truncated metadata[{orig_meta_len}] → [{n_pairs}] column-pairs")
        # ────────────────────────────────────────────────────────────────────

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

        # ----------------------------------------------------------------------
        # (A) helper : FWHM on *any* 1-D signal (baseline already removed)
        # ------------------------------------------------------------------
        def fwhm(x, y, peak_idx):
            """
            Return FWHM in x-units for the peak at peak_idx.
            y should already be baseline-corrected (so half-max = ½·y[p]).
            Uses linear interpolation between the two samples that straddle
            the half-height on each side.
            """
            half = 0.5 * y[peak_idx]

            # --- left side -------------------------------------------------
            i_left = peak_idx
            while i_left > 0 and y[i_left] > half:
                i_left -= 1
            if i_left == 0:
                return np.nan                  # clipped at start
            # linear interpolate between (i_left, i_left+1)
            x_left = np.interp(
                half,
                [y[i_left], y[i_left+1]],
                [x[i_left], x[i_left+1]]
            )

            # --- right side ------------------------------------------------
            i_right = peak_idx
            while i_right < len(y)-1 and y[i_right] > half:
                i_right += 1
            if i_right == len(y)-1:
                return np.nan                  # clipped at end
            x_right = np.interp(
                half,
                [y[i_right-1], y[i_right]],
                [x[i_right-1], x[i_right]]
            )

            return abs(x_right - x_left)
        # ---------- end helper --------------------------------------------

        results = []
        # ------------------------------------------------------------------
        # iterate one voltammogram / metadata entry at a time
        # ------------------------------------------------------------------
        for index, row in enumerate(parsed_metadata):

            # --- 1. grab X/Y arrays -------------------------------------------------
            voltage_data = data_array[:, voltage_columns[index]]
            current_data = data_array[:, current_columns[index]]

            # drop NaNs
            keep = ~np.isnan(voltage_data) & ~np.isnan(current_data)
            voltage_data = voltage_data[keep]
            current_data = current_data[keep]

            # --- 2. peak finding & baseline corr. -----------------------------------
            peaks, smoothed_current = find_peaks_in_data(voltage_data, current_data)

            baseline         = baseline_als(smoothed_current)
            adjusted_current = smoothed_current - baseline

            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            # KEEP A COPY **before** any window-mask so FWHM has full data
            orig_V    = voltage_data.copy()
            orig_adj  = adjusted_current.copy()
            # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            # ------------------------------------------------------------------
            # apply user-defined peak window(s) if provided
            # ------------------------------------------------------------------
            my_windows = None
            if peak_regions:
                for key, win in peak_regions.items():        # let "UA" match "HX,Xan,UA"
                    if key.lower() in row["Analytes"].lower():
                        my_windows = win
                        break

            if my_windows:
                peaks = np.array([
                    p for p in peaks
                    if any(lo <= orig_V[p] <= hi for lo, hi in my_windows)
                ])

            # ------------------------------------------------------------------
            # write out results  (including NEW FWHM column)
            # ------------------------------------------------------------------
            for peak_idx in peaks:
                peak_voltage         = float(orig_V[peak_idx])
                peak_current         = float(current_data[peak_idx])
                adj_peak_current     = float(orig_adj[peak_idx])

                # ---- NEW: FWHM computed on full (uncropped) trace --------------
                this_fwhm = fwhm(orig_V, orig_adj, peak_idx)

                # ───────── DEBUG HELP ───────────────────────────────────────────
                if "UA" in row["Analytes"]:
                    if np.isnan(this_fwhm):
                        print(f"[DEBUG-UA]  peak@{peak_voltage:.3f} V   FWHM = NaN  "
                            f"(left or right half-height was clipped)")
                    else:
                        print(f"[DEBUG-UA]  peak@{peak_voltage:.3f} V   FWHM = {this_fwhm:.4f} V")
                # ────────────────────────────────────────────────────────────────

                results.append({
                    'Entry'                : index,
                    'Electrode'            : row['Electrode'],
                    'Analytes'             : row['Analytes'],
                    'Method'               : row['Method'],
                    'Concentration'        : row['Concentration'],
                    'Peak Voltage'         : peak_voltage,
                    'Peak Current'         : peak_current,
                    'Adjusted Peak Current': adj_peak_current,
                    'Y Offset'             : peak_current - adj_peak_current,
                    'Baseline at Peak'     : float(baseline[peak_idx]),
                    'FWHM'                 : this_fwhm,        # ← new column
                })
        # ----------------------------------------------------------------------
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
        # 1) Build analyte → { concentration_value : { 'mean':…, 'std':… } }
        analytes_dict = defaultdict(lambda: defaultdict(dict))
        number_re     = re.compile(r"([0-9]*\.?[0-9]+)")

        # populate 'mean'
        for (conc_str, analytes), means in mean_peak_currents.items():
            m = number_re.match(conc_str)
            if not m:
                continue                        # skip empty or non-numeric
            conc = float(m.group(1))
            for i, analyte in enumerate(analytes.split(',')):
                if i >= len(means):
                    break
                name = analyte.strip()
                # if there were multi-analyte entries you prefixed them with '3'
                if ',' in analytes:
                    name = '3' + name
                analytes_dict[name][conc]['mean'] = means[i]

        # populate 'std'
        for (conc_str, analytes), stds in std_peak_currents.items():
            m = number_re.match(conc_str)
            if not m:
                continue
            conc = float(m.group(1))
            for i, analyte in enumerate(analytes.split(',')):
                if i >= len(stds):
                    break
                name = analyte.strip()
                if ',' in analytes:
                    name = '3' + name
                analytes_dict[name][conc]['std'] = stds[i]

        # 2) Now do pairwise t-tests
        t_test_results = []
        for analyte, conc_dict in analytes_dict.items():
            # need at least two concentrations to compare
            concentrations = sorted(conc_dict.keys())
            if len(concentrations) < 2:
                continue

            n1 = n2 = sample_size or 3
            for c1, c2 in combinations(concentrations, 2):
                mean1 = conc_dict[c1].get('mean', 0)
                mean2 = conc_dict[c2].get('mean', 0)
                std1  = conc_dict[c1].get('std',  0)
                std2  = conc_dict[c2].get('std',  0)

                # pooled standard deviation
                pooled = np.sqrt(((n1 - 1)*std1**2 + (n2 - 1)*std2**2) / (n1 + n2 - 2)) \
                        if (n1 + n2 - 2) > 0 else 0

                t_stat = 0.0
                if pooled > 0:
                    t_stat = (mean1 - mean2) / (pooled * np.sqrt(1/n1 + 1/n2))

                df = n1 + n2 - 2
                crit = t.ppf(1 - 0.025, df) if df > 0 else np.nan

                t_test_results.append({
                    'Analyte'       : analyte,      # use analyte_name
                    'Concentration1': c1,
                    'Concentration2': c2,
                    'T-statistic'   : t_stat,
                    'Critical value': crit,
                    'Significant'   : abs(t_stat) > crit if not np.isnan(crit) else False
                })

        return t_test_results

    # Main execution using the new loader:
    headers, data_array, parsed_metadata = load_data_and_metadata(file_path, metadata_fields=metadata_fields)

    if headers is None or data_array is None or parsed_metadata is None:
        return None
    else:
        results = find_peaks_and_process_data(headers, data_array, parsed_metadata, peak_regions)
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