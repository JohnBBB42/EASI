##############################
# Part 1: Imports & Basic Setup
##############################
import streamlit as st
import pandas as pd
import numpy as np
import os
import sys
import tempfile
from collections import defaultdict
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.stats import linregress
from PIL import Image
import matplotlib.pyplot as plt
import streamlit as st

from electrochemical_methods.Analysis_DPV import Analysis_DPV
from electrochemical_methods.Analysis_CV import Analysis_CV
from electrochemical_methods.Analysis_EIS import Analysis_EIS

def run_cv_analysis(df):
    """
    Mimics your old CVApp run_analysis logic, but uses Streamlit
    instead of Tkinter prompts. 
    """
    st.write("## CV Analysis")

    # 1) Let user pick how many lines to skip, columns, etc.
    st.write("### Step 1: Configure Analysis")
    values_row_start = st.number_input("values_row_start (skip lines)", value=2)
    potential_column = st.number_input("potential_column (1-based)", value=1)
    current_column   = st.number_input("current_column (1-based)", value=2)
    scan_column      = st.number_input("scan_column (1-based, 0=none)", value=0)
    scan_number      = st.number_input("scan_number", value=1)
    linreg_start_idx = st.number_input("linreg_start_index", value=15)
    r2_value         = st.number_input("r2_accept_value", value=0.90)
    pot_unit         = st.text_input("potential_unit", value="V")
    cur_unit         = st.text_input("current_unit", value="A")
    num_decimals     = st.number_input("num_decimals", value=3)

    # Let user pick a folder path or name for saving (but with Streamlit, 
    # you might do local saving or skip. We'll do a text box as example).
    saving_folder = st.text_input("Saving Folder Path", value=".")

    if st.button("Run CV Analysis"):
        # 2) Call Analysis_CV with these parameters
        try:
            Analysis_CV(
                df=df,
                values_row_start=values_row_start,
                potential_column=potential_column,
                current_column=current_column,
                scan_column=scan_column,
                scan_number=scan_number,
                linreg_start_index=linreg_start_idx,
                r2_accept_value=r2_value,
                potential_unit=pot_unit,
                current_unit=cur_unit,
                num_decimals=num_decimals,
                saving_folder=saving_folder
            )
            st.success("CV analysis completed successfully!")
            st.write(f"Check your output files in: {saving_folder}")
        except Exception as e:
            st.error(f"Error during CV Analysis: {e}")

def run_eis_analysis(df):
    """
    Mimics your old EISApp run_analysis logic, but uses Streamlit.
    """
    st.write("## EIS Analysis")

    # Basic parameters, as your old code uses:
    values_row_start = st.number_input("values_row_start", value=1)
    real_col         = st.number_input("real_col (1-based)", value=3)
    imag_col         = st.number_input("imag_col (1-based)", value=4)
    x_start          = st.text_input("x_start", "")
    x_end            = st.text_input("x_end", "")
    y_start          = st.text_input("y_start", "")
    y_end            = st.text_input("y_end", "")
    unit             = st.text_input("impedance_unit", value="Ω")
    circle_pt1       = st.number_input("circle_point1_index", value=0)
    circle_pt2       = st.number_input("circle_point2_index", value=0)

    saving_folder = st.text_input("Saving Folder Path", value=".")

    if st.button("Run EIS Analysis"):
        # Convert x_start, x_end, y_start, y_end to floats or None
        def parse_or_none(s):
            s = s.strip()
            if not s:
                return None
            try:
                return float(s)
            except:
                return None

        x_s = parse_or_none(x_start)
        x_e = parse_or_none(x_end)
        y_s = parse_or_none(y_start)
        y_e = parse_or_none(y_end)

        try:
            result = Analysis_EIS(
                df=df,
                values_row_start=values_row_start,
                real_col=real_col,
                imag_col=imag_col,
                x_start=x_s,
                x_end=x_e,
                y_start=y_s,
                y_end=y_e,
                unit=unit,
                circle_pt1_index=circle_pt1,
                circle_pt2_index=circle_pt2,
                saving_folder=saving_folder
            )
            st.success("EIS Analysis completed successfully!")
            if isinstance(result, dict) and "plot_path" in result:
                st.write(f"See your results in: {result['plot_path']}")
            else:
                st.write("Check your output folder for saved PDF plots.")
        except Exception as e:
            st.error(f"Error during EIS Analysis: {e}")


##############################################
# streamlit_dpv.py
# A complete refactor of your DPVApp code
# for a Streamlit-based interface
##############################################

import os
import io
import tempfile
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.stats import linregress
from PIL import Image

##############################################
# 1) If your old code used: 
# from electrochemical_methods.Analysis_DPV import Analysis_DPV
# ensure you can import that function now. 
# This must be updated to accept in-memory data, not file dialogs.
##############################################
from electrochemical_methods.Analysis_DPV import Analysis_DPV

##############################################
# 2) Helper Functions for Plots
#    (Directly from your DPVApp code, adapted for in-memory usage)
##############################################

def baseline_als(y, lam=1e5, p=0.01, niter=10):
    """
    Used by 'plot_data_array_with_corrected_baseline'
    to do baseline correction with an ALS approach.
    """
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for _ in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = np.linalg.solve(Z.toarray(), (w * y))
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def plot_data_array_with_corrected_baseline(data_array, headers, parsed_metadata, selected_indices):
    """
    Equivalent to 'plot_data_array_with_corrected_baseline' from DPVApp, 
    but returns a matplotlib figure we can display via Streamlit.
    """
    fig, ax = plt.subplots(figsize=(10,6))

    if selected_indices is None:
        selected_indices = range(0, len(headers), 2)

    for i in selected_indices:
        potential = data_array[:, i]
        current   = data_array[:, i+1]
        valid_idx = ~np.isnan(potential) & ~np.isnan(current)
        potential = potential[valid_idx]
        current   = current[valid_idx]
        if len(potential) == 0 or len(current) == 0:
            continue

        base = baseline_als(current)
        corrected = current - base

        # Suppose each pair of columns has metadata...
        meta        = parsed_metadata[i//2]
        analyte     = meta.get('Analytes', 'Unknown')
        concentration = meta.get('Concentration', '?')
        label       = f"{analyte} - {concentration}"

        ax.plot(potential, corrected, label=label)

    ax.set_xlabel("Potential (V)")
    ax.set_ylabel("Current (µA)")
    ax.set_title("Corrected DPV Data")
    ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')
    fig.tight_layout()
    return fig


def analyze_peak_currents(mean_peak_currents, std_peak_currents, selected_indices=None):
    """
    Equivalent to 'analyze_peak_currents' from DPVApp.
    Returns a list of matplotlib figures if multiple analytes.
    """
    figures = []
    analytes_dict = defaultdict(list)
    for key, mean_peaks in mean_peak_currents.items():
        # key is (concentration, analyte)
        concentration, analytes = key
        # we assume std_peak_currents is a dict with the same keys
        stds = std_peak_currents.get(key, [0]*len(mean_peaks))
        analytes_dict[analytes].append((float(concentration[:-2]), mean_peaks, stds))

    for analyte, data in analytes_dict.items():
        data.sort(key=lambda x: x[0])  # sort by concentration
        concentrations, mean_peaks_list, stds_list = zip(*data) # each is a tuple
        # might have multiple peaks
        # build a figure
        fig, ax = plt.subplots(figsize=(10,6))
        # we do a log scale in x
        for peak_idx in range(len(mean_peaks_list[0])):
            # gather y-values for this peak_idx across all data
            y_vals = [mp[peak_idx] for mp in mean_peaks_list]
            errs   = [sd[peak_idx] for sd in stds_list]
            # optional selected_indices check
            if (selected_indices is not None) and (peak_idx*2 not in selected_indices):
                continue

            ax.errorbar(concentrations, y_vals, yerr=errs, fmt='o', label=f"Peak {peak_idx+1}")

            # do linear regression on log10(conc)
            if len(concentrations) >= 2:
                logc = np.log10(concentrations)
                slope, intercept, r_val, _, _ = linregress(logc, y_vals)
                fit_line = slope * logc + intercept
                # rebuild x for plotting
                ax.plot(concentrations, fit_line, label=f"Fit {peak_idx+1} (R²={r_val**2:.2f})")

        ax.set_xscale('log')
        ax.set_xlabel("Concentration (uM)")
        ax.set_ylabel("Mean Peak Current")
        ax.set_title(f"Analyte: {analyte}")
        ax.legend()
        fig.tight_layout()
        figures.append(fig)
    return figures


def plot_observed_vs_expected(mean_peak_currents, selected_indices=None):
    """
    Similar to 'plot_observed_vs_expected_concentration'.
    Returns a list of figures, one per analyte grouping.
    """
    figures = []
    analytes_dict = defaultdict(list)
    for key, peaks in mean_peak_currents.items():
        concentration_str, analyte = key
        try:
            c_val = float(concentration_str[:-2])
        except:
            continue
        analytes_dict[analyte].append((c_val, peaks))

    for analyte, data in analytes_dict.items():
        data.sort(key=lambda x: x[0])
        max_peaks = max(len(x[1]) for x in data)
        fig, ax = plt.subplots(figsize=(10,6))
        for peak_idx in range(max_peaks):
            x_vals, y_vals = [], []
            for (conc, arr) in data:
                if peak_idx < len(arr):
                    x_vals.append(conc)
                    y_vals.append(arr[peak_idx])
            if len(x_vals)<2:
                continue
            # regression
            slope, intercept, r_val, _, _ = linregress(x_vals, y_vals)
            c_obs = [(y - intercept)/slope for y in y_vals]
            label_str = f"Peak {peak_idx+1} (R²={r_val**2:.3f})"
            ax.scatter(x_vals, c_obs, label=label_str)

        if data:
            all_concs = [x[0] for x in data]
            min_val, max_val = min(all_concs), max(all_concs)
            ax.plot([min_val, max_val], [min_val, max_val], 'k--')

        ax.set_xlabel("Expected Conc (uM)")
        ax.set_ylabel("Observed Conc (uM)")
        ax.set_title(f"Analyte: {analyte}")
        ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')
        fig.tight_layout()
        figures.append(fig)
    return figures

def perform_lod_analysis(analysis_results: dict):
    """
    For 'LOD Analysis' from your code. 
    We'll just read results['lod_results'] or so and display it.
    """
    lod_info = analysis_results.get('lod_results', {})
    return lod_info

def perform_t_test_analysis(analysis_results: dict):
    """
    For 'T-test Analysis'.
    """
    t_tests = analysis_results.get('t_test_results', [])
    return t_tests

##############################################
# 3) Main Streamlit Flow 
#    Replacing the entire DPVApp logic
##############################################

def main():
    st.title("DPV Analysis (Streamlit Refactor)")

    st.write("## Step 1: Upload Main DPV File")
    dpv_file = st.file_uploader("DPV File (CSV/Excel)", type=["csv","xlsx","xls"])
    main_df = None

    if dpv_file:
        # Let user choose if first row is header
        header_choice = st.checkbox("First row is header?", value=False)
        try:
            if dpv_file.name.lower().endswith(".csv"):
                main_df = pd.read_csv(dpv_file, header=0 if header_choice else None)
            else:
                main_df = pd.read_excel(dpv_file, header=0 if header_choice else None)
            st.write("Preview of DPV data:")
            st.dataframe(main_df.head())
        except Exception as e:
            st.error(f"Could not read file: {e}")
            main_df = None

    st.write("## Step 2: (Optional) Upload Blank File for LOD")
    blank_file = st.file_uploader("Blank File (CSV/Excel)", type=["csv","xlsx","xls"])
    blank_df   = None
    if blank_file:
        # Let user choose if first row is header
        blank_header = st.checkbox("Blank file has header row?", value=False, key="blank_header_box")
        try:
            if blank_file.name.lower().endswith(".csv"):
                blank_df = pd.read_csv(blank_file, header=0 if blank_header else None)
            else:
                blank_df = pd.read_excel(blank_file, header=0 if blank_header else None)
            st.write("Preview of Blank data:")
            st.dataframe(blank_df.head())
        except Exception as e:
            st.error(f"Could not read blank file: {e}")
            blank_df = None

    st.write("## Step 3: Concentration Conversion (Optional)")
    unit_choices = ["", "nM", "µM", "mM", "M", "g/L", "mg/L"]
    from_unit = st.selectbox("From Unit:", unit_choices)
    to_unit   = st.selectbox("To Unit:",   unit_choices)

    st.write("## Step 4: Choose Analysis Option")
    analysis_options = [
        "DPV analysis",
        "Plot Data Array with Corrected Baseline",
        "Analyze Peak Currents",
        "Observed vs Expected Concentration",
        "LOD Analysis",
        "T-test Analysis"
    ]
    choice = st.selectbox("Analysis Type", analysis_options)

    # We'll store the analysis results dict here
    analysis_results = None

    if st.button("Run DPV Analysis"):
        if main_df is None or main_df.empty:
            st.error("Please upload a valid main DPV file first.")
            return

        # 1) Actually run the DPV analysis logic
        st.info("Running 'Analysis_DPV' on the uploaded DataFrame ...")
        # We assume 'Analysis_DPV' can handle a DataFrame. If not, adapt it.
        # or do Analysis_DPV(filepath=some_tempfile, blank_responses=some_tempfile).
        # For demonstration, let's do:
        analysis_results = Analysis_DPV(df=main_df, blank_responses=blank_df)

        if not analysis_results:
            st.warning("No analysis results returned.")
            return

        # 2) Possibly do concentration conversion
        if from_unit and to_unit and from_unit != to_unit:
            st.info(f"Converting concentrations from {from_unit} to {to_unit} ...")
            # A direct adaptation of your old convert_concentrations
            convert_concentrations_in_results(analysis_results, from_unit, to_unit)
            st.write("Concentrations updated in analysis_results['parsed_metadata']")

        # 3) Based on the 'choice', do sub-analyses or plotting
        st.write(f"### Chosen Analysis: {choice}")
        if choice == "DPV analysis":
            st.write("Raw DPV Analysis Results (JSON):")
            st.json(analysis_results)

        elif choice == "Plot Data Array with Corrected Baseline":
            # Suppose analysis_results has 'data_array', 'headers', 'parsed_metadata'
            data_array = analysis_results.get('data_array', None)
            headers    = analysis_results.get('headers', [])
            parsed_md  = analysis_results.get('parsed_metadata', [])
            if data_array is None:
                st.error("No 'data_array' found in results.")
            else:
                # let user pick which columns to plot
                num_plots = len(headers)//2
                # checkboxes
                st.write("Select which plots to correct (in pairs of potential/current columns):")
                checks = []
                for i in range(num_plots):
                    c = st.checkbox(f"Plot Pair {i+1}", value=True)
                    checks.append(c)
                selected_indices = [i*2 for i, c in enumerate(checks) if c]
                fig = plot_data_array_with_corrected_baseline(data_array, headers, parsed_md, selected_indices)
                st.pyplot(fig)

        elif choice == "Analyze Peak Currents":
            # We read 'mean_peak_currents' and 'std_peak_currents' from results
            mpc = analysis_results.get('mean_peak_currents', {})
            spc = analysis_results.get('std_peak_currents', {})
            if not mpc:
                st.error("No mean_peak_currents found in results.")
            else:
                # Let user pick which peak indexes
                st.write("Analyze which peaks? We'll just do them all for now.")
                figs = analyze_peak_currents(mpc, spc)
                for f in figs:
                    st.pyplot(f)

        elif choice == "Observed vs Expected Concentration":
            mpc = analysis_results.get('mean_peak_currents', {})
            if not mpc:
                st.error("No mean_peak_currents found in results.")
            else:
                figs = plot_observed_vs_expected(mpc)
                for f in figs:
                    st.pyplot(f)

        elif choice == "LOD Analysis":
            lod = perform_lod_analysis(analysis_results)
            st.write("LOD Analysis Results:")
            st.json(lod)

        elif choice == "T-test Analysis":
            t_tests = perform_t_test_analysis(analysis_results)
            st.write("T-Test Analysis Results:")
            st.json(t_tests)

        st.success("DPV analysis step completed.")

    # Optionally let user "Save Results" 
    if st.button("Save Results"):
        st.warning("Implement your logic to save 'analysis_results' to disk or let user download. "
                   "Because we're in Streamlit, you might use 'st.download_button'.")


##############################################
# 4) Launch if run as main
##############################################
if __name__ == "__main__":
    main()

