##############################
# Part 1: Imports & Basic Setup
##############################
import streamlit as st
import pandas as pd
import numpy as np
import os
import re
import io
import csv
import sys
import tempfile
import matplotlib.ticker as mtick
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from collections import defaultdict
from itertools import combinations
from scipy.stats import t, linregress, ttest_ind
from scipy.signal import find_peaks, savgol_filter
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.ndimage import gaussian_filter1d
from PIL import Image
import matplotlib.pyplot as plt
import zipfile


# Import your analysis functions:
from electrochemical_methods.Analysis_DPV import Analysis_DPV
from electrochemical_methods.Analysis_CV import Analysis_CV
from electrochemical_methods.Analysis_EIS import Analysis_EIS


def open_plot_selection_modal(analysis_option, headers, data_array):
    """
    Opens a modal dialog that mimics the Tkinter plot selection window.
    Displays a "Select All" checkbox, individual plot checkboxes,
    and a button to generate the plots.
    """
    with st.expander("Select Plots"):
        st.write("Select the plots you want to display:")
        # Number of plots (assume one plot per two header columns)
        num_plots = len(range(0, len(headers), 2))
        # "Select All" checkbox
        select_all = st.checkbox("Select All", value=True, key="select_all_plots_modal")
        # Create checkboxes for each plot and record in a dict
        plot_selection = {}
        for i in range(num_plots):
            plot_selection[i] = st.checkbox(f"Plot {i+1}", value=select_all, key=f"plot_select_modal_{i}")
        # When user clicks "Plot Selected" within the modal:
        if st.button("Plot Selected", key="modal_plot_button"):
            selected_indices = [i * 2 for i, checked in plot_selection.items() if checked]
            st.write(f"Selected plot indices: {selected_indices}")
            # Call the appropriate plotting function based on analysis_option.
            if analysis_option == "Plot Data Array with Corrected Baseline":
                fig = plot_data_array_with_corrected_baseline(data_array, headers, selected_indices)
            elif analysis_option == "Analyze Peak Currents":
                fig = plot_peak_currents_modal(data_array, headers, selected_indices)
            else:
                st.write("Plot selection for this option not implemented.")
                return
            # Display the figure in a modal:
            with st.expander("Plot Output"):
                st.pyplot(fig)
            # Also, save the figure in session_state for exporting.
            buf = io.BytesIO()
            fig.savefig(buf, format="png")
            buf.seek(0)
            st.session_state.plot_images = st.session_state.get("plot_images", []) + [("Plot_output.png", buf.getvalue())]

# Example plotting function for baseline-corrected data.
def plot_data_array_with_corrected_baseline(data_array, headers, selected_indices=None):
    if selected_indices is None:
        selected_indices = range(0, len(headers), 2)
    fig, ax = plt.subplots(figsize=(8,6))
    for i in selected_indices:
        pot = data_array[:, i]
        cur = data_array[:, i+1]
        # Filter NaNs
        mask = ~np.isnan(pot) & ~np.isnan(cur)
        pot = pot[mask]
        cur = cur[mask]
        if pot.size == 0 or cur.size == 0:
            continue
        # Baseline correction using ALS (as in DPVApp)
        L = len(cur)
        D = sparse.csc_matrix(np.diff(np.eye(L), 2))
        w = np.ones(L)
        lam, p, niter = 1e5, 0.01, 10
        for _ in range(niter):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.T)
            baseline = spsolve(Z, w * cur)
            w = p * (cur > baseline) + (1 - p) * (cur < baseline)
        corrected = cur - baseline
        label = f"Entry {i//2 + 1}"
        ax.plot(pot, corrected, label=label)
    ax.set_xlabel("Potential (V)")
    ax.set_ylabel("Corrected Current (µA)")
    ax.set_title("Baseline-Corrected DPV Data")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    return fig

# Dummy example for peak currents plot; you would replace with your actual logic.
def plot_peak_currents_modal(data_array, headers, selected_indices=None):
    fig, ax = plt.subplots(figsize=(8,6))
    # Dummy implementation: plot random bar chart for demonstration.
    if selected_indices is None or not selected_indices:
        st.warning("No plot indices selected.")
        return fig
    values = [np.random.uniform(0.1, 10) for _ in selected_indices]
    labels = [f"Peak {i//2+1}" for i in selected_indices]
    ax.bar(labels, values, color="cyan")
    ax.set_ylabel("Peak Current (µA)")
    ax.set_title("Peak Current Analysis")
    return fig


##############################################
# CV Analysis
##############################################
def run_cv_analysis(df):
    st.write("## CV Analysis")
    st.write("### Step 1: Configure Analysis")
    values_row_start = st.number_input("values_row_start (skip lines)", value=2)
    potential_column = st.number_input("potential_column (1-based)", value=1)
    current_column   = st.number_input("current_column (1-based)", value=3)
    scan_column      = st.number_input("scan_column (1-based, 0=none)", value=5)
    scan_number      = st.number_input("scan_number", value=1)
    linreg_start_idx = st.number_input("linreg_start_index", value=15)
    r2_value         = st.number_input("r2_accept_value", value=0.90)
    pot_unit         = st.text_input("potential_unit", value="V")
    cur_unit         = st.text_input("current_unit", value="A")
    num_decimals     = st.number_input("num_decimals", value=3)
    saving_folder    = st.text_input("Saving Folder Path", value=".")
    
    if st.button("Run CV Analysis"):
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

##############################################
# EIS Analysis
##############################################
def run_eis_analysis(df):
    st.write("## EIS Analysis")
    
    # We assume that the Real (Z') and Imag (-Z'') data are in the 3rd and 4th columns (1-based indexing)
    real_idx = 3 - 1  # zero-based index for column 3
    imag_idx = 4 - 1  # zero-based index for column 4
    
    try:
        # Check if the first entry in the real column is non-numeric (likely a header)
        first_val_real = df.iloc[0, real_idx]
        # Replace commas with dots and remove at most one dot for a proper float test.
        first_val_str = str(first_val_real).replace(',', '.')
        # isdigit() is not enough (for floats), so we try to convert it
        try:
            float(first_val_str)
        except ValueError:
            st.info("Detected header text in data. Dropping the first row.")
            df = df.iloc[1:].reset_index(drop=True)
        
        # Now convert the designated columns to numeric
        df.iloc[:, real_idx] = pd.to_numeric(df.iloc[:, real_idx].astype(str).str.replace(',', '.'), errors='raise')
        df.iloc[:, imag_idx] = pd.to_numeric(df.iloc[:, imag_idx].astype(str).str.replace(',', '.'), errors='raise')
    except Exception as conv_err:
        st.error(f"Error converting columns to numeric: {conv_err}")
        return

    # Let the user specify a folder path for saving results (optional)
    saving_folder = st.text_input("Saving Folder Path", value=".")

    # When the user clicks the button, run the EIS Analysis using your defaults
    if st.button("Run EIS Analysis"):
        try:
            result = Analysis_EIS(
                df=df,
                values_row_start=1,   # same as your Tkinter default
                real_col=3,           # 3rd column (1-based)
                imag_col=4,           # 4th column (1-based)
                x_start=None,
                x_end=None,
                y_start=None,
                y_end=None,
                unit="Ω",
                circle_pt1_index=0,
                circle_pt2_index=0,
                saving_folder=saving_folder
            )
            st.success("EIS Analysis completed successfully!")
            if isinstance(result, dict) and "plot_path" in result:
                st.write(f"Plot saved at: {result['plot_path']}")
            else:
                st.write("Check your output folder for saved PDF plots.")
        except Exception as e:
            st.error(f"Error during EIS Analysis: {e}")

##############################################
# DPV Analysis
##############################################
def run_dpv_analysis():
    st.write("## DPV Analysis")
    # --- File Upload Section (always visible) ---
    uploaded_file = st.file_uploader("Upload DPV data file", type=["csv", "xlsx", "xls", "txt"])
    blank_file = st.file_uploader("Upload Blank Responses File (Optional)", type=["csv", "xlsx", "xls", "txt"])
    
    if uploaded_file:
        tfile = tempfile.NamedTemporaryFile(delete=False, suffix="." + uploaded_file.name.split('.')[-1])
        tfile.write(uploaded_file.getvalue())
        tfile.close()
        file_path = tfile.name
        
        if blank_file:
            tbfile = tempfile.NamedTemporaryFile(delete=False, suffix="." + blank_file.name.split('.')[-1])
            tbfile.write(blank_file.getvalue())
            tbfile.close()
            blank_path = tbfile.name
        else:
            blank_path = None
        
        if st.button("Run DPV Analysis"):
            try:
                dpv_result = Analysis_DPV(file_path, blank_responses=blank_path) if blank_path else Analysis_DPV(file_path)
                st.success("DPV analysis completed successfully!")
                st.write("### DPV Analysis Results")
                if not dpv_result or not isinstance(dpv_result, dict):
                    st.warning("No structured results found. Showing raw output:")
                    st.write(dpv_result)
                    return
                # Display basic results
                if 'mean_peak_currents' in dpv_result:
                    st.markdown("**Mean Peak Currents:**")
                    for key, val in dpv_result['mean_peak_currents'].items():
                        st.write(f"{key}: {val}")
                if 'std_peak_currents' in dpv_result:
                    st.markdown("**Standard Deviation of Peak Currents:**")
                    for key, val in dpv_result['std_peak_currents'].items():
                        st.write(f"{key}: {val}")
                if 'cov_peak_currents' in dpv_result:
                    st.markdown("**Coefficient of Variation of Peak Currents:**")
                    for key, val in dpv_result['cov_peak_currents'].items():
                        st.write(f"{key}: {val}")
                # Ensure entry names exist for extended options
                if 'parsed_metadata' in dpv_result:
                    entry_names = []
                    for meta in dpv_result['parsed_metadata']:
                        name = f"{meta.get('Analytes', '')} {meta.get('Concentration', '')}".strip()
                        entry_names.append(name if name else "Unknown")
                    dpv_result['entry_names'] = entry_names
                st.session_state.analysis = dpv_result
            except Exception as e:
                st.error(f"Error during DPV Analysis: {e}")
    else:
        st.info("Please upload a DPV data file.")
    
    # --- Sidebar Extended Options (always visible) ---
    analysis_options = [
        "Plot Data Array with Corrected Baseline",
        "Analyze Peak Currents",
        "Observed vs Expected Concentration",
        "LOD Analysis",
        "T-test Analysis",
        "Concentration Conversion"
    ]
    choice = st.sidebar.selectbox("Choose extended DPV analysis:", analysis_options)
    
    # Sidebar: Entry selection checkboxes
    selected_entries = []
    if st.session_state.get("analysis") and "entry_names" in st.session_state.analysis:
        st.sidebar.markdown("**Select data entries to include:**")
        for i, entry in enumerate(st.session_state.analysis["entry_names"]):
            if st.sidebar.checkbox(entry, value=True, key=f"entry_{i}"):
                selected_entries.append(entry)
    else:
        st.sidebar.info("Upload and run base DPV analysis to see entry selection options.")
    
    # Sidebar: Button to run extended analysis
    run_ext = st.sidebar.button("Run Extended Analysis")
    
    if run_ext:
        if not st.session_state.get("analysis"):
            st.warning("Please run the base DPV analysis first!")
        else:
            if choice == "Plot Data Array with Corrected Baseline":
                # Open modal for plot selection (pop-up)
                open_plot_selection_modal(
                    analysis_option=choice,
                    headers=st.session_state.analysis.get("headers", []),
                    data_array=st.session_state.analysis.get("data_array", None)
                )
            elif choice == "Analyze Peak Currents":
                st.subheader("Peak Current Analysis")
                if selected_entries:
                    peak_data = []
                    for entry in selected_entries:
                        peak = st.session_state.analysis.get_peak_current(entry)
                        peak_data.append({"Entry": entry, "Peak Current (µA)": peak})
                    st.table(pd.DataFrame(peak_data))
                    # Additionally, show a bar chart
                    values = [d["Peak Current (µA)"] for d in peak_data]
                    labels = [d["Entry"] for d in peak_data]
                    fig, ax = plt.subplots()
                    ax.bar(labels, values, color="cyan")
                    ax.set_ylabel("Peak Current (µA)")
                    ax.set_xlabel("Entry")
                    st.pyplot(fig)
                    # Save plot image to session for export
                    buf = io.BytesIO()
                    fig.savefig(buf, format="png")
                    buf.seek(0)
                    st.session_state.plot_images = st.session_state.get("plot_images", []) + [("Peak_Currents.png", buf.getvalue())]
                else:
                    st.write("No entries selected for analysis.")
            elif choice == "Observed vs Expected Concentration":
                st.subheader("Observed vs Expected Concentration")
                data = st.session_state.analysis.get_concentration_comparison() if hasattr(st.session_state.analysis, "get_concentration_comparison") else None
                if data is not None:
                    plot_df = data[data["Entry"].isin(selected_entries)] if selected_entries else data
                    if plot_df.empty:
                        st.write("No data to display for the selected entries.")
                    else:
                        import altair as alt
                        scatter = alt.Chart(plot_df).mark_point(size=100, color='blue').encode(
                            x=alt.X('Expected', title='Expected Concentration'),
                            y=alt.Y('Observed', title='Observed Concentration'),
                            tooltip=['Entry', 'Expected', 'Observed']
                        )
                        line = alt.Chart(pd.DataFrame({'x': [plot_df["Expected"].min(), plot_df["Expected"].max()]})).mark_line(color='red').encode(x='x', y='x')
                        st.altair_chart(scatter + line, use_container_width=True)
                else:
                    st.write("Concentration comparison not available.")
            elif choice == "LOD Analysis":
                st.subheader("Limit of Detection Analysis")
                lod_value = st.session_state.analysis.get_LOD() if hasattr(st.session_state.analysis, "get_LOD") else None
                loq_value = st.session_state.analysis.get_LOQ() if hasattr(st.session_state.analysis, "get_LOQ") else None
                if lod_value is not None:
                    st.write(f"**Calculated LOD:** {lod_value:.3f} (in same units as input)")
                    if loq_value is not None:
                        st.write(f"**Calculated LOQ:** {loq_value:.3f}")
                    calib_df = st.session_state.analysis.get_calibration_curve() if hasattr(st.session_state.analysis, "get_calibration_curve") else None
                    if calib_df is not None:
                        fig, ax = plt.subplots()
                        ax.scatter(calib_df["Concentration"], calib_df["Signal"], label="Calibration data")
                        ax.axvline(lod_value, color='red', linestyle='--', label=f"LOD ~ {lod_value:.3f}")
                        ax.set_xlabel("Concentration")
                        ax.set_ylabel("Signal (peak current)")
                        ax.legend()
                        st.pyplot(fig)
                        buf = io.BytesIO()
                        fig.savefig(buf, format="png")
                        buf.seek(0)
                        st.session_state.plot_images = st.session_state.get("plot_images", []) + [("LOD_Calibration.png", buf.getvalue())]
                else:
                    st.write("LOD analysis could not be performed (insufficient data or blank not provided).")
            elif choice == "T-test Analysis":
                st.subheader("T-test Statistical Comparison")
                entries = st.session_state.analysis.get("entry_names", [])
                if len(entries) < 2:
                    st.write("Need at least two data entries for a T-test.")
                else:
                    col1, col2 = st.columns(2)
                    group1 = col1.selectbox("Select first group", entries, index=0)
                    group2 = col2.selectbox("Select second group", entries, index=1)
                    if group1 == group2:
                        st.warning("Please select two different entries.")
                    else:
                        if st.button("Run T-test"):
                            try:
                                data1 = st.session_state.analysis.get_data(group1)
                                data2 = st.session_state.analysis.get_data(group2)
                            except AttributeError:
                                st.error("Analysis object does not support direct data retrieval for t-test.")
                                return
                            from scipy import stats
                            t_stat, p_val = stats.ttest_ind(data1, data2, equal_var=False, nan_policy='omit')
                            st.write(f"t-statistic = {t_stat:.3f}, p-value = {p_val:.3f}")
                            alpha = 0.05
                            if p_val < alpha:
                                st.write(f"Result: The difference between **{group1}** and **{group2}** is statistically significant (p < {alpha}).")
                            else:
                                st.write(f"Result: No significant difference between **{group1}** and **{group2}** (p ≥ {alpha}).")
            elif choice == "Concentration Conversion":
                st.subheader("Concentration Unit Conversion")
                units = ["M", "mM", "µM", "nM", "mg/L", "µg/mL"]
                from_unit = st.selectbox("From unit:", units, index=units.index("µM"))
                to_unit = st.selectbox("To unit:", units, index=units.index("nM"))
                value = st.number_input(f"Value in {from_unit}:", value=1.0)
                if st.button("Convert"):
                    factors = {
                        "M": 1,
                        "mM": 1e-3,
                        "µM": 1e-6,
                        "nM": 1e-9,
                        "mg/L": 1e-3,
                        "µg/mL": 1e-3
                    }
                    if from_unit in factors and to_unit in factors:
                        value_in_M = value * factors[from_unit]
                        converted_value = value_in_M / factors[to_unit]
                        st.write(f"{value} {from_unit} = **{converted_value:.4g} {to_unit}**")
                    else:
                        st.write("Conversion for selected units is not defined.")

    # --- Saving Results Section ---
    if st.session_state.get("analysis"):
        with st.expander("Save / Export Results"):
            # Prepare CSV of base analysis results
            if 'mean_peak_currents' in st.session_state.analysis:
                csv_buffer = io.StringIO()
                writer = csv.writer(csv_buffer)
                writer.writerow(["Entry", "Mean Peak Current", "Std Dev", "CoV"])
                for key, mean_val in st.session_state.analysis['mean_peak_currents'].items():
                    std_val = st.session_state.analysis.get('std_peak_currents', {}).get(key, "")
                    cov_val = st.session_state.analysis.get('cov_peak_currents', {}).get(key, "")
                    writer.writerow([key, mean_val, std_val, cov_val])
                csv_data = csv_buffer.getvalue().encode('utf-8')
                st.download_button("Download Base Results (CSV)", data=csv_data, file_name="DPV_base_results.csv", mime="text/csv")
            # If extended analysis plots were saved, offer them as a ZIP
            if st.session_state.get("plot_images"):
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zipf:
                    for img_name, img_bytes in st.session_state.plot_images:
                        zipf.writestr(img_name, img_bytes)
                zip_data = zip_buffer.getvalue()
                st.download_button("Download All Extended Plots (ZIP)", data=zip_data, file_name="DPV_extended_plots.zip", mime="application/zip")



##############################################
# Main Application Interface
##############################################
def main():
    st.title("Electrochemical Analysis Software Interface (EASI)")
    st.markdown("""
    Welcome to the Electrochemical Analysis Software Interface (EASI).
    
    This tool provides various electrochemical analysis methods:
    1. Cyclic Voltammetry (CV) Analysis
    2. Differential Pulse Voltammetry (DPV) Analysis
    3. Electrochemical Impedance Spectroscopy (EIS) Analysis
    
    Please select an analysis from the dropdown below to proceed.
    """)
    
    analysis_choice = st.selectbox("Select Analysis", ["(None)", "CV Analysis", "DPV Analysis", "EIS Analysis"])
    if analysis_choice == "(None)":
        st.info("Please choose an analysis to begin.")
        return
    
    if analysis_choice == "CV Analysis":
        st.header("Cyclic Voltammetry (CV) Analysis")
        uploaded_file = st.file_uploader("Upload a CV data file", type=["csv","xlsx","xls"])
        df = None
        if uploaded_file:
            has_header = st.checkbox("Does the first row contain headers?", value=False)
            try:
                if uploaded_file.name.lower().endswith(".csv"):
                    df = pd.read_csv(uploaded_file, header=(0 if has_header else None))
                else:
                    df = pd.read_excel(uploaded_file, header=(0 if has_header else None))
                st.write("Preview of CV data:")
                st.dataframe(df.head(10))
            except Exception as e:
                st.error(f"Could not read CV file: {e}")
                df = None
        if df is not None and not df.empty:
            st.subheader("Configure & Run CV Analysis")
            run_cv_analysis(df)
        else:
            st.warning("Please upload a valid CV file.")
            
    elif analysis_choice == "DPV Analysis":
        st.header("Differential Pulse Voltammetry (DPV) Analysis")
        run_dpv_analysis()
    
    elif analysis_choice == "EIS Analysis":
        st.header("Electrochemical Impedance Spectroscopy (EIS) Analysis")
        # Upload once here
        uploaded_file = st.file_uploader("Upload an EIS data file", type=["csv","xlsx","xls"])
        df = None
        if uploaded_file:
            has_header = st.checkbox("Does the first row contain headers?", value=False)
            try:
                if uploaded_file.name.lower().endswith(".csv"):
                    df = pd.read_csv(uploaded_file, header=(0 if has_header else None))
                else:
                    df = pd.read_excel(uploaded_file, header=(0 if has_header else None))
                st.write("Preview of EIS data:")
                st.dataframe(df.head(10))
            except Exception as e:
                st.error(f"Could not read EIS file: {e}")
                df = None

        if df is not None and not df.empty:
            st.subheader("Configure & Run EIS Analysis")
            run_eis_analysis(df)
        else:
            st.warning("Please upload a valid EIS file.")

if __name__ == "__main__":
    main()


