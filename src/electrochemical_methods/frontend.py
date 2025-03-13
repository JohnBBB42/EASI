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
from scipy.stats import t, linregress
from scipy.signal import find_peaks, savgol_filter
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.ndimage import gaussian_filter1d
from PIL import Image
import matplotlib.pyplot as plt

# from electrochemical_methods.Analysis_DPV import Analysis_DPV
from electrochemical_methods.Analysis_DPV import analysis_dpv_streamlit
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
# 4) The main Streamlit UI 
#    So we can replicate the old DPV usage
##############################################
def main():
    st.title("DPV Analysis - Refactored from Analysis_DPV")

    st.markdown("""
    This app refactors your original 'Analysis_DPV' to use data
    from Streamlit uploads instead of reading from a file path.
    """)

    # Step A: user uploads main DPV file
    dpv_file = st.file_uploader("Upload your DPV file (CSV or Excel)", type=["csv","xlsx","xls"])
    dpv_df   = None
    if dpv_file:
        # user picks whether first row is a real header or not
        has_header = st.checkbox("First row is header line in DPV file?", value=False)
        try:
            if dpv_file.name.lower().endswith(".csv"):
                dpv_df = pd.read_csv(dpv_file, header=0 if has_header else None, encoding='utf-16', sep=None, engine='python')
            else:
                dpv_df = pd.read_excel(dpv_file, header=0 if has_header else None)
            st.write("Preview of DPV data:")
            st.dataframe(dpv_df.head(15))
        except Exception as e:
            st.error(f"Could not parse DPV file: {e}")
            dpv_df=None

    # Step B: optional blank file
    blank_file=st.file_uploader("Optionally upload a blank file for LOD analysis:", type=["csv","xlsx","xls"])
    blank_df = None
    if blank_file:
        blank_header = st.checkbox("Blank file has header?", value=False)
        try:
            if blank_file.name.lower().endswith(".csv"):
                blank_df = pd.read_csv(blank_file, header=0 if blank_header else None, encoding='utf-16', sep=None, engine='python')
            else:
                blank_df = pd.read_excel(blank_file, header=0 if blank_header else None)
            st.write("Preview of Blank data:")
            st.dataframe(blank_df.head(15))
        except Exception as e:
            st.error(f"Could not parse blank file: {e}")
            blank_df=None

    # Step C: run analysis
    if st.button("Run DPV Analysis"):
        if dpv_df is None or dpv_df.empty:
            st.error("Please upload a valid DPV file before running analysis.")
            return
        with st.spinner("Running analysis..."):
            analysis_results = analysis_dpv_streamlit(dpv_df=dpv_df, blank_df=blank_df)

        if not analysis_results:
            st.warning("No results returned from analysis.")
            return

        st.success("DPV analysis complete!")
        # We can show them
        st.write("**Mean Peak Currents**:")
        st.json(analysis_results['mean_peak_currents'])

        st.write("**STD of Peak Currents**:")
        st.json(analysis_results['std_peak_currents'])

        st.write("**COV of Peak Currents**:")
        st.json(analysis_results['cov_peak_currents'])

        if analysis_results['lod_results'] is not None:
            st.write("**LOD Results**:")
            st.json(analysis_results['lod_results'])

        st.write("**T-Test Results**:")
        st.json(analysis_results['t_test_results'])

        st.write("**Parsed Metadata**:")
        st.json(analysis_results['parsed_metadata'])

        st.write("**Headers**:")
        st.json(analysis_results['headers'])

        # Optionally show the peak detection "results"
        st.write("**Collected Peak Results**:")
        st.dataframe(pd.DataFrame(analysis_results['results']))

        # If you want to replicate the old dropdown logic 
        # (Plot Data Array with Corrected Baseline, etc.),
        # you can do that. Example:
        sub_options = ["None","Plot Data Array with Corrected Baseline","Analyze Peak Currents","Observed vs Expected Concentration"]
        sub_choice  = st.selectbox("Additional Analysis Step", sub_options)
        if sub_choice == "Plot Data Array with Corrected Baseline":
            data_array    = analysis_results['data_array']
            headers       = analysis_results['headers']
            parsed_md     = analysis_results['parsed_metadata']
            if data_array is None:
                st.error("No data_array found in results.")
            else:
                # Let user pick how many pairs to plot
                # For simplicity, we do all
                # We'll build a selected_indices by default
                n_pairs = len(headers)//2
                checks=[]
                for i in range(n_pairs):
                    c=st.checkbox(f"Plot Pair {i+1} (cols {i*2} & {i*2+1})", value=True)
                    checks.append(c)
                selected_indices=[i*2 for i,c in enumerate(checks) if c]
                fig = plot_data_array_with_corrected_baseline(data_array, headers, parsed_md, selected_indices)
                st.pyplot(fig)

        elif sub_choice=="Analyze Peak Currents":
            # already done in main analysis, but we can do more 
            st.write("Peak Currents analysis was done. Shown above.")

        elif sub_choice=="Observed vs Expected Concentration":
            from math import isnan
            mean_peaks=analysis_results['mean_peak_currents']
            figs = []
            def plot_obs_vs_exp(mean_peaks):
                # replicate your old 'plot_observed_vs_expected_concentration'
                from math import isnan
                figs = []
                analytes_dict=defaultdict(list)
                for key, arr in mean_peaks.items():
                    conc_str, analyte=key
                    try:
                        c_val=float(conc_str[:-2])
                    except:
                        c_val=0
                    analytes_dict[analyte].append((c_val, arr))
                for analyte, data in analytes_dict.items():
                    data.sort(key=lambda x: x[0])
                    max_peaks=max(len(p[1]) for p in data)
                    fig, ax=plt.subplots()
                    for peak_idx in range(max_peaks):
                        x_vals=[]
                        y_vals=[]
                        for (c_val, arr) in data:
                            if peak_idx<len(arr):
                                x_vals.append(c_val)
                                y_vals.append(arr[peak_idx])
                        if len(x_vals)<2: 
                            continue
                        slope, intercept, r_val, p_val, std_err = linregress(x_vals, y_vals)
                        c_obs=[(y - intercept)/slope for y in y_vals]
                        label_str=f"Peak {peak_idx+1} (R²={r_val**2:.3f})"
                        ax.scatter(x_vals, c_obs, label=label_str)
                    if data:
                        all_concs=[x[0] for x in data]
                        min_val, max_val=min(all_concs), max(all_concs)
                        ax.plot([min_val,max_val],[min_val,max_val],'k--',alpha=0.4)
                    ax.set_xlabel("Expected (uM)")
                    ax.set_ylabel("Observed (uM)")
                    ax.set_title(f"Observed vs Expected: {analyte}")
                    ax.legend()
                    figs.append(fig)
                return figs

            figs=plot_obs_vs_exp(mean_peaks)
            for f in figs:
                st.pyplot(f)

        st.info("Done with sub-steps.")
