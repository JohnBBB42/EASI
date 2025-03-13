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
from electrochemical_methods.Analysis_DPV_streamlit import analysis_dpv_streamlit
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
def run_dpv_analysis():
    def parse_complex_dpv_file(dpv_file):
        """
        Reads all lines from dpv_file (Streamlit upload).
        Finds the line that starts with 'V,µA' (the real header),
        parses from that line down as CSV, returns (metadata_lines, dpv_df).
        """
        import io
        content = dpv_file.getvalue()  # raw bytes
        text = content.decode('utf-16', errors='replace')  # your file is often utf-16
        lines = text.split('\n')
        
        # find line that starts with V,µA
        header_line_index = None
        for i, line in enumerate(lines):
            if line.strip().startswith("V,µA"):
                header_line_index = i
                break
        
        if header_line_index is None:
            # none found => return lines as metadata, no DF
            return lines, None
        
        metadata_lines = lines[:header_line_index]
        csv_text       = '\n'.join(lines[header_line_index:])
        import pandas as pd
        from io import StringIO
        df = pd.read_csv(
            StringIO(csv_text),
            sep=',',  # if truly comma-delimited
            engine='python'
        )
        return metadata_lines, df

    ######################################################
    # The main body of run_dpv_analysis
    ######################################################
    st.title("DPV Analysis - Refactored from Analysis_DPV")

    st.markdown("""
    This app refactors your original 'Analysis_DPV' to use data
    from Streamlit uploads instead of reading from a file path.
    """)

    # Step A: user uploads main DPV file
    dpv_file = st.file_uploader("Upload your DPV file (CSV or Excel)", type=["csv","xlsx","xls"])
    dpv_df   = None
    metadata_lines = []

    if dpv_file:
        # user used to pick whether first row is header, but we ignore that now 
        # because we do line-based search
        try:
            # *** Instead of read_csv(...), use parse_complex_dpv_file ***
            metadata_lines, dpv_df = parse_complex_dpv_file(dpv_file)
            if dpv_df is None:
                # no line started with "V,µA"
                st.error("Could not find a line starting with 'V,µA'. Please check your CSV format.")
                st.write("Below are lines we treated as metadata or leftover text:")
                for mline in metadata_lines:
                    st.text(mline)
                dpv_df = None  # just to be safe
            else:
                st.subheader("Metadata lines above the real header:")
                for line in metadata_lines:
                    st.text(line)

                st.write("Preview of DPV data after parsing:")
                st.dataframe(dpv_df.head(15))

        except Exception as e:
            st.error(f"Could not parse DPV file with line-based approach: {e}")
            dpv_df = None

    # Step B: optional blank file
    blank_file=st.file_uploader("Optionally upload a blank file for LOD analysis:", type=["csv","xlsx","xls"])
    blank_df = None
    if blank_file:
        try:
            # This blank file might be simpler, so we can do normal read 
            # or do parse_complex_dpv_file again if it also has weird lines
            blank_df = pd.read_csv(blank_file, engine='python', sep=None, encoding='utf-16')
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
            # Now pass dpv_df to your "analysis_dpv_streamlit" logic
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

        # The rest of your sub-options (Plot Data Array, etc.) unchanged
        sub_options = ["None","Plot Data Array with Corrected Baseline","Analyze Peak Currents","Observed vs Expected Concentration"]
        sub_choice  = st.selectbox("Additional Analysis Step", sub_options)
        if sub_choice == "Plot Data Array with Corrected Baseline":
            data_array    = analysis_results['data_array']
            headers       = analysis_results['headers']
            parsed_md     = analysis_results['parsed_metadata']
            if data_array is None:
                st.error("No data_array found in results.")
            else:
                n_pairs = len(headers)//2
                checks=[]
                for i in range(n_pairs):
                    c=st.checkbox(f"Plot Pair {i+1} (cols {i*2} & {i*2+1})", value=True)
                    checks.append(c)
                selected_indices=[i*2 for i,c in enumerate(checks) if c]
                fig = plot_data_array_with_corrected_baseline(data_array, headers, parsed_md, selected_indices)
                st.pyplot(fig)

        elif sub_choice=="Analyze Peak Currents":
            st.write("Peak Currents analysis was done above.")
        elif sub_choice=="Observed vs Expected Concentration":
            from math import isnan
            mean_peaks=analysis_results['mean_peak_currents']
            # ... do your old plotting ...
            st.write("Plot observed vs expected, etc.")
        st.info("Done with sub-steps.")


def main():
    st.title("Electrochemical Analysis Software Interface (EASI)")
    st.markdown("""
    Welcome to the Electrochemical Analysis Software Interface (EASI).
    
    This tool provides various electrochemical analysis methods:
    1. Cyclic Voltammetry (CV) Analysis
    2. Differential Pulse Voltammetry (DPV) Analysis
    3. Electrochemical Impedance Spectroscopy (EIS) Analysis
    
    Please select an option from the dropdown below to proceed.
    """)

    analysis_choice = st.selectbox(
        "Select Analysis",
        ["(None)", "CV Analysis", "DPV Analysis", "EIS Analysis"]
    )

    if analysis_choice == "(None)":
        st.info("Please choose an analysis to begin.")
        return

    ############################################
    # CV Analysis
    ############################################
    if analysis_choice == "CV Analysis":
        st.header("Cyclic Voltammetry (CV) Analysis")

        # Let user upload a CSV/Excel
        uploaded_file = st.file_uploader(
            "Upload a CV data file", 
            type=["csv","xlsx","xls"]
        )
        df = None
        if uploaded_file:
            # Ask if there's a header
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
            # Now call your run_cv_analysis(df)
            st.subheader("Configure & Run CV Analysis")
            run_cv_analysis(df)
        else:
            st.warning("Please upload a valid CV file.")

    ############################################
    # DPV Analysis
    ############################################
    elif analysis_choice == "DPV Analysis":
        st.header("Differential Pulse Voltammetry (DPV) Analysis")

        # Your DPV logic is already a self-contained function 
        # that does the file uploads inside it. 
        # So we just call it directly. 
        # (It will show its own UI for uploading main + blank file.)
        run_dpv_analysis()

    ############################################
    # EIS Analysis
    ############################################
    elif analysis_choice == "EIS Analysis":
        st.header("Electrochemical Impedance Spectroscopy (EIS) Analysis")

        # Let user upload a CSV/Excel
        uploaded_file = st.file_uploader(
            "Upload an EIS data file", 
            type=["csv","xlsx","xls"]
        )
        df = None
        if uploaded_file:
            # Ask if there's a header
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
            # Now call your run_eis_analysis(df)
            st.subheader("Configure & Run EIS Analysis")
            run_eis_analysis(df)
        else:
            st.warning("Please upload a valid EIS file.")


#############################################
# The usual Streamlit entry point
#############################################
if __name__ == "__main__":
    main()

