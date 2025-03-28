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

# Import your analysis functions:
from electrochemical_methods.Analysis_DPV import Analysis_DPV
from electrochemical_methods.Analysis_CV import Analysis_CV
from electrochemical_methods.Analysis_EIS import Analysis_EIS

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
    uploaded_file = st.file_uploader("Upload a DPV data file", type=["csv","xlsx","xls"])
    blank_file = st.file_uploader("Upload a Blank Responses File (Optional)", type=["csv","xlsx","xls"])
    
    if uploaded_file is not None:
        # 1) Save the uploaded DPV file to a temp file (Analysis_DPV expects a file path)
        tfile = tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(uploaded_file.name)[-1])
        tfile.write(uploaded_file.getvalue())
        tfile.close()
        file_path = tfile.name
        
        # 2) Handle optional blank file
        if blank_file is not None:
            tbfile = tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(blank_file.name)[-1])
            tbfile.write(blank_file.getvalue())
            tbfile.close()
            blank_path = tbfile.name
        else:
            blank_path = None
        
        # 3) When the user clicks the button, run the DPV analysis
        if st.button("Run DPV Analysis"):
            try:
                # Call Analysis_DPV with or without the blank responses
                if blank_path:
                    dpv_result = Analysis_DPV(file_path, blank_responses=blank_path)
                else:
                    dpv_result = Analysis_DPV(file_path)
                
                st.success("DPV analysis completed successfully!")
                
                # 4) Display results in the same style as your old Tkinter app
                st.write("### DPV Analysis Results")

                if not dpv_result or not isinstance(dpv_result, dict):
                    # If there's no structured dictionary, just show raw output
                    st.warning("No structured results found. Showing raw output:")
                    st.write(dpv_result)
                    return
                
                # If the analysis returned a dictionary, let's format the main sections:
                # (mean_peak_currents, std_peak_currents, cov_peak_currents, etc.)
                
                # Mean Peak Currents
                if 'mean_peak_currents' in dpv_result:
                    st.write("**Mean Peak Currents:**")
                    for key, val in dpv_result['mean_peak_currents'].items():
                        # key might be something like ('10uM', 'HX,Xan,UA'), val is a list of floats
                        st.write(f"{key}: {val}")
                
                # Standard Deviation of Peak Currents
                if 'std_peak_currents' in dpv_result:
                    st.write("\n**Standard Deviation of Peak Currents:**")
                    for key, val in dpv_result['std_peak_currents'].items():
                        st.write(f"{key}: {val}")
                
                # Coefficient of Variation
                if 'cov_peak_currents' in dpv_result:
                    st.write("\n**Coefficient of Variation of Peak Currents:**")
                    for key, val in dpv_result['cov_peak_currents'].items():
                        st.write(f"{key}: {val}")
                
                # You can add more sections here if `Analysis_DPV` returns other data
                # e.g. 'lod_results', 't_test_results', etc.

            except Exception as e:
                st.error(f"Error during DPV Analysis: {e}")
    else:
        st.warning("Please upload a valid DPV data file.")


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


