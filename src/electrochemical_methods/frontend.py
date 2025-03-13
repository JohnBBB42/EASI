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
    # Updated defaults to match tkinter app
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
    saving_folder    = st.text_input("Saving Folder Path", value=".")
    
    if st.button("Run EIS Analysis"):
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
            # Convert the chosen columns to numeric before analysis
            real_idx = real_col - 1  # zero-based index for df
            imag_idx = imag_col - 1
            df.iloc[:, real_idx] = pd.to_numeric(df.iloc[:, real_idx], errors='coerce')
            df.iloc[:, imag_idx] = pd.to_numeric(df.iloc[:, imag_idx], errors='coerce')

            # Now call your EIS analysis
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
# DPV Analysis
##############################################
def run_dpv_analysis():
    st.write("## DPV Analysis")
    uploaded_file = st.file_uploader("Upload a DPV data file", type=["csv","xlsx","xls"])
    blank_file = st.file_uploader("Upload a Blank Responses File (Optional)", type=["csv","xlsx","xls"])
    
    if uploaded_file is not None:
        # Save the uploaded DPV file to a temporary file (Analysis_DPV expects a file path)
        tfile = tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(uploaded_file.name)[-1])
        tfile.write(uploaded_file.getvalue())
        tfile.close()
        file_path = tfile.name
        
        if blank_file is not None:
            tbfile = tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(blank_file.name)[-1])
            tbfile.write(blank_file.getvalue())
            tbfile.close()
            blank_path = tbfile.name
        else:
            blank_path = None
        
        if st.button("Run DPV Analysis"):
            try:
                # Call Analysis_DPV with or without the blank responses file
                if blank_path:
                    result = Analysis_DPV(file_path, blank_responses=blank_path)
                else:
                    result = Analysis_DPV(file_path)
                st.success("DPV analysis completed successfully!")
                st.write("Analysis Results:")
                st.write(result)
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


