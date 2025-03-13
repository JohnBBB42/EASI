import streamlit as st
import pandas as pd
import numpy as np
import os

# If these are local modules (like your old code):
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
