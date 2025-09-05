"""
Created on Fri Aug 21 14:11:08 2024
@author: jonatanemilsvendsen
"""
#_____________________________
# Presetting
#_____________________________

# Tweaks and ideas
# Have to import files after choosing interval
# Choose range after import
# Choose range before baseline correction
# Results for specific graphs - DPV Analysis, LOD and T-test

# Methods - Check
# Baseline correction methods: Gaussian, polynomial, linear, Savitzky-Golay
# Peak fitting Gaussian fit or rarely Lorentzian fit
# Calibration curve - Linear curve!
# LOQ
# Levenberg Marquardt and gaussian fitting for iteration

# Important
# Peak deconvolution - working on it - DONE
# Residual vs potential(peak location) plot (points should be close to 0) - Bland-Altman plot - DONE
# Baseline correction should be for each analyte - DONE
# Take multiple pstrace files

# Easy to do
# Report Peak height and center of gravity (peak location) - DONE
# Say how many peaks are in the data
# peak maxima, position, FWHM, Relative standard deviation (RSD)
# Peak placement vs concentration plot 
# In analyze peak current show ax+b formula

# Reporting
# Make a data sheet with all the results - Done
# Look at which Statistics to include
# Report of the fit analysis
# Look at Neetis notebook - Done
# Excel file output with all the results

# User friendliness


import os, sys, re
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
sys.path.append(dname)

# Import packages
import csv
import json
import pandas as pd
import tempfile
import numpy as np
import tkinter as tk
from scipy.signal import find_peaks
from collections import defaultdict
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.stats import linregress
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from datetime import datetime
from tkinter import ttk
import pandas as pd
from tkinter import Button, END, Entry, Label, LabelFrame, OptionMenu, StringVar, messagebox, Toplevel, filedialog
import hashlib
from PIL import Image
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from electrochemical_methods.Analysis_CV import Analysis_CV
from electrochemical_methods.Analysis_DPV import Analysis_DPV
from electrochemical_methods.Analysis_EIS import Analysis_EIS
from electrochemical_methods.Plotting import Plots_all_imported
# from electrochemical_methods.EASI_Plotting import (
#     plot_data_array,
#     plot_data_array_with_corrected_baseline,
#     plot_concentrations_vs_mean_peak_currents,
#     analyze_peak_currents,
#     plot_observed_vs_expected_concentration,
# )
# from electrochemical_methods.EASI_Utils import (
#     open_plot_selection_window,
#     toggle_all,
# )

class CVApp:
    def __init__(self, master):
        self.master = master
        master.title("CV Analysis")
        master.geometry("500x300")

        self.df = None                # We'll store the loaded DataFrame here
        self.filepath = None

        # Dropdown for future expansions
        self.options = ["CV Analysis"]
        self.clicked = StringVar(value=self.options[0])

        frame = LabelFrame(master, text="Select Analysis")
        frame.pack(padx=10, pady=10, fill="x")
        OptionMenu(frame, self.clicked, *self.options).pack(fill="x")

        Button(master, text="Import File", command=self.import_file).pack(pady=5)
        self.file_label = Label(master, text="No file selected")
        self.file_label.pack(pady=5)

        Button(master, text="Run Analysis", command=self.run_analysis).pack(pady=5)

    def import_file(self):
        """
        Prompt user for a single CSV or Excel, load into self.df exactly once.
        No headers are assumed. If your file has a header row, adjust as needed.
        """
        self.filepath = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx;*.xls")]
        )
        if not self.filepath:
            return

        self.file_label.config(text=self.filepath)

        # Attempt to load
        try:
            if self.filepath.lower().endswith(".csv"):
                # If the file actually has a header row, remove 'header=None'
                self.df = pd.read_csv(self.filepath, header=None)
            else:
                self.df = pd.read_excel(self.filepath, header=None)
        except Exception as e:
            messagebox.showerror("File Error", f"Could not read file: {e}")
            self.df = None

    def run_analysis(self):
        """
        Called when user clicks 'Run Analysis':
          1) Check if df is loaded
          2) Auto-detect columns
          3) Ask user for a folder to save all results
          4) Call analysis
        """
        if self.df is None or self.df.empty:
            messagebox.showwarning("No Data", "Please import a valid file first.")
            return

        # Step 1: Auto-detect columns
        # We'll do a simple approach:
        #  - If df has at least 2 columns: col0=Potential, col1=Current
        #  - If df has ≥3 columns, we assume col2=Scan
        #  - Otherwise, error
        num_cols = self.df.shape[1]
        if num_cols < 2:
            messagebox.showerror("Column Error", "Data must have at least 2 columns.")
            return

        pot_col = 1  # 1-based index for Analysis_CV
        cur_col = 2
        scan_col= 0  # default = 0 means 'no scan column used'

        if num_cols >= 3:
            scan_col = 3  # 1-based index if we see there's a 3rd column

        # Step 2: Ask user for a saving folder
        saving_folder = filedialog.askdirectory(title="Select folder to save results")
        if not saving_folder:
            messagebox.showwarning("No Folder", "No folder selected, analysis canceled.")
            return

        # Step 3: Call the Analysis_CV function
        try:
            # We'll pass the loaded df, plus we tell it row_start=1 if there's no header
            # or row_start=2 if we want to skip the first row, etc.
            # Adjust as needed. For now, we'll do row_start=1 if your data starts on row 1.
            Analysis_CV(
                df=self.df,             # DataFrame you already read
                values_row_start=2,     # skip the header row
                potential_column=1,     # A
                current_column=3,       # C
                scan_column=5,          # E
                scan_number=1,          # must be an actual scan number in col E
                linreg_start_index=15,
                r2_accept_value=0.90,
                potential_unit="V",
                current_unit="A",
                num_decimals=3,
                saving_folder="."
            )

            messagebox.showinfo("Success", "CV analysis completed successfully!")
        except Exception as e:
            messagebox.showerror("Analysis Error", f"An error occurred: {e}")

class EISApp:
    def __init__(self, master):
        self.master = master
        master.title("EIS Analysis")
        master.geometry("500x300")

        self.df = None
        self.filepath = None

        self.options = ["EIS Analysis"]
        self.clicked = StringVar(value=self.options[0])

        frame = LabelFrame(master, text="Select Analysis")
        frame.pack(padx=10, pady=10, fill="x")
        OptionMenu(frame, self.clicked, *self.options).pack(fill="x")

        Button(master, text="Import File", command=self.import_file).pack(pady=5)
        self.file_label = Label(master, text="No file selected")
        self.file_label.pack(pady=5)

        Button(master, text="Run Analysis", command=self.run_analysis).pack(pady=5)

    def import_file(self):
        """
        User picks a single CSV or Excel. We assume the first row has headers,
        so we do header=0. If your data uses commas as decimals, add decimal=','.
        """
        self.filepath = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx;*.xls")]
        )
        if not self.filepath:
            return

        self.file_label.config(text=self.filepath)
        try:
            if self.filepath.lower().endswith(".csv"):
                # If your decimals are normal dots, no change needed.
                # If your decimals are commas, you might do: decimal=','
                self.df = pd.read_csv(self.filepath, header=0)
            else:
                self.df = pd.read_excel(self.filepath, header=0)
        except Exception as e:
            messagebox.showerror("File Error", f"Could not read file: {e}")
            self.df = None

    def run_analysis(self):
        if self.df is None or self.df.empty:
            messagebox.showwarning("No Data", "Please import a valid file first.")
            return

        selected_option = self.clicked.get()
        if selected_option == "EIS Analysis":
            try:
                saving_folder = filedialog.askdirectory(title="Select folder to save EIS results")
                if not saving_folder:
                    messagebox.showwarning("No folder", "Analysis canceled.")
                    return

                # Suppose from your data screenshot:
                #   Column A: Index
                #   Column B: Frequency
                #   Column C: Z' (Ω)
                #   Column D: -Z'' (Ω)
                #   ...
                # That means: real_col=3, imag_col=4 (1-based indices).
                # If your actual data starts on row 2 (since row 1 was headers),
                # often you do values_row_start=1 or 2 depending on how you want to skip rows.
                # Try values_row_start=1 if your data is now row 0-based after header=0.

                result = Analysis_EIS(
                    df=self.df,
                    values_row_start=1,   # Don't skip any additional lines now that header=0
                    real_col=3,           # If Z' is the 3rd column in 1-based indexing
                    imag_col=4,           # If -Z'' is the 4th column
                    x_start=None,
                    x_end=None,
                    y_start=None,
                    y_end=None,
                    unit="Ω",
                    circle_pt1_index=0,
                    circle_pt2_index=0,
                    saving_folder=saving_folder
                )
                messagebox.showinfo("Success", f"EIS analysis done. Plot at {result['plot_path']}")
            except Exception as e:
                messagebox.showerror("Error", str(e))

class Custom_Plot_App:
    def __init__(self,master):
        self.master = master
        master.title('Customized Plotting')
        master.geometry("500x400")
        
        # Drop down boxes
        self.options = ["Plot imported files"]
        self.clicked = StringVar()
        self.clicked.set(self.options[0])
        
        # Make frame
        self.frame = LabelFrame(master,text="Select analysis", height=5, width=5)
        self.frame.grid(row=0,column=1)
        #_____________________________
        # Dropdown menu
        #_____________________________
        dropAnalys = OptionMenu(self.frame, self.clicked, *self.options)
        dropAnalys.grid(row=0, column=1)
        
        self.clicked.trace("w",self.change) #track if something happens to variable
        
        #_____________________________
        # Initial Entry boxes
        #_____________________________
        self.variable1 = self.entry_box('values_row_start',2,1,2)
        self.variable2 = self.entry_box('x_column',3,1,2)
        self.variable3 = self.entry_box('y_column',4,1,5)
        self.variable4 = self.entry_box('x_start',5,1,"")
        self.variable5 = self.entry_box('x_end',6,1,"")
        self.variable6 = self.entry_box('y_start',7,1,"")
        self.variable7 = self.entry_box('y_end',8,1,"")
        self.variable8 = self.entry_box('x_label',9,1,"")
        self.variable9 = self.entry_box('y_label',10,1,"")
        self.variable10 = self.entry_box('plot_title',11,1,"")
        
        #_____________________________
        # Info boxes
        #_____________________________
        self.button1 = Button(master, text = "info",command=self.info_box1).grid(row=2, column=2)
        self.button2 = Button(master, text = "info", command=self.info_box2).grid(row=3, column=2)
        self.button3 = Button(master, text = "info", command=self.info_box3).grid(row=4, column=2)
        self.button4 = Button(master, text = "info", command=self.info_box4).grid(row=5, column=2)
        self.button5 = Button(master, text = "info", command=self.info_box5).grid(row=6, column=2)
        self.button6 = Button(master, text = "info", command=self.info_box6).grid(row=7, column=2)
        self.button7 = Button(master, text = "info", command=self.info_box7).grid(row=8, column=2)
        self.button8 = Button(master, text = "info", command=self.info_box8).grid(row=9, column=2)
        self.button9 = Button(master, text = "info", command=self.info_box9).grid(row=10, column=2)
        self.button10 = Button(master, text = "info", command=self.info_box10).grid(row=11, column=2)

                
        #_____________________________
        # Run button
        #_____________________________
        self.run_button = Button(master,text="Run!",command = self.onclick)
        self.run_button.grid(row=12,column=1)
        
    #_____________________________
    # Info boxes - Text definitions
    #_____________________________ 
    def info_box1(self):
        Text = "The row for which the First value appears in the excel/csv file."
        tk.messagebox.showinfo("Info",str(Text))
    
    def info_box2(self):
        if self.clicked.get() == self.options[0]:
            Text = "CCheck your excel/csv file and list the column number of where the x-values appear."
        tk.messagebox.showinfo("Info",str(Text))
    
    def info_box3(self):
        if self.clicked.get() == self.options[0]:
            Text = "Check your excel/csv file and list the column number of where the y-values appear."
        tk.messagebox.showinfo("Info",str(Text))
    
    def info_box4(self):
        if self.clicked.get() == self.options[0]:
            Text = "Insert value for Start interval of x-values. NB: both x_start and x_end must be inserted for the interval range to be created."
        tk.messagebox.showinfo("Info",str(Text))
    
    def info_box5(self):
        if self.clicked.get() == self.options[0]:
            Text = "Insert value for End interval of x-values. NB: both x_start and x_end must be inserted for the interval range to be created."
        tk.messagebox.showinfo("Info",str(Text))
    
    def info_box6(self):
        if self.clicked.get() == self.options[0]:
            Text = "Insert value for Start interval of y-values. NB: both y_start and y_end must be inserted for the interval range to be created."
        tk.messagebox.showinfo("Info",str(Text))
    
    def info_box7(self):
        if self.clicked.get() == self.options[0]:
            Text = "Insert value for End interval of y-values. NB: both y_start and y_end must be inserted for the interval range to be created."
        tk.messagebox.showinfo("Info",str(Text))        

    def info_box8(self):
        if self.clicked.get() == self.options[0]:
            Text = "Type the x-label."
        tk.messagebox.showinfo("Info",str(Text))  
    
    def info_box9(self):
        if self.clicked.get() == self.options[0]:
            Text = "Type the y-label."
        tk.messagebox.showinfo("Info",str(Text))  
    
    def info_box10(self):
        if self.clicked.get() == self.options[0]:
            Text = "Type the plot-title."
        tk.messagebox.showinfo("Info",str(Text))  
    #_____________________________
    # Entry boxes
    #_____________________________
    def entry_box(self,Text, row, column, default_input, width = 20):
        Text = str(Text)
        Label(self.master, text=str(Text),width = width).grid(row=row)
        E = Entry(self.master)
        E.insert(END,str(default_input))
        E.grid(row=row, column=column) 
        return E

    def change(self, *args):
        if self.clicked.get() == self.options[0]:
            self.variable1 = self.entry_box('values_row_start',2,1,2)
            self.variable2 = self.entry_box('x_column',3,1,2)
            self.variable3 = self.entry_box('y_column',4,1,5)
            self.variable4 = self.entry_box('x_start',5,1,"")
            self.variable5 = self.entry_box('x_end',6,1,"")
            self.variable6 = self.entry_box('y_start',7,1,"")
            self.variable7 = self.entry_box('y_end',8,1,"")
            self.variable8 = self.entry_box('x_label',9,1,"")
            self.variable9 = self.entry_box('y_label',10,1,"")
            self.variable10 = self.entry_box('plot_title',11,1,"")
        else:
            pass

    #_____________________________
    # Save input in entry boxes and function for run-command
    #_____________________________
    
    def onclick(self,*args): 
        self.variable1_get = self.variable1.get()
        self.variable2_get = self.variable2.get()
        self.variable3_get = self.variable3.get()
        self.variable4_get = self.variable4.get()
        self.variable5_get = self.variable5.get()
        self.variable6_get = self.variable6.get()
        self.variable7_get = self.variable7.get()
        self.variable8_get = self.variable8.get()
        self.variable9_get = self.variable9.get()
        self.variable10_get = self.variable10.get()
        values = [int(self.variable1_get),int(self.variable2_get),int(self.variable3_get),str(self.variable4_get),str(self.variable5_get),
                  str(self.variable6_get),str(self.variable7_get),str(self.variable8_get),str(self.variable9_get),str(self.variable10_get)]
        
        #option to close tkinter window
        #self.close_window

        if self.clicked.get() == self.options[0]:
            Plots_all_imported(values[0], values[1],values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9])
        else:
            pass
        
    def close_window(self): 
        self.master.destroy()

class DPVApp:
    def __init__(self, master, main_window):
        self.master = master
        self.main_window = main_window
        master.title('DPV analysis')
        # Positioning the DPV window next to the main window
        master.geometry("420x640+{}+{}".format(main_window.winfo_x() + main_window.winfo_width(), main_window.winfo_y()))

        # Create a main frame to center everything
        self.main_frame = tk.Frame(master)
        self.main_frame.pack(expand=True, fill='both')

        # Variable to store the directory path for saving results
        self.save_directory = None
        self.analysis_results = None  # Store analysis results
        self.filepath = None
        self.blank_filepath = None

        # Drop down boxes
        self.options = [
            "DPV analysis", 
            #"Plot Data Array", 
            "Plot Data Array with Corrected Baseline", 
            #"Plot Concentrations vs Mean Peak Currents", 
            "Analyze Peak Currents",
            "Observed vs Expected Concentration",
            "LOD Analysis",
            "T-test Analysis",
            "Generate Peak Report",    # ← NEW
            "Peak Deconvolution",
            "Residual vs Potential Plot",
        ]
        self.clicked = StringVar()
        self.clicked.set(self.options[0])

        # Make frame
        self.frame = LabelFrame(self.main_frame, text="Select analysis", height=5, width=5)
        self.frame.pack(padx=10, pady=10, fill="x")

        # Dropdown menu
        dropAnalys = OptionMenu(self.frame, self.clicked, *self.options)
        dropAnalys.pack(padx=10, pady=10, fill="x")

        self.clicked.trace("w", self.change)  # track if something happens to variable
        
        # ------------------------
        # Add Metadata Fields Entry:
        # ------------------------
        self.metadata_fields_label = Label(self.main_frame, text="Metadata Fields (comma separated):")
        self.metadata_fields_label.pack(padx=10, pady=5)
        self.metadata_fields_entry = Entry(self.main_frame, width=50)
        # Set a default string that represents the default field names
        self.metadata_fields_entry.insert(END, "Electrode,Analytes,Concentration,Method")
        self.metadata_fields_entry.pack(padx=10, pady=5)

        # -------------------------------------------------
        # Peak-region selector (one line, free-text)
        # -------------------------------------------------
        self.peak_label = tk.Label(self.main_frame,
                                text="Peak regions e.g.  HX:0.78-0.88; UA:0.32-0.40,0.82-0.92")
        self.peak_label.pack(padx=10, pady=(15, 5))

        self.peak_entry = tk.Entry(self.main_frame, width=45)
        self.peak_entry.pack(padx=10, pady=5)

        # keep a parsed copy here
        self.peak_regions = None

        #-------------------------
        # Import file(s)
        #-------------------------

        # Import file button
        self.import_button = Button(self.main_frame, text="Import File", command=self.import_file)
        self.import_button.pack(pady=10)

        # Display file path
        self.file_path_label = Label(self.main_frame, text="No file selected", wraplength=300)
        self.file_path_label.pack(pady=10)

        # Import blank responses file button
        self.import_blank_button = Button(self.main_frame, text="Import Blank Responses File (Optional)", command=self.import_blank_file)
        self.import_blank_button.pack(pady=10)

        # Display blank responses file path
        self.blank_file_path_label = Label(self.main_frame, text="No blank responses file selected", wraplength=300)
        self.blank_file_path_label.pack(pady=10)

        # -------------------------
        # Concentration Conversion
        # -------------------------

        self.unit_frame = LabelFrame(self.main_frame, text="Concentration Settings")
        self.unit_frame.pack(padx=10, pady=5, fill="x")

        # For demonstration, let's define possible units:
        self.unit_choices = [
            "",        # blank → no conversion
            "nM",
            "µM",
            "mM",
            "M",
            "g/L",
            "mg/L",
            # ... add others as needed
        ]

        Label(self.unit_frame, text="From Unit:").pack(anchor='w')
        self.from_unit_var = StringVar()
        self.from_unit_var.set("")  # default is blank
        self.from_unit_menu = OptionMenu(self.unit_frame, self.from_unit_var, *self.unit_choices)
        self.from_unit_menu.pack(padx=10, fill='x')

        Label(self.unit_frame, text="To Unit:").pack(anchor='w')
        self.to_unit_var = StringVar()
        self.to_unit_var.set("")  # default is blank (no conversion)
        self.to_unit_menu = OptionMenu(self.unit_frame, self.to_unit_var, *self.unit_choices)
        self.to_unit_menu.pack(padx=10, fill='x')

        # Button to do the conversion
        self.convert_button = Button(
            self.unit_frame, 
            text="Convert Concentrations", 
            command=self.convert_concentrations
        )
        self.convert_button.pack(pady=5)

        # -------------------------
        # Run and Save
        # -------------------------

        # Run button
        self.run_button = Button(self.main_frame, text="Run Analysis", command=self.run_analysis)
        self.run_button.pack(pady=10)

        # Save results button
        self.save_as_button = Button(self.main_frame, text="Save Results", command=self.save_results, state=tk.DISABLED)
        self.save_as_button.pack(pady=10)


    # ---------------------------------------------------------
    # helper: show what goes into Analysis_DPV
    # ---------------------------------------------------------
    def _debug_dump(self, csv_path: str, n_lines: int = 12) -> None:
        print("\n=== DEBUG — first lines of temp CSV ===")
        with open(csv_path, "r", encoding="utf-16") as f:
            for i in range(n_lines):
                try:
                    print(f"{i:2}: {next(f).rstrip()}")
                except StopIteration:
                    break
        print("=== end preview =====================================\n")

    def _parse_peak_region_string(self, txt: str) -> dict | None:
        """
        Convert the free-text field into the dict expected by Analysis_DPV.
        Returns None if the field is blank.
        """
        import re
        txt = txt.strip()
        if not txt:
            return None                   # user left it blank → analyze everything

        regions = {}
        for chunk in txt.split(';'):
            if ':' not in chunk:
                continue
            name, ranges = chunk.split(':', 1)
            windows = []
            for r in re.split(r'[,\s]+', ranges.strip()):
                if '-' in r:
                    try:
                        lo, hi = map(float, r.split('-', 1))
                        windows.append((lo, hi))
                    except ValueError:
                        pass
            if windows:
                regions[name.strip()] = windows
        return regions or None

    # ---------------------------------------------------------
    # robust .pssession loader  (measurement‑level metadata)
    # ---------------------------------------------------------
    def load_pssession(self, path: str) -> tuple[pd.DataFrame, list[str]]:
        """
        df      – numeric columns  V, µA, V_2, µA_2 …  (ready for Analysis_DPV)
        tokens  – one metadata token per curve, e.g.
                ['pLM_HX_D29_60µM_DPV', 'pLM_HX_D29_60µM_DPV', …]
        """
        import json, re, pandas as pd

        # ── pull the JSON header out of the .pssession file ──────────────────────
        with open(path, "r", encoding="utf-16") as fh:
            txt = fh.read()

        l = txt.find("{")
        depth = 0
        for i, ch in enumerate(txt[l:], l):
            depth += (ch == "{") - (ch == "}")
            if depth == 0:
                r = i
                break
        sess = json.loads(txt[l : r + 1])

        # ── iterate over all DPV measurements ────────────────────────────────────
        cols, tokens = {}, []
        col_idx = 1                                           # V_2 / µA_2 suffixes

        for meas in sess["Measurements"]:
            if "METHOD_ID=DPV" not in meas.get("Method", "").upper():
                continue

            title = meas.get("Title") or meas.get("Name") or ""
            print(f"[DBG] meas title: {title!r}")
            body  = title.split(". ", 1)[-1]                  # drop leading “1. ”
            elec  = body.split("_", 1)[0] if "_" in body else "UNK"

            # ── 1) first try the special  D<day>.<conc><unit>*  pattern ─────────
            #     e.g.  “…D29.60uM* HX‑repeat…”
            m_special = re.search(
                r"D(?P<day>\d+)[.,](?P<val>\d+(?:[.,]\d+)?)(?P<unit>[munpµ]M)\*?",
                body, flags=re.I,
            )

            if m_special:
                day      = m_special.group("day")
                val_txt  = m_special.group("val").replace(".", ",")      # keep comma
                raw_unit = m_special.group("unit")
                unit     = "µM" if raw_unit.lower().startswith(("u", "µ")) else "mM"
                conc_str = f"{val_txt}{unit}"

                # analyte = first word after the concentration
                tail  = body[m_special.end():]
                ama   = re.search(r"([A-Za-z]{1,15})", tail)
                ana   = (ama.group(1) if ama else "UNK").upper()

                analyte_label = f"{ana}_D{day}"               # ← HX_D29

            else:
                # ── 2) generic patterns  HX_50µM  or  50µM_HX ─────────────────
                m_generic = (
                    re.search(
                        r"(?P<ana>[A-Za-z]+)[\s*_\-]*(?P<val>\d+(?:[.,]\d+)?)(?P<unit>[munpµ]M)\*?",
                        body, flags=re.I,
                    )
                    or re.search(
                        r"(?P<val>\d+(?:[.,]\d+)?)(?P<unit>[munpµ]M)\*?[\s*_\-]*(?P<ana>[A-Za-z]+)",
                        body, flags=re.I,
                    )
                )
                if m_generic:
                    val_txt = m_generic.group("val").replace(".", ",")
                    raw_unit= m_generic.group("unit")
                    unit    = "µM" if raw_unit.lower().startswith(("u", "µ")) else "mM"
                    conc_str= f"{val_txt}{unit}"
                    analyte_label = m_generic.group("ana").upper()
                else:
                    conc_str      = "0µM"
                    analyte_label = "UNK"

            token = f"{elec}_{analyte_label}_{conc_str}_DPV"
            print(f"[DBG] token produced: {token}")

            # ── add every curve belonging to this measurement ───────────────────
            for curve in meas.get("Curves", []):
                raw = curve.get("RawDataArray") or {}
                X   = raw.get("Potential") or raw.get("X")
                Y   = raw.get("Current")   or raw.get("Y")

                if X is None or Y is None:                   # PalmSens 3 fallback
                    xb = curve.get("XAxisDataArray") or {}
                    yb = curve.get("YAxisDataArray") or {}
                    X  = xb.get("DataValues");  Y = yb.get("DataValues")
                    if X and isinstance(X[0], dict):         # dict‑style arrays
                        X = [float(d["V"]) for d in X]
                        Y = [float(d["V"]) for d in Y]

                if X is None or Y is None or len(X) != len(Y):
                    raise ValueError("curve arrays missing or length mismatch")

                suf = "" if col_idx == 1 else f"_{col_idx}"
                cols[f"V{suf}"]  = X
                cols[f"µA{suf}"] = Y
                tokens.append(token)
                col_idx += 1

        if not cols:
            raise ValueError("No DPV curves found in session")

        df = pd.DataFrame(cols)
        print(f"→ Loaded {len(tokens)} curve(s), {len(df)} points each")
        return df, tokens

    
    def import_file(self):
        """
        Import a data file (CSV, Excel, or PalmSens session) and integrate metadata extraction for DPV measurements.
        """
        # Prompt user for file selection (CSV, Excel, or .pssession)
        self.filepath = filedialog.askopenfilename(
            filetypes=[
                ("All files", "*.*"),
                ("PalmSens session", "*.pssession"),
                ("CSV files", "*.csv"),
                ("Excel files", "*.xlsx;*.xls"),
            ]
        )
        if not self.filepath:
            # No file selected, exit early
            self.file_path_label.config(text="No file selected")
            self.analysis_results = None
            return

        # Update the GUI label to show the selected file path
        self.file_path_label.config(text=self.filepath)

        # 1) figure out which metadata fields to extract
        custom_fields = [
            f.strip() for f in 
            self.metadata_fields_entry.get().split(',')
            if f.strip()
        ]
        if not custom_fields:
            custom_fields = ["Electrode","Analytes","Concentration","Method"]

        ext = os.path.splitext(self.filepath)[1].lower()

        # ─────────────────────────── *.pssession* ────────────────────────────
        if ext == ".pssession":
            try:
                df, tokens = self.load_pssession(self.filepath)   # uses code above
                meta_entry = ",,".join(tokens)

                with tempfile.NamedTemporaryFile(
                    suffix=".csv", delete=False,
                    mode="w", encoding="utf-16", newline=""
                ) as tmp:
                    tmp.write("\n" * 3)
                    tmp.write(f"Metadata row: {meta_entry}\n")
                    tmp.write("Date and time measurement:\n")
                    df.to_csv(tmp, index=False, header=True)
                    file_to_pass = tmp.name

            except Exception as e:
                messagebox.showerror("Load Error", f"Could not parse session file:\n{e}")
                self.analysis_results = None
                return


        elif ext in ('.xlsx', '.xls'):
            # 0) read the sheet WITHOUT skipping rows (so we keep titles + V/µA row)
            df_raw = pd.read_excel(self.filepath,
                                engine='openpyxl',
                                header=None)

            # 1) split into:  (row-0) titles  (row-1) V/µA header  (row-2…) numeric data
            title_row   = df_raw.iloc[0]
            header_row  = df_raw.iloc[1].fillna(method='ffill')       # 'V' / 'µA' repeats
            df          = df_raw.iloc[2:].reset_index(drop=True)
            df.columns  = [str(c) for c in header_row]                # keep 'V' 'µA' …

            # 2) make metadata token for every V/µA column-pair
            def parse_token(title: str) -> str:
                parts = str(title).split('.')
                elec  = parts[2] if len(parts) > 2 else ""            # → 'pLM'

                ana_match = re.search(
                    r'_(?:\d+x)?([\d.,]+uM|[\d.,]+mM)\s+([A-Za-z]+)',
                    str(title)
                )
                conc   = ana_match.group(1).replace(',', '.') if ana_match else ""
                analy  = ana_match.group(2) if ana_match else ""

                return f"{elec}_{analy}_{conc}_DPV"

            tokens       = [parse_token(c) for c in title_row[::2]]   # every 2nd cell
            meta_entry   = ",,".join(tokens)                          #  ,, → blank cell
            n_reps       = len(tokens)

            # 3) custom field list from GUI
            custom_fields = [
                f.strip() for f in self.metadata_fields_entry.get().split(',')
                if f.strip()
            ] or ["Electrode", "Analytes", "Concentration", "Method"]

            # 4) write a temporary UTF-16 CSV identical to the .pssession layout
            with tempfile.NamedTemporaryFile(
                    suffix='.csv',
                    delete=False,
                    mode='w',
                    encoding='utf-16',
                    newline=''
            ) as tmp:
                tmp.write("\n"*3)
                tmp.write(f"Metadata row: {meta_entry}\n")
                tmp.write("Date and time measurement:\n")
                df.to_csv(tmp, sep=',', decimal='.', index=False, header=True)
                tmp_excel_csv = tmp.name

                # ------------ EXCEL DEBUG ---------------------------------------
                print("\n=== EXCEL IMPORT DEBUG ===")
                print("Original Excel file :", self.filepath)
                print("title_row (first 6) :", list(title_row[:6]))
                print("header_row (first 6):", list(header_row[:6]))
                print("DataFrame shape     :", df.shape)
                print("custom_fields       :", custom_fields)
                print("n_reps detected     :", n_reps)
                print("meta_entry written  :", meta_entry)
                print("temporary CSV path  :", tmp_excel_csv)
                print("=== END EXCEL DEBUG ===\n")
                # -----------------------------------------------------------------

            file_to_pass = tmp_excel_csv
        # ────────────────────────────────────────────────────────────────────────────────

        else:
            # every other file type goes through unchanged
            file_to_pass = self.filepath

        # ------------------------------------------------------------------
        # build peak_regions from the GUI text-field
        # ------------------------------------------------------------------
        self.peak_regions = self._parse_peak_region_string(self.peak_entry.get())

        # ------------- preview right before handing to Analysis_DPV ------------
        print("\n>>> Hand-off to Analysis_DPV")
        print("peak_regions :", self.peak_regions)
        print("Preview of CSV passed to Analysis_DPV (first 10 lines):")
        with open(file_to_pass, encoding='utf-16') as fh:
            for _ in range(10):
                line = fh.readline()
                if not line:
                    break
                print(line.rstrip())
        print(">>> END PREVIEW\n")
        # -----------------------------------------------------------------------

        # finally call the unchanged DPV loader
        self.analysis_results = Analysis_DPV(
            file_to_pass,
            metadata_fields=custom_fields,
            peak_regions=self.peak_regions
        )

        if not self.analysis_results:
            messagebox.showwarning(
                "Load Failure",
                "Could not analyze the selected file.\n"
                "Please check the file format and metadata."
            )
        else:
            self.save_as_button.config(state=tk.NORMAL)
        # ─────────────────────────────────  END OF BLOCK  ─────────────────────────────

    def import_blank_file(self):
        """
        Prompt user for a blank responses file, load it,
        then re-run Analysis_DPV with blank_responses so LOD is updated.
        """
        self.blank_filepath = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx;*.xls")]
        )
        if self.blank_filepath:
            self.blank_file_path_label.config(text=self.blank_filepath)

            # If we already have a main data file, let's re-run analysis with blank
            if self.filepath:
                blank_responses = self.load_blank_responses(self.blank_filepath)
                self.analysis_results = Analysis_DPV(
                    self.filepath,
                    blank_responses=blank_responses
                )
                if not self.analysis_results:
                    messagebox.showwarning(
                        "Blank Analysis Failed",
                        "Could not analyze data with the provided blank. Check your file."
                    )
                else:
                    messagebox.showinfo(
                        "Blank Loaded",
                        "Blank responses file imported and analysis updated for LOD."
                    )
                    # Optionally enable "Save Results" if you'd like
                    self.save_as_button.config(state=tk.NORMAL)
            else:
                messagebox.showinfo(
                    "No Main File",
                    "You selected a blank file, but haven't imported a main file yet."
                )
        else:
            self.blank_file_path_label.config(text="No blank responses file selected")


    def change(self, *args):
        # Logic to handle changes in dropdown, if needed
        pass

    def convert_concentrations(self):
        """
        Convert the concentration in parsed_metadata from the chosen
        'from_unit' to 'to_unit', if the user picked something.
        """

        def unify_micro_symbol(unit_str):
            """
            If unit_str is 'uM', rewrite it to 'µM' so that
            everything in the code sees 'µM' as the standard micro symbol.
            """
            if unit_str == 'uM':
                return 'µM'
            return unit_str

        def convert_concentration(value, from_unit, to_unit):
            """
            Convert 'value' (float) from 'from_unit' to 'to_unit'.
            Return the converted float.
            """

            # First, unify the from_unit and to_unit (so 'uM' becomes 'µM')
            from_unit = unify_micro_symbol(from_unit)
            to_unit   = unify_micro_symbol(to_unit)

            # If either from_unit or to_unit is blank, or they are the same, do nothing
            if not from_unit or not to_unit or (from_unit == to_unit):
                return value

            # Dictionary for factor from each unit to M
            to_m_factor = {
                'nM': 1e-9,
                'µM': 1e-6,   # We can remove 'uM' if we unify them anyway
                'mM': 1e-3,
                'M':  1.0
            }

            # If either unit is missing from the dictionary, skip
            if from_unit not in to_m_factor or to_unit not in to_m_factor:
                return value

            # Convert from from_unit => M
            in_molar = value * to_m_factor[from_unit]

            # Convert from M => to_unit
            result = in_molar / to_m_factor[to_unit]
            return result

        from_unit = self.from_unit_var.get().strip()
        to_unit   = self.to_unit_var.get().strip()

        if not self.analysis_results:
            messagebox.showwarning(
                "No Data", 
                "Please run or load data first before converting concentrations."
            )
            return

        parsed_metadata = self.analysis_results.get('parsed_metadata', None)
        if parsed_metadata is None:
            messagebox.showwarning(
                "No Metadata", 
                "No parsed_metadata found in analysis results."
            )
            return

        # If user didn't pick or if same unit, do nothing
        if (not from_unit) or (not to_unit) or (from_unit == to_unit):
            messagebox.showinfo("Info", "No valid conversion was selected.")
            return

        import re

        count_converted = 0
        for meta in parsed_metadata:
            # e.g. meta['Concentration'] = "50uM", "100µM", "10.0"
            conc_str = meta['Concentration']

            # Extract numeric portion
            match = re.match(r"(\d+\.?\d*)", conc_str)
            if not match:
                continue  # skip if we can't parse

            numeric_val = float(match.group(1))

            # Now do the conversion
            new_val = convert_concentration(numeric_val, from_unit, to_unit)

            # Store as e.g. "5.0µM"
            meta['Concentration'] = f"{new_val:.3g}{to_unit}"
            count_converted += 1

        messagebox.showinfo(
            "Conversion Complete",
            f"Converted {count_converted} metadata entries from {from_unit} to {to_unit}."
        )
        
    def generate_peak_report_page(self, save_path=None, selected_indices=None):
        """
        Generate a single PDF or PNG page like the OriginLab report
        with:
        - title and metadata
        - corrected DPV plot
        - table of peaks
        """

        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        from matplotlib import image as mpimg
        import matplotlib.patches as patches

        # ───────────────────────────────────────────────
        # STEP 1 - Generate plot
        # ───────────────────────────────────────────────

        temp_dir = tempfile.mkdtemp()
        plot_path = os.path.join(temp_dir, "corrected_plot.png")

        # Create the corrected baseline plot as a separate image file
        self.plot_data_array_with_corrected_baseline(
            self.analysis_results['data_array'],
            self.analysis_results['headers'],
            self.analysis_results['parsed_metadata'],
            plot_path,
            selected_indices
        )

        # Read that plot back in as an image
        plot_img = mpimg.imread(plot_path)

        # ───────────────────────────────────────────────
        # STEP 2 - Create new figure for the whole report
        # ───────────────────────────────────────────────

        fig = plt.figure(figsize=(8.5, 11))  # US letter size

        # Three rows: Title, Plot, Table
        gs = gridspec.GridSpec(
            3, 1,
            height_ratios=[0.15, 0.55, 0.30],
            hspace=0.5
        )

        # ───────── TITLE PANEL ─────────
        ax_title = plt.subplot(gs[0])
        ax_title.axis("off")

        title_text = (
            "Peak Analysis Report\n\n"
            f"File: {os.path.basename(self.filepath)}\n"
            f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        )
        ax_title.text(
            0.0, 1.0,
            title_text,
            fontsize=10,
            va='top',
            ha='left',
            family='monospace'
        )

        # ───────── PLOT PANEL ─────────
        ax_img = plt.subplot(gs[1])
        ax_img.imshow(plot_img)
        ax_img.axis("off")

        # Add a border around the image for clarity
        ax_img.add_patch(
            patches.Rectangle(
                (0, 0),
                plot_img.shape[1],
                plot_img.shape[0],
                linewidth=1,
                edgecolor='black',
                facecolor='none'
            )
        )

        # ───────── TABLE PANEL ─────────
        ax_table = plt.subplot(gs[2])
        ax_table.axis("off")

        # If user picked plots, filter peaks accordingly
        peak_table_all = self.analysis_results.get("results", [])

        if selected_indices is None:
            peak_table = peak_table_all
        else:
            # Only keep peaks belonging to selected curves
            selected_curve_indices = set(idx // 2 for idx in selected_indices)
            peak_table = [
                row for row in peak_table_all
                if row['Entry'] in selected_curve_indices
            ]


        if not peak_table:
            table_text = "No peaks found."
            ax_table.text(
                0, 1,
                table_text,
                va='top',
                ha='left',
                fontsize=10
            )
        else:
            # Construct table rows
            headers = [
                "Index", "Analyte", "Conc.",
                "Peak Voltage (V)", "Peak Height (µA)", "FWHM (V)"
            ]

            rows = []
            for idx, row in enumerate(peak_table, start=1):
                analyte = row.get("Analytes", "")
                conc = row.get("Concentration", "")
                voltage = f"{row.get('Peak Voltage', float('nan')):.4g}"
                height = f"{row.get('Adjusted Peak Current', float('nan')):.4g}"
                fwhm = row.get("FWHM")
                fwhm_str = f"{fwhm:.5f}" if fwhm is not None and not np.isnan(fwhm) else "-"

                rows.append([
                    str(idx),
                    analyte,
                    conc,
                    voltage,
                    height,
                    fwhm_str
                ])

            # Build table
            table = ax_table.table(
                cellText=rows,
                colLabels=headers,
                loc='center',
                cellLoc='center',
                colLoc='center'
            )
            table.auto_set_font_size(False)
            table.set_fontsize(9)

            # Adjust column widths
            for k, cell in table.get_celld().items():
                cell.set_edgecolor('black')
                cell.set_linewidth(0.3)

            table.scale(1.2, 1.5)

        # ───────── Save the entire page ─────────
        if save_path is None:
            save_path = os.path.join(tempfile.gettempdir(), "peak_report.pdf")

        fig.tight_layout()
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

        print(f"Peak report saved to: {save_path}")
        messagebox.showinfo("Report Saved", f"Peak report saved to:\n{save_path}")


    def perform_peak_report_analysis(self, filepath, selected_indices=None):
        """
        Generates a single-page PDF/PNG report.
        """
        if not self.analysis_results:
            messagebox.showwarning(
                "No Data",
                "Please import a file first. No analysis results available."
            )
            return

        # Let user choose where to save
        filetypes = [("PDF files", "*.pdf"), ("PNG files", "*.png")]
        save_path = filedialog.asksaveasfilename(
            title="Save Peak Report",
            defaultextension=".pdf",
            filetypes=filetypes
        )
        if not save_path:
            return

        # Pass selected_indices to the report generator
        self.generate_peak_report_page(save_path, selected_indices)

    def run_analysis(self):
        if not self.filepath:
            messagebox.showwarning("Warning", "Please select a file first.")
            return

        analysis_option = self.clicked.get()

        if analysis_option == "DPV analysis":
            # Use existing data or do a specialized step
            self.perform_dpv_analysis(self.filepath)

        elif analysis_option in [
            "Plot Data Array",
            "Plot Data Array with Corrected Baseline",
            "Plot Concentrations vs Mean Peak Currents",
            "Analyze Peak Currents",
            "Observed vs Expected Concentration",
            "Residual vs Potential Plot",
        ]:
            # Now, do not re-run. Just ensure we have something in self.analysis_results
            if not self.analysis_results:
                messagebox.showwarning("No Data", "Please import a file first.")
                return

            # open the plot selection window based on the analysis option
            self.open_plot_selection_window(
                analysis_option,
                self.analysis_results['headers'],
                self.analysis_results['data_array']
            )

        elif analysis_option == "LOD Analysis":
            self.perform_lod_analysis(self.filepath)

        elif analysis_option == "T-test Analysis":
            self.perform_t_test_analysis(self.filepath)

        elif analysis_option == "Peak Deconvolution":
            # open the same plot-picker UI
            if not self.analysis_results:
                messagebox.showwarning("No Data", "Please import a file first.")
                return
            self.open_plot_selection_window(
                analysis_option,
                self.analysis_results['headers'],
                self.analysis_results['data_array']
            )

        elif analysis_option == "Generate Peak Report":
            # Open the plot selection window for the report
            self.open_plot_selection_window(
                analysis_option,
                self.analysis_results['headers'],
                self.analysis_results['data_array']
            )

        elif analysis_option == "Residual vs Potential Plot":
            # Open the plot selection window for the report
            self.open_plot_selection_window(
                analysis_option,
                self.analysis_results['headers'],
                self.analysis_results['data_array']
            )

        else:
            print(f"Running {analysis_option} on file {self.filepath}")
            # Implement call to other analysis functions here

        # Enable "Save Results" if we have data
        if self.analysis_results:
            self.save_as_button.config(state=tk.NORMAL)


    def perform_dpv_analysis(self, filepath):
        """
        Formerly ran Analysis_DPV; now we just do something special with the
        already-loaded data, e.g. display results.
        """
        if not self.analysis_results:
            messagebox.showwarning(
                "No Data",
                "Please import a file first. No analysis results available."
            )
            return

        # We already have self.analysis_results from import_file
        self.display_results(self.analysis_results)
        return self.analysis_results

    def perform_lod_analysis(self, filepath):
        """
        If we re-run inside import_blank_file, we don't need to call
        Analysis_DPV again. Just do anything special needed for LOD display.
        """
        if not self.analysis_results:
            messagebox.showwarning("No Data", "No analysis results to display.")
            return

        # Now just display the LOD portion of self.analysis_results
        self.display_lod_results(self.analysis_results)
        return self.analysis_results

    def perform_t_test_analysis(self, filepath):
        """
        Similarly, do not re-run Analysis_DPV. Just use self.analysis_results.
        """
        if not self.analysis_results:
            messagebox.showwarning(
                "No Data",
                "Please import a file first. No analysis results available."
            )
            return

        # Possibly do extra T-test steps if needed
        self.display_t_test_results(self.analysis_results)
        return self.analysis_results
    
    def _windows_for_analyte(self, analyte_label: str):
        """
        Return the window‐list for this analyte, or for the special key 'all'
        (case‐insensitive) if present.  This lets you write ALL:0.65-0.90
        and have it apply to every curve.
        """
        if not self.peak_regions:
            return None

        # 1) case‐insensitive 'all' catch-all
        for key, winlist in self.peak_regions.items():
            if key.lower() == "all":
                return winlist

        # 2) otherwise, look for a matching analyte key
        for key, winlist in self.peak_regions.items():
            if key.lower() in analyte_label.lower():
                return winlist

        return None
    
    def _show_deconv_results_table(self, peak_dicts):
        """
        Render a ttk.Treeview table:
            Entry | Analyte | Peak | Center (V) | Amplitude (µA) | Width (V)

        Each row corresponds to ONE Gaussian (so every curve produces 2 rows).
        """
        import tkinter as tk
        from tkinter import ttk

        if not peak_dicts:
            tk.messagebox.showinfo("No data", "Nothing to show.")
            return

        # ─────── create window & basic widgets ──────────────────────────────────
        win = tk.Toplevel(self.master)
        win.title("Peak-deconvolution table")

        frame = ttk.Frame(win, padding=4)
        frame.pack(fill="both", expand=True)

        # Treeview with scrollbars
        cols = ("entry", "analyte", "peak", "center", "amp", "width")
        tree = ttk.Treeview(frame, columns=cols, show="headings", height=12)

        # prettier column headers
        hdr_text = {
            "entry": "Entry",
            "analyte": "Analyte",
            "peak": "Peak #",              # thin space avoids wrapping
            "center": "Center (V)",
            "amp": "Amplitude (µA)",
            "width": "Width (V)",
        }
        for c in cols:
            tree.heading(c, text=hdr_text[c])
            tree.column(c, anchor="center", width=110 if c in ("center", "amp", "width") else 80)

        vsb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")

        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)

        # ─────── populate rows: one per Gaussian ────────────────────────────────
        for rec in peak_dicts:
            tree.insert(
                "", "end",
                values=(
                    rec["Entry"],
                    rec["Analytes"],
                    1,
                    f"{rec['Center1']:.5f}",
                    f"{rec['Amp1']:.5g}",
                    f"{rec['Width1']:.5f}",
                )
            )
            tree.insert(
                "", "end",
                values=(
                    rec["Entry"],
                    rec["Analytes"],
                    2,
                    f"{rec['Center2']:.5f}",
                    f"{rec['Amp2']:.5g}",
                    f"{rec['Width2']:.5f}",
                )
            )


    def perform_peak_deconvolution(self, filepath, selected_indices):
        """
        Peak deconvolution using lmfit's two-component Gaussian model.
        This automatically handles two peaks even if the raw peak finder only sees one.
        """
        import numpy as np
        import matplotlib.pyplot as plt
        import tkinter as tk
        from lmfit import Model
        from scipy.signal import savgol_filter
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from tkinter import Toplevel, messagebox
        from scipy import sparse
        from scipy.sparse.linalg import spsolve

        data = self.analysis_results['data_array']
        meta = self.analysis_results['parsed_metadata']
        regions = self.peak_regions

        # Baseline correction using asymmetric least squares (ALS) method
        def baseline_als(y, lam=1e5, p=0.01, niter=10):
            L = len(y)
            D = sparse.csc_matrix(np.diff(np.eye(L), 2))
            w = np.ones(L)
            for _ in range(niter):
                W = sparse.spdiags(w, 0, L, L)
                z = spsolve(W + lam * D.dot(D.T), w * y)
                w = p * (y > z) + (1 - p) * (y < z)
            return z

        # Build a combined two-Gaussian model using lmfit
        gauss1 = Model(lambda x, amp1, cen1, wid1: amp1 * np.exp(-(x-cen1)**2/(2*wid1**2)), prefix='g1_')
        gauss2 = Model(lambda x, amp2, cen2, wid2: amp2 * np.exp(-(x-cen2)**2/(2*wid2**2)), prefix='g2_')
        model = gauss1 + gauss2

        out = []
        for idx in selected_indices:
            V = data[:, idx]
            I = data[:, idx + 1]
            mask = ~np.isnan(V) & ~np.isnan(I)
            V, I = V[mask], I[mask]

            # Apply optional window cropping for the analyte
            wins = self._windows_for_analyte(meta[idx//2]['Analytes']) if regions else None
            if regions and wins is None:
                # Skip this curve if windows are defined but none for this analyte
                continue
            if wins:
                sel_mask = np.logical_or.reduce([(V >= lo) & (V <= hi) for lo, hi in wins])
                V, I = V[sel_mask], I[sel_mask]
                if V.size == 0:
                    continue  # no data in the specified window

            if V.size < 5:
                messagebox.showwarning("Data Error", "Not enough points in region.")
                continue

            # Baseline correction using ALS
            baseline = baseline_als(I)
            y = I - baseline
            y -= y.min()  # shift baseline to zero

            # Smooth the baseline-corrected data for initial peak estimate
            ys = savgol_filter(y, 7, 3)
            # Determine a midpoint to initialze two peak centers on left/right halves
            mid = 0.5 * (V.max() + V.min())
            left_mask = V < mid
            right_mask = V >= mid
            # Initial guesses for peak centers: max position in each half (fallback to edges)
            c1_init = V[left_mask][np.argmax(ys[left_mask])] if left_mask.any() else V[0]
            c2_init = V[right_mask][np.argmax(ys[right_mask])] if right_mask.any() else V[-1]
            # Initial amplitude guesses at those centers (or half max as fallback)
            amp1_init = ys[V == c1_init][0] if c1_init in V else ys.max() / 2
            amp2_init = ys[V == c2_init][0] if c2_init in V else ys.max() / 2
            # Initial width guess: one-sixth of the total span
            wid_init = (V.max() - V.min()) / 6

            # Set up parameters with bounds
            params = model.make_params(
                g1_amp1=amp1_init, g1_cen1=c1_init, g1_wid1=wid_init,
                g2_amp2=amp2_init, g2_cen2=c2_init, g2_wid2=wid_init
            )
            # Enforce ordering: g1 center on left half, g2 center on right half
            params['g1_cen1'].min = V.min();  params['g1_cen1'].max = mid
            params['g2_cen2'].min = mid;     params['g2_cen2'].max = V.max()
            # Ensure widths are positive
            params['g1_wid1'].min = 1e-6;    params['g2_wid2'].min = 1e-6

            # Perform the two-Gaussian fit
            try:
                result = model.fit(y, params, x=V, method='leastsq')
            except Exception as e:
                messagebox.showwarning("Fit Error", str(e))
                continue

            # Extract fitted parameter values
            r = result.params
            out.append({
                'Entry': idx // 2,
                'Analytes': meta[idx // 2]['Analytes'],
                'Center1': r['g1_cen1'].value, 'Amp1': r['g1_amp1'].value, 'Width1': r['g1_wid1'].value,
                'Center2': r['g2_cen2'].value, 'Amp2': r['g2_amp2'].value, 'Width2': r['g2_wid2'].value,
            })

            # Plot the baseline-corrected data and the two Gaussian components
            fig, ax = plt.subplots(figsize=(5, 3))
            analyte_name = meta[idx // 2]['Analytes']
            concentration = meta[idx // 2].get('Concentration', '')
            label_data = f"{analyte_name} – {concentration}" if concentration else analyte_name
            ax.plot(V, y, 'b-', label=label_data)
            ax.plot(V, result.eval_components(x=V)['g1_'], 'r--', label='Gaussian 1')
            ax.plot(V, result.eval_components(x=V)['g2_'], 'g--', label='Gaussian 2')
            ax.set_xlabel('Potential (V)')
            ax.set_ylabel('Current (µA)')
            ax.legend()
            # Open a new window for this deconvolution plot
            win = Toplevel(self.master)
            win.title(f"Peak Deconvolution: {analyte_name} – {concentration}")
            FigureCanvasTkAgg(fig, master=win).get_tk_widget().pack(fill='both', expand=True)

        # # If any results were collected, show them in a summary window
        # if out:
        #     hdr = list(out[0].keys())
        #     summary = 'Peak Deconvolution Results\n\n' + '  '.join(hdr) + '\n' + '-' * 60 + '\n'
        #     for row in out:
        #         summary += '  '.join(str(row[k]) for k in hdr) + '\n'
        #     result_win = Toplevel(self.master)
        #     result_win.title('Deconvolution Results')
        #     txt = tk.Text(result_win, width=80, height=10, wrap=tk.NONE)
        #     txt.insert('1.0', summary)
        #     txt.pack(fill='both', expand=True)

        if out:
            self._show_deconv_results_table(out)

    def save_results(self):
            if not self.analysis_results:
                messagebox.showerror("No Results", "No results to save.")
                return

            # Step 1: Prompt the user to select a directory to save the file
            save_directory = filedialog.askdirectory()
            if not save_directory:
                messagebox.showwarning("No Directory Selected", "Please select a directory to save the results.")
                return

            # Step 2: Prompt the user to enter just the file name
            file_name = tk.simpledialog.askstring("File Name", "Enter the file name:")
            if not file_name:
                messagebox.showwarning("No Filename", "Please enter a file name to save the results.")
                return

            # Construct the full file path by combining directory, file name, and .csv extension
            save_path = os.path.join(save_directory, f"{file_name}.csv")

            # Step 3: Save the results in the selected directory with the constructed file name
            if 'mean_peak_currents' in self.analysis_results:
                with open(save_path, 'w') as f:
                    f.write("Mean Peak Currents,Standard Deviation, Coefficient of Variation\n")
                    for key in self.analysis_results['mean_peak_currents'].keys():
                        mean_current = self.analysis_results['mean_peak_currents'][key]
                        std_dev = self.analysis_results['std_peak_currents'].get(key, "N/A")
                        cov = self.analysis_results['cov_peak_currents'].get(key, "N/A")
                        f.write(f"{key},{mean_current},{std_dev},{cov}\n")
                messagebox.showinfo("Results Saved", f"Results saved to {save_path}")

            # Step 4: Save plots if there are any
            if 'plot_filenames' in self.analysis_results:
                for key, plot_files in self.analysis_results['plot_filenames'].items():
                    if isinstance(plot_files, list):  # Handle multiple files
                        for i, plot_file in enumerate(plot_files):
                            plot_save_path = os.path.join(save_directory, f"{file_name}_{key}_{i+1}.png")
                            os.rename(plot_file, plot_save_path)
                    else:  # Single file case
                        plot_save_path = os.path.join(save_directory, f"{file_name}_{key}.png")
                        os.rename(plot_files, plot_save_path)
                messagebox.showinfo("Plots Saved", "Plots have been saved to the selected directory.")


    def load_blank_responses(self, filepath):
        # Placeholder: Load blank responses from file (implement as needed)
        # Here we assume it loads a list or array of blank responses
        # For example, assuming the file is a CSV with one column of blank responses
        return np.loadtxt(filepath, delimiter=',')  # Example placeholder

    def open_plot_selection_window(self, analysis_option, headers, data_array):
        plot_selection_window = Toplevel(self.master)
        plot_selection_window.title("Select Plots")
        plot_selection_window.geometry("350x500+{}+{}".format(self.master.winfo_x() - 350, self.master.winfo_y()))

        tk.Label(plot_selection_window, text="Select the plots you want to display:").pack(anchor='w')

        checkbox_frame = tk.Frame(plot_selection_window)
        checkbox_frame.pack(fill='both', expand=True)

        self.plot_vars = []
        num_plots = len(range(0, len(headers), 2))

        # Initialize select_all_var and ensure it's set to 1 (checked) by default
        self.select_all_var = tk.IntVar(value=1)
        select_all_checkbox = tk.Checkbutton(checkbox_frame, text="Select All", variable=self.select_all_var, command=self.toggle_all)
        select_all_checkbox.grid(row=0, column=0, sticky="w")

        # Set all individual checkboxes to checked by default
        for i in range(num_plots):
            var = tk.IntVar(value=1)  # Initialize as checked
            self.plot_vars.append(var)
            checkbox = tk.Checkbutton(checkbox_frame, text=f'Plot {i+1}', variable=var)
            checkbox.grid(row=i % 15 + 1, column=i // 15, sticky="w", padx=5)

        tk.Button(plot_selection_window, text="Plot Selected", command=lambda: self.plot_selected(analysis_option)).pack(pady=10)

    def toggle_all(self):
        # Check if select_all_var is checked
        if self.select_all_var.get() == 1:
            # Check all plot checkboxes
            for var in self.plot_vars:
                var.set(1)
        else:
            # Uncheck all plot checkboxes
            for var in self.plot_vars:
                var.set(0)

    def display_results(self, results):
        """
        Summarise each distinguishable peak (Peak # 1, 2, …) for every
        (Analyte, Concentration) pair across replicates.

        Columns:
            Analyte | Conc. | Peak # | Peak Max (µA) | Peak Pos (V) |
            FWHM (V) | σ (µA) | RSD %
        """

        if not results or 'results' not in results:
            messagebox.showinfo("No data", "Nothing to show.")
            return

        # ────────── helper: aggregate replicate data by peak rank ────────────
        # Build { (conc, analyte) : [ [peak1 dict, peak2 …]  ← replicate A,
        #                             [peak1 dict, …]        ← replicate B, … ] }
        grouped = defaultdict(list)
        for rec in results['results']:
            key = (rec['Concentration'], rec['Analytes'])
            grouped[key].append(rec)

        # We need the replicate structure (rows grouped by Entry) so we re-nest:
        per_entry = defaultdict(lambda: defaultdict(list))
        for rec in results['results']:
            key = (rec['Concentration'], rec['Analytes'])
            per_entry[key][rec['Entry']].append(rec)

        # Sort each replicate’s peaks left→right (ascending voltage)
        for key in per_entry:
            for entry in per_entry[key]:
                per_entry[key][entry].sort(key=lambda r: r['Peak Voltage'])

        # ────────── build rows for the table ─────────────────────────────────
        rows = []   # each is (Analyte, Conc, idx, mean_amp, mean_pos, mean_w, σ, RSD)
        for key, entry_dict in per_entry.items():
            conc, analyte = key
            replicates = list(entry_dict.values())
            max_peaks = max(len(lst) for lst in replicates)

            for idx in range(max_peaks):        # Peak # idx+1
                amps  = [lst[idx]['Adjusted Peak Current']
                        for lst in replicates if idx < len(lst)]
                volts = [lst[idx]['Peak Voltage']
                        for lst in replicates if idx < len(lst)]
                widths = [lst[idx]['FWHM'] for lst in replicates
                        if idx < len(lst) and lst[idx]['FWHM'] is not None]

                if not amps:          # no data for this rank in any replicate
                    continue

                mean_amp = np.mean(amps)
                std_amp  = np.std (amps, ddof=0)
                rsd_amp  = 100 * std_amp / mean_amp if mean_amp else np.nan

                mean_pos = np.mean(volts)
                mean_w   = np.mean(widths) if widths else np.nan

                rows.append((
                    analyte, conc, idx+1,
                    mean_amp, mean_pos, mean_w,
                    std_amp, rsd_amp
                ))

        # ────────── Tk window & Treeview set-up ──────────────────────────────
        win = tk.Toplevel(self.master)
        win.title("DPV Peak Summary")

        frame = ttk.Frame(win, padding=6)
        frame.pack(fill="both", expand=True)

        cols = ("analyte", "conc", "peak",
                "max", "pos", "fwhm", "std", "rsd")
        hdr = {
            "analyte": "Analyte",
            "conc"   : "Conc.",
            "peak"   : "Peak #",
            "max"    : "Peak Max (µA)",
            "pos"    : "Peak Pos (V)",
            "fwhm"   : "FWHM (V)",
            "std"    : "σ (µA)",
            "rsd"    : "RSD (%)",
        }

        tree = ttk.Treeview(frame, columns=cols, show="headings", height=14)
        for c in cols:
            tree.heading(c, text=hdr[c])
            w = 120 if c in ("max", "pos", "fwhm", "std", "rsd") else 90
            tree.column(c, anchor="center", width=w)

        vsb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")

        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)

        # ────────── insert rows ──────────────────────────────────────────────
        fmt = lambda x, digits=5: ("—" if x is None or np.isnan(x)
                                else f"{x:.{digits}g}")
        for a, c, p, m, v, w, s, r in rows:
            tree.insert(
                "", "end",
                values=(a, c, p,
                        fmt(m, 5), fmt(v, 5), fmt(w, 5),
                        fmt(s, 5), fmt(r, 4))
            )



    def display_lod_results(self, results):
        # Display the LOD results in a new window or a suitable UI component
        result_window = Toplevel(self.master)
        result_window.title("LOD Analysis Results")
        
        text_area = tk.Text(result_window, wrap=tk.WORD, width=80, height=20)
        text_area.pack(padx=10, pady=10)

        lod_results = results.get('lod_results', {})
        result_text = "LOD Analysis Results\n\n"

        for key, value in lod_results.items():
            result_text += f"{key}: {value}\n"
        
        text_area.insert(tk.END, result_text)

    def display_t_test_results(self, results):
        # Display the t-test results in a new window or a suitable UI component
        result_window = Toplevel(self.master)
        result_window.title("T-test Analysis Results")
        
        text_area = tk.Text(result_window, wrap=tk.WORD, width=80, height=20)
        text_area.pack(padx=10, pady=10)

        t_test_results = results.get('t_test_results', [])
        result_text = "T-test Analysis Results\n\n"

        if isinstance(t_test_results, list):
            for i, result in enumerate(t_test_results, 1):
                result_text += f"Test {i}: {result}\n"
        else:
            result_text += "No t-test results found."

        text_area.insert(tk.END, result_text)

    def plot_selected(self, analysis_option):
        selected_indices = [i * 2 for i, var in enumerate(self.plot_vars) if var.get()]
        print(f"Selected plots for {analysis_option}: {selected_indices}")

        if (analysis_option == "Plot Data Array with Corrected Baseline"
                and self.peak_regions):                 # user entered at least one window
            all_regions = self.peak_regions.copy()      # keep full dict for later restore

            for key in all_regions:                     # e.g. 'HX', then 'UA', …
                print(f"→ plotting only '{key}' window")
                self.peak_regions = {key: all_regions[key]}   # TEMP: 1-key dict
                # call the usual single-plot helper – it now sees only that window
                self.display_and_save_plot(
                    self.plot_data_array_with_corrected_baseline,
                    selected_indices                       # keep user’s curve selection
                )

            self.peak_regions = all_regions              # restore original dict
            return                                       # done – skip generic section

        if analysis_option == "Plot Data Array":
            self.display_and_save_plot(self.plot_data_array, selected_indices)
        elif analysis_option == "Plot Data Array with Corrected Baseline":
            self.display_and_save_plot(self.plot_data_array_with_corrected_baseline, selected_indices)
        elif analysis_option == "Plot Concentrations vs Mean Peak Currents":
            self.display_and_save_multiple_plots(self.plot_concentrations_vs_mean_peak_currents, selected_indices)
        elif analysis_option == "Analyze Peak Currents":
            self.display_and_save_multiple_plots(
                self.analyze_peak_currents,
                selected_indices,
                self.analysis_results['mean_peak_currents'],
                self.analysis_results['std_peak_currents']
            )
        elif analysis_option == "Observed vs Expected Concentration":
            # Force None so we always plot all peaks
            self.display_and_save_multiple_plots(
                self.plot_observed_vs_expected_concentration,
                None,
                self.analysis_results['mean_peak_currents']
            )
        elif analysis_option == "Peak Deconvolution":
            # calls your new deconvolution routine, passing the user’s curve picks:
            self.perform_peak_deconvolution(self.filepath, selected_indices)
        elif analysis_option == "Generate Peak Report":
            # For the report, store the selected indices and continue
            self.perform_peak_report_analysis(
                self.filepath,
                selected_indices
            )
        elif analysis_option == "Residual vs Potential Plot":
            # call the routine directly (or wrap with display_and_save_plot
            # if you want it auto-saved in the same folder)
            self.plot_residual_vs_peak_location(selected_indices)
        else:
            print(f"Selected plots not handled for: {analysis_option}")

    def display_and_save_plot(self, plot_func, selected_indices):
        temp_dir = tempfile.mkdtemp()
        save_path = os.path.join(temp_dir, "plot.png")
        
        # Call the plot function
        plot_func(self.analysis_results['data_array'], self.analysis_results['headers'], self.analysis_results['parsed_metadata'], save_path, selected_indices)
        
        if isinstance(save_path, str):
            try:
                img = Image.open(save_path)
                img.show()
            except FileNotFoundError as e:
                print(f"File not found: {e}")
        else:
            print(f"Error: save_path is not a string but a {type(save_path)}")

    def display_and_save_multiple_plots(self, plot_func, selected_indices, *args):
        # Create a temporary directory to save the plot images
        temp_dir = tempfile.mkdtemp()

        # Call the plot function to generate and save the plots
        plot_filenames = plot_func(*args, temp_dir, selected_indices)

        # Display each plot
        for plot_file in plot_filenames:
            img = Image.open(plot_file)
            img.show()

    # --------------------------------------------------------------------------
    # 1)  RAW DATA PLOT
    # --------------------------------------------------------------------------
    def plot_data_array(
            self, data_array, headers, parsed_metadata,
            save_path, selected_indices=None):

        import matplotlib.pyplot as plt          # local import is fine
        plt.figure(figsize=(10, 6))

        if selected_indices is None:
            selected_indices = range(0, len(headers), 2)

        for i in selected_indices:
            potential = data_array[:, i]
            current   = data_array[:, i + 1]

            # drop NaNs ----------------------------------------------------------
            keep = ~np.isnan(potential) & ~np.isnan(current)
            potential = potential[keep]
            current   = current[keep]

            # -------- NEW: crop to user-defined window(s) -----------------------
            meta     = parsed_metadata[i // 2]
            analyte  = meta['Analytes']
            #winlist  = self.peak_regions.get(analyte) if getattr(self, "peak_regions", None) else None
            winlist = self._windows_for_analyte(analyte)

            # ----------------- DEBUG -----------------
            print(f"[{i//2:02}] {analyte:<15}  winlist={winlist}")
            # -----------------------------------------

            if self.peak_regions and winlist is None:
                continue                 # ← jump to next curve

            if winlist:
                mask = np.logical_or.reduce(
                    [(potential >= lo) & (potential <= hi) for lo, hi in winlist]
                )

                # ----- DEBUG #2 -----
                print(f"      points before={potential.size:4}  after={mask.sum():4}")
                # --------------------

                potential = potential[mask]
                current   = current[mask]
                if potential.size == 0:          # nothing falls inside → skip
                    continue
            # -------------------------------------------------------------------

            label = f"{analyte} – {meta['Concentration']}"
            plt.plot(potential, current, label=label)

        plt.xlabel('Potential (V)')
        plt.ylabel('Current (µA)')
        plt.title('Overlapping Plots of Potential vs Current')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close()


    # --------------------------------------------------------------------------
    # 2)  BASELINE-CORRECTED PLOT
    # --------------------------------------------------------------------------
    def plot_data_array_with_corrected_baseline(
            self, data_array, headers, parsed_metadata,
            save_path, selected_indices=None):

        import matplotlib.pyplot as plt
        from scipy import sparse
        from scipy.sparse.linalg import spsolve
        plt.figure(figsize=(10, 6))

        if selected_indices is None:
            selected_indices = range(0, len(headers), 2)

        def baseline_als(y, lam=1e5, p=0.01, niter=10):
            L = len(y)
            D = sparse.csc_matrix(np.diff(np.eye(L), 2))
            w = np.ones(L)
            for _ in range(niter):
                W = sparse.spdiags(w, 0, L, L)
                z = spsolve(W + lam * D.dot(D.T), w * y)
                w = p * (y > z) + (1 - p) * (y < z)
            return z

        for i in selected_indices:
            potential = data_array[:, i]
            current   = data_array[:, i + 1]

            keep = ~np.isnan(potential) & ~np.isnan(current)
            potential = potential[keep]
            current   = current[keep]
            if potential.size == 0:
                continue

            # -------- NEW: crop to user-defined window(s) -----------------------
            meta     = parsed_metadata[i // 2]
            analyte  = meta['Analytes']
            winlist = self._windows_for_analyte(analyte)

            # ----------------- DEBUG -----------------
            print(f"[{i//2:02}] {analyte:<15}  winlist={winlist}")
            # -----------------------------------------

            if self.peak_regions and winlist is None:
                continue      # skip unrelated analytes

            if winlist:
                mask = np.logical_or.reduce(
                    [(potential >= lo) & (potential <= hi) for lo, hi in winlist]
                )

                # ----- DEBUG #2 -----
                print(f"      points before={potential.size:4}  after={mask.sum():4}")
                # --------------------

                potential = potential[mask]
                current   = current[mask]
                if potential.size == 0:
                    continue
            # -------------------------------------------------------------------

            baseline          = baseline_als(current)
            corrected_current = current - baseline

            label = f"{analyte} – {meta['Concentration']}"
            plt.plot(potential, corrected_current, label=label)

        plt.xlabel('Potential (V)')
        plt.ylabel('Current (µA)')
        plt.title('Corrected Plots of Potential vs Current')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
                fontsize='small', ncol=2)
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

            save_path = os.path.join(temp_dir, f"{base_save_path}_{analytes}.png")
            plt.savefig(save_path)
            plt.close()
            plot_filenames.append(save_path)

        return plot_filenames

    def analyze_peak_currents(self,
                            mean_peak_currents: dict,
                            std_peak_currents: dict,
                            base_save_path: str,
                            selected_indices=None):

        print("↳  Keys arriving in mean_peak_currents:")
        for k in list(mean_peak_currents)[:10]:
            print("   ", k)
        print("   …")

        def numeric_conc(s: str | None) -> float | None:
            """
            Return the first *decimal* number in the string, ignoring the day part in
            labels such as 'HX_D29_60µM'.
            """
            if not s:
                return None
            # strip a leading 'Dxx_' if present e.g. D29_60µM → 60µM
            s = re.sub(r'^D\d+[._]', '', s, flags=re.I)      # NEW
            s = s.replace(',', '.')                          # allow comma decimals
            m  = re.search(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', s)
            return float(m.group(0)) if m else None


        analytes_dict: dict[str, list] = defaultdict(list)

        for key, mean_pk in mean_peak_currents.items():
            conc_str, analyte = key
            c_val = numeric_conc(conc_str)
            # new guard: ignore 0 or negative concentrations
            if c_val is None or c_val <= 0:
                print(f"⚠  skipped row – invalid conc '{conc_str}' → {c_val}")
                continue

            std_err = std_peak_currents.get(key, [0] * len(mean_pk))
            analytes_dict[analyte].append((c_val, mean_pk, std_err))

        # ---- 3) fit & plot -----------------------------------------------
        calib_models = {}
        plot_filenames = []
        temp_dir = base_save_path                    # folder that came in

        for analyte, rows in analytes_dict.items():
            # make sure we have data
            rows.sort(key=lambda x: x[0])            # sort by concentration
            if not rows:
                continue

            concentrations, mean_pk_list, std_list = zip(*rows)
            plt.figure(figsize=(10, 6))

            for i, mean_pk in enumerate(zip(*mean_pk_list)):
                if selected_indices is not None and i not in selected_indices:
                    continue            # user did not select this peak

                if len(concentrations) < 2:
                    print(f"⚠  not enough points to fit {analyte} – peak {i+1}")
                    continue

                # optional error bars
                errors = [s[i] for s in std_list] if std_list else None
                if errors:
                    plt.errorbar(concentrations, mean_pk, yerr=errors,
                                fmt='o', label=f"{analyte} – Peak {i+1}")
                else:
                    plt.scatter(concentrations, mean_pk,
                                label=f"{analyte} – Peak {i+1}")

                # ----- fit log10(C) vs current ----------------------------
                log_c = np.log10(concentrations)
                slope, intercept, r2, *_ = linregress(log_c, mean_pk)
                calib_models[(analyte, i)] = (slope, intercept)   # ❹ store fit
                plt.plot(concentrations, slope*log_c + intercept,
                        label=(f"Fit Peak {i+1} "
                                f"(R²={r2**2:.2f}, m={slope:.2f}, b={intercept:.2f})"))

            plt.xscale('log')
            plt.xlabel('Concentration (µM)')
            plt.ylabel('Mean Peak Current (µA)')
            plt.title(f'Mean Peak Currents – {analyte}')
            plt.legend()
            plt.tight_layout()

            img_path = os.path.join(temp_dir, f"{base_save_path}_{analyte}.png")
            plt.savefig(img_path)
            plt.close()
            plot_filenames.append(img_path)

        # ---- 4) hand the fits to the rest of the app ---------------------
        self._cal_models = calib_models
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

    def plot_residual_vs_peak_location(self, selected_indices=None):
        """
        Bland–Altman scatter: residual (µA) vs peak potential (V).
        """
        import numpy as np, matplotlib.pyplot as plt
        from collections import defaultdict
        if not getattr(self, "_cal_models", None) or not self._cal_models:
            tk.messagebox.showwarning(
                "No calibration",
                "Run “Analyze Peak Currents” first (with at least two concentrations per peak)."
            )
            return
        # filter curves if user pre-selected some
        chosen_entries = None
        if selected_indices is not None:
            chosen_entries = {i // 2 for i in selected_indices}

        # group residuals per (Analyte, peak#) so we can colour-code
        data = defaultdict(list)   #  {(analyte, k): [(V,res), …]}

        for row in self.analysis_results['results']:
            if chosen_entries and row['Entry'] not in chosen_entries:
                continue

            analyte = row['Analytes']
            conc    = row['Concentration']
            peak_V  = row['Peak Voltage']
            Iobs    = row['Adjusted Peak Current']

            # Determine which peak number this is within its Entry
            k = sum(r['Entry']==row['Entry'] and r['Peak Voltage']<peak_V
                    for r in self.analysis_results['results']
                    if r['Analytes']==analyte)      # 0 for 1st peak, 1 for 2nd …

            model_key = (analyte, k)
            if model_key not in self._cal_models:
                # calibration not available (e.g. peak outside fit); skip
                continue

            slope, intercept = self._cal_models[model_key]

            # need the concentration in *numeric* form for prediction
            import re
            #c_val = float(re.match(r"([0-9]*\.?[0-9]+)", conc).group(1))
            match = re.match(r"([0-9]*\.?[0-9]+)", conc)
            if not match:
                continue                           # weird label → skip
            c_val = float(match.group(1))
            if c_val <= 0:
                continue                           # cannot take log10(0) → skip
            Ipred = slope * np.log10(c_val) + intercept
            resid = Iobs - Ipred
            data[model_key].append((peak_V, resid))
            Ipred = slope * np.log10(c_val) + intercept
            resid = Iobs - Ipred
            data[model_key].append((peak_V, resid))

        # ── draw ──────────────────────────────────────────────────────────────
        plt.figure(figsize=(7,5))
        for (analyte, k), pts in data.items():
            V, R = zip(*pts)
            lbl  = f"{analyte} – Peak {k+1}"
            plt.scatter(V, R, label=lbl, alpha=0.8)

        plt.axhline(0, color='k', lw=0.8, ls='--')
        plt.xlabel("Peak potential  /  V")
        plt.ylabel("Residual (Iₒᵦₛ – I_fitted)  /  µA")
        plt.title("Residual vs Peak Potential")
        plt.legend(bbox_to_anchor=(1.04,1), loc='upper left', fontsize='small')
        plt.tight_layout()
        plt.show()


#_____________________________
# Create Main Application
#_____________________________

from tkinter import font

# Main Application Window (For completeness)
class MainApp:
    def __init__(self, master):
        self.master = master
        master.title("Electrochemical Analysis Software Interface (EASI)")
        master.geometry("550x300")

        headline_font = tk.font.Font(family="Helvetica", size=18, weight="bold")

        self.header_frame = tk.Frame(master)
        self.header_frame.grid(row=0, column=0, columnspan=2, padx=20, pady=10, sticky="w")

        self.headline = Label(self.header_frame, text="EASI", font=headline_font)
        self.headline.pack(side="right")

#        logo_path = r"C:\Users\Jonat\Desktop\Electrochemical Analysis Software Interface(EASI)\LOGO.gif"
#        self.logo_image = tk.PhotoImage(file=logo_path)
#        self.logo_image = self.logo_image.subsample(15, 15)
#        self.logo_label = Label(self.header_frame, image=self.logo_image)
#        self.logo_label.pack(side="left", padx=10)

        self.content_frame = tk.Frame(master)
        self.content_frame.grid(row=1, column=0, columnspan=2, sticky="nw")

        self.description = Label(self.content_frame, text="Welcome to the Electrochemical Analysis Software Interface (EASI).\n"
                                                          "This tool provides various electrochemical analysis methods:\n"
                                                          "1. Cyclic Voltammetry(CV) Analysis\n"
                                                          "2. Differential Pulse Voltammetry(DPV) Analysis\n"
                                                          "3. Electrochemical Impedance Spectroscopy(EIS) Analysis\n"
                                                          "4. Customised Plotting Analysis\n"
                                                          "\n"
                                                          "Please select an option on the right to proceed.",
                                 justify="left", padx=10, pady=10)
        self.description.grid(row=0, column=0, sticky="nw")

        self.button_frame = tk.Frame(self.content_frame)
        self.button_frame.grid(row=0, column=1, padx=10, pady=10, sticky="ne")

        self.open_cv_button = Button(self.button_frame, text=" Open CV Analysis ", command=self.open_cv)
        self.open_cv_button.grid(row=0, column=0, pady=5)

        self.open_dpv_button = Button(self.button_frame, text="Open DPV Analysis", command=lambda: self.open_dpv(master))
        self.open_dpv_button.grid(row=1, column=0, pady=5)

        self.open_eis_button = Button(self.button_frame, text=" Open EIS Analysis ", command=self.open_eis)
        self.open_eis_button.grid(row=2, column=0, pady=5)

        self.open_plot_button = Button(self.button_frame, text="  Custom Plotting  ", command=self.open_plot)
        self.open_plot_button.grid(row=3, column=0, pady=5)

    def open_cv(self):
        new_window = Toplevel(self.master)
        CVApp(new_window)

    def open_dpv(self, master):
        new_window = Toplevel(master)
        DPVApp(new_window, master)

    def open_eis(self):
        new_window = Toplevel(self.master)
        EISApp(new_window)

    def open_plot(self):
        new_window = Toplevel(self.master)
        Custom_Plot_App(new_window)

if __name__ == "__main__":
    root = tk.Tk()
    main_app = MainApp(root)
    root.mainloop()
