
import os, sys
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
from datetime import datetime
from tkinter import ttk
import pandas as pd
from tkinter import Button, END, Entry, Label, LabelFrame, OptionMenu, StringVar, messagebox, Toplevel, filedialog
import hashlib
from PIL import Image
from matplotlib import pyplot as plt

from electrochemical_methods.Analysis_CV import Analysis_CV
from electrochemical_methods.Analysis_DPV import Analysis_DPV
from electrochemical_methods.Analysis_EIS import Analysis_EIS
from electrochemical_methods.Plotting import Plots_all_imported

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
# robust .pssession loader  (works for your 1. S0_Donor Q file)
# ---------------------------------------------------------
def load_pssession(self, path: str) -> tuple[pd.DataFrame, int]:
    """
    Return (dataframe, n_replicates).

    The dataframe columns are:
        Potential (V)           – replicate 1
        Current (µA)            – replicate 1
        Potential (V)_2         – replicate 2
        Current (µA)_2          – replicate 2
        …
    """
    import json

    # ---------- read & isolate JSON header ----------------
    with open(path, "r", encoding="utf-16") as f:
        txt = f.read()

    start = txt.find("{")
    depth = 0
    for i, ch in enumerate(txt[start:], start):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                end = i
                break
    sess = json.loads(txt[start : end + 1])

    # ---------- gather *all* DPV curves -------------------
    curves = []
    for m in sess["Measurements"]:
        if "METHOD_ID=DPV" not in m.get("Method", "").upper():
            continue
        curves.extend(m.get("Curves", []))

    if not curves:
        raise ValueError("No DPV curves in session")

    columns = {}
    for idx, curve in enumerate(curves, start=1):
        # Try RawDataArray first --------------------------
        raw = curve.get("RawDataArray") or {}
        X = raw.get("Potential") or raw.get("X")
        Y = raw.get("Current")   or raw.get("Y")

        # Fallback: XAxisDataArray / YAxisDataArray -------
        if X is None or Y is None:
            xb = curve.get("XAxisDataArray") or {}
            yb = curve.get("YAxisDataArray") or {}
            X = xb.get("DataValues")
            Y = yb.get("DataValues")
            if X and isinstance(X[0], dict):
                X = [float(d["V"]) for d in X]
                Y = [float(d["V"]) for d in Y]

        if X is None or Y is None or len(X) != len(Y):
            raise ValueError(f"Curve #{idx} has inconsistent arrays")

        # Column names: first pair has no suffix, the rest “_2”, “_3”…
        suf = "" if idx == 1 else f"_{idx}"
        columns[f"Potential (V){suf}"] = X
        columns[f"Current (µA){suf}"]  = Y

    df = pd.DataFrame(columns)
    print(f"→ Loaded {len(curves)} curve(s), {len(df)} points each")
    return df, len(curves)

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
    if ext == '.pssession':
        try:
            # --- load the DPV curve itself ---
            df, n_reps = self.load_pssession(self.filepath)

            # --- grab the Title from the JSON for metadata ---
            with open(self.filepath, 'r', encoding='utf-16') as f:
                txt = f.read()
            start = txt.find('{'); brace = 0; end = None
            for i,ch in enumerate(txt[start:], start):
                if   ch=='{': brace += 1
                elif ch=='}': brace -= 1
                if brace==0:
                    end = i
                    break
            sess    = json.loads(txt[start:end+1])
            dpv_meas= next(m for m in sess["Measurements"]
                            if "METHOD_ID=DPV" in m.get("Method","").upper())
            raw_title = dpv_meas.get("Title","").strip()
            print("DEBUG raw Title:", raw_title)

            # --- strip the leading "1. ", then cut off at "_DPV" ---
            body     = raw_title.split('. ',1)[-1]
            name_only= body.split('_DPV',1)[0]

            # --- split into exactly (N-1) pieces, N = len(custom_fields) ---
            n_meta   = len(custom_fields)
            parts    = name_only.split('_', n_meta-1)

            # --- pad to (N-1) tokens if needed ---
            while len(parts) < n_meta-1:
                parts.append("")

            # --- finally tack on "DPV" as the Method field ---
            parts.append("DPV")

            # --- join back to one Metadata row entry ---

            #meta_entry = "_".join(parts)
            meta_entry_single = "_".join(parts)
            meta_entry = ",,".join([meta_entry_single] * n_reps)
            print("DEBUG meta_entry:", meta_entry)

            
            with tempfile.NamedTemporaryFile(
                    suffix='.csv', delete=False,
                    mode='w', encoding='utf-16', newline=''
            ) as tmp:
                tmp.write("\n"*3)                   # 0-3: blank lines
                tmp.write(f"Metadata row: {meta_entry}\n")   # 4: metadata
                tmp.write("Date and time measurement: \n") 
                df.to_csv(tmp, index=False, header=True)      # 5-…: header + data
                tmp_path = tmp.name

            # (optional) quick preview
            print("=== Temp CSV Preview ===")
            with open(tmp_path, 'r', encoding='utf-16') as dbg:
                for _ in range(5):
                    line = dbg.readline()
                    if not line: break
                    print(line.strip())
            print("=== End Preview ===")

            file_to_pass = tmp_path

        except Exception as e:
            messagebox.showerror(
                "Load Error",
                f"Couldn’t parse session file:\n{e}"
            )
            self._debug_dump(file_to_pass)

            self.analysis_results = None
            return

    else:
        file_to_pass = self.filepath

    # -----------------------------------------------
    # 2) build peak_regions from the text-field
    # -----------------------------------------------
    self.peak_regions = self._parse_peak_region_string(
        self.peak_entry.get()
    )

    # 2) hand it off to your unchanged DPV loader
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
        "Observed vs Expected Concentration"
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

    elif analysis_option == "Generate Peak Report":
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

def display_results(self, results):
    """Show all numeric summaries (means, stdev, COV, FWHM …) in a pop-up."""
    # ── basic Tk window set-up ───────────────────────────────────────────────
    result_window = Toplevel(self.master)
    result_window.title("DPV Analysis Results")

    text_area = tk.Text(result_window, wrap=tk.WORD, width=80, height=20)
    text_area.pack(padx=10, pady=10)

    # ── build text report ───────────────────────────────────────────────────
    import numpy as np
    from collections import defaultdict

    result_text = "DPV Analysis Results\n\n"

    # -------- mean peak currents -------------------------------------------
    result_text += "Mean Peak Currents:\n"
    for key, value in results['mean_peak_currents'].items():
        result_text += f"{key}: {value}\n"

    # -------- std dev -------------------------------------------------------
    result_text += "\nStandard Deviation of Peak Currents:\n"
    for key, value in results['std_peak_currents'].items():
        result_text += f"{key}: {value}\n"

    # -------- coefficient of variation -------------------------------------
    result_text += "\nCoefficient of Variation of Peak Currents:\n"
    for key, value in results['cov_peak_currents'].items():
        result_text += f"{key}: {value}\n"

    # ────────────────────────────────────────────────────────────────────────
    # NEW: F W H M  per analyte  + overall stats
    # ────────────────────────────────────────────────────────────────────────
    per_analyte = defaultdict(list)        # {'HX': [ … ], 'UA': [ … ], …}

    for row in results['results']:         # one dict per detected peak
        w = row.get('FWHM')
        if w is not None and not np.isnan(w):
            per_analyte[row['Analytes']].append(w)

    if per_analyte:                         # only print if we have at least one
        result_text += "\nFWHM per analyte (mean ± SD):\n"
        all_vals = []

        for analyte, lst in per_analyte.items():
            all_vals.extend(lst)
            μ = np.mean(lst)
            σ = np.std(lst, ddof=0)
            n = len(lst)
            nums = ", ".join(f"{v:.3f}" for v in lst)
            result_text += (
                f"  {analyte}: {μ:.4g} ± {σ:.4g} V  (n={n})   "
                f"values=[{nums}]\n"
            )

        # result_text += (
        #     "\nFWHM summary (all peaks):\n"
        #     f"  mean   = {np.mean(all_vals):.4g} V\n"
        #     f"  median = {np.median(all_vals):.4g} V\n"
        #     f"  min / max = {min(all_vals):.4g} / {max(all_vals):.4g} V\n"
        # )

    # ────────────────────────────────────────────────────────────────────────

    text_area.insert(tk.END, result_text)


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


