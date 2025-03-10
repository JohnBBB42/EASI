"""
Created on Fri Aug 21 14:11:08 2024

@author: jonatanemilsvendsen
"""
#_____________________________
# Presetting
#_____________________________
#Change directory to location of Python script
import os, sys
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
sys.path.append(dname)

# Import packages
import numpy as np
import tkinter as tk
from tkinter import Button, END, Entry, Label, LabelFrame, OptionMenu, StringVar, messagebox, Toplevel, filedialog
import hashlib
from PIL import Image
from matplotlib import pyplot as plt

from electrochemical_methods.Analysis_CV import Analysis_CV
from electrochemical_methods.Analysis_DPV import Analysis_DPV
from electrochemical_methods.Analysis_EIS import Analysis_EIS
from electrochemical_methods.Plotting import Plots_all_imported


#_____________________________
# GUI setup
#_____________________________
class CVApp:
    def __init__(self, master):
        self.master = master
        master.title('CV Analysis')
        master.geometry("400x600")

        self.filepath = None
        self.analysis_results = None

        # Dropdown options
        self.options = ["CV Analysis", "Custom CV Plotting"]
        self.clicked = StringVar(value=self.options[0])

        # Dropdown Menu
        self.frame = LabelFrame(master, text="Select Analysis")
        self.frame.pack(padx=10, pady=10, fill="x")
        OptionMenu(self.frame, self.clicked, *self.options).pack(fill="x")

        # File import
        Button(master, text="Import File", command=self.import_file).pack(pady=10)
        self.file_label = Label(master, text="No file selected")
        self.file_label.pack(pady=5)

        # Run analysis
        Button(master, text="Run Analysis", command=self.run_analysis).pack(pady=10)

        # Save results
        Button(master, text="Save Results", command=self.save_results).pack(pady=10)

    def import_file(self):
        self.filepath = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx;*.xls")]
        )
        if self.filepath:
            self.file_label.config(text=self.filepath)

    def run_analysis(self):
        if not self.filepath:
            messagebox.showwarning("No File", "Please select a file first.")
            return

        df = pd.read_csv(self.filepath) if self.filepath.endswith('.csv') else pd.read_excel(self.filepath)

        if self.clicked.get() == "CV Analysis":
            try:
                saving_folder = filedialog.askdirectory(title="Select folder to save results")
                if not saving_folder:
                    messagebox.showwarning("No Folder", "Please select a folder to save results.")
                    return

                self.analysis_results = Analysis_CV(
                    df,
                    values_row_start=2,
                    potential_column=1,
                    current_column=2,
                    scan_column=0,
                    scan_number=1,
                    linreg_start_index=15,
                    r2_accept_value=0.90,
                    potential_unit="V",
                    current_unit="A",
                    num_decimals=3,
                    saving_folder=saving_folder
                )
                messagebox.showinfo("Success", "CV analysis completed successfully!")
            except Exception as e:
                messagebox.showerror("Analysis Error", f"An error occurred: {e}")

        elif self.clicked.get() == "Custom CV Plotting":
            messagebox.showinfo("Info", "Custom plotting feature not yet implemented.")

    def save_results(self):
        if self.analysis_results:
            save_dir = filedialog.askdirectory()
            if save_dir:
                plot_path = self.analysis_results['plot_path']
                new_plot_path = os.path.join(save_dir, os.path.basename(plot_path))
                os.rename(plot_path, new_plot_path)
                messagebox.showinfo("Saved", f"Results saved to {new_plot_path}")
        else:
            messagebox.showwarning("No Results", "Run analysis first.")



class EISApp:
    def __init__(self,master):
        self.master = master
        master.title('EIS analysis')
        master.geometry("500x400")
        
        # Drop down boxes
        self.options = ["EIS-analysis"]
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
        self.variable2 = self.entry_box('Z_R_column',3,1,3)
        self.variable3 = self.entry_box('Z_Im_column',4,1,4)
        self.variable4 = self.entry_box('Z_R_start',5,1,"")
        self.variable5 = self.entry_box('Z_R_end',6,1,"")
        self.variable6 = self.entry_box('Z_Im_start',7,1,"")
        self.variable7 = self.entry_box('Z_Im_end',8,1,"")
        self.variable8 = self.entry_box('impedance_unit',9,1,"Ω")
        self.variable9 = self.entry_box('circle_point1_index',10,1,0)
        self.variable10 = self.entry_box('circle_point2_index',11,1,0)
        
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
            Text = "Check your excel/csv file and list the column number of where the Real part of impedance values appear."
        tk.messagebox.showinfo("Info",str(Text))
    
    def info_box3(self):
        if self.clicked.get() == self.options[0]:
            Text = "Check your excel/csv file and list the column number of where the Imaginary part of impedance values appear."
        tk.messagebox.showinfo("Info",str(Text))
    
    def info_box4(self):
        if self.clicked.get() == self.options[0]:
            Text = "Insert value for Start interval of Real part of impedance values. NB: both Z_R_start and Z_R_end must be inserted for the interval range to be created."
        tk.messagebox.showinfo("Info",str(Text))
    
    def info_box5(self):
        if self.clicked.get() == self.options[0]:
            Text = "Insert value for End interval of Real part of impedance values. NB: both Z_R_start and Z_R_end must be inserted for the interval range to be created."
        tk.messagebox.showinfo("Info",str(Text))
    
    def info_box6(self):
        if self.clicked.get() == self.options[0]:
            Text = "Insert value for Start interval of Imaginary part of impedance values. NB: both Z_Im_start and Z_Im_end must be inserted for the interval range to be created."
        tk.messagebox.showinfo("Info",str(Text))
    
    def info_box7(self):
        if self.clicked.get() == self.options[0]:
            Text = "Insert value for End interval of Imaginary part of impedance values. NB: both Z_Im_start and Z_Im_end must be inserted for the interval range to be created."
        tk.messagebox.showinfo("Info",str(Text))        

    def info_box8(self):
        if self.clicked.get() == self.options[0]:
            Text = "Type the Impedance unit."
        tk.messagebox.showinfo("Info",str(Text))  
    
    def info_box9(self):
        if self.clicked.get() == self.options[0]:
            Text = "Type the starting index number to be used in the semicircle fit."
        tk.messagebox.showinfo("Info",str(Text))  
    
    def info_box10(self):
        if self.clicked.get() == self.options[0]:
            Text = "Type the number of indexes to the left of the semicircle end that defines the final interval of the semicircle fit."
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
            self.variable2 = self.entry_box('Z_R_column',3,1,3)
            self.variable3 = self.entry_box('Z_Im_column',4,1,4)
            self.variable4 = self.entry_box('Z_R_start',5,1,"")
            self.variable5 = self.entry_box('Z_R_end',6,1,"")
            self.variable6 = self.entry_box('Z_Im_start',7,1,"")
            self.variable7 = self.entry_box('Z_Im_end',8,1,"")
            self.variable8 = self.entry_box('impedance_unit',9,1,"Ω")
            self.variable9 = self.entry_box('circle_point1_index',10,1,0)
            self.variable10 = self.entry_box('circle_point2_index',11,1,0)

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
            Analysis_EIS(values[0], values[1],values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9])          
        else:
            pass


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

import os
import csv
import tempfile
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from collections import defaultdict
from tkinter import filedialog, messagebox, Toplevel, StringVar, Label, LabelFrame, Button, OptionMenu
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.stats import linregress



class DPVApp:
    def __init__(self, master, main_window):
        self.master = master
        self.main_window = main_window
        master.title('DPV analysis')
        # Positioning the DPV window next to the main window
        master.geometry("350x600+{}+{}".format(main_window.winfo_x() + main_window.winfo_width(), main_window.winfo_y()))

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
            "T-test Analysis"
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

    def import_file(self):
        self.filepath = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx;*.xls")]
        )
        if self.filepath:
            self.file_path_label.config(text=self.filepath)

            # Automatically run Analysis_DPV and store its result
            self.analysis_results = Analysis_DPV(self.filepath)
            if not self.analysis_results:
                messagebox.showwarning("Load Failure", "Could not analyze the selected file.")
            else:
                # Optionally, enable the Save Results button now
                self.save_as_button.config(state=tk.NORMAL)

        else:
            self.file_path_label.config(text="No file selected")
            self.analysis_results = None  # Clear out old results if user canceled


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
        # Display the results in a new window or a suitable UI component
        result_window = Toplevel(self.master)
        result_window.title("DPV Analysis Results")
        
        text_area = tk.Text(result_window, wrap=tk.WORD, width=80, height=20)
        text_area.pack(padx=10, pady=10)

        result_text = "DPV Analysis Results\n\n"

        result_text += "Mean Peak Currents:\n"
        for key, value in results['mean_peak_currents'].items():
            result_text += f"{key}: {value}\n"
        
        result_text += "\nStandard Deviation of Peak Currents:\n"
        for key, value in results['std_peak_currents'].items():
            result_text += f"{key}: {value}\n"
        
        result_text += "\nCoefficient of Variation of Peak Currents:\n"
        for key, value in results['cov_peak_currents'].items():
            result_text += f"{key}: {value}\n"
        
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



    def plot_selected(self, analysis_option):
        selected_indices = [i * 2 for i, var in enumerate(self.plot_vars) if var.get()]
        print(f"Selected plots for {analysis_option}: {selected_indices}")

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

    def plot_data_array(self, data_array, headers, parsed_metadata, save_path, selected_indices=None):
        plt.figure(figsize=(10, 6))

        if selected_indices is None:
            selected_indices = range(0, len(headers), 2)

        for i in selected_indices:
            potential = data_array[:, i]
            current = data_array[:, i + 1]
            
            meta = parsed_metadata[i//2]  # if each pair of columns matches one metadata entry
            analyte = meta['Analytes']
            concentration = meta['Concentration']
            label = f"{analyte} - {concentration}"
            plt.plot(potential, current, label=label)

        plt.xlabel('Potential (V)')
        plt.ylabel('Current (µA)')
        plt.title('Overlapping Plots of Potential vs Current')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(save_path)
        plt.close()

    def plot_data_array_with_corrected_baseline(self, data_array, headers, parsed_metadata, save_path, selected_indices=None):
        plt.figure(figsize=(10, 6))

        if selected_indices is None:
            selected_indices = range(0, len(headers), 2)

        def baseline_als(y, lam=1e5, p=0.01, niter=10):
            L = len(y)
            D = sparse.csc_matrix(np.diff(np.eye(L), 2))
            w = np.ones(L)
            for i in range(niter):
                W = sparse.spdiags(w, 0, L, L)
                Z = W + lam * D.dot(D.transpose())
                z = spsolve(Z, w * y)
                w = p * (y > z) + (1 - p) * (y < z)
            return z

        for i in selected_indices:
            potential = data_array[:, i]
            current = data_array[:, i + 1]

            # Filter out NaN values
            valid_idx = ~np.isnan(potential) & ~np.isnan(current)
            potential = potential[valid_idx]
            current = current[valid_idx]

            if len(potential) == 0 or len(current) == 0:
                print(f"Skipping plot for index {i//2} due to insufficient data")
                continue

            # Apply baseline correction
            baseline = baseline_als(current)
            corrected_current = current - baseline

            meta = parsed_metadata[i//2]  # if each pair of columns matches one metadata entry
            analyte = meta['Analytes']
            concentration = meta['Concentration']
            label = f"{analyte} - {concentration}"

            plt.plot(potential, corrected_current, label=label)

        plt.xlabel('Potential (V)')
        plt.ylabel('Current (µA)')
        plt.title('Corrected Plots of Potential vs Current')
        plt.legend(
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        fontsize='small',        # or 'small', 'medium', etc.
        ncol=2                      
        ) 
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

    def analyze_peak_currents(self, mean_peak_currents, std_peak_currents, base_save_path, selected_indices=None):
        # Ensure std_peak_currents is a dictionary
        if not isinstance(std_peak_currents, dict):
            raise TypeError(f"Expected std_peak_currents to be a dictionary, but got {type(std_peak_currents)} instead.")

        temp_dir = base_save_path
        plot_filenames = []

        analytes_dict = defaultdict(list)
        for key, mean_peak_current in mean_peak_currents.items():
            concentration, analytes = key  # Key is a tuple (Concentration, Analytes)

            # Ensure std_peak_currents is accessed with the correct key
            if key not in std_peak_currents:
                std_error = [0] * len(mean_peak_current)  # Default to zeros if key is missing
            else:
                std_error = std_peak_currents[key]

            analytes_dict[analytes].append((float(concentration[:-2]), mean_peak_current, std_error))

        for analytes, data in analytes_dict.items():
            data.sort(key=lambda x: x[0])  # Sort by concentration
            concentrations, mean_peak_currents, std_errors = zip(*data)
            plt.figure(figsize=(10, 6))

            for i, mean_peak_current in enumerate(zip(*mean_peak_currents)):
                if selected_indices is None or i in selected_indices:
                    error = [std[i] for std in std_errors] if std_errors else None
                    # Scatter or errorbar
                    if error:
                        plt.errorbar(
                            concentrations,
                            mean_peak_current,
                            yerr=error,
                            fmt='o',
                            label=f"{analytes} - Peak {i+1}"
                        )
                    else:
                        plt.scatter(concentrations, mean_peak_col, label=f"{analytes} - Peak {i+1}")
                    
                    log_concentrations = np.log10(concentrations)
                    slope, intercept, r_value, p_value, std_err = linregress(log_concentrations, mean_peak_current)
                    fit_line = slope * log_concentrations + intercept
                    plt.plot(concentrations, fit_line, label=f'Fit Peak {i+1} (R²={r_value ** 2:.2f}, Slope={slope:.2f}, Intercept={intercept:.2f})')

            plt.xscale('log')
            plt.xlabel('Concentration (uM)')
            plt.ylabel('Mean Peak Current')
            plt.title(f'Mean Peak Currents and Linear Fit for Analytes: {analytes}')
            plt.legend()
            plt.tight_layout()

            save_path = os.path.join(temp_dir, f"{base_save_path}_{analytes}.png")
            plt.savefig(save_path)
            plt.close()
            plot_filenames.append(save_path)

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
