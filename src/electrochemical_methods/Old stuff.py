def Analysis_DPV(file_path, blank_responses=None):
    #_________________________________________
    # DEFINE FUNCTIONS
    #_________________________________________

    def load_data_and_metadata(file_path):
        # Function to read lines from the file
        def read_file_lines(file_path, encoding='utf-16'):
            with open(file_path, 'r', encoding=encoding) as file:
                return file.readlines()

        def extract_metadata_row(lines, index=3):
            # Get the fourth row, remove leading 'Metadata row: ' if present
            metadata_row = lines[index].strip()
            metadata_row = re.sub(r'^Metadata row: ', '', metadata_row)
            metadata_entries = metadata_row.split(',,')

            # Clean each metadata entry by removing everything after 'DPV'
            cleaned_metadata_entries = []
            missing_parts = []
            for entry in metadata_entries:
                entry = entry.strip()  # Remove leading/trailing whitespace

                # If 'DPV' exists, remove everything after it
                if 'DPV' in entry:
                    entry = entry.split('DPV')[0] + 'DPV'
                else:
                    missing_parts.append("Method")
                cleaned_metadata_entries.append(entry)
            # Filter out any empty strings that might exist due to trailing delimiters
            cleaned_metadata_entries = [entry for entry in cleaned_metadata_entries if entry]
            return cleaned_metadata_entries, missing_parts
        
        def validate_metadata_entry(entry, entry_number, missing_parts):
            for entry_number, entry in enumerate(metadata_entries, start=1):
                entry_to_validate = re.sub(r'PP_', '', entry)
                parts_to_validate = entry_to_validate.split('_')

                # Check if parts_to_validate has at least 1 part to access index [0]
                if len(parts_to_validate) > 0:
                    if '-' in parts_to_validate[0]:
                        missing_parts.append("Electrode")

                # Check if parts_to_validate has at least 2 parts to access index [1]
                if len(parts_to_validate) > 1:
                    if '-' not in parts_to_validate[1]:
                        if any(char.isdigit() for char in parts_to_validate[1]):
                            missing_parts.append("Analytes")
                        else:
                            missing_parts.append("Concentration")
                
                if len(parts_to_validate) > 2:
                    if 'DPV' not in entry:
                        missing_parts.append("Method")

                # Raise an error if there are missing parts
                if missing_parts:
                    raise CustomValueError(f"Error in Entry #{entry_number}: Entry is missing required components: {', '.join(missing_parts)}. Expected format: 'Electrode_Analytes-Concentration_Method'.")

        def clean_metadata_entry(entry):
            entry = re.sub(r'PP_', '', entry)
            entry = entry.replace('-', '_')
            parts = entry.split('_')
            entry = '_'.join(parts[:4])
            return entry

        def parse_metadata_entry(entry):
            parts = entry.split('_')
            return {
                'Electrode': parts[0],
                'Analytes': parts[1],
                'Concentration': parts[2],
                'Method': parts[3]
            }

        def load_and_clean_data(file_path, encoding='utf-16', header_index=4, metadata_identifiers=['Date and time measurement:']):
            with open(file_path, 'r', encoding=encoding) as file:
                reader = csv.reader(file)
                data = list(reader)[header_index + 1:]

            headers = data[0]
            data_cleaned = [row for row in data if not any(identifier in row[0] for identifier in metadata_identifiers)]
            valid_data = [row for row in data_cleaned[1:] if len(row) == len(headers)]

            return headers, valid_data

        def convert_data_to_array(valid_data):
            def safe_float_convert(value):
                try:
                    return float(value)
                except ValueError:
                    return np.nan
            return np.array([[safe_float_convert(cell) for cell in row] for row in valid_data])

        # Read lines from the file
        lines = read_file_lines(file_path)

        # Extract and clean metadata
        metadata_entries, missing_parts = extract_metadata_row(lines)
        metadata_entries_val, missing_parts = extract_metadata_row(lines)
        
        # Validate Metadata Entries - Loop through all entries and stop if all are valid
        entry_number = []
        error_occurred = False

        try:
            for entry_number, entry in enumerate(metadata_entries_val, start=1):
                validate_metadata_entry(entry, entry_number, missing_parts)
        except CustomValueError as e:
            print(e)
            error_occurred = True

        if not error_occurred:
            # Clean metadata entries and parse metadata
            cleaned_metadata_entries = [clean_metadata_entry(entry) for entry in metadata_entries]
            parsed_metadata = [parse_metadata_entry(entry) for entry in cleaned_metadata_entries]

            # Load and clean the data
            headers, valid_data = load_and_clean_data(file_path)
            data_array = convert_data_to_array(valid_data)

            return headers, data_array, parsed_metadata
        else:
            print("Stopping due to validation error. No further processing will be performed, check your dataset.")
            return None, None, None
        
    # Main execution
    headers, data_array, parsed_metadata = load_data_and_metadata(file_path)

    # Check if any of the returned values are None
    if headers is None or data_array is None or parsed_metadata is None:
        results = mean_peak_currents = std_peak_currents = cov_peak_currents = lod_results = t_test_results = None
        return results
    else:
        results = find_peaks_and_process_data(headers, data_array, parsed_metadata)
        cov_peak_currents, mean_peak_currents, std_peak_currents = calculate_cov_peak_currents(results)
        lod_results = calculate_lod(mean_peak_currents, blank_responses)
        t_test_results = perform_t_tests(mean_peak_currents, std_peak_currents)
        return {
            'results': results,
            'mean_peak_currents': mean_peak_currents,
            'std_peak_currents': std_peak_currents,
            'cov_peak_currents': cov_peak_currents,
            'lod_results': lod_results,
            't_test_results': t_test_results,
            'headers': headers,
            'data_array': data_array,
            'parsed_metadata': parsed_metadata
        }
    

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


                def calculate_lod(mean_peak_currents, blank_responses=None):
        if blank_responses is None:
            blank_responses = [0.1, 0.15, 0.05, 0.1, 0.08]
        
        lod_results = {}
        std_blank = np.std(blank_responses)
        analytes_dict = defaultdict(list)
        for key, mean_peak_current in mean_peak_currents.items():
            concentration, analytes = key
            analytes_dict[analytes].append((float(concentration[:-2]), mean_peak_current))
        
        for analytes, data in analytes_dict.items():
            data.sort(key=lambda x: x[0])
            concentrations, mean_peak_currents = zip(*data)
            for i, mean_peak_current in enumerate(zip(*mean_peak_currents)):
                epsilon = 1e-10
                log_concentrations = np.log10(np.array(concentrations) + epsilon)
                slope, intercept, r_value, p_value, std_err = linregress(log_concentrations, mean_peak_current)
                lod = (3.3 * std_blank) / slope
                lod_results[(analytes, f'Peak {i+1}')] = lod
        
        return lod_results
    

        def import_file(self):
        self.filepath = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx;*.xls")]
        )
        if self.filepath:
            self.file_path_label.config(text=self.filepath)
            # Retrieve custom metadata fields entered by the user:
            custom_fields_str = self.metadata_fields_entry.get()  # e.g., "Sensor,Analyte,Dose,Technique"
            # Convert the comma-separated string into a list and strip any extra whitespace:
            custom_fields = [field.strip() for field in custom_fields_str.split(',')]
            # Pass the custom fields into Analysis_DPV:
            self.analysis_results = Analysis_DPV(self.filepath, metadata_fields=custom_fields)
            if not self.analysis_results:
                messagebox.showwarning("Load Failure", "Could not analyze the selected file.")
            else:
                # Optionally, enable the Save Results button now
                self.save_as_button.config(state=tk.NORMAL)
        else:
            self.file_path_label.config(text="No file selected")
            self.analysis_results = None  # Clear out old results if user canceled