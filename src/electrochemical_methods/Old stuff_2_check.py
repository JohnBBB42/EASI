def Analysis_DPV(file_path, blank_responses=None, metadata_fields=None):
    import os, re, csv
    import numpy as np
    from collections import defaultdict
    from scipy import sparse
    from scipy.sparse.linalg import spsolve
    from scipy.signal import find_peaks, savgol_filter
    from scipy.stats import linregress, t



    # ---------------------------
    # Flexible metadata loading function
    # ---------------------------
    def load_data_and_metadata(file_path, 
                               metadata_fields=metadata_fields, 
                               encoding='utf-16', 
                               metadata_header_index=3, 
                               metadata_delimiter=",,"):
        """
        Load the file, extract and validate metadata using a flexible naming scheme,
        and load the accompanying numerical data.
        
        Parameters:
            file_path (str): Path to the file.
            metadata_fields (list, optional): Expected metadata field names.
                Default is ["Electrode", "Analytes", "Concentration", "Method"].
            encoding (str): File encoding. Default 'utf-16'.
            metadata_header_index (int): Row index for the metadata header (default is 3, i.e. the fourth line).
            metadata_delimiter (str): Delimiter between metadata entries (default is ",,").
        
        Returns:
            tuple: (headers, data_array, parsed_metadata) on success; (None, None, None) if errors occur.
        """

        if metadata_fields is None:
            metadata_fields = ["Electrode", "Analytes", "Concentration", "Method"]

        # Read all lines from the file.
        def read_file_lines(fp, encoding=encoding):
            with open(fp, 'r', encoding=encoding) as file:
                return file.readlines()
        
        # Extract the metadata row and split its entries.
        def extract_metadata_row(lines, index=metadata_header_index, delimiter=metadata_delimiter):
            metadata_row = lines[index].strip()
            metadata_row = re.sub(r'^Metadata row: ', '', metadata_row)
            metadata_entries = [entry.strip() for entry in metadata_row.split(delimiter) if entry.strip()]
            # For each entry, if "DPV" is present, keep content up to (and including) "DPV"
            cleaned_entries = []
            missing_parts = []
            for entry in metadata_entries:
                entry = entry.strip()
                if 'DPV' in entry:
                    entry = entry.split('DPV')[0] + 'DPV'
                else:
                    missing_parts.append("Method")
                cleaned_entries.append(entry)
            return cleaned_entries, missing_parts

        # New separate cleaning function (exactly as in your original preprocessing)
        def clean_metadata_entry(entry, num_fields):
            entry = re.sub(r'^PP_', '', entry)
            entry = entry.replace('-', '_')
            parts = entry.split('_')
            # Only keep the first num_fields parts
            entry = '_'.join(parts[:num_fields])
            return entry

        # Validate and parse each metadata entry into a dictionary.
        def validate_and_parse_metadata_entries(entries, fields):
            parsed_entries = []
            errors = []
            for i, entry in enumerate(entries, start=1):
                cleaned = clean_metadata_entry(entry, len(fields))
                parts = cleaned.split('_')
                if len(parts) < len(fields):
                    missing = fields[len(parts):]
                    errors.append(f"Error in Entry #{i}: Missing field(s): {', '.join(missing)}. Expected format: {'_'.join(fields)}.")
                    continue
                # Use exactly the first len(fields) parts
                parts = parts[:len(fields)]
                entry_dict = {field: part for field, part in zip(fields, parts)}
                parsed_entries.append(entry_dict)
            return parsed_entries, errors

        # Load and clean the numerical data starting after the metadata.
        def load_and_clean_data(fp, encoding=encoding, values_header_index=metadata_header_index+1,
                                metadata_identifiers=['Date and time measurement:']):
            with open(fp, 'r', encoding=encoding) as file:
                reader = csv.reader(file)
                data = list(reader)[values_header_index:]
            headers = data[0]
            # Exclude rows containing metadata identifiers in the first column.
            clean_data = [row for row in data[1:] if not any(identifier in row[0] for identifier in metadata_identifiers)]
            # Only keep rows that have the same length as headers.
            valid_data = [row for row in clean_data if len(row) == len(headers)]
            return headers, valid_data

        def safe_float_convert(value):
            try:
                return float(value)
            except ValueError:
                return np.nan

        def convert_data_to_array(data_rows):
            return np.array([[safe_float_convert(cell) for cell in row] for row in data_rows])
        
        # --- Execute the metadata and data loading ---
        lines = read_file_lines(file_path)
        raw_metadata_entries, missing_parts = extract_metadata_row(lines)
        parsed_metadata, errors = validate_and_parse_metadata_entries(raw_metadata_entries, metadata_fields)
        if errors:
            for error in errors:
                print(error)
            print("Validation errors occurred. Stopping further processing.")
            return None, None, None

        headers, valid_data = load_and_clean_data(file_path)
        data_array = convert_data_to_array(valid_data)
        return headers, data_array, parsed_metadata