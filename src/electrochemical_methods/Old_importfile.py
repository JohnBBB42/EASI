    def load_pssession(self, path):
        # 1) Read file as UTF‑16 so the BOM is stripped automatically
        with open(path, 'r', encoding='utf-16') as f:
            txt = f.read()

        # 2) Extract only the JSON header by counting braces
        start = txt.find('{')
        if start < 0:
            raise ValueError("No JSON object in session file")
        brace = 0
        for i, ch in enumerate(txt[start:], start):
            if ch == '{':
                brace += 1
            elif ch == '}':
                brace -= 1
                if brace == 0:
                    end = i
                    break
        else:
            raise ValueError("Could not find end of JSON object")
        json_text = txt[start:end+1]

        # 3) Parse JSON
        sess = json.loads(json_text)

        # 4) DEBUG: see your measurements
        print("SESSION KEYS:", sess.keys())
        for mi, meas in enumerate(sess["Measurements"]):
            print(f"\nMeasurement #{mi}:")
            # This long 'Method' contains your METHOD_ID=dpv line
            print(" Method:", meas.get("Method", "")[:60].replace('\r\n','\\n') + "…")
            curves = meas.get("Curves", [])
            print("  → # of curves:", len(curves))

        # 5) Pick the DPV measurement by looking for 'METHOD_ID=dpv'
        dpv_meas = None
        for meas in sess["Measurements"]:
            method_lower = meas.get("Method", "").lower()
            if "method_id=dpv" in method_lower:
                dpv_meas = meas
                break

        if dpv_meas is None:
            raise ValueError("No measurement with METHOD_ID=dpv found")

        # 6) Pull out its curve list (your file has exactly one)
        curves = dpv_meas.get("Curves", [])
        if not curves:
            raise ValueError("DPV measurement has no 'Curves' block")

        # 7) Take the first (and only) curve object
        curve = curves[0]

        # 8) Read the X/Y arrays
        x = curve.get("XAxisDataArray")
        y = curve.get("YAxisDataArray")
        if x is None or y is None:
            raise ValueError("Curve missing XAxisDataArray or YAxisDataArray")

        # 9) Build a DataFrame and return
        df = pd.DataFrame({"V": x, "C": y})
        print(f"→ Loaded {len(df)} DPV points from session")
        return df

def import_file(self):
        self.filepath = filedialog.askopenfilename(
            filetypes=[
                ("All files",       "*.*"),
                ("PalmSens session","*.pssession"),
                ("CSV files",       "*.csv"),
                ("Excel files",     "*.xlsx;*.xls"),
            ]
        )
        if not self.filepath:
            self.file_path_label.config(text="No file selected")
            self.analysis_results = None
            return

        self.file_path_label.config(text=self.filepath)

        # parse user’s metadata‐fields string
        custom_fields = [
            fld.strip()
            for fld in self.metadata_fields_entry.get().split(',')
        ]

        # if it’s a .pssession, turn it into a temp CSV first
        ext = os.path.splitext(self.filepath)[1].lower()
        if ext == '.pssession':
            try:
                df = self.load_pssession(self.filepath)
            except Exception as e:
                messagebox.showerror("Load Error",
                                     f"Couldn’t parse session file:\n{e}")
                return

            # dump to a temp CSV encoded as UTF-16 (with BOM)
            tmp = tempfile.NamedTemporaryFile(
                suffix='.csv',
                delete=False,
                mode='w',
                newline='',
                encoding='utf-16'        # <<< add this
            )
            # to_csv also needs the same encoding
            df.to_csv(tmp.name, index=False, encoding='utf-16')
            file_to_pass = tmp.name
        else:
            file_to_pass = self.filepath

        # now call your existing DPV loader
        self.analysis_results = Analysis_DPV(
            file_to_pass,
            metadata_fields=custom_fields
        )

        if not self.analysis_results:
            messagebox.showwarning("Load Failure",
                                   "Could not analyze the selected file.")
        else:
            self.save_as_button.config(state=tk.NORMAL)


            
        # --- PALMSENS CSV GUARD: no "Metadata row:" → skip metadata parsing ---
        if metadata_header_index >= len(lines) \
           or not lines[metadata_header_index].startswith("Metadata row:"):

            # 1) Load the numeric table exactly as you do in load_and_clean_data:
            with open(file_path, 'r', encoding=encoding) as f:
                reader = csv.reader(f)
                all_rows = list(reader)

            headers   = all_rows[0]
            data_rows = all_rows[1:]
            data_array = convert_data_to_array(data_rows)

            # 2) Infer metadata: fill only the Method field with "DPV"
            n_curves = data_array.shape[1] // 2
            parsed_metadata = []
            for _ in range(n_curves):
                md = {field: "" for field in metadata_fields}
                md["Method"] = "DPV"
                parsed_metadata.append(md)

            return headers, data_array, parsed_metadata