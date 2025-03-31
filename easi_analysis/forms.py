from django import forms

class CVAnalysisForm(forms.Form):
    file = forms.FileField(label="Upload your CV data (CSV or Excel)")


class DPVAnalysisForm(forms.Form):
    file = forms.FileField(required=False, label="Upload your DPV CSV or Excel file:")
    blank_file = forms.FileField(required=False, label="Upload Blank Responses File (Optional):")
    selected_analysis = forms.ChoiceField(
        choices=[
            ("DPV analysis", "DPV analysis"),
            ("LOD Analysis", "LOD Analysis"),
            ("T-test Analysis", "T-test Analysis"),
            ("Plot Data Array with Corrected Baseline", "Plot Data Array with Corrected Baseline"),
            ("Analyze Peak Currents", "Analyze Peak Currents"),
            ("Observed vs Expected Concentration", "Observed vs Expected Concentration")
        ],
        label="Select Analysis Type:"
    )

class EISAnalysisForm(forms.Form):
    file = forms.FileField(label="Upload your EIS data (CSV or Excel)")
