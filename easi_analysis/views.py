from django.shortcuts import render
from .forms import CVAnalysisForm, DPVAnalysisForm, EISAnalysisForm
from .utils.analysis_cv import Analysis_CV
from .utils.analysis_eis import Analysis_EIS
from .utils.analysis_dpv import Analysis_DPV
from .utils.plotting_dpv import DPVPlotting
from .utils.basics import Load_data
import pandas as pd
import tempfile
import os
import shutil
import pandas as pd
import numpy as np



# Create your views here.
def home(request):
    return render(request, 'easi_analysis/home.html')

def cv_analysis_view(request):
    result_message = None
    if request.method == 'POST':
        form = CVAnalysisForm(request.POST, request.FILES)
        if form.is_valid():
            uploaded_file = request.FILES['file']
            ext = uploaded_file.name.split('.')[-1].lower()

            try:
                if ext == 'csv':
                    df = pd.read_csv(uploaded_file, header=None)
                elif ext in ['xls', 'xlsx']:
                    df = pd.read_excel(uploaded_file, header=None)
                else:
                    result_message = "Unsupported file format."
                    df = None

                if df is not None:
                    temp_dir = tempfile.mkdtemp()
                    Analysis_CV(
                        df=df,
                        values_row_start=2,
                        potential_column=1,
                        current_column=3,
                        scan_column=5,
                        scan_number=1,
                        linreg_start_index=15,
                        r2_accept_value=0.90,
                        potential_unit="V",
                        current_unit="A",
                        num_decimals=3,
                        saving_folder=temp_dir
                    )
                    result_message = "✅ CV Analysis completed and saved to server (temp folder)."
            except Exception as e:
                result_message = f"❌ Error: {e}"
    else:
        form = CVAnalysisForm()

    return render(request, 'easi_analysis/cv.html', {'form': form, 'result': result_message})

def dpv_analysis_view(request):
    result_message = None
    results = None
    plot_filenames = []
    plot_choices = [] 
    analysis_results = None

    # If a file has been uploaded previously, retrieve its path from session.
    file_path = request.session.get('dpv_file_path')
    blank_path = request.session.get('dpv_blank_path')  # optional blank file path
    file_name = request.session.get('dpv_file_name')

    # If uploaded previously, run initial analysis and populate plot choices
    if file_path and os.path.exists(file_path):
        try:
            blank_array = np.loadtxt(blank_path, delimiter=',') if blank_path else None
            analysis_results = Analysis_DPV(file_path, blank_responses=blank_array)
            headers = analysis_results.get('headers', [])
            plot_choices = [(str(i), f"Plot {i // 2 + 1}") for i in range(0, len(headers), 2)]
        except Exception as e:
            result_message = f"⚠️ Initial analysis error: {e}"

    # Dynamic form with plot choices
    class DynamicDPVForm(DPVAnalysisForm): pass
    DynamicDPVForm.base_fields['selected_plots'].choices = plot_choices
    form = DynamicDPVForm(request.POST or None, request.FILES or None)

    if request.method == 'POST' and form.is_valid():
        selected_analysis = form.cleaned_data['selected_analysis']
        plot_dir = os.path.join('easi_analysis', 'static', 'easi_analysis', 'plots')
        os.makedirs(plot_dir, exist_ok=True)
        [os.remove(os.path.join(plot_dir, f)) for f in os.listdir(plot_dir)]

        file = request.FILES.get('file')
        blank_file = request.FILES.get('blank_file')

        if file:
            file_path = os.path.join(plot_dir, file.name)
            with open(file_path, 'wb+') as dest:
                for chunk in file.chunks():
                    dest.write(chunk)
            file_name = file.name
            request.session['dpv_file_path'] = file_path
            request.session['dpv_file_name'] = file.name

        if blank_file:
            blank_path = os.path.join(plot_dir, blank_file.name)
            with open(blank_path, 'wb+') as dest:
                for chunk in blank_file.chunks():
                    dest.write(chunk)
            request.session['dpv_blank_path'] = blank_path

        try:
            blank_array = np.loadtxt(blank_path, delimiter=',') if blank_path else None
            analysis_results = Analysis_DPV(file_path, blank_responses=blank_array)
            plotter = DPVPlotting()
            selected_indices = [int(i) for i in form.cleaned_data.get('selected_plots', []) or []]

            # Process analysis depending on the selected method
            if selected_analysis == "DPV analysis":
                results = {
                    'mean_peak_currents': analysis_results.get('mean_peak_currents', {}),
                    'std_peak_currents': analysis_results.get('std_peak_currents', {}),
                    'cov_peak_currents': analysis_results.get('cov_peak_currents', {})
                }
                result_message = "✅ DPV Analysis completed!"

            elif selected_analysis == "LOD Analysis":
                results = {
                    'lod_results': analysis_results.get('lod_results', {})
                }
                result_message = "✅ LOD Analysis completed!"

            elif selected_analysis == "T-test Analysis":
                results = {
                    't_test_results': analysis_results.get('t_test_results', [])
                }
                result_message = "✅ T-test Analysis completed!"
            elif selected_analysis == "Plot Data Array with Corrected Baseline":
                plot_path = os.path.join(plot_dir, "plot_corrected.png")
                plotter.plot_data_array_with_corrected_baseline(
                    analysis_results['data_array'],
                    analysis_results['headers'],
                    analysis_results['parsed_metadata'],
                    plot_path,
                    selected_indices=selected_indices
                )
                plot_filenames.append(f"easi_analysis/plots/{os.path.basename(plot_path)}")
                result_message = "✅ Plot created!"

            elif selected_analysis == "Analyze Peak Currents":
                peak_paths = plotter.analyze_peak_currents(
                    analysis_results['mean_peak_currents'],
                    analysis_results['std_peak_currents'],
                    plot_dir,
                    selected_indices=selected_indices
                )
                plot_filenames.extend([f"easi_analysis/plots/{os.path.basename(p)}" for p in peak_paths])
                result_message = "✅ Peak plots created!"

            elif selected_analysis == "Observed vs Expected Concentration":
                obs_paths = plotter.plot_observed_vs_expected_concentration(
                    analysis_results['mean_peak_currents'],
                    plot_dir,
                    selected_indices=selected_indices
                )
                plot_filenames.extend([f"easi_analysis/plots/{os.path.basename(p)}" for p in obs_paths])
                result_message = "✅ Observed vs Expected plots created!"

        except Exception as e:
            result_message = f"❌ Analysis failed: {e}"

    return render(request, 'easi_analysis/dpv.html', {
        'form': form,
        'result': result_message,
        'metrics': results,
        'plot_files': plot_filenames,
        'file_name': file_name or request.session.get('dpv_file_name')
    })


def eis_analysis_view(request):
    result_message = None

    if request.method == 'POST':
        form = EISAnalysisForm(request.POST, request.FILES)
        if form.is_valid():
            uploaded_file = request.FILES['file']
            ext = uploaded_file.name.split('.')[-1].lower()

            try:
                if ext == 'csv':
                    df = pd.read_csv(uploaded_file, header=0)
                elif ext in ['xls', 'xlsx']:
                    df = pd.read_excel(uploaded_file, header=0)
                else:
                    result_message = "Unsupported file format."
                    df = None

                if df is not None:
                    temp_dir = tempfile.mkdtemp()

                    result = Analysis_EIS(
                        df=df,
                        values_row_start=1,
                        real_col=3,
                        imag_col=4,
                        x_start=None,
                        x_end=None,
                        y_start=None,
                        y_end=None,
                        unit="Ω",
                        circle_pt1_index=0,
                        circle_pt2_index=0,
                        saving_folder=temp_dir
                    )

                    plot_path = result.get('plot_path')
                    result_message = f"✅ EIS Analysis complete. Plot saved to: {plot_path}"
            except Exception as e:
                result_message = f"❌ Error: {e}"

    else:
        form = EISAnalysisForm()

    return render(request, 'easi_analysis/eis.html', {'form': form, 'result': result_message})



