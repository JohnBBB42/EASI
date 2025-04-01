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

    if request.method == 'POST':
        form = DPVAnalysisForm(request.POST, request.FILES)

        if form.is_valid():
            selected_analysis = form.cleaned_data['selected_analysis']
            plot_dir = os.path.join('easi_analysis', 'static', 'easi_analysis', 'plots')
            os.makedirs(plot_dir, exist_ok=True)
            for f in os.listdir(plot_dir):
                os.remove(os.path.join(plot_dir, f))  # Clear old plots

            file = request.FILES.get('file')  # get() won't crash if not present
            blank_file = request.FILES.get('blank_file')

            # ✅ Handle file upload or fallback to session
            if file:
                file_path = os.path.join(plot_dir, file.name)

                with open(file_path, 'wb+') as dest:
                    for chunk in file.chunks():
                        dest.write(chunk)
                file_name = file.name
                request.session['dpv_file_path'] = file_path
                request.session['dpv_file_name'] = file.name

                # ✅ Now the file exists → it's safe to load it
                try:
                    temp_array, headers, metadata = Load_data(file_path)
                    plot_choices = [(str(i), f"Plot {i//2 + 1}") for i in range(0, len(headers), 2)]
                    form.fields['selected_plots'].choices = plot_choices
                except:
                    pass
            else:
                file_path = request.session.get('dpv_file_path')
                file_name = request.session.get('dpv_file_name')
                if not file_path or not os.path.exists(file_path):
                    result_message = "❌ No file selected and no previous file found."
                    return render(request, 'easi_analysis/dpv.html', {
                        'form': form,
                        'result': result_message,
                        'metrics': None,
                        'plot_files': [],
                        'file_name': None
                    })

            # ✅ Handle blank file as usual
            blank_array = None
            if blank_file:
                blank_path = os.path.join(plot_dir, blank_file.name)
                with open(blank_path, 'wb+') as dest:
                    for chunk in blank_file.chunks():
                        dest.write(chunk)
                try:
                    blank_array = np.loadtxt(blank_path, delimiter=',')
                except Exception as e:
                    result_message = f"Blank file failed to load: {e}"

            # ✅ Analysis logic (same as before)
            try:
                analysis_results = Analysis_DPV(file_path, blank_responses=blank_array)
                results = {}
                plotter = DPVPlotting()

                selected_plot_indices = form.cleaned_data.get('selected_plots', [])
                selected_indices = [int(idx) for idx in selected_plot_indices] if selected_plot_indices else None


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
                    plot_filenames.append("easi_analysis/plots/plot_corrected.png")
                    result_message = "✅ Plot Data Array with Corrected Baseline created!"

                elif selected_analysis == "Analyze Peak Currents":
                    peak_paths = plotter.analyze_peak_currents(
                        analysis_results['mean_peak_currents'],
                        analysis_results['std_peak_currents'],
                        plot_dir,
                        selected_indices=selected_indices
                        
                    )

                    for p in peak_paths:
                        # ✅ No need to copy again — it's already saved in plot_dir
                        plot_filenames.append("easi_analysis/plots/" + os.path.basename(p))
                    result_message = "✅ Peak currents plot created!"

                elif selected_analysis == "Observed vs Expected Concentration":
                    obs_paths = plotter.plot_observed_vs_expected_concentration(
                        analysis_results['mean_peak_currents'],
                        plot_dir,
                        selected_indices=selected_indices
                    )

                    for p in obs_paths:
                        plot_filenames.append("easi_analysis/plots/" + os.path.basename(p))
                    result_message = "✅ Observed vs Expected plot created!"



                # 🖼 Plot copying (unchanged)
                if 'plot_filenames' in analysis_results:
                    for key, file_or_list in analysis_results['plot_filenames'].items():
                        if isinstance(file_or_list, list):
                            for i, file_path in enumerate(file_or_list):
                                if os.path.exists(file_path):
                                    filename = f"{key}_{i+1}.png"
                                    new_path = os.path.join(plot_dir, filename)
                                    shutil.copy(file_path, new_path)
                                    plot_filenames.append(f"easi_analysis/plots/{filename}")
                        else:
                            if os.path.exists(file_or_list):
                                filename = f"{key}.png"
                                new_path = os.path.join(plot_dir, filename)
                                shutil.copy(file_or_list, new_path)
                                plot_filenames.append(f"easi_analysis/plots/{filename}")

            except Exception as e:
                result_message = f"❌ DPV Analysis failed: {e}"

    else:
        form = DPVAnalysisForm()

    return render(request, 'easi_analysis/dpv.html', {
        'form': form,
        'result': result_message,
        'metrics': results,
        'plot_files': plot_filenames,
        'file_name': file_name if 'file_name' in locals() else request.session.get('dpv_file_name')
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



