import os
import sys
import pylab as p
import numpy as np
from electrochemical_methods.basics import Load_data, Table, getFilepath
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression

def Analysis_CV(
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
    saving_folder="."
):
    """
    Single-File Approach:
      - 'df' is a user-loaded DataFrame (headerless).
      - 'values_row_start' determines which row to start from. 2 => skip 1 line
      - 'potential_column', 'current_column', 'scan_column' are 1-based columns
         (0 means 'no scan column used').
      - If 'scan_column'>0 => filter df by (df[scan_column-1] == scan_number).
      - Plot baseline, table, and 'all in one' plot, saving them as PDF to saving_folder.
    """

    import pylab as p
    font_size = 14

    # Trim rows if user wants to skip the first few lines
    df = df.iloc[values_row_start-1:].copy()
    df.columns = range(df.shape[1])  # rename columns to 0,1,2,...

    # If scan_column > 0, filter on scan_number
    if scan_column > 0:
        sc_idx = scan_column - 1
        df = df[df.iloc[:, sc_idx] == scan_number].reset_index(drop=True)

    if df.empty:
        raise ValueError("No data left after applying row_start/scan filter.")

    #----------------------------------
    # 1) We'll define the 'Peak_finder' function
    #----------------------------------
    def Peak_finder(Data, LinReg_start_index, R2_accept_value, Scan_number):
        upperPeak_index = Data.iloc[:,1].idxmax()
        lowerPeak_index = Data.iloc[:,1].idxmin()

        x_upperPeak = Data.iloc[upperPeak_index, 0]
        y_upperPeak = Data.iloc[upperPeak_index, 1]
        x_lowerPeak = Data.iloc[lowerPeak_index, 0]
        y_lowerPeak = Data.iloc[lowerPeak_index, 1]

        max_potential = Data.iloc[:,0].idxmax()
        min_potential = Data.iloc[:,0].idxmin()

        def safe_baseline(idx_peak, idx_extreme):
            # clamp indices
            nrows = len(Data)
            idx_peak    = max(0, min(idx_peak, nrows-1))
            idx_extreme = max(0, min(idx_extreme, nrows-1))

            start = min(idx_peak, idx_extreme) + LinReg_start_index
            end   = max(idx_peak, idx_extreme) + 100
            if start < 0: start = 0
            if end   > nrows: end = nrows
            if start >= end:
                return np.zeros(nrows)  # fallback

            x_lin = Data.iloc[start:end, 0].values.reshape(-1,1)
            y_lin = Data.iloc[start:end, 1].values.reshape(-1,1)

            fit = LinearRegression().fit(x_lin, y_lin)
            y_pred = fit.intercept_ + fit.coef_[0]*Data.iloc[:,0]
            return y_pred

        y_pred1 = safe_baseline(upperPeak_index, max_potential)
        y_pred2 = safe_baseline(lowerPeak_index, min_potential)

        y_upperPeak_baseline = y_pred1[upperPeak_index]
        y_lowerPeak_baseline = y_pred2[lowerPeak_index]

        return (
            y_pred1, y_pred2,
            x_upperPeak, y_upperPeak,
            x_lowerPeak, y_lowerPeak,
            y_upperPeak_baseline, y_lowerPeak_baseline
        )

    #----------------------------------
    # 2) Single-file approach => we do everything in one pass
    #----------------------------------
    # columns are potential_column-1 for x, current_column-1 for y
    xcol = potential_column - 1
    ycol = current_column - 1

    # We'll rename them so it's consistent with the old code: 0 => potential, 1 => current
    # if xcol != 0 or ycol != 1, we reorder columns
    if xcol != 0 or ycol != 1:
        df = df[[xcol, ycol] + [c for c in df.columns if c not in [xcol, ycol]]]

    # Now the potential is in df.iloc[:,0], current in df.iloc[:,1]
    # Let's do the analysis
    dec_fmt = '.' + '%s' % +num_decimals + 'g'

    # 2a) Baseline + peaks
    y_pred1, y_pred2, x_up, y_up, x_low, y_low, up_base, low_base = Peak_finder(
        df, linreg_start_index, r2_accept_value, scan_number
    )

    # 2b) Plot (single-file baseline)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(df.iloc[:,0], df.iloc[:,1], label=f"Scan: {scan_number}")
    ax.plot(df.iloc[:,0], y_pred1, color='green')
    ax.plot(df.iloc[:,0], y_pred2, color='grey')

    # arrows
    ax.arrow(
        x_up, y_up,
        0, (up_base - y_up),
        color="green", head_width=0.02,
        length_includes_head=True
    )
    ax.arrow(
        x_low, y_low,
        0, (low_base - y_low),
        color="grey", head_width=0.02,
        length_includes_head=True
    )

    ax.set_xlabel(f"Potential ({potential_unit})", fontsize=font_size)
    ax.set_ylabel(f"Current ({current_unit})", fontsize=font_size)
    ax.set_title(
        "Cyclic Voltammetry\n"
        + f"Oxidation peak: E_pa={x_up:.{num_decimals}g}{potential_unit}, "
          f"i_pa={(y_up - up_base):.{num_decimals}g}{current_unit}\n"
        + f"Reduction peak: E_pc={x_low:.{num_decimals}g}{potential_unit}, "
          f"i_pc={(y_low - low_base):.{num_decimals}g}{current_unit}",
        fontsize=font_size
    )
    ax.legend()
    plt.tight_layout()

    baseline_plot = os.path.join(saving_folder, f"CV_analysis_scan_{scan_number}.pdf")
    plt.savefig(baseline_plot)
    plt.show()

    # 2c) Print summary
    print("=== CV Analysis Summary ===")
    dE = abs(x_up - x_low)
    print(f"E_pa={x_up}, E_pc={x_low}, ΔE={dE}")
    print(f"i_pa={y_up-up_base}, i_pc={y_low-low_base}")

    # 2d) If you want a table or all-in-one plot, add them here.
    # We'll skip since this is single-file. (Add your Table code if needed.)

    # Return a dictionary with the path (so the GUI can optionally re-save it).
    return {
        "plot_path": baseline_plot,
        "E_pa": x_up,
        "E_pc": x_low,
        "E_diff": dE,
        "i_pa": (y_up - up_base),
        "i_pc": (y_low - low_base)
    }
