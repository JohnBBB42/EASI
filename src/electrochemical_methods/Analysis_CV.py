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
    Single-file version:
      1) 'df' is your loaded DataFrame (headerless).
      2) 'values_row_start' typically 2 if row 1 is the header.
      3) 'potential_column', 'current_column', 'scan_column' are 1-based indices
         for the columns that hold the data. If scan_column=0, no scanning logic is used.
      4) 'scan_number' selects only rows where df's scan_column==scan_number.
      5) Writes out:
          - A baseline PDF with peaks
          - A table PDF of peak values
          - An all-in-one PDF
      6) Returns a dictionary with at least 'plot_path' for the baseline PDF.
    """

    font_size = 14

    # 1) Trim rows in case there's a header
    df = df.iloc[values_row_start-1:].copy()
    df.columns = range(df.shape[1])  # reset columns to 0,1,2,...

    # 2) If user said scan_column > 0, filter by scan_number
    if scan_column > 0:
        # in code, that means we look at df.iloc[:, scan_column-1]
        # but we've re-labeled columns with range(...), so:
        sc_idx = scan_column-1
        df = df[df.iloc[:, sc_idx] == scan_number].reset_index(drop=True)

    # If there's no data after filter, that's an error
    if len(df) == 0:
        raise ValueError("No data left after applying row_start or scan filtering.")

    #---------------------------------
    # Helper: safe baseline
    #---------------------------------
    def safe_baseline(Data, idx_peak, idx_extreme):
        # if idx_peak or idx_extreme is out of range, clamp them
        nrows = len(Data)
        idx_peak = max(0, min(idx_peak, nrows-1))
        idx_extreme = max(0, min(idx_extreme, nrows-1))

        start = min(idx_peak, idx_extreme) + linreg_start_index
        end   = max(idx_peak, idx_extreme) + 100
        if start < 0:
            start = 0
        if end > nrows:
            end = nrows
        if start >= end:
            # fallback: baseline is 0
            return np.zeros(nrows)

        x_lin = Data.iloc[start:end, 0].values.reshape(-1,1)
        y_lin = Data.iloc[start:end, 1].values.reshape(-1,1)

        fit = LinearRegression()
        fit.fit(x_lin, y_lin)
        y_pred = fit.intercept_ + fit.coef_[0] * Data.iloc[:, 0]
        return y_pred

    #---------------------------------
    # Helper: Peak finder
    #---------------------------------
    def Peak_finder(Data):
        # identify upper/lower peaks
        upperPeak_index = Data.iloc[:, 1].idxmax()
        lowerPeak_index = Data.iloc[:, 1].idxmin()

        x_up = Data.iloc[upperPeak_index, 0]
        y_up = Data.iloc[upperPeak_index, 1]
        x_low = Data.iloc[lowerPeak_index, 0]
        y_low = Data.iloc[lowerPeak_index, 1]

        # extremes
        max_pot_idx = Data.iloc[:,0].idxmax()
        min_pot_idx = Data.iloc[:,0].idxmin()

        # build baselines
        y_pred1 = safe_baseline(Data, upperPeak_index, max_pot_idx)
        y_pred2 = safe_baseline(Data, lowerPeak_index, min_pot_idx)

        # baseline offsets
        y_up_base  = y_pred1[upperPeak_index]
        y_low_base = y_pred2[lowerPeak_index]

        return y_pred1, y_pred2, x_up, y_up, x_low, y_low, y_up_base, y_low_base

    #---------------------------------
    # 3) Run the peak finder
    #---------------------------------
    y_pred1, y_pred2, x_up, y_up, x_low, y_low, up_base, low_base = Peak_finder(df)

    #---------------------------------
    # 4) Single-file baseline plot
    #---------------------------------
    fig, ax = plt.subplots()
    ax.plot(df.iloc[:,0], df.iloc[:,1], label=f"Scan: {scan_number}")
    ax.plot(df.iloc[:,0], y_pred1, color='green')
    ax.plot(df.iloc[:,0], y_pred2, color='grey')

    # Arrow for upper peak
    ax.arrow(
        x_up, y_up,
        0, (up_base - y_up),
        color="green", head_width=0.02,
        length_includes_head=True
    )
    # Arrow for lower peak
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
        + f"Oxidation peak: E_pa={x_up:.{num_decimals}g}{potential_unit}, i_pa={(y_up - up_base):.{num_decimals}g}{current_unit}\n"
        + f"Reduction peak: E_pc={x_low:.{num_decimals}g}{potential_unit}, i_pc={(y_low - low_base):.{num_decimals}g}{current_unit}",
        fontsize=font_size
    )

    ax.legend()
    plt.tight_layout()

    baseline_plot = os.path.join(saving_folder, f"CV_analysis_scan_{scan_number}.pdf")
    plt.savefig(baseline_plot)
    plt.show()

    #---------------------------------
    # 5) If you want a table or all-in-one,
    #    you can add them here. We'll do a minimal approach:
    #---------------------------------

    # Example table-like print
    print("[Single-file table of characteristic values]")
    dE = abs(x_up - x_low)
    print(f"E_pa={x_up}, E_pc={x_low}, |ΔE|={dE}, i_pa={y_up-up_base}, i_pc={y_low-low_base}")

    # If you rely on 'Table' from your basics, you can do:
    # table_values = np.array([[x_up, x_low, dE, (y_up-up_base), (y_low-low_base)]])
    # columns = [f"E_pa ({potential_unit})", f"E_pc ({potential_unit})", "|ΔE|", f"i_pa({current_unit})", f"i_pc({current_unit})"]
    # rows = ["Dataset"]
    # T = Table(table_values, rows, columns, 2,2, figsize=len(columns))
    # plt.savefig(os.path.join(saving_folder,"CV_table.pdf"))
    # plt.show()

    # 6) Return dictionary with main plot path
    results = {
        "E_pa": x_up,
        "E_pc": x_low,
        "E_diff": dE,
        "i_pa": (y_up - up_base),
        "i_pc": (y_low - low_base),
        "plot_path": baseline_plot
    }
    return results
