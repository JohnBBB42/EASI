#_________________________________________
# EIS ANALYSIS
#_________________________________________
def Analysis_EIS(
    df,
    values_row_start=2,
    real_col=1,
    imag_col=2,
    x_start=None,
    x_end=None,
    y_start=None,
    y_end=None,
    unit="Ω",
    circle_pt1_index=0,
    circle_pt2_index=0,
    saving_folder="."
):
    """
    Single-file EIS analysis:
      1) 'df' is already loaded externally (no repeated Load_data).
      2) 'values_row_start' to skip header rows if needed.
      3) 'real_col','imag_col' => 1-based column indices (like your old code).
      4) 'x_start','x_end','y_start','y_end' => optionally limit x/y axis.
      5) 'circle_pt1_index','circle_pt2_index' => for circle fitting.
      6) 'unit' => label for your impedance (Ω or something).
      7) 'saving_folder' => a directory path where we save PDFs.

    The function:
      - Trims rows,
      - Sorts df by real axis,
      - Does circle fitting,
      - Plots the single-file Nyquist,
      - Creates the same style diameter bar chart (but for 1 file it won't show much).
    """

    import os, sys
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    from sklearn.metrics import r2_score
    from sklearn.linear_model import LinearRegression
    from scipy.optimize import curve_fit
    # If you rely on your local Table function:
    # from electrochemical_methods.basics import Table

    # 1) Decide the analysis type (impedance or capacitance)
    Analysis = ["impedance", str(unit)]  # default

    # 2) Convert your old variables
    values_row_start = int(values_row_start)
    Real_column      = int(real_col)
    Imag_column      = int(imag_col)

    # 3) Convert user-limits for x, y
    def parse_limit(val_start, val_end):
        if val_start is None or val_end is None or val_start == "" or val_end == "":
            return []
        return [float(val_start), float(val_end)]

    x_lim_plot = parse_limit(x_start, x_end)
    y_lim_plot = parse_limit(y_start, y_end)

    # 4) Indices for circle fitting
    cir_pt1_index = int(circle_pt1_index)
    cir_pt2_index = int(circle_pt2_index)

    num_decimals = 3
    font_size    = 14

    # 5) Trim the DataFrame, skip the first (values_row_start - 1) lines
    df = df.iloc[values_row_start - 1:].copy()
    df.columns = range(df.shape[1])  # rename columns to 0,1,2,..

    # Real/Imag are now:
    real_idx = Real_column - 1  # zero-based index
    imag_idx = Imag_column - 1  # zero-based

    # 6) Sort the DataFrame by Real axis
    df.sort_values(by=df.columns[real_idx], inplace=True)

    # 7) Optionally limit the data by x-limits
    if len(x_lim_plot) == 2:
        df = df[df.iloc[:, real_idx] > x_lim_plot[0]]
        df = df[df.iloc[:, real_idx] < x_lim_plot[1]]

    # 8) Circle + linear fit logic
    #    We'll do it all for a single dataset now
    sorted_data = df.reset_index(drop=True)
    x_lin1 = sorted_data.iloc[:, real_idx]
    y_lin1 = sorted_data.iloc[:, imag_idx]

    # Because your old logic flips y to detect slope:
    y_lin1_flipped = y_lin1[::-1]
    slope_2 = np.gradient(y_lin1_flipped)
    slope_2_index = np.where(np.round(slope_2,0) > -1)

    print("Single File EIS Analysis")

    if len(slope_2_index[0]) == 0:
        # no slope found => entire circle?
        circle_end_index = len(y_lin1_flipped)
        print("No linear curve was found.")
    else:
        slope_2_value = slope_2_index[0][0]
        circle_end_index = len(y_lin1_flipped) - slope_2_value
        # linear regression from the 'end' of circle to last point:
        x_reg = np.array(sorted_data.iloc[circle_end_index:, real_idx]).reshape(-1,1)
        y_reg = np.array(sorted_data.iloc[circle_end_index:, imag_idx]).reshape(-1,1)
        fit1  = LinearRegression().fit(x_reg, y_reg)
        y_pred1 = fit1.intercept_ + fit1.coef_[0]*sorted_data.iloc[:,real_idx]
        R2_1_score = r2_score(y_reg.ravel(), y_pred1[circle_end_index:].to_numpy())

        num_dec = '.' + f'{num_decimals}g'
        print("Linear function:",
              f"{fit1.coef_[0][0]:{num_dec}} x + {fit1.intercept_[0]:{num_dec}}")
        if circle_end_index < len(x_lin1):
            print("Linear interval: [",
                  f"{sorted_data.iloc[circle_end_index,real_idx]:{num_dec}},",
                  f"{sorted_data.iloc[len(x_lin1)-1,real_idx]:{num_dec}} ]")

    # 9) Circle fit
    def func_circle(x, X_c, R):
        # Y^2 = R^2 - (x - X_c)^2
        return R**2 - (x - X_c)**2

    # pick data from cir_pt1_index to circle_end_index-cir_pt2_index
    xdata = sorted_data.iloc[cir_pt1_index: circle_end_index - cir_pt2_index, real_idx]
    ydata = sorted_data.iloc[cir_pt1_index: circle_end_index - cir_pt2_index, imag_idx]

    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(func_circle, xdata.values, (ydata.values**2))
    # diameter:
    D = 2 * popt[1]
    print("Diameter:", D)

    # 10) For single-file, we can do a bar chart or table, but it’s just 1 file
    # We'll skip the multi-file difference table logic.
    # We'll show a single "Resist Charge Transfer" line:
    print(f"Resist Charge Transfer: {D:.3f} {unit}")

    # 11) Plot the all-in-one
    symbol = "Z"
    Title  = "Impedance Nyquist"
    # If you had 'capacitance', you'd do something else

    fig = plt.figure()
    plt.rcParams["font.family"] = "georgia"
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal', adjustable='box')

    # single-file
    ax.plot(sorted_data.iloc[:,real_idx],
            sorted_data.iloc[:,imag_idx],
            'o-', label='Data')

    ax.legend(loc='best')
    ax.set_xlabel(f"{symbol}'  ; Re({symbol}) / {unit}", fontsize=font_size)
    ax.set_ylabel(f"-{symbol}''  ; -Im({symbol}) / {unit}", fontsize=font_size)
    ax.set_title(Title, fontsize=font_size)

    if len(x_lim_plot) == 2:
        plt.xlim(x_lim_plot[0], x_lim_plot[1])
    if len(y_lim_plot) == 2:
        plt.ylim(y_lim_plot[0], y_lim_plot[1])

    # If near zero, use scientific notation
    if round(sorted_data.iloc[:,real_idx].max()) == 0:
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    if round(sorted_data.iloc[:,imag_idx].max()) == 0:
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

    plot_name = os.path.join(saving_folder, "EIS_single_file_plot.pdf")
    plt.savefig(plot_name, bbox_inches='tight')
    plt.show()

    print(f"Plot saved: {plot_name}")

    # Return dictionary with relevant data
    return {
      "diameter": D,
      "plot_path": plot_name
    }
