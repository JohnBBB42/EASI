import os
import sys
import pylab as p
import numpy as np
from electrochemical_methods.basics import Load_data, Table, getFilepath
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression

def Analysis_CV(df, values_row_start, potential_column, current_column, scan_column,
                scan_number, linreg_start_index, r2_accept_value, potential_unit,
                current_unit, num_decimals, saving_folder):
    
    df = df.iloc[values_row_start-1:].copy()
    df.columns = range(df.shape[1])  # reset columns to integers for simplicity

    if scan_column > 0:
        df = df[df[scan_column-1] == scan_number].reset_index(drop=True)

    font_size = 14

    def Peak_finder(Data):
        upperPeak_index = Data.iloc[:, 1].idxmax()
        lowerPeak_index = Data.iloc[:, 1].idxmin()
    
        x_upperPeak = Data.iloc[upperPeak_index, 0]
        y_upperPeak = Data.iloc[upperPeak_index, 1]
        x_lowerPeak = Data.iloc[lowerPeak_index, 0]
        y_lowerPeak = Data.iloc[lowerPeak_index, 1]
    
        min_potential = Data.iloc[:, 0].idxmin()
        max_potential = Data.iloc[:, 0].idxmax()
    
        def safe_baseline(idx_peak, idx_extreme):
            start_idx = min(idx_peak, idx_extreme) + linreg_start_index
            end_idx = min(max(idx_peak, idx_extreme), len(Data)-1)
    
            # Ensure safe indices
            start = min(start_index, Data.shape[0]-1)
            end = min(idx_extreme + 100, Data.shape[0]-1)
    
            if start >= end:
                start, end = idx_extreme, Data.shape[0]-1
    
            x_lin = Data.iloc[start:end, 0].values.reshape(-1,1)
            y_lin = Data.iloc[start:end, 1].values.reshape(-1,1)
    
            fit = LinearRegression()
            fit.fit(x_lin, y_lin)
            y_pred = fit.intercept_ + fit.coef_[0] * Data.iloc[:, 0]
    
            return y_pred
    
        y_pred1 = safe_baseline(upperPeak_index, max_potential)
        y_pred2 = safe_baseline(lowerPeak_index, min_potential)
    
        return y_pred1, y_pred2, x_upperPeak, y_upperPeak, x_lowerPeak, y_lowerPeak, y_pred1[upperPeak_index], y_pred2[lowerPeak_index]
    
        y_pred1, y_pred2, x_upperPeak, y_upperPeak, x_lowerPeak, y_lowerPeak, y_upper_baseline, y_lower_baseline = Peak_finder(df)
    
        fig, ax = plt.subplots()
        ax.plot(df.iloc[:, 0], df.iloc[:, 1], label=f'Scan: {scan_number}')
        ax.plot(df.iloc[:, 0], y_pred1, color="green")
        ax.plot(df.iloc[:, 0], y_pred2, color="grey")
    
        ax.arrow(x_upperPeak, y_upperPeak, 0, y_upper_baseline-y_upperPeak, color="green",
                 head_width=0.02, head_length=0.02)
        ax.arrow(x_lowerPeak, y_lowerPeak, 0, y_lower_baseline-y_lowerPeak, color="grey",
                 head_width=0.02, head_length=0.02)
    
        ax.set_xlabel(f"Potential ({potential_unit})", fontsize=font_size)
        ax.set_ylabel(f"Current ({current_unit})", fontsize=font_size)
        ax.set_title(f"Cyclic Voltammetry\nOxidation peak: E_pa={x_upperPeak:.{num_decimals}g}{potential_unit}, i_pa={(y_upperPeak - y_upper_baseline):.{num_decimals}g}{current_unit}\n"
                     f"Reduction peak: E_pc={x_lowerPeak:.{num_decimals}g}{potential_unit}, i_pc={(y_lowerPeak - y_lower_baseline):.{num_decimals}g}{current_unit}", fontsize=font_size)
    
        plt.legend()
        plt.tight_layout()
        plot_path = os.path.join(saving_folder, f'CV_analysis_scan_{scan_number}.pdf')
        plt.savefig(plot_path)
        plt.show()
    
        results = {
            "E_pa": x_upperPeak,
            "E_pc": x_lowerPeak,
            "E_diff": abs(x_upperPeak - x_lowerPeak),
            "i_pa": y_upperPeak - y_upper_baseline,
            "i_pc": y_lowerPeak - y_lower_baseline,
            "plot_path": plot_path
        }
    
        return results

        
    def CV_plot(filename, Data, number_decimals, y_pred1, y_pred2, x_upperPeak, y_upperPeak, x_lowerPeak, y_lowerPeak, y_upperPeak_baseline, y_lowerPeak_baseline):
        num_dec = '.'+'%s'% +num_decimals+'g'
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1)
        #Plot baselines
        plt.plot(Data.iloc[:,0], y_pred1, color = "green")
        plt.plot(Data.iloc[:,0], y_pred2, color = "grey")
        
        #Draw arrows to indicate upper- and lower peak
        p.arrow( x_upperPeak, y_upperPeak, 0.0, y_upperPeak_baseline-y_upperPeak, fc="green", ec="green",head_width=0.001*y_upperPeak_baseline, head_length=0.001*y_upperPeak_baseline)
        p.arrow( x_lowerPeak, y_lowerPeak, 0.0, y_lowerPeak_baseline-y_lowerPeak, fc="grey", ec="grey",head_width=0.001*y_lowerPeak_baseline, head_length=0.001*y_lowerPeak_baseline)
        #Make the CV plot
        plt.plot(Data.iloc[:,0], Data.iloc[:,1], label='%s' % filename + '   Scan: ' + str(Scan_number))
        plt.rcParams["font.family"] = "georgia"
        #Add legend
        plt.legend(loc = 'lower right')
        #Add labels
        plt.xlabel(str("Potential applied / ")+str(potential_unit),fontsize=font_size)
        plt.ylabel(str("WE Current / ")+str(current_unit),fontsize=font_size)
        plt.title("Cyclic voltammetry, IUPAC conv.\n Oxidation peak: "+"E_pa="+ str(format(x_upperPeak,num_dec))+str(potential_unit)+",  i_pa="+
                  str(format(y_upperPeak-y_upperPeak_baseline,num_dec))+str(current_unit)+"\n" + "Reduction peak: "+"E_pc="+ 
                  str(format(x_lowerPeak,num_dec))+str(potential_unit)+",  i_pc="+str(format(y_lowerPeak-y_lowerPeak_baseline,num_dec))+str(current_unit),fontsize=font_size)
        ax1.yaxis.set_label_position("left")
        print(str("File: "),filename)
        print("|\u0394E|=", str(format(np.abs(x_upperPeak-x_lowerPeak),num_dec)),str(potential_unit)) #Find potential difference
        print("Oxidation peak: E_pa=", str(format(x_upperPeak,num_dec)),str(potential_unit)+",  i_pa=",str(format(y_upperPeak-y_upperPeak_baseline,num_dec)),str(current_unit))
        print("Reduction peak: E_pc=", str(format(x_lowerPeak,num_dec)),str(potential_unit)+",  i_pc=",str(format(y_lowerPeak-y_lowerPeak_baseline,num_dec)),str(current_unit))
        return ax1
    
    #_________________________________________
    #_________________________________________
    # CV analysis
    #_________________________________________
    try:
        #Load data
        if scan_column > 0: #check whether scan_column is defined or not
            use_cols = [potential_column-1,current_column-1, scan_column-1] #columns to extract from imported data               
        else:
            use_cols = [potential_column-1,current_column-1] #columns to extract from imported data 
        message = "Select your CV file(s)"
        Loaded_data, filename, name = Load_data(message, use_cols,header=None, skiprows=values_row_start-1)
        
        #Choose where to save plots
        savingFolder = str(getFilepath("Choose saving folder"))
        
        E_pa = []
        E_pc = []
        E_dif = []
        i_pa = []
        i_pc = []
        
        try:
            for j in range(0,len(filename)):   
                num_dec = '.'+'%s'% +num_decimals+'g'
                #Define only data by selected scan
                if scan_column > 0:
                    Scan_select = Loaded_data[j][Loaded_data[j].iloc[:,2] == Scan_number] 
                else: 
                    Scan_select = Loaded_data[j]
                Scan_select = Scan_select.reset_index(drop=True)
                
                #Find the peaks
                y_pred1, y_pred2, x_upperPeak, y_upperPeak, x_lowerPeak, y_lowerPeak, y_upperPeak_baseline, y_lowerPeak_baseline = Peak_finder(Scan_select, LinReg_start_index, R2_accept_value, Scan_number)
        
                E_pa.append(float(format(x_upperPeak,num_dec)))
                E_pc.append(float(format(x_lowerPeak,num_dec)))
                E_dif.append(float(format(np.abs(x_upperPeak-x_lowerPeak),num_dec)))
                i_pa.append(float(format(y_upperPeak-y_upperPeak_baseline,num_dec)))
                i_pc.append(float(format(y_lowerPeak-y_lowerPeak_baseline,num_dec)))
                
                #Make a plot of the baselines and found peaks
                ax = CV_plot(filename[j], Scan_select, num_decimals, y_pred1, y_pred2, x_upperPeak, y_upperPeak, x_lowerPeak, y_lowerPeak, y_upperPeak_baseline, y_lowerPeak_baseline)
                #Display ticks as scientific if they are rounded to 0
                if round(max(Loaded_data[j].iloc[:,0])) == 0: 
                    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
                if round(max(Loaded_data[j].iloc[:,1])) == 0:
                    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
                plt.savefig(os.path.join(savingFolder, 'CV_with_baselines-' + str(filename[j]) + '.pdf'),bbox_inches='tight')
                plt.show()
                
            #_________________________________________
            #_________________________________________
            # OVERVIEW TABLE
            #_________________________________________
    
            #Table of characteristic values
            #Define column labels
            columns = ["E_pa / "+str(potential_unit), "E_pc / "+str(potential_unit),"|\u0394E| / "+str(potential_unit), 
                       "i_pa / "+str(current_unit), "i_pc / "+str(current_unit)]
            #Define row labels
            rows = []
            for i in filename:
                rows.append(str(i))
            #Make value matrix
            value_matrix = np.array(([E_pa,E_pc,E_dif,i_pa,i_pc])).T
               
            #Table - Main table
            print(str("Overview of imported data:"))
            table = Table(value_matrix,rows,columns,2,2,figsize=len(columns))
            plt.tight_layout()
            #Save plot
            plt.savefig(os.path.join(savingFolder, 'CV_Table-' + filename[0] + '.pdf'),bbox_inches='tight')
            plt.show()
        except:
            pass
        
        #Make a plot of all selected files, so redo same procedure as CV analysis step
        print("All in one plot: ")
        fig = plt.figure()
        ax3 = fig.add_subplot(1,1,1)
        markers = ['.',',','o','v','^','>','<','1','2','3','4','s','p','*','h','H','+','x','D','d','|','_']
        for j in range(0,len(filename)):    
            #Define only data by selected scan
            if scan_column > 0:
                Scan_select = Loaded_data[j][Loaded_data[j].iloc[:,2] == Scan_number]
                legend_label = '%s' % filename[j] + '   Scan: ' + str(Scan_number)
            else: 
                Scan_select = Loaded_data[j]
                legend_label = '%s' % filename[j] 
            Scan_select = Scan_select.reset_index(drop=True)
            plt.plot(Scan_select.iloc[:,0], Scan_select.iloc[:,1], label=legend_label, marker=markers[j])
            plt.rcParams["font.family"] = "georgia"
            #Display ticks as scientific if they are rounded to 0
            if round(max(Loaded_data[j].iloc[:,0])) == 0:
                ax3.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
            if round(max(Loaded_data[j].iloc[:,1])) == 0:
                ax3.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
            
        # Make legend to the right of the plot
        box = ax.get_position()
        ax3.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        #ax3.legend(loc='lower right')    
        plt.xlabel(str("Potential applied / ")+str(potential_unit),fontsize=font_size)
        plt.ylabel(str("WE current / ")+str(current_unit),fontsize=font_size)
        plt.title("Cyclic voltammetry, IUPAC conv.",fontsize=font_size)
        ax3.yaxis.set_label_position("left")
        
        #Save plot
        plt.savefig(os.path.join(savingFolder, 'CV_Plot-' + filename[0] + '.pdf'),bbox_inches='tight')
        plt.show()
    
    except:
        sys.exit("Error in loading data. Please check you have selected a cyclic voltammetry excel/csv file and written the correct columns to be included.")
