#_________________________________________
# EIS ANALYSIS
#_________________________________________
def Analysis_EIS(values_row_start_get, x_column_get,y_column_get,x_start_get, x_end_get,
                 y_start_get, y_end_get, unit,cir_pt1_index, cir_pt2_index):
    #_________________________________________
    # IMPORT PACKAGES
    #_________________________________________
    import os
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    from sklearn.metrics import r2_score
    from sklearn.linear_model import LinearRegression
    from electrochemical_methods.basics import getFilepath, Load_data, Table
    
    #_________________________________________
    #_________________________________________
    # VARIABLES
    #_________________________________________
    
    Analysis = ["impedance",str(unit)]  #[method,unit] 
                                  #method is "impedance" or "capacitance"
                                  #unit is only used as label
    #Analysis = ["capacitance","F"]
    
    values_row_start = int(values_row_start_get) #Check EIS/ECS excel file,
                         #if values start from 2nd row
                         #then values_row_start = 2
    Real_column = int(x_column_get) #Check EIS/ECS excel file,
                    #if Re(Z) or Re(C) values in third column, 
                    #then Real_column = 3
    Imag_column = int(y_column_get) #Check EIS/ECS excel file,
                    #if Im(Z) or Im(C) values in fourth column, 
                    #then Imag_column = 4 (assumed to be column of -Im(Z) or -Im(C))                              
    if len(x_start_get)==0 or len(x_end_get)==0:
        x_lim_plot = [ ]
    else:
        x_lim_plot = [float(x_start_get),float(x_end_get)] # if empty [ ] => no limit on plot x-axis is set! Otherwise, e.g [0,200]
    
    if len(y_start_get)==0 or len(y_end_get)==0:
        y_lim_plot = [ ]
    else:
        y_lim_plot = [float(y_start_get),float(y_end_get)] # if empty [ ] => no limit on plot y-axis is set! Otherwise, e.g [0,200]
    
    cir_pt1_index = int(cir_pt1_index) #to circle fit, nr. index away from first index
    cir_pt2_index = int(cir_pt2_index) #to circle fit, nr. index away to left of end of circle!!
    
    num_decimals = 3 #number of decimals to round numbers
    font_size = 14 #fontsize in labels
    
    
    #_________________________________________
    #_________________________________________
    # IMPORT DATA
    #_________________________________________
    try:
        #Load data
        use_cols = [Real_column-1, Imag_column-1] #columns to extract from imported data
        message = "Select your EIS file(s)"
        Loaded_data, filename, name = Load_data(message, use_cols, header=None, skiprows=values_row_start-1)
        
        #Only use data in selected x_limit
        for j in range(0,len(filename)):
            if len(x_lim_plot) == 2: 
                Loaded_data[j] = Loaded_data[j].sort_values(by=Loaded_data[j].columns[0])
                Loaded_data[j] = Loaded_data[j][Loaded_data[j].iloc[:,0]>x_lim_plot[0]]
                Loaded_data[j] = Loaded_data[j][Loaded_data[j].iloc[:,0]<x_lim_plot[1]]
    
        #Choose where to save plots
        savingFolder = str(getFilepath("Choose saving folder"))
        
        #_________________________________________
        #_________________________________________
        # Curve fit
        #_________________________________________
        
        num_dec = '.'+'%s'% +num_decimals+'g'
        D_list = []
        Lin_func_list = []
        Lin_interval = []
            
        try:
            for j in range(0,len(filename)):
                #_________________________________________
                #_________________________________________
                # Linear fit
                #_________________________________________ 
                
                #sort data from smallest to greatest x-value
                sorted_data = Loaded_data[j].sort_values(by=Loaded_data[j].columns[0])
                x_lin1 = sorted_data.iloc[:,0]
                y_lin1 = sorted_data.iloc[:,1]
    
                y_lin1 = y_lin1[::-1] # flips y-coordinates, so the characteristic linear curve gets negative slope
                slope_2 = np.gradient(y_lin1) # calculate slope
                slope_2_index = np.where(np.round(slope_2,0) > -1) # find index where linear curce stops (end of semicircle)
    
                print(str("File: "),filename[j])
                if slope_2_index[0][0] == 0: #if no positive slope was found then slope_2_index[0][0]=0
                    circle_end_index = len(y_lin1)
                    print("No linear curve was found.")
                else:
                    slope_2_value = slope_2_index[0][0] #index where slope > -1
                    circle_end_index = len(y_lin1)-slope_2_value 
                    x_reg = np.array(sorted_data.iloc[circle_end_index:,0]).reshape(-1,1) #x_data from end of circle to last point
                    y_reg = np.array(sorted_data.iloc[circle_end_index:,1]).reshape(-1,1) #y_data from end of circle to last point
                    fit1 = LinearRegression() 
                    fit1.fit(x_reg, y_reg) #make linear fit
                    y_pred1 = fit1.intercept_ + fit1.coef_[0]*sorted_data.iloc[:,0] #define regression on form y=b+ax
                    R2_1_score = r2_score(y_reg.reshape(1,-1)[0],y_pred1[circle_end_index:].to_numpy()) #compute R^2 value
                    
                    print(str("Linear function:"),format(fit1.coef_[0][0],num_dec),str("x +"), format(fit1.intercept_[0],num_dec))
                    print(str("Linear interval: ["),format(sorted_data.iloc[circle_end_index,0],num_dec),str(","),format(sorted_data.iloc[len(x_lin1)-1,0],num_dec),str("]"))
                #_________________________________________
                #_________________________________________
                # Circle fit
                #_________________________________________
                def func(x, X_c, R): #define function for "non-linear least squares" fit
                    return R**2-(x-X_c)**2 #Y^2=R^2-(x-X_c)^2, rewritten form of (x-X_c)^2+Y^2=R^2
                #only include data from selected start point to end of semicircle minus second selected point
                xdata = sorted_data.iloc[cir_pt1_index:circle_end_index-cir_pt2_index,0]
                #only include data from selected start point to end of semicircle minus second selected point
                ydata = sorted_data.iloc[cir_pt1_index:circle_end_index-cir_pt2_index,1]
                
                #solve non linear system
                popt, pcov = curve_fit(func, np.array(xdata), np.array(ydata**2))
                D=2*popt[1] #diameter of semicircle
                print("Diameter:",D)
                print(str("Resist Charge Transfer:"),format(D,num_dec))
                print("\n")
                D_list.append(float(format(D,num_dec)))
            #_________________________________________
            #_________________________________________
            # DIAMETER TABLE
            #_________________________________________
            # Make Diameter table
            if len(D_list)>1: #if only 1 selected file, no need to compare to any other file
                columns_1 = ["Resist Charge Transfer / "+str(Analysis[1])]
                rows_1 = []
                for i in filename: #filename is list of all selected files with their filename
                    rows_1.append(str(i))
                value_matrix = np.array((D_list)).reshape(1,len(D_list)).T #define value matrix
                        
                #Table - Main table
                print(str("Diameter overview:"))
                Table(value_matrix,rows_1,columns_1,scale_width=2,scale_height=2, figsize=1.5)
                plt.tight_layout()
                #Save plot
                plt.savefig(os.path.join(savingFolder, 'Diameter-' + filename[0] + '.pdf'), bbox_inches='tight')
                plt.show() 
            
            #_________________________________________
            #_________________________________________
            # DIAMETER DIFFERENCE TABLE
            #_________________________________________
            # Make Diameter difference table
            if len(D_list)>1: #if only 1 selected file, no need to compare to any other file
                columns = []
                
                rows = []
                for i in range(0,len(filename)): #create labels to top column and left row
                    columns.append(str(filename[i])+str("\n Resist Charge Transfer: ")+str(D_list[i])+str(" ")+str(Analysis[1]))
                    rows.append(str(filename[i])+str("\n Resist Charge Transfer: ")+str(D_list[i])+str(" ")+str(Analysis[1]))
    
                D_difference = []
                for j in range(0,len(D_list)):
                    for i in range(0,len(D_list)):
                        D_difference.append(D_list[i]-D_list[j])
                D_matrix = np.array(np.round(D_difference,5)).reshape(len(D_list),len(D_list))
                        
                #Table - Main table
                print(str("Resist Charge Transfer difference matrix, \u0394Diameter:"))
                Table(D_matrix,rows,columns,scale_width=2.5, scale_height=2,figsize=len(columns)) #this scale works best
                #Save plot
                plt.savefig(os.path.join(savingFolder, 'Dif_diameter' + filename[0] + '.pdf'),bbox_inches='tight')
                plt.show()
    
            #_________________________________________
            #_________________________________________
            # DIAMETER BAR CHART
            #_________________________________________
            # Make Diameter bar chart 
            if len(D_list)>1: #if only 1 selected file, no need to compare to any other file
                fig, ax = plt.subplots()  
                y_pos = np.arange(len(filename)) #for ticks on y-axis
                plt.barh(y_pos, D_list, align='center', alpha=0.8, color="cornflowerblue")
                plt.yticks(y_pos, filename) #place the ticks on y-axis
                plt.xlabel("Resist Charge Transfer / "+str(Analysis[1]))
                plt.title("Overview of Resist Charge Transfer for each file")
                plt.xlim((0,1.23*max(D_list))) #make xlim great enough for text to be next to bar
                #Position diameter value next to each bar
                for i, v in enumerate(D_list):
                    ax.text(v+(0.01*max(D_list)), i, str(v), color='black', fontweight='bold')
                #Display ticks as scientific if they are rounded to 0
                if round(max(Loaded_data[i].iloc[:,0])) == 0:
                    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
                if round(max(Loaded_data[i].iloc[:,1])) == 0:
                    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
                plt.savefig(os.path.join(savingFolder, 'Bar_diameter-' + filename[0] + '.pdf'),bbox_inches='tight')
                print("Diameter overview (Bar chart):")
                plt.show()
        except:
            pass
        #_________________________________________
        #_________________________________________
        # EIS/ECS plot
        #_________________________________________
        # Make the EIS/ECS plot
        # Based on input determine whether to label as EIS or ECS
        symbol = ["Z","C"]
        if str.lower(Analysis[0]) == "impedance":
            symbol = "Z"
            Title = "Impedance"
        elif str.lower(Analysis[0]) == "capacitance":
            symbol = "C"
            Title = "Capacitance"
        else:
            symbol = "unit"
            Title = "Undefined"
        
        print("All in one plot: ")
        fig = plt.figure()
        plt.rcParams["font.family"] = "georgia"
        
        #Constrain axes
        plt.gca().set_aspect('equal', adjustable='box')
        ax = fig.add_subplot(1, 1, 1)
        
        #Define list of different plotting markers
        markers = ['.',',','o','v','^','>','<','1','2','3','4','s','p','*','h','H','+','x','D','d','|','_']
        for i in range(0,len(filename)): 
            plt.plot(Loaded_data[i].iloc[:,0], Loaded_data[i].iloc[:,1], label='%s' % filename[i], marker=markers[i])
    
        # Put a legend to the right of the current axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        #Make labels 
        plt.xlabel(str(symbol)+"'  ; Re("+str(symbol)+") / "+str(Analysis[1]),fontsize=font_size)
        plt.ylabel("-"+str(symbol)+"''  ; -Im("+str(symbol)+") / "+str(Analysis[1]),fontsize=font_size)
        plt.title(str(Title)+" Nyquist", fontsize=font_size)
        ax.yaxis.set_label_position("left")
        #Set limits on x- and y-range
        if len(x_lim_plot) == 2: 
            plt.xlim((x_lim_plot[0], x_lim_plot[1]))
        if len(y_lim_plot) == 2:
            plt.ylim((y_lim_plot[0], y_lim_plot[1]))
        
        #Display ticks as scientific if they are rounded to 0
        if round(max(Loaded_data[i].iloc[:,0])) == 0:
            ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        if round(max(Loaded_data[i].iloc[:,1])) == 0:
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
            
        #Save plot
        plt.savefig(os.path.join(savingFolder, 'Plot-' + filename[0] + '.pdf'),bbox_inches='tight')
        plt.show()
    except:
        sys.exit("Error in loading data. Please check you have selected an EIS/ECS excel/csv file and written the correct columns to be included.")
