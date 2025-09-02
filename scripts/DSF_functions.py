#!/usr/bin/env python3

import warnings, sys, os, gc, multiprocessing
from scipy import stats
from scipy.interpolate import splev, splrep
from scipy.optimize import OptimizeWarning
from scipy.signal import find_peaks
from functools import partial, reduce
import numpy as np
import pandas as pd
import scipy.optimize as optimization
from scipy.optimize import OptimizeWarning

def read_default_metadata(metadata_file, delimiter,dr_choice):
    print("Reading in metadata...")
    try:
        map_raw_df = pd.read_csv(metadata_file, sep=delimiter, header=0)
        if dr_choice == "yes":
            map_raw_df = map_raw_df.rename(columns={'ASSAY_PLATE': 'Assay_Plate', 'SOURCE_PLATE': 'Source_Plate', 'WELL': 'Well', 'COMPOUND': 'Compound', 'FRACTION': 'Fraction', 'REPLICATE':'Replicate','CONCENTRATION':'Concentration'})
        else:
            map_raw_df = map_raw_df.rename(columns={'ASSAY_PLATE': 'Assay_Plate', 'SOURCE_PLATE': 'Source_Plate', 'WELL': 'Well', 'COMPOUND': 'Compound', 'FRACTION': 'Fraction'})
        map_raw_df['Row'] = map_raw_df['Well'].str[0]
        map_raw_df['Column'] = map_raw_df['Well'].str[1:].str.replace(r'^(0+)', '', regex=True)
        map_raw_df['Unique_key'] = map_raw_df['Assay_Plate'] + '_' + map_raw_df['Row'] + map_raw_df['Column']
        map_raw_df['Column'] = map_raw_df['Column'].astype(int)
        map_raw_df.Fraction = map_raw_df.Fraction.fillna('NA')
        print("Metadata successfully loaded.")
        return map_raw_df
    except Exception as e:
        print(f"A problem was encountered while parsing metadata: {e}")
        sys.exit()

def roche_import(file_list,delim):
    if file_list: # Check if there are any files in the specified input folder
        dfs = []  # Initialize list of all files (converted to DataFrames)
        for file in file_list:
            try:
                data = pd.read_csv(file, sep=delim, header=0)
                data['Origin of data'] = file.stem  # Get filename without extension
                dfs.append(data)
            except Exception as e:
                print(f"There was a problem concatenating file: {file}. Please confirm the file is tab-delimited and has correct headers. Error: {e}")
                sys.exit()
    else:
        print("No .txt or .csv files found in specified input folder.")
        sys.exit()
    try:
        raw_input_df = pd.concat(dfs, ignore_index=True)#Smush them all together
        raw_input_df.columns = raw_input_df.columns.str.replace("X.", "X (", regex = True)
        raw_input_df.columns = [x+')' if 'X (' in x else x for x in raw_input_df.columns]
        parsed_df = raw_input_df.T.drop_duplicates().T #Remove duplicated temp columns
        parsed_df  = parsed_df.rename(columns={'X': 'Temp', 'Origin of data':'Origin'}) #Rename temp column
        parsed_df = parsed_df.drop(columns=[col for col in parsed_df.columns if 'X (' in col])
        melted_df = pd.melt(parsed_df,id_vars=['Temp','Origin']) #Melt dataframe
        melted_df['Well'] = melted_df.variable.str.split(':', expand=True)[0] #create well column
        melted_df  = melted_df.rename(columns={'value': 'Fluorescence'}) #Rename columns
        melted_df['Assay_Plate'] = melted_df.Origin.str.split('.', expand=True)[0] #Create plate column from origin data
        melted_df['Row'] = melted_df['Well'].str[0] #Take first character from Well string and make it Row value 
        melted_df['Column'] = melted_df['Well'].str[1:]
        melted_df = melted_df.astype({'Column':'int'})
        semifinal_df = melted_df.drop(['Origin','variable'],axis = 1) #Get rid of useless/redundant columns
        semifinal_df['Unique_key'] = semifinal_df['Assay_Plate']+'_'+semifinal_df['Well'] #Add in new column with generated unique key
    except Exception as e: 
        print("There was an issue parsing concatenated data. Please confirm input data format is correct. Error: ", e) 
        sys.exit()
    return semifinal_df

def roche_import_platewise(file,delim):
    try:
        data = pd.read_csv(file, sep=delim, header=0)
        data['Origin of data'] = file.stem  # Get filename without extension
    except Exception as e:
        print(f"There was a problem concatenating file: {file}. Please confirm the file is tab-delimited and has correct headers. Error: {e}")
        sys.exit()
    try:
        raw_input_df = data
        raw_input_df.columns = raw_input_df.columns.str.replace("X.", "X (", regex = True)
        raw_input_df.columns = [x+')' if 'X (' in x else x for x in raw_input_df.columns]
        parsed_df = raw_input_df.T.drop_duplicates().T #Remove duplicated temp columns
        parsed_df  = parsed_df.rename(columns={'X': 'Temp', 'Origin of data':'Origin'}) #Rename temp column
        parsed_df = parsed_df.drop(columns=[col for col in parsed_df.columns if 'X (' in col])
        melted_df = pd.melt(parsed_df,id_vars=['Temp','Origin']) #Melt dataframe
        melted_df['Well'] = melted_df.variable.str.split(':', expand=True)[0] #create well column
        melted_df  = melted_df.rename(columns={'value': 'Fluorescence'}) #Rename columns
        melted_df['Assay_Plate'] = melted_df.Origin.str.split('.', expand=True)[0] #Create plate column from origin data
        melted_df['Row'] = melted_df['Well'].str[0] #Take first character from Well string and make it Row value 
        melted_df['Column'] = melted_df['Well'].str[1:]
        melted_df = melted_df.astype({'Column':'int'})
        semifinal_df = melted_df.drop(['Origin','variable'],axis = 1) #Get rid of useless/redundant columns
        semifinal_df['Unique_key'] = semifinal_df['Assay_Plate']+'_'+semifinal_df['Well'] #Add in new column with generated unique key
    except Exception as e: 
        print("There was an issue parsing concatenated data. Please confirm input data format is correct. Error: ", e) 
        sys.exit()
    return semifinal_df

def biorad_import(file_list,delim):
    if file_list: # Check if there are any files in the specified input folder
        dfs = []  # Initialize list of all files (converted to DataFrames)
        for file in file_list:
            try:
                data = pd.read_csv(file, sep=delim,skiprows=4, header=0) #Skip the first 4 rows
                data = data.drop(['Step','Cycle','Dye'],axis = 1) #Drop the first 3 columns (Step, Cycle, Dye)
                data['Origin of data'] = file.stem # Get filename without extension and add as new column
                data = data.dropna(subset=['Temp.'])
                dfs.append(data)
            except Exception as e:
                print(f"There was a problem concatenating file: {file}. Please confirm the file is tab-delimited and has correct headers. Error: {e}")
                sys.exit()
    else:
        print("No .txt or .csv files found in specified input folder.")
        sys.exit()
    try:
        raw_input_df = pd.concat(dfs, ignore_index=True)#Smush them all together
        raw_input_df  = raw_input_df.rename(columns={'Temp.': 'Temp', 'Origin of data':'Origin'}) #Rename temp column
        melted_df = pd.melt(raw_input_df,id_vars=['Temp','Origin']) #Melt dataframe
        melted_df = melted_df.rename(columns={'variable': 'Well','value': 'Fluorescence'})
        melted_df['Assay_Plate'] = melted_df.Origin.str.split('.', expand=True)[0] #Create plate column from origin data
        melted_df['Row'] = melted_df['Well'].str[0] #Take first character from Well string and make it Row value 
        melted_df['Column'] = melted_df['Well'].str[1:]
        melted_df = melted_df.astype({'Column':'int'})
        semifinal_df = melted_df.drop(['Origin'],axis = 1) #Get rid of useless/redundant columns
        semifinal_df['Unique_key'] = semifinal_df['Assay_Plate']+'_'+semifinal_df['Well'] #Add in new column with generated unique key
    except Exception as e: 
        print("There was an issue parsing concatenated data. Please confirm input data format is correct. Error: ", e) 
        sys.exit()
    return semifinal_df

def biorad_import_platewise(file,delim):
    try:
        data = pd.read_csv(file, sep=delim,skiprows=4, header=0) #Skip the first 4 rows
        data = data.drop(['Step','Cycle','Dye'],axis = 1) #Drop the first 3 columns (Step, Cycle, Dye)
        data['Origin of data'] = file.stem # Get filename without extension and add as new column
        data = data.dropna(subset=['Temp.'])
        dfs.append(data)
    except Exception as e:
        print(f"There was a problem concatenating file: {file}. Please confirm the file is tab-delimited and has correct headers. Error: {e}")
        sys.exit()
    try:
        raw_input_df = pd.concat(dfs, ignore_index=True)#Smush them all together
        raw_input_df  = raw_input_df.rename(columns={'Temp.': 'Temp', 'Origin of data':'Origin'}) #Rename temp column
        melted_df = pd.melt(raw_input_df,id_vars=['Temp','Origin']) #Melt dataframe
        melted_df = melted_df.rename(columns={'variable': 'Well','value': 'Fluorescence'})
        melted_df['Assay_Plate'] = melted_df.Origin.str.split('.', expand=True)[0] #Create plate column from origin data
        melted_df['Row'] = melted_df['Well'].str[0] #Take first character from Well string and make it Row value 
        melted_df['Column'] = melted_df['Well'].str[1:]
        melted_df = melted_df.astype({'Column':'int'})
        semifinal_df = melted_df.drop(['Origin'],axis = 1) #Get rid of useless/redundant columns
        semifinal_df['Unique_key'] = semifinal_df['Assay_Plate']+'_'+semifinal_df['Well'] #Add in new column with generated unique key
    except Exception as e: 
        print("There was an issue parsing concatenated data. Please confirm input data format is correct. Error: ", e) 
        sys.exit()
    return semifinal_df

def concatenate_files(files, output_dir_string, output_filename): 
    dfs = [] 
    for file in files: 
        try: 
            file_path = os.path.join(output_dir_string, file) 
            data = pd.read_csv(file_path, sep='\t', header=0) 
            dfs.append(data) 
        except Exception as e: 
            print(f"There was a problem reading the file: {file}. Error: {e}") 
            continue 
    try: 
        dataframe = pd.concat(dfs, ignore_index=True) 
        output_path = os.path.join(output_dir_string, output_filename) 
        dataframe.to_csv(output_path, sep="\t", index=False) 
    except Exception as e: 
        print(f"There was a problem writing the output file: {output_filename}. Error: {e}")

def find_inflection(x,y):
    deriv1 = np.gradient(y,x) #Find first deriv
    deriv2 = list(np.gradient(deriv1,x)) #Find second deriv
    for i in range(len(deriv2) - 1): #Start looking for a sign change
        if deriv2[i] * deriv2[i+1] < 0:
            x_zero = x[i] - deriv2[i] * (x[i+1] - x[i]) / (deriv2[i+1] - deriv2[i]) #Find x-coord when sign flips = inflection
            return x_zero

def slice_list(positions,input_list):
    returned_list = []
    for i in range(len(positions)-1):
        start = positions[i]
        end = positions[i + 1]
        returned_list.append(input_list[start:end]) #Generate new list of coord lists for each slice position pair
    return returned_list

def clean_curve(x,spl):
    pred_y = splev(x, spl)
    pred_x = x
    minima, _n = find_peaks(-pred_y) #Find local minima of curve
    maxima, _x = find_peaks(pred_y, height=0) #Find local maxima of curve
    index_positions = np.concatenate((minima, maxima)) #Create an array of local minima and maxima
    index_positions = np.sort(np.append(index_positions, [0,len(pred_y)])) #Add in 0 and last position of coordinates
    x_slices = slice_list(index_positions, pred_x) #Slice x-coordinates at minima and maxima
    y_slices = slice_list(index_positions, pred_y) #Slice y-coordinates at minima and maxima
    x_slices = [sublist for sublist in x_slices if len(sublist) >= 5] #Remove any coord list with less than 5 data points
    y_slices = [sublist for sublist in y_slices if len(sublist) >= 5] #Remove any coord list with less than 5 data points
    clean_x_list = []
    clean_y_list = []
    for i in range(len(x_slices)):
        average_slope = np.nanmean(np.gradient(y_slices[i], x_slices[i])) #Ignore Nan values when calculating mean gradient
        #Next, let's calculate the gradient of the second derivative slop for each coord set. 
        #We would expect a mix of positive and negative values if the curve is actually sigmoidal (i.e first deriv curve actually has a peak of sorts)
        first_deriv = list(np.gradient(y_slices[i], x_slices[i]))
        second_deriv = list(np.gradient(first_deriv, x_slices[i]))
        neg_count = len(list(filter(lambda x: (x < 0), second_deriv)))
        pos_count = len(list(filter(lambda x: (x > 0), second_deriv)))
        total_count = len(second_deriv)
        neg_perc = neg_count/total_count*100
        pos_perc = pos_count/total_count*100
        if average_slope < 0:
            continue #If average slope is negative do not add it to final list
        elif (neg_perc < 10)|(pos_perc <10): #If less than 10% of values of the 2nd deriv gradient are either positive or negative, the segment isn't sigmoidal and we shouldn't save it as a legit curve
            continue
        else:
            clean_x_list.append(x_slices[i])
            clean_y_list.append(y_slices[i])
    return clean_x_list, clean_y_list

def split_curves(list_of_x_coords, list_of_y_coords):
    new_list_of_x_coords = []
    new_list_of_y_coords = []
    for i in range(len(list_of_x_coords)): #For each set of coordinates in the list...
        first_deriv_y_coordinates = np.gradient(list_of_x_coords[i], list_of_y_coords[i]) #...Find the first derivative
        maxima, _x = find_peaks(first_deriv_y_coordinates) #Find local maxima of first derivative curve
        maxima = np.sort(np.append(maxima, [0,len(first_deriv_y_coordinates)])) #Add in 0 and last position of coordinates
        if len(maxima) == 2: #If there is no maxima, we have a single curve and can return the original coordinate list
            new_list_of_x_coords.append(list_of_x_coords[i])
            new_list_of_y_coords.append(list_of_y_coords[i])
        else: #If there is a maxima, we must split the graph
            sliced_x = slice_list(maxima, list_of_x_coords[i])
            sliced_y = slice_list(maxima, list_of_y_coords[i])
            sliced_x = [sublist for sublist in sliced_x if len(sublist) >= 5]
            sliced_y = [sublist for sublist in sliced_y if len(sublist) >= 5]
            new_list_of_x_coords.extend(sliced_x)
            new_list_of_y_coords.extend(sliced_y)
    return new_list_of_x_coords,new_list_of_y_coords

def add_curve_data(unique_identifier, sub_plot, temps_list, smooth_y_list, boltzmann_y_list):
    keys = [unique_identifier for _ in range(len(temps_list))]
    sub_plots = [sub_plot for _ in range(len(temps_list))]
    new_row = {
        "Unique_key": keys, 
        "Subplot": sub_plots,
        "Temps": temps_list,
        "Smooth Fluorescence": smooth_y_list,
        "Boltzmann": boltzmann_y_list
    }
    return new_row

def boltzmann_sigmoid(x, A, B, C, D):
    return A + (D - A) / (1 + np.exp(-B * (x - C)))

def initial_params(smoothed_y_coords, smoothed_x_coords):
    y_grad = np.gradient(smoothed_y_coords, smoothed_x_coords)
    max_grad_index = int(np.where((y_grad == np.nanmax(y_grad)))[0][0])
    inflection_point = smoothed_x_coords[max_grad_index]
    exact_inflection_point = find_inflection(smoothed_x_coords,smoothed_y_coords)
    infl_slope = np.nanmax(y_grad) #Find slope at inflection point
    low_asm = np.nanmin(smoothed_y_coords) #Get y coord of lower asymptotes
    high_asm = np.nanmax(smoothed_y_coords) #Get y coord of upper asymptote
    return low_asm,infl_slope,inflection_point,high_asm #return all four values for 4 parameter logistic regression model

def Model_data(model, predx, predy, init_params, maxfevopt):
    try:
        popt, pcov = optimization.curve_fit(model, predx, predy, init_params, maxfev = maxfevopt,method='dogbox') #Run scipy curve_fit to fit data to model curve. Dogbox: https://stackoverflow.com/questions/55725139/fit-sigmoid-function-s-shape-curve-to-data-using-python
        model_y_list_loop = model(predx, *popt) #Get model y-coordinates
        residuals = predy - model_y_list_loop #Calculate residuals (difference between data points and model)
        RSS = np.sum(residuals**2) #Get residual sum of squares
        MSE = np.mean(residuals**2) #Get residual mean of squares
        y_grad = np.gradient(model_y_list_loop, predx)
        Tm = predx[int(np.where((y_grad==np.nanmax(y_grad)))[0][0])]
        message = ''
    except TypeError:       
        try:
            popt, pcov = optimization.curve_fit(model, predx, predy, init_params, maxfev = maxfevopt,method='lm') #Try alternate method which can deal with slightly more misbehaved curves
            model_y_list_loop = model(predx, *popt) #Get model y-coordinates
            residuals = predy - model_y_list_loop #Calculate residuals (difference between data points and model)
            RSS = np.sum(residuals**2) #Get residual sum of squares
            MSE = np.mean(residuals**2) #Get residual mean of squares
            y_grad = np.gradient(model_y_list_loop, predx)
            Tm = predx[int(np.where((y_grad==np.nanmax(y_grad)))[0][0])]
            message = ''
        except Exception:
            raise
    except RuntimeError as errorstring:
        if "Optimal parameters not found" in str(errorstring):
            try:
                popt, pcov = optimization.curve_fit(model, predx, predy, init_params, maxfev = maxfevopt,method='lm') #Try alternate method which can deal with slightly more misbehaved curves
                model_y_list_loop = model(predx, *popt) #Get model y-coordinates
                residuals = predy - model_y_list_loop #Calculate residuals (difference between data points and model)
                RSS = np.sum(residuals**2) #Get residual sum of squares
                MSE = np.mean(residuals**2) #Get residual mean of squares
                y_grad = np.gradient(model_y_list_loop, predx)
                Tm = predx[int(np.where((y_grad==np.nanmax(y_grad)))[0][0])]
                message = ''
            except RuntimeError as errorstring2:
                if "Optimal parameters not found" in str(errorstring2):
                    model_y_list_loop = np.array([None for _ in range(len(predx))])
                    Tm = None
                    MSE = None
                    RSS = None
                    message = "Parameter optimization failed for both methods"
                else:
                    raise
            except TypeError:
                    model_y_list_loop = np.array([None for _ in range(len(predx))])
                    Tm = None
                    MSE = None
                    RSS = None
                    message = "Modelling failed for both methods"
        else:
            raise
    return Tm, model_y_list_loop, MSE,RSS, message

def process_well(sub_df,smoothing_factor,normalize):
    original_curve_rows = []
    sub_curve_rows = []
    temps = []
    smooth_fluorescence = []
    boltzmann_y_list = []
    new_Tm_rows = []
    key = sub_df['Unique_key'].unique()[0]
    if sub_df.Well_type.unique()[0] == 'Blank': #Skip the loop if the well is classified as blank. Total waste of time to analyze rubbish data
        new_Tm_rows = [{'Assay_Plate': sub_df.Assay_Plate.unique()[0], 'Source_Plate': sub_df.Source_Plate.unique()[0],'Well': sub_df.Well.unique()[0],'Unique_key': sub_df.Unique_key.unique()[0],'Subplot':None,'Well_type': sub_df.Well_type.unique()[0], 'Compound': None,'Fraction': None,'Smooth_Tm': None,'Boltzmann_Tm':None,'Boltzmann_RSS': None,'Amplitude':None,'Curve_height': None, 'Error':'Blank well','Warning':''}]
        return original_curve_rows, sub_curve_rows, new_Tm_rows
    else:
        x = list(sub_df["Temp"].values.tolist()) #Convert Temp and fluorescence columns to lists of values
        if normalize == "y":
            y = list(sub_df["Normalized_Fluorescence"].values.tolist())
        elif normalize == "n":
            y = list(sub_df["Fluorescence"].values.tolist())
        else:
            print("Unknown argument for normalization provided. Argument for normalization can only be y or n.")
        try:
            spl = splrep(x,y, s = smoothing_factor) #Smooth data. This used to be 0.05 when working with raw data. This changes drastically when dealing with normalized data
        except RuntimeWarning:
            list_length = len(x[0::4])
            final_message = "Failed: Smoothing failed"   
            temps.append([None for _ in range(list_length)])
            smooth_fluorescence.append([None for _ in range(list_length)])
            boltzmann_y_list = [None for _ in range(len(x[0::4]))]
            original_curve_rows.append(add_curve_data(key, None,temps[0], smooth_fluorescence[0], boltzmann_y_list))
            new_Tm_rows = [{'Assay_Plate': sub_df.Assay_Plate.unique()[0], 'Source_Plate': sub_df.Source_Plate.unique()[0],'Well': sub_df.Well.unique()[0],'Unique_key': sub_df.Unique_key.unique()[0],'Subplot':None,'Well_type': sub_df.Well_type.unique()[0], 'Compound': sub_df.Compound.unique()[0],'Fraction': sub_df.Fraction.unique()[0],'Smooth_Tm': None,'Boltzmann_Tm':None,'Boltzmann_RSS': None,'Amplitude':None,'Curve_height': None, 'Error':final_message,'Warning':''}] 
            return original_curve_rows, sub_curve_rows, new_Tm_rows
        #Store original curve coordinates of smoothed data for well (i.e. curve isn't split yet)
        temps.append(x[0::4])
        smooth_fluorescence.append(list((splev(x, spl))[0::4]))
        #Check for multiple peaks in first deriv graph
        list_of_x_coords, list_of_y_coords = clean_curve(x,spl) # return list of lists where each index pair is a new "sub-curve"
        curve_count = len(list_of_x_coords)
        boltzmann_y_list = [None for _ in range(len(x[0::4]))]
        original_curve_rows.append(add_curve_data(key, None,temps[0], smooth_fluorescence[0],boltzmann_y_list))
        if curve_count == 0:
            final_message = "Well failed: No region of positive slope or shape not sigmoidal"
            new_Tm_rows = [{'Assay_Plate': sub_df.Assay_Plate.unique()[0], 'Source_Plate': sub_df.Source_Plate.unique()[0],'Well': sub_df.Well.unique()[0],'Unique_key': sub_df.Unique_key.unique()[0],'Subplot':None,'Well_type': sub_df.Well_type.unique()[0], 'Compound': sub_df.Compound.unique()[0],'Fraction': sub_df.Fraction.unique()[0],'Smooth_Tm': None,'Boltzmann_Tm':None,'Boltzmann_RSS': None,'Amplitude':None,'Curve_height': None, 'Error':final_message,'Warning':''}]
        elif curve_count > 4:
            final_message = "Well failed: Too many subcurves identified - messy curve assumed"
            new_Tm_rows = [{'Assay_Plate': sub_df.Assay_Plate.unique()[0], 'Source_Plate': sub_df.Source_Plate.unique()[0],'Well': sub_df.Well.unique()[0],'Unique_key': sub_df.Unique_key.unique()[0],'Subplot':None,'Well_type': sub_df.Well_type.unique()[0], 'Compound': sub_df.Compound.unique()[0],'Fraction': sub_df.Fraction.unique()[0],'Smooth_Tm': None,'Boltzmann_Tm':None,'Boltzmann_RSS': None,'Amplitude':None,'Curve_height': None, 'Error':final_message,'Warning':''}]
        else:
            sub_plot_count = 0
            final_list_of_x_coords, final_list_of_y_coords = split_curves(list_of_x_coords, list_of_y_coords) #Submit clean curve for second round of peak identification
            final_list_of_y_coords = [list(arr) for arr in final_list_of_y_coords]
            for k in range(len(final_list_of_x_coords)):
                sub_plot_count +=1 #Get a count for the subplots
                p0 = initial_params(final_list_of_y_coords[k],final_list_of_x_coords[k]) #Find initial parameters for model
                maxfev_opt = 100*(len(final_list_of_y_coords[k])+1) #Optimize maxfev parameter for each dataset https://www.physics.utoronto.ca/~phy326/python/curve_fit_to_data.py
                exact_inflection_point = find_inflection(final_list_of_x_coords[k],final_list_of_y_coords[k])
                if exact_inflection_point is None:
                    boltzmann_y_coordinates = [None for _ in range(len(final_list_of_x_coords[k]))]
                    sub_curve_rows.append(add_curve_data(key, sub_plot_count, final_list_of_x_coords[k], final_list_of_y_coords[k],boltzmann_y_coordinates))
                    new_Tm_rows.extend([{'Assay_Plate': sub_df.Assay_Plate.unique()[0], 'Source_Plate': sub_df.Source_Plate.unique()[0],'Well': sub_df.Well.unique()[0],'Unique_key': sub_df.Unique_key.unique()[0],'Subplot':sub_plot_count,'Well_type': sub_df.Well_type.unique()[0], 'Compound': sub_df.Compound.unique()[0],'Fraction': sub_df.Fraction.unique()[0],'Smooth_Tm': exact_inflection_point,'Boltzmann_Tm':None,'Boltzmann_RSS': None,'Amplitude':p0[3] - p0[0],'Curve_height': p0[3], 'Error':'Failed: Inflection point could not be determined','Warning':''}])
                else:
                    boltzmann_Tm, boltzmann_y_coordinates, MSE, RSS, final_message = Model_data(boltzmann_sigmoid, final_list_of_x_coords[k],final_list_of_y_coords[k], p0, maxfev_opt) #Fit Boltzmann
                    boltzmann_y_coordinates = (boltzmann_y_coordinates.tolist())
                    sub_curve_rows.append(add_curve_data(key, sub_plot_count, final_list_of_x_coords[k], final_list_of_y_coords[k],boltzmann_y_coordinates))
                    new_Tm_rows.extend([{'Assay_Plate': sub_df.Assay_Plate.unique()[0], 'Source_Plate': sub_df.Source_Plate.unique()[0],'Well': sub_df.Well.unique()[0],'Unique_key': sub_df.Unique_key.unique()[0],'Subplot':sub_plot_count,'Well_type': sub_df.Well_type.unique()[0], 'Compound': sub_df.Compound.unique()[0],'Fraction': sub_df.Fraction.unique()[0],'Smooth_Tm': exact_inflection_point,'Boltzmann_Tm':boltzmann_Tm,'Boltzmann_RSS': RSS,'Amplitude':p0[3] - p0[0],'Curve_height': p0[3], 'Error':final_message,'Warning':''}])
        return original_curve_rows, sub_curve_rows, new_Tm_rows

def default_analysis(final_df,control_cols,delimiter,normalization,smoothing_factor,processors,only_tm,dose_response,failed_control_wells,ctrl_tm_cutoff,ctrl_amp_cutoff,min_amp_cutoff,max_amp_cutoff):
    # Assigning Well Types and Normalization
    print("Assigning control columns and well types...")
    if control_cols.endswith((".txt",".tab")):
        well_mapping_df = pd.read_csv(control_cols, sep=delimiter, header=0)
        well_mapping_df['Row'] = well_mapping_df['Well'].str[0]
        well_mapping_df['Column'] = well_mapping_df['Well'].str[1:].str.replace(r'^(0+)', '', regex=True)
        well_mapping_df['Unique_key'] = well_mapping_df['Assay_Plate'] + '_' + well_mapping_df['Row'] + well_mapping_df['Column']
        well_mapping_df = well_mapping_df.drop(['Assay_Plate','Well'], axis=1)
        final_df = pd.merge(final_df, well_mapping_df, on='Unique_key', how='left')
    else:
        control_list = [int(control) for control in control_cols.split(",")]
        if not control_list:
            print("ERROR: No control columns specified")
            sys.exit()
        is_control = final_df['Column'].isin(control_list)
        is_blank = final_df['Compound'].isnull()
        final_df['Well_type'] = np.select([is_control, is_blank], ['Control', 'Blank'], default='Experimental')
    
    # Normalize Fluorescence data
    if normalization == "y":
        print("Normalizing data...")
        final_df['Normalized_Fluorescence'] = final_df.groupby('Assay_Plate')['Fluorescence'].transform(lambda s: s / s.max())
    elif normalization == "n":
        print("Skipping normalization.")
        final_df['Normalized_Fluorescence'] = final_df.groupby('Assay_Plate')['Fluorescence'].transform(lambda s: s / s.max())
    else:
        print("ERROR: Invalid normalization option. Please select either y or n.")
        sys.exit()

    # PART 4: Processing wells using multiprocessing
    print("Processing wells via multiprocessing...")
    list_of_sub_dfs = [group.reset_index(drop=True) for _, group in final_df.groupby('Unique_key')]
    
    all_original_curve_rows = []
    all_subcurves = []
    all_tm_rows = []
    
    with multiprocessing.Pool(processes=int(processors)) as pool:
        partial_process_wrapper = partial(process_well, smoothing_factor=smoothing_factor, normalize=normalization)
        results = pool.map(partial_process_wrapper, list_of_sub_dfs)
        for original_curve_rows, sub_curve_rows, Tm_rows in results:
            all_original_curve_rows.extend(original_curve_rows)
            all_subcurves.extend(sub_curve_rows)
            all_tm_rows.extend(Tm_rows)
        del results
        gc.collect()
    
    print("Well processing complete.")

    # Consolidate results
    original_curve_df = pd.concat([pd.DataFrame(d) for d in all_original_curve_rows], ignore_index=True)
    all_subcurves_df = pd.concat([pd.DataFrame(d) for d in all_subcurves], ignore_index=True)

    semifinal_curves = pd.concat([original_curve_df, all_subcurves_df], ignore_index=True)
    semifinal_curves['Subplot'].fillna('Original', inplace=True)
    semifinal_curves['Subplot'] = semifinal_curves['Subplot'].astype(str)
    semifinal_curves[['Assay_Plate', 'Well']] = semifinal_curves['Unique_key'].str.rsplit('_', 1, expand=True)
    
    Tm_df = pd.DataFrame(all_tm_rows)
    print(f"Tm DataFrame created with {len(Tm_df)} rows.")
    
    if only_tm:
        return {'Tm_df': Tm_df, 'curves_df': semifinal_curves, 'plate_report_df': None, 'problems_df': None}

    # PART 5 & 6: Flagging bad results and calculating final values
    print("Flagging bad wells and calculating final results...")
    
    control_z_scores = Tm_df[Tm_df['Well_type'] == 'Control'].groupby('Assay_Plate', as_index=False)[['Smooth_Tm', 'Amplitude']].transform(stats.zscore).abs()
    control_z_scores = control_z_scores.rename(columns={'Smooth_Tm': 'Ctrl_Tm_z-score', 'Amplitude': 'Ctrl_Amplitude_z-score'}).fillna(0)
    Tm_df = Tm_df.join(control_z_scores)

    Tm_df['Error'] = Tm_df['Error'].fillna('')
    Tm_df['Error'] = np.where(Tm_df['Ctrl_Tm_z-score'] > ctrl_tm_cutoff, Tm_df['Error'] + f"Failed: Tm > {ctrl_tm_cutoff} std away from mean. ", Tm_df['Error'])
    Tm_df['Error'] = np.where(Tm_df['Ctrl_Amplitude_z-score'] > ctrl_amp_cutoff, Tm_df['Error'] + f"Failed: Amplitude > {ctrl_amp_cutoff} std away from mean. ", Tm_df['Error'])

    #Find differences between smoothed and modeled data
    Tm_df['Tm_SE'] = Tm_df[['Smooth_Tm','Boltzmann_Tm']].sem(axis=1, skipna = False) #calculate std error between smoothed data and Boltzmann
    Tm_df['Tm_difference'] = (Tm_df['Boltzmann_Tm']-Tm_df['Smooth_Tm']).abs() #Find absolute difference between smoothed and boltzmann Tm

    avg_ctrl_amplitude_df = Tm_df[(Tm_df['Well_type'] == 'Control') & (~Tm_df['Error'].str.contains('ailed'))].groupby('Assay_Plate', as_index=False)['Amplitude'].mean()
    avg_ctrl_amplitude_df = avg_ctrl_amplitude_df.rename(columns={'Amplitude': 'Ctrl_avg_amplitude'})
    Tm_df = pd.merge(Tm_df, avg_ctrl_amplitude_df, on='Assay_Plate', how='left')

    Tm_df['Amplitude_lower_cutoff'] = Tm_df['Ctrl_avg_amplitude'] * min_amp_cutoff
    Tm_df['Error'] = np.where(Tm_df['Amplitude'] < Tm_df['Amplitude_lower_cutoff'], Tm_df['Error'] + " Failed: Amplitude too small. ", Tm_df['Error'])
    Tm_df['Amplitude_upper_cutoff'] = Tm_df['Ctrl_avg_amplitude'] * max_amp_cutoff
    Tm_df['Error'] = np.where(Tm_df['Amplitude'] > Tm_df['Amplitude_upper_cutoff'], Tm_df['Error'] + " Failed: Amplitude too large. ", Tm_df['Error'])

    def final_decision(group):
        if group['Error'].eq('').any():
            group['Final_decision'] = np.where(group['Error'] == '', 'Pass', 'Removed')
        else:
            group['Final_decision'] = 'Failed'
        return group
    
    grouped = Tm_df.groupby('Unique_key')
    Tm_df = pd.concat([final_decision(group) for _, group in grouped], ignore_index=True)
    Tm_df['Final_decision'] = np.where(Tm_df['Well_type'] == 'Blank', 'Removed', Tm_df['Final_decision'])
    Tm_df['Final_Tm'] = np.where(Tm_df['Final_decision'] == 'Pass', Tm_df['Smooth_Tm'], np.NaN)

    well_type_df = Tm_df.loc[:, ['Unique_key', 'Well_type']]
    decision_df = Tm_df.loc[:, ['Unique_key', 'Subplot', 'Final_decision', 'Error', 'Ctrl_Tm_z-score']]
    decision_df['Subplot'] = decision_df['Subplot'].astype(pd.Int64Dtype()).astype(str).replace('<NA>', 'Original')
    
    final_curves = pd.merge(semifinal_curves, well_type_df.drop_duplicates(), on='Unique_key', how='inner')
    final_curves = pd.merge(final_curves, Tm_df[['Unique_key', 'Source_Plate']].drop_duplicates(), on='Unique_key', how='left')
    final_curves = pd.merge(final_curves, decision_df, on=['Unique_key', 'Subplot'], how='left')

    if dose_response:
        dr_metadata = final_df[['Unique_key','Concentration','Replicate']].drop_duplicates()
        Tm_df = pd.merge(Tm_df, dr_metadata, on='Unique_key', how='inner')
        Tm_df['Unique_key_subplot'] = Tm_df['Unique_key'] + "_" + Tm_df['Subplot'].astype(str)
        Tm_df['Relative_amplitude'] = Tm_df['Amplitude']/Tm_df['Ctrl_avg_amplitude']
        return {'Tm_df': Tm_df, 'curves_df': final_curves, 'plate_report_df': None, 'problems_df': None}
    else: 
        #Get averages for control wells (DOES NOT INCLUDE FAILED WELLS)
        ctrl_df = Tm_df[Tm_df['Well_type'] == 'Control'].groupby(['Assay_Plate'],as_index=False)[['Final_Tm']].mean() #Get average ctrl Tm
        ctrl_df = ctrl_df.rename(columns={'Final_Tm': 'Avg_ctrl_melting_temp'})
        Tm_df = pd.merge(Tm_df, ctrl_df, on = 'Assay_Plate') #Merge that info in for easy calcs

        Tm_df = Tm_df.sort_values(by=["Unique_key", "Subplot"])
        Tm_df = Tm_df.reset_index(drop = True)

        Tm_df['Row'] = Tm_df['Well'].str[0] #Generate Row column (needed to build visualizations)
        Tm_df['Column'] = (Tm_df['Well'].str[1:]).astype(int) #Generate Column column (needed to build visualizations)
        Tm_df['Relative_amplitude'] = Tm_df['Amplitude']/Tm_df['Ctrl_avg_amplitude']

        #Calculating final z-scores per well
        plate_zscores = Tm_df[Tm_df['Final_Tm'].notna()].groupby(['Assay_Plate'],as_index=False)['Final_Tm'].transform(stats.zscore) #Find z-score for all wells regardless of well type
        plate_zscores.rename(columns={'Final_Tm': 'Well_zscore'}, inplace=True)
        Tm_df = Tm_df.join(plate_zscores) #merge on index

        #Finding max ctrl z-score per plate
        ctrl_max_zscore = Tm_df[Tm_df["Well_type"] == 'Control'].groupby(['Assay_Plate'],as_index=False)['Well_zscore'].max()
        ctrl_max_zscore.rename(columns={'Well_zscore': 'Max_ctrl_zscore_for_plate'}, inplace=True)
        Tm_df = pd.merge(Tm_df, ctrl_max_zscore, on = 'Assay_Plate')

        #Finding min ctrl z-score per plate
        ctrl_min_zscore = Tm_df[Tm_df["Well_type"] == 'Control'].groupby(['Assay_Plate'],as_index=False)['Well_zscore'].min()
        ctrl_min_zscore.rename(columns={'Well_zscore': 'Min_ctrl_zscore_for_plate'}, inplace=True)
        Tm_df = pd.merge(Tm_df, ctrl_min_zscore, on = 'Assay_Plate')

        #Finding difference between well and ctrl average
        Tm_df['Diff from ctrl avg'] = Tm_df['Final_Tm'] - Tm_df['Avg_ctrl_melting_temp']

        #Finding no. std devs (z-score) - i.e. the traditional route for finding hits
        plate_std_dev = Tm_df[Tm_df['Final_Tm'].notna()].groupby('Assay_Plate')['Final_Tm'].std().reset_index()
        plate_std_dev.rename(columns={'Final_Tm': 'Plate_Std_Dev'}, inplace=True)
        Tm_df = pd.merge(Tm_df, plate_std_dev, on = 'Assay_Plate')
        Tm_df['Std. devs from ctrl mean'] = (Tm_df['Final_Tm']- Tm_df['Avg_ctrl_melting_temp'])/Tm_df['Plate_Std_Dev']
        
        print("Calculations complete!")

        # Generating reports
        from functools import reduce

        print("Generating reports...")

        failed_df_raw = Tm_df[Tm_df['Final_decision'] == 'Failed']
        failed_df = failed_df_raw.drop_duplicates(subset='Unique_key', keep='first')
        all_plates = Tm_df['Assay_Plate'].unique()

        #Count total failed controls per assay plate
        total_failures_series = failed_df[(failed_df['Well_type'] == 'Control')].groupby(['Assay_Plate']).size()
        total_failures_series = total_failures_series.reindex(all_plates, fill_value = 0)
        total_failures_df = total_failures_series.reset_index(name='Total_failures')
        total_failures_df['Plate_status'] = np.where(total_failures_df['Total_failures'] > failed_control_wells, "Failed","OK")

        #Count how many of the errors were down to the Tm being outside the desired z-score range
        melting_temp_failures_series = failed_df[(failed_df['Well_type'] == 'Control')&(failed_df['Error'].str.contains('Tm >'))].groupby('Assay_Plate').size()
        melting_temp_failures_series = melting_temp_failures_series.reindex(all_plates, fill_value = 0)
        melting_temp_failures_df = melting_temp_failures_series.reset_index(name='Ctrl_Tm_Z-score_failures')

        #Count ctrl wells that failed modelling
        model_failure_series = failed_df[(failed_df['Well_type'] == 'Control')&((failed_df['Error'].str.contains('No region of positive slope')) | \
                                                                                (failed_df['Error'].str.contains('Parameter optimization')) | \
                                                                                (failed_df['Error'].str.contains('Modelling failed for both methods')))].groupby('Assay_Plate').size()
        model_failure_series = model_failure_series.reindex(all_plates, fill_value = 0)
        model_failure_df = model_failure_series.reset_index(name='Modelling_failures')

        #Count ctrl wells that failed amplitude restrictions
        amp_failures_series = failed_df[(failed_df['Well_type'] == 'Control')&(failed_df['Error'].str.contains('Failed: Amplitude'))].groupby('Assay_Plate').size()
        amp_failures_series = amp_failures_series.reindex(all_plates, fill_value = 0)
        amp_failures_df = amp_failures_series.reset_index(name='Amp_Z-score_failures')

        #List out failed wells
        failed_wells = failed_df[failed_df['Well_type'] == 'Control'][failed_df[failed_df['Well_type'] == 'Control']['Final_decision'].str.contains('Failed', case=False)].groupby(['Assay_Plate'])['Well'].agg(lambda x: ','.join(x)).reset_index()
        failed_wells= failed_wells.rename(columns={'Well':'Failed ctrl wells'})
        
        #Put it all together
        plate_status_df = Tm_df[['Assay_Plate','Source_Plate']].drop_duplicates()
        list_of_dfs = [plate_status_df,total_failures_df,failed_wells,melting_temp_failures_df,model_failure_df,amp_failures_df]
        plate_report = reduce(lambda  left,right: pd.merge(left,right,on=['Assay_Plate'],how='outer'), list_of_dfs)

        plate_status_only = plate_report.loc[:, ['Assay_Plate', 'Plate_status']]
        Tm_df = pd.merge(Tm_df, plate_status_only, on = 'Assay_Plate') #merge plate status back in for visualization purposes
        Tm_df['Platename'] = np.where(Tm_df['Plate_status'] =='Failed', Tm_df['Source_Plate']+" ("+Tm_df['Plate_status']+")", Tm_df['Source_Plate']) #Add "Failed" to plate name to make it obvious in visulizations

        #tmp = failed_df[failed_df['Well_type'] == 'Control'].groupby(['Well'],as_index=False)['Error'].apply(lambda x: x[x.str.contains('Failed')].count())
        tmp = (failed_df[failed_df['Well_type'] == 'Control'].groupby(['Well'], as_index=False).apply(lambda x: x['Error'].str.contains('Failed').sum()).rename(columns={None: 'FailedCount'}))
        if tmp.empty:
            well_error_count_df = pd.DataFrame(columns=['Well','Wells that repeatedly failed'])
        else:
            well_error_count_df = tmp[tmp['FailedCount'] >= 3]
            well_error_count_df = well_error_count_df.rename(columns={'FailedCount':'Wells that repeatedly failed'})

        Tm_df['Unique_key_subplot'] = Tm_df['Unique_key']+"_"+Tm_df['Subplot'].astype(str)

        print("Analysis complete!")

        return {'Tm_df': Tm_df, 'curves_df': final_curves, 'plate_report_df': plate_report, 'problems_df': well_error_count_df}
