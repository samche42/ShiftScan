#!/usr/bin/env python3

import scipy.optimize as optimization
from scipy.interpolate import splev, splrep
from scipy.optimize import OptimizeWarning
import warnings, sys, os
from scipy.signal import find_peaks
import numpy as np
import pandas as pd

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


def find_inflection(x,y):
    deriv1 = np.gradient(y,x) #Find first deriv
    deriv2 = list(np.gradient(deriv1,x)) #Find second deriv
    for i in range(len(deriv2) - 1): #Start looking for a sign change
        if deriv2[i] * deriv2[i+1] < 0:
            x_zero = x[i] - deriv2[i] * (x[i+1] - x[i]) / (deriv2[i+1] - deriv2[i]) #Find x-coord when sign flips = inflection
            return x_zero

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

#Function to model Boltzmann curve
def boltzmann_sigmoid(x, A, B, C, D):
    return A + (D - A) / (1 + np.exp(-B * (x - C)))

#Function calculating initial parameters for 4PL model
def initial_params(smoothed_y_coords, smoothed_x_coords):
	y_grad = np.gradient(smoothed_y_coords, smoothed_x_coords)
	max_grad_index = int(np.where((y_grad == np.nanmax(y_grad)))[0][0])
	inflection_point = smoothed_x_coords[max_grad_index]
	exact_inflection_point = find_inflection(smoothed_x_coords,smoothed_y_coords)
	infl_slope = np.nanmax(y_grad) #Find slope at inflection point
	low_asm = np.nanmin(smoothed_y_coords) #Get y coord of lower asymptotes
	high_asm = np.nanmax(smoothed_y_coords) #Get y coord of upper asymptote
	return low_asm,infl_slope,inflection_point,high_asm #return all four values for 4 parameter logistic regression model

#Function to return melting temp, residual sum of squares and y-coodinate list for each model
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
