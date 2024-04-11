import scipy.optimize as optimization
from scipy.interpolate import splev, splrep
from scipy.optimize import OptimizeWarning
import warnings
from scipy.signal import find_peaks
import numpy as np

def slice_list(positions,input_list):
	returned_list = []
	for i in range(len(positions)-1):
		start = positions[i]
		end = positions[i + 1]
		returned_list.append(input_list[start:end])
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
	x_slices = [sublist for sublist in x_slices if len(sublist) >= 20] #Remove any coord list with less than 20 data points
	y_slices = [sublist for sublist in y_slices if len(sublist) >= 20] #Remove any coord list with less than 20 data points
	clean_x_list = []
	clean_y_list = []
	for i in range(len(x_slices)): #If average slope is negative do not add it to final list
		average_slope = (np.gradient(y_slices[i], x_slices[i])).mean()
		if average_slope < 0:
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
			sliced_x = [sublist for sublist in sliced_x if len(sublist) >= 20]
			sliced_y = [sublist for sublist in sliced_y if len(sublist) >= 20]
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
	inflection_point = smoothed_x_coords[int(np.where((y_grad==max(y_grad)))[0])]
	infl_slope = max(y_grad) #Find slope at inflection point
	low_asm = min(smoothed_y_coords) #Get y coord of lower asymptotes
	high_asm = max(smoothed_y_coords) #Get y coord of upper asymptote
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
		Tm = predx[int(np.where((y_grad==max(y_grad)))[0])]
		message = ''
	except TypeError:		
		try:
			popt, pcov = optimization.curve_fit(model, predx, predy, init_params, maxfev = maxfevopt,method='lm') #Try alternate method which can deal with slightly more misbehaved curves
			model_y_list_loop = model(predx, *popt) #Get model y-coordinates
			residuals = predy - model_y_list_loop #Calculate residuals (difference between data points and model)
			RSS = np.sum(residuals**2) #Get residual sum of squares
			MSE = np.mean(residuals**2) #Get residual mean of squares
			y_grad = np.gradient(model_y_list_loop, predx)
			Tm = predx[int(np.where((y_grad==max(y_grad)))[0])]
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
				Tm = predx[int(np.where((y_grad==max(y_grad)))[0])]
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

def process_well(key, sub_df):
	original_curve_rows = []
	sub_curve_rows = []
	temps = []
	smooth_fluorescence = []
	boltzmann_y_list = []
	new_Tm_rows = []
	if sub_df.Well_type.unique()[0] == 'Blank': #Skip the loop if the well is classified as blank. Total waste of time to analyze rubbish data
		new_Tm_rows = [{'Assay_Plate': sub_df.Assay_Plate.unique()[0], 'Source_Plate': sub_df.Source_Plate.unique()[0],'Well': sub_df.Well.unique()[0],'Unique_key': sub_df.Unique_key.unique()[0],'Subplot':None,'Well_type': sub_df.Well_type.unique()[0], 'Compound': None,'Fraction': None,'Smooth_Tm': None,'Boltzmann_Tm':None,'Boltzmann_RSS': None,'Amplitude':None,'Curve_height': None, 'Error':'','Warning':''}]
	else:
		x = list(sub_df["Temp"].values.tolist()) #Convert Temp and fluorescence columns to lists of values
		y = list(sub_df["Normalized_Fluorescence"].values.tolist())
		try:
			spl = splrep(x,y, s = 0.0005) #Smooth data. This used to be 0.05 when working with raw data. This changes drastically when dealing with normalized data
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
			final_message = "Well failed: No region of positive slope"
			new_Tm_rows = [{'Assay_Plate': sub_df.Assay_Plate.unique()[0], 'Source_Plate': sub_df.Source_Plate.unique()[0],'Well': sub_df.Well.unique()[0],'Unique_key': sub_df.Unique_key.unique()[0],'Subplot':None,'Well_type': sub_df.Well_type.unique()[0], 'Compound': sub_df.Compound.unique()[0],'Fraction': sub_df.Fraction.unique()[0],'Smooth_Tm': None,'Boltzmann_Tm':None,'Boltzmann_RSS': None,'Amplitude':None,'Curve_height': None, 'Error':final_message,'Warning':''}]
		else:
			sub_plot_count = 0
			final_list_of_x_coords, final_list_of_y_coords = split_curves(list_of_x_coords, list_of_y_coords) #Submit clean curve for second round of peak identification
			final_list_of_y_coords = [list(arr) for arr in final_list_of_y_coords]
			for k in range(len(final_list_of_x_coords)):
				sub_plot_count +=1 #Get a count for the subplots
				p0 = initial_params(final_list_of_y_coords[k],final_list_of_x_coords[k]) #Find initial parameters for model
				maxfev_opt = 100*(len(final_list_of_y_coords[k])+1) #Optimize maxfev parameter for each dataset https://www.physics.utoronto.ca/~phy326/python/curve_fit_to_data.py
				boltzmann_Tm, boltzmann_y_coordinates, MSE, RSS, final_message = Model_data(boltzmann_sigmoid, final_list_of_x_coords[k],final_list_of_y_coords[k], p0, maxfev_opt) #Fit Boltzmann
				boltzmann_y_coordinates = (boltzmann_y_coordinates.tolist())
				sub_curve_rows.append(add_curve_data(key, sub_plot_count, final_list_of_x_coords[k], final_list_of_y_coords[k],boltzmann_y_coordinates))
				new_Tm_rows.extend([{'Assay_Plate': sub_df.Assay_Plate.unique()[0], 'Source_Plate': sub_df.Source_Plate.unique()[0],'Well': sub_df.Well.unique()[0],'Unique_key': sub_df.Unique_key.unique()[0],'Subplot':sub_plot_count,'Well_type': sub_df.Well_type.unique()[0], 'Compound': sub_df.Compound.unique()[0],'Fraction': sub_df.Fraction.unique()[0],'Smooth_Tm': p0[2],'Boltzmann_Tm':boltzmann_Tm,'Boltzmann_RSS': RSS,'Amplitude':p0[3] - p0[0],'Curve_height': p0[3], 'Error':final_message,'Warning':''}])
		return original_curve_rows, sub_curve_rows, new_Tm_rows

def process_wrapper(key,final_df):
    sub_df = final_df.loc[final_df['Unique_key'] == key]
    return process_well(key, sub_df)

