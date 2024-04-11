################################
#
# PART 1: Reading in and parsing the the experimental data 
#
################################

import pandas as pd
import numpy as np

#raw_input_df = pd.read_csv("trial_concatenated.txt",sep='\t',header=0)
parsed_df = raw_input_df.T.drop_duplicates().T #Remove duplicated temp columns
parsed_df  = parsed_df.rename(columns={'X': 'Temp', 'Origin of data':'Origin'}) #Rename temp column
melted_df = pd.melt(parsed_df,id_vars=['Temp','Origin']) #Melt dataframe
melted_df['Well'] = melted_df.variable.str.split(':', expand=True)[0] #create well column
melted_df['Assay_Plate'] = melted_df.Origin.str.split('.', expand=True)[0] #Create plate column from origin data
melted_df['Row'] = melted_df['Well'].str[0] #Take first character from Well string and make it Row value 
melted_df['Column'] = melted_df['Well'].str[1:]
melted_df = melted_df.astype({'Column':'int'})
melted_df  = melted_df.rename(columns={'value': 'Fluorescence'}) #Rename columns
semifinal_df = melted_df.drop(['Origin','variable'],axis = 1) #Get rid of useless/redundant columns
semifinal_df['Unique_key'] = semifinal_df['Assay_Plate']+'_'+semifinal_df['Well'] #Add in new column with generated unique key

################################
#
# PART 2: Reading in and parsing the the mapping data 
#
################################

#map_raw_df = pd.read_csv("RBD_Dataset.csv",sep='\t',header=0)
map_raw_df = map_raw_df.rename(columns={'ASSAY_PLATE': 'Assay_Plate','SOURCE_PLATE':'Source_Plate','WELL':'Well','COMPOUND':'Compound','FRACTION':'Fraction'}) 
map_raw_df['Row'] = map_raw_df['Well'].str[0] #Take first character from Well string and make it Row value 
map_raw_df['Column'] = map_raw_df['Well'].str[1:] #Take everything beyond first character as Column (includes padded zeros, as temporarily kept as string)
map_raw_df['Column'] = map_raw_df['Column'].str.replace(r'^(0+)', '',regex=True) #Strip padding zeros
map_raw_df['Unique_key'] = map_raw_df['Assay_Plate']+'_'+map_raw_df['Row']+map_raw_df['Column'] #Add in new column with generated unique key
map_raw_df = map_raw_df.astype({'Column':'int'}) #Convert Column column to integer type
map_raw_df.Fraction = map_raw_df.Fraction.fillna('NA')

#Now let's subset and do an inner merge
tmp = map_raw_df[['Assay_Plate','Source_Plate']].drop_duplicates()
semifinal_df2 = pd.merge(semifinal_df,tmp,on= 'Assay_Plate', how = 'inner')
final_df = pd.merge(semifinal_df2,map_raw_df[['Unique_key','Compound','Fraction']].drop_duplicates(), on = 'Unique_key', how = 'left') #Left keeps all wells (i.e controls with no compounds)

#Working in tsinghua substring check into this monstrous super nested, glorified if then statement

final_df['Well_type'] = (np.where(final_df['Source_Plate'].str.contains('singhua'), \
	np.where((final_df['Column'] == 1) | (final_df['Column'] == 24), 'Control', (np.where((final_df['Compound'].isnull()), 'Blank', 'Experimental'))), \
	np.where((final_df['Column'] == 1) | (final_df['Column'] == 2), 'Control', (np.where((final_df['Compound'].isnull()), 'Blank', 'Experimental')))))

#Finally, we're going to try and bypass the noise of early high fluorescence by removing all data points with Temps lower than 35 degrees. 
final_df = final_df[final_df['Temp'] > 35]
final_df = final_df[final_df['Temp'] < 70]

#Let's delete a few dataframes taking up memory
del parsed_df, melted_df, semifinal_df,semifinal_df2,tmp

# Function to normalize values between 0 and 1
def normalize(series):
    return series / series.max()

# Applying normalization to the 'Fluorescence' column grouped by 'Assay_Plate'
final_df['Normalized_Fluorescence'] = final_df.groupby('Assay_Plate')['Fluorescence'].apply(normalize)

################################
#
# PART 3: Getting everything ready: Generating blank lists/dataframes, defining functions etc
#
################################

import scipy.optimize as optimization
from scipy.interpolate import splev, splrep
from scipy.optimize import OptimizeWarning
import warnings
from scipy.signal import find_peaks

#Create empty dataframe for melting temp summary output
Tm_df = pd.DataFrame(columns = ['Assay_Plate','Source_Plate','Well','Unique_key','Subplot','Well_type','Compound','Fraction','Smooth_Tm','Boltzmann_Tm','Boltzmann_RSS','Amplitude','Curve_height','Error','Warning'])

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

################################
#
# PART 4: Processing wells
#
################################

all_original_curve_rows = []
all_tm_rows = []
all_subcurves = []

#Process each well
unique_keys = list(final_df.Unique_key.unique()) #Create list of all unique keys for subsetting
for i in unique_keys:
	sub_df = final_df.loc[final_df['Unique_key'] == i]
	original_curve_rows, sub_curve_rows, Tm_rows = process_well(i, sub_df)
	all_original_curve_rows.extend(original_curve_rows) #Add new curve rows to overall list of lists
	all_subcurves.extend(sub_curve_rows)
	all_tm_rows.extend(Tm_rows) #Add new rows of Tm data

#Clean up raw curve coordinate tables
original_curve_dataframes = [pd.DataFrame(d) for d in all_original_curve_rows]
original_curve_df = pd.concat(original_curve_dataframes, ignore_index=True)

#Get subcurve table set up
all_subcurves_dataframes = [pd.DataFrame(d) for d in all_subcurves]
all_subcurves_df = pd.concat(all_subcurves_dataframes, ignore_index=True)

#Merge curves and subcurves
semifinal_curves = pd.concat([original_curve_df,all_subcurves_df], ignore_index=True)
semifinal_curves['Subplot'].fillna('Original', inplace=True)
semifinal_curves['Subplot'] = semifinal_curves['Subplot'].astype(str) 
semifinal_curves[['Assay_Plate', 'Well']] = semifinal_curves['Unique_key'].str.rsplit('_', 1, expand=True) #Add in plate and well info by breaking up Unique_key

#Generate Tm df
Tm_df = pd.DataFrame(all_tm_rows)

################################
#
# PART 5: Finding and flagging bad results
#
################################
from scipy import stats

#Count how many wells failed in modelling stage 
ctrl_error_df = Tm_df[Tm_df['Well_type'] == 'Control'].groupby(['Assay_Plate'],as_index=False)['Error'].apply(lambda x: x[x.str.contains('ailed')].count())
ctrl_error_df = ctrl_error_df.rename(columns={'Error': 'Modelling failures'}) 

control_z_scores = Tm_df[Tm_df['Well_type'] == 'Control'].groupby(['Assay_Plate'],as_index=False)[['Smooth_Tm','Amplitude']].transform(stats.zscore).abs() #Generate z-scores for control wells
control_z_scores = control_z_scores.rename(columns={'Smooth_Tm': 'Ctrl_Tm_z-score', 'Amplitude':'Ctrl_Amplitude_z-score'}) #Rename columns
Tm_df = Tm_df.join(control_z_scores) #merge on index

#Find differences between smootehd and modeled data
Tm_df['Tm_SE'] = Tm_df[['Smooth_Tm','Boltzmann_Tm']].sem(axis=1, skipna = False) #calculate std error between smoothed data and Boltzmann
Tm_df['Tm_difference'] = (Tm_df['Boltzmann_Tm']-Tm_df['Smooth_Tm']).abs() #Find absolute difference between smoothed and boltzmann Tm

#Fail ctrls where the z-scores are beyond the designated cutoffs
Tm_cutoff = 2
Amp_cutoff = 3
Tm_df['Error'] = np.where((Tm_df['Ctrl_Tm_z-score'] > Tm_cutoff), Tm_df['Error']+ "Failed: Tm > 2 std away from mean. ", Tm_df['Error'])
Tm_df['Error'] = np.where((Tm_df['Ctrl_Amplitude_z-score'] > Amp_cutoff), Tm_df['Error']+ "Failed: Amplitude > 3 std away from mean. ", Tm_df['Error'])

#Count up ctrls that fail due to deviation from Tm mean
ctrl_tm_failures = Tm_df[Tm_df['Well_type'] == 'Control'].groupby(['Assay_Plate'], as_index=False)['Ctrl_Tm_z-score'].apply(lambda x: (x > Tm_cutoff).sum())
ctrl_tm_failures = ctrl_tm_failures.rename(columns={'Ctrl_Tm_z-score': 'Tm_Z-score_failures'}) #Rename column
ctrl_error_df = pd.merge(ctrl_error_df, ctrl_tm_failures, on=['Assay_Plate'])

#Count up ctrls that fail due to deviation from amplitude or height means
ctrl_amplitude_failures = Tm_df[Tm_df['Well_type'] == 'Control'].groupby(['Assay_Plate'],as_index=False)['Ctrl_Amplitude_z-score'].apply(lambda x: (x > Amp_cutoff).sum())
ctrl_amplitude_failures = ctrl_amplitude_failures.rename(columns={'Ctrl_Amplitude_z-score': 'Amp_Z-score_failures'}) #Rename column
ctrl_error_df = pd.merge(ctrl_error_df, ctrl_amplitude_failures, on=['Assay_Plate'])

################################
#
# PART 6: Calulating results with bad wells removed
#
################################
import math
import numpy as np

#First, let's QC our experimentals
avg_ctrl_amplitude_df = Tm_df[(Tm_df['Well_type'] == 'Control') & (~Tm_df['Error'].str.contains('Failed'))].groupby('Assay_Plate', as_index=False)['Amplitude'].mean() #Get average amplitude of all non-failed ctrl wells per plate
avg_ctrl_amplitude_df = avg_ctrl_amplitude_df.rename(columns={'Amplitude': 'Ctrl_avg_amplitude'}) #Rename column
Tm_df = pd.merge(Tm_df, avg_ctrl_amplitude_df, on = 'Assay_Plate') #Merge info back in to Tm_df

#Fail wells that have a amplitude less than 1/4 of the ctrl amplitude
amp_lower_cutoff = 0.25 
Tm_df['Amplitude_lower_cutoff'] = Tm_df['Ctrl_avg_amplitude']*amp_lower_cutoff #Just add the cutoff value to Tm_df for easy calculation
Tm_df['Error'] = np.where((Tm_df['Amplitude'] < Tm_df['Amplitude_lower_cutoff']), Tm_df['Error']+ "Failed: Amplitude too small. ", Tm_df['Error'])
amp_upper_cutoff = 6
Tm_df['Amplitude_upper_cutoff'] = Tm_df['Ctrl_avg_amplitude']*amp_upper_cutoff #Just add the cutoff value to Tm_df for easy calculation
Tm_df['Error'] = np.where((Tm_df['Amplitude'] > Tm_df['Amplitude_upper_cutoff']), Tm_df['Error']+ "Failed: Amplitude too large. ", Tm_df['Error'])

#IF TWO SUBPLOTS EXISTS AND ONLY ONE FAILS, REPLACE THE STRING "FAILED" WITH "REMOVED". 
#THIS WAY IT WILL NOT COUNT TOWARD THE FINAL FAILURE COUNT PER PLATE
#ADD IN NEW COLUMN TO Tm_df called "Final_decision": Can have "Pass, Failed, or Removed"
decision_df = pd.DataFrame()
for key in unique_keys:
    test_df = (Tm_df[Tm_df['Unique_key'] == key]).loc[:, ['Unique_key', 'Subplot','Error']]
    if len(test_df[test_df['Error'] == '']) > 0:
        test_df['Final_decision'] = np.where(test_df['Error'] == '', 'Pass','Removed')
    else:
        test_df['Final_decision'] = "Failed"
    decision_df = pd.concat([decision_df,test_df], ignore_index=True)

decision_df = decision_df.drop(['Error'],axis = 1)
Tm_df =pd.merge(Tm_df, decision_df, on = ["Unique_key","Subplot"], how = 'inner')

#Summarize and count errors
failed_wells = Tm_df[Tm_df['Well_type'] == 'Control'][Tm_df[Tm_df['Well_type'] == 'Control']['Final_decision'].str.contains('Failed', case=False)].groupby(['Assay_Plate'])['Well'].agg(lambda x: ','.join(x)).reset_index() #Get list of failed wells
failed_wells['Total'] = failed_wells['Well'].apply(lambda x: len(x.split(','))) #Count failed wells listed
ctrl_error_df = pd.merge(ctrl_error_df, failed_wells, on=['Assay_Plate'], how = 'outer') #Merge that back in to the ctrl df
ctrl_error_df['Total'].fillna(0, inplace=True) #Fill in nans with 0s
ctrl_error_df['Total'] = ctrl_error_df['Total'].astype(int) #Convert column values to integers instead of float (just looks neater to me)
plate_status_df = ctrl_error_df.groupby('Assay_Plate')['Total'].sum().reset_index() #Count failures per plate
plate_status_df['Plate_status'] = np.where(plate_status_df['Total'] > 8, 'Failed','OK') #Assign plate status according to ctrl failure cutoff
Tm_df = pd.merge(Tm_df, plate_status_df, on = 'Assay_Plate') #merge that all back in
Tm_df['Spotfire_Platename'] = np.where(Tm_df['Total'] >8, Tm_df['Source_Plate']+" ("+Tm_df['Plate_status']+")", Tm_df['Source_Plate']) #Add "Failed" to plate name to make it obvious in visulizations

#next, we'll get a final Tm for all wells that have passed (ctrl and experi)
Tm_df['Final_Tm'] = np.where(Tm_df['Error'].str.contains('Failed'), np.NaN, Tm_df['Smooth_Tm']) #Report Smooth Tm as final Tm unless well failed

#Get averages for control wells (DOES NOT INCLUDE FAILED WELLS)
ctrl_df = Tm_df[Tm_df['Well_type'] == 'Control'].groupby(['Assay_Plate'],as_index=False)[['Final_Tm']].mean() #Get average ctrl Tm
ctrl_df = ctrl_df.rename(columns={'Final_Tm': 'Avg_ctrl_melting_temp'})
Tm_df = pd.merge(Tm_df, ctrl_df, on = 'Assay_Plate') #Merge that info in for easy calcs

#Calculate no. of std dev away from control Tm mean (z-score) for all experimental wells
#Based on: https://www.scribbr.com/statistics/standard-deviation/ 
Tm_df['Squared deviation from ctrl mean'] = (Tm_df['Final_Tm']-Tm_df['Avg_ctrl_melting_temp'])**2
sum_squares_df = Tm_df.groupby(['Assay_Plate'],as_index=False)[['Squared deviation from ctrl mean']].sum() #Find sum of squares for each plate
sum_squares_df = sum_squares_df.rename(columns={'Squared deviation from ctrl mean': 'Sum of squares'})
#successful_samples_per_plate = Tm_df[Tm_df['Final_decision'] == 'Pass'].groupby(['Assay_Plate'], as_index=False)['Final_decision'].value_counts() #Count all passed wells per plate
#successful_samples_per_plate =successful_samples_per_plate.drop(['Final_decision'],axis = 1)
successful_samples_per_plate = Tm_df[Tm_df['Final_decision'] == 'Pass'].groupby(['Assay_Plate']).size().reset_index(name='count')
std_dev_df = pd.merge(sum_squares_df,successful_samples_per_plate, on = 'Assay_Plate')
std_dev_df['Variance'] = std_dev_df['Sum of squares']/(std_dev_df['count']-1)
std_dev_df['Std_dev_for_plate'] = std_dev_df['Variance'].apply(lambda x: np.sqrt(x))
Tm_df = pd.merge(Tm_df, std_dev_df, on = 'Assay_Plate')
Tm_df['No. Std Dev'] = (Tm_df['Final_Tm']-Tm_df['Avg_ctrl_melting_temp'])/Tm_df['Std_dev_for_plate']
Tm_df['No. Std Dev'] = np.where(Tm_df['Well_type'] == 'Control', np.NaN, Tm_df['No. Std Dev']) #Nullify control std devs
Tm_df['Row'] = Tm_df['Well'].str[0] #Generate Row column (needed to build visualizations)
Tm_df['Column'] = (Tm_df['Well'].str[1:]).astype(int) #Generate Column column (needed to build visualizations)

well_type_df = Tm_df.loc[:, ['Unique_key', 'Well_type']]
decision_df = Tm_df.loc[:, ['Unique_key','Subplot','Final_decision','Error']]
decision_df['Subplot'] = decision_df['Subplot'].astype(pd.Int64Dtype()).astype(str)
decision_df['Subplot'] = decision_df['Subplot'].replace('<NA>', 'Original')
final_curves0 = pd.merge(semifinal_curves,well_type_df, on = 'Unique_key', how = 'inner')
final_curves = pd.merge(final_curves0,decision_df, on = ['Unique_key','Subplot'], how = 'outer')

################################
#
# PART 7: Generating reports
#
################################

plate_report_tmp = pd.merge(ctrl_error_df, plate_status_df, on = 'Assay_Plate')
plate_report_tmp = plate_report_tmp.drop('Total_x',axis = 1)
plate_report = plate_report_tmp.rename(columns={'Total_y':'Total_failures'})

tmp = Tm_df[Tm_df['Well_type'] == 'Control'].groupby(['Well'],as_index=False)['Error'].apply(lambda x: x[x.str.contains('Failed')].count())
well_error_count_df = tmp[tmp['Error'] >= 3]
well_error_count_df = well_error_count_df.rename(columns={'Error':'No. failed'})
