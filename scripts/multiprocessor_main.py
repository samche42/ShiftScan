#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
import os
import pandas as pd
import gc 

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_dir", help="Full file path to directory with input files")
parser.add_argument("-m", "--metadata", help="Full file path to metadata file")
parser.add_argument("-p", "--processors", help="No. of processors you want to use", default = 4)
parser.add_argument("-c", "--control_cols", help="The column numbers of your controls", default = "1,2")
parser.add_argument("-o", "--output_dir", help="Full path to desired output directory", default = "./")
parser.add_argument("-x", "--failed_control_wells", help="The number of controls allowed to fail before a plate is failed", default = 8)
parser.add_argument("-t", "--ctrl_tm_cutoff", help="Maximum z-score of ctrl curve melting temps", default = 1.5)
parser.add_argument("-a", "--ctrl_amp_cutoff", help="Maximum z-score of ctrl curve amplitude", default = 3)
parser.add_argument("-u", "--max_amp_cutoff", help="Maximum relative amplitude of curves allowed", default = 6)
parser.add_argument("-l", "--min_amp_cutoff", help="Minimum relative amplitude of curves allowed", default = 0.25)
parser.add_argument("-s", "--smoothing_factor", help="Desired smoothing factor", default = 0.0005)
parser.add_argument("-n", "--normalization", help="Should data be normalized, y or n", default = "y")

args = parser.parse_args()

#Concatenate all input files to main file
print("Reading in raw data...")

files = Path(args.input_dir).glob('*.txt')  #generate list of all txt files
#Confirm there are at least some files in specified input folder
if any(Path(args.input_dir).glob('*.txt')):
    dfs = list() #Empty list of all files (converted to dfs)
    for file in files:
        try:
            data = pd.read_csv(file, sep='\t',header=0)
            data['Origin of data'] = file.stem #Stem: Get filename without extension
            dfs.append(data)
        except:
            print("There was a problem concatenating file: "+str(file)+". Please confirm the file is tab delimited and has correct headers")
else:
    print("No text files found in specified input folder")
    sys.exit() 

raw_input_df = pd.concat(dfs, ignore_index=True)
raw_input_df.columns = raw_input_df.columns.str.replace("X.", "X (", regex = True)
raw_input_df.columns = [x+')' if 'X (' in x else x for x in raw_input_df.columns]

#Prep output directory string
output_dir_string = args.output_dir
if output_dir_string.endswith('/'):
    output_dir_string = output_dir_string[:-1]

#Confirm that output directory exists, if not, create it:
if not os.path.exists(output_dir_string):
    # If it doesn't exist, create it
    os.makedirs(output_dir_string)
    print(f"Output directory not found. Creating directory at '{output_dir_string}'")

#Define control columns
control_string = args.control_cols
control_list = [int(control) for control in control_string.split(",")]
if len(control_list) == 0:
    print("ERROR: No control columns specified")
    sys.exit() 


if __name__ == '__main__':

    ################################
    #
    # PART 1: Reading in and parsing the the experimental data 
    #
    ################################
    print("Reading in concatenated data")

    import numpy as np
    try:
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
        print("Concatenated data successfully loaded. Moving on...")
    except:
        print("There was an issue parsing concatenated data. Please confirm input data format is correct")

    ################################
    #
    # PART 2: Reading in and parsing the the mapping data 
    #
    ################################
    print("Reading in metadata")
    try:
        map_raw_df = pd.read_csv(args.metadata,sep='\t',header=0)
    except:
        print("Metadata could not be read in. Please check that data is tab delimited and has a header in the first row")
        sys.exit()
    try:
        map_raw_df = map_raw_df.rename(columns={'ASSAY_PLATE': 'Assay_Plate','SOURCE_PLATE':'Source_Plate','WELL':'Well','COMPOUND':'Compound','FRACTION':'Fraction'}) 
    except:
        print("Metadata could not be read in. Please confirm that headers are as expected")
        sys.exit()
    try:
        map_raw_df['Row'] = map_raw_df['Well'].str[0] #Take first character from Well string and make it Row value 
        map_raw_df['Column'] = map_raw_df['Well'].str[1:] #Take everything beyond first character as Column (includes padded zeros, as temporarily kept as string)
        map_raw_df['Column'] = map_raw_df['Column'].str.replace(r'^(0+)', '',regex=True) #Strip padding zeros
        map_raw_df['Unique_key'] = map_raw_df['Assay_Plate']+'_'+map_raw_df['Row']+map_raw_df['Column'] #Add in new column with generated unique key
        map_raw_df = map_raw_df.astype({'Column':'int'}) #Convert Column column to integer type
        map_raw_df.Fraction = map_raw_df.Fraction.fillna('NA')
        print("Metadata successfully loaded. Onwards...")
    except: 
        print("A problem was encountered while parsing metadata. Please confirm that format is correct")
        sys.exit()

    #Now let's subset and do an inner merge
    tmp = map_raw_df[['Assay_Plate','Source_Plate']].drop_duplicates()
    semifinal_df2 = pd.merge(semifinal_df,tmp,on= 'Assay_Plate', how = 'inner')
    final_df = pd.merge(semifinal_df2,map_raw_df[['Unique_key','Compound','Fraction']].drop_duplicates(), on = 'Unique_key', how = 'left') #Left keeps all wells (i.e controls with no compounds)

    #Assign well types
    is_control = final_df['Column'].isin(control_list)
    is_blank = final_df['Compound'].isnull()
    final_df['Well_type'] = np.select([is_control, is_blank],['Control', 'Blank'],default='Experimental')
    print("Control columns assigned")

    #Let's delete a few dataframes taking up memory
    del parsed_df, melted_df, semifinal_df,semifinal_df2,tmp, map_raw_df,raw_input_df

    # Function to normalize values between 0 and 1
    def normalize(series):
        return series / series.max()

    # Applying normalization to the 'Fluorescence' column grouped by 'Assay_Plate'
    smoothing_fact = float(args.smoothing_factor)
    normalize_data = args.normalization

    if normalize_data == "y":
        final_df['Normalized_Fluorescence'] = final_df.groupby('Assay_Plate')['Fluorescence'].apply(normalize)
        print("Data normalized")
    elif normalize_data == "n":
        final_df['Normalized_Fluorescence'] = final_df.groupby('Assay_Plate')['Fluorescence'].apply(normalize)
        print("Normalization skipped")
    else:
        print("ERROR: Invalid normalization option. Please select either y or n for this parameter")
        sys.exit()

    ################################
    #
    # PART 4: Processing wells
    #
    ################################
    import multiprocessing
    from multiprocessing import cpu_count
    from functools import partial
    from DSF_functions_fiddling import slice_list, clean_curve, split_curves,add_curve_data, boltzmann_sigmoid, initial_params, Model_data, process_well

    #Create empty dataframe for melting temp summary output
    Tm_df = pd.DataFrame(columns = ['Assay_Plate','Source_Plate','Well','Unique_key','Subplot','Well_type','Compound','Fraction','Smooth_Tm','Boltzmann_Tm','Boltzmann_RSS','Amplitude','Curve_height','Error','Warning'])

    all_original_curve_rows = []
    all_tm_rows = []
    all_subcurves = []

    # Group by 'Unique_key' and convert each group to a DataFrame
    list_of_sub_dfs = [group.reset_index(drop=True) for _, group in final_df.groupby('Unique_key')]

    print("Multiprocessing of wells started")

    with multiprocessing.Pool(processes=int(args.processors)) as pool:
        partial_process_wrapper = partial(process_well, smoothing_factor=smoothing_fact,normalize=normalize_data)
        results = pool.map(partial_process_wrapper, list_of_sub_dfs)
        for original_curve_rows, sub_curve_rows, Tm_rows in results:
            all_original_curve_rows.extend(original_curve_rows)
            all_subcurves.extend(sub_curve_rows)
            all_tm_rows.extend(Tm_rows)
        del results
        gc.collect()

    print("Multiprocessing of wells ended")

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
    print("Tm_df created with "+str(len(Tm_df))+" rows of data")

    ################################
    #
    # PART 5: Finding and flagging bad results
    #
    ################################
    from scipy import stats
    print("Looking for bad control wells")

    control_z_scores = Tm_df[Tm_df['Well_type'] == 'Control'].groupby(['Assay_Plate'],as_index=False)[['Smooth_Tm','Amplitude']].transform(stats.zscore).abs() #Generate z-scores for control wells
    control_z_scores = control_z_scores.rename(columns={'Smooth_Tm': 'Ctrl_Tm_z-score', 'Amplitude':'Ctrl_Amplitude_z-score'}) #Rename columns
    if control_z_scores.isnull().values.any() == True:
        print("Wow! At least one plate with perfect controls! Nice!")
        control_z_scores = control_z_scores.fillna(0)
    Tm_df = Tm_df.join(control_z_scores) #merge on index

    #Find differences between smootehd and modeled data
    Tm_df['Tm_SE'] = Tm_df[['Smooth_Tm','Boltzmann_Tm']].sem(axis=1, skipna = False) #calculate std error between smoothed data and Boltzmann
    Tm_df['Tm_difference'] = (Tm_df['Boltzmann_Tm']-Tm_df['Smooth_Tm']).abs() #Find absolute difference between smoothed and boltzmann Tm

    #Fail ctrls where the z-scores are beyond the designated cutoffs
    Tm_cutoff = float(args.ctrl_tm_cutoff)
    Amp_cutoff = float(args.ctrl_amp_cutoff)
    Tm_df['Error'] = Tm_df['Error'].fillna('')
    Tm_df['Error'] = np.where((Tm_df['Ctrl_Tm_z-score'] > Tm_cutoff), Tm_df['Error']+ "Failed: Tm > "+str(Tm_cutoff)+" std away from mean. ", Tm_df['Error'])
    Tm_df['Error'] = np.where((Tm_df['Ctrl_Amplitude_z-score'] > Amp_cutoff), Tm_df['Error']+ "Failed: Amplitude > "+str(Amp_cutoff)+" std away from mean. ", Tm_df['Error'])

    print("Bad controls identified and removed from downstream analyses")

    ################################
    #
    # PART 6: Calulating results with bad wells removed
    #
    ################################
    import math
    print("Flagging bad experimental wells and calculating z-scores of good wells..")

    #First, let's QC our experimentals
    avg_ctrl_amplitude_df = Tm_df[(Tm_df['Well_type'] == 'Control') & (~Tm_df['Error'].str.contains('ailed'))].groupby('Assay_Plate', as_index=False)['Amplitude'].mean() #Get average amplitude of all non-failed ctrl wells per plate
    avg_ctrl_amplitude_df = avg_ctrl_amplitude_df.rename(columns={'Amplitude': 'Ctrl_avg_amplitude'}) #Rename column
    Tm_df = pd.merge(Tm_df, avg_ctrl_amplitude_df, on = 'Assay_Plate') #Merge info back in to Tm_df

    #Fail wells that have a amplitude less than 1/4 of the ctrl amplitude
    amp_lower_cutoff = float(args.min_amp_cutoff)
    Tm_df['Amplitude_lower_cutoff'] = Tm_df['Ctrl_avg_amplitude']*amp_lower_cutoff #Just add the cutoff value to Tm_df for easy calculation
    Tm_df['Error'] = np.where((Tm_df['Amplitude'] < Tm_df['Amplitude_lower_cutoff']), Tm_df['Error']+ " Failed: Amplitude too small. ", Tm_df['Error'])
    amp_upper_cutoff = float(args.max_amp_cutoff)
    Tm_df['Amplitude_upper_cutoff'] = Tm_df['Ctrl_avg_amplitude']*amp_upper_cutoff #Just add the cutoff value to Tm_df for easy calculation
    Tm_df['Error'] = np.where((Tm_df['Amplitude'] > Tm_df['Amplitude_upper_cutoff']), Tm_df['Error']+ " Failed: Amplitude too large. ", Tm_df['Error'])

    #IF TWO SUBPLOTS EXISTS AND ONLY ONE FAILS, REPLACE THE STRING "FAILED" WITH "REMOVED". 
    #THIS WAY IT WILL NOT COUNT TOWARD THE FINAL FAILURE COUNT PER PLATE
    #ADD IN NEW COLUMN TO Tm_df called "Final_decision": Can have "Pass, Failed, or Removed"

    # Group by unique_key
    grouped = Tm_df.groupby('Unique_key')

    # Define a function to apply to each group
    def final_decision(group):
        if group['Error'].eq('').any():
            group['Final_decision'] = np.where(group['Error'] == '', 'Pass', 'Removed')
        else:
            group['Final_decision'] = 'Failed'
        return group

    # Apply the function to each group and concatenate the results
    Tm_df = pd.concat([final_decision(group) for _, group in grouped], ignore_index=True)

    #next, we'll get a final Tm for all wells that have passed (ctrl and experi)
    Tm_df['Final_Tm'] = np.where(Tm_df['Error'].str.contains('ailed'), np.NaN, Tm_df['Smooth_Tm']) #Report Smooth Tm as final Tm unless well failed

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

    well_type_df = Tm_df.loc[:, ['Unique_key', 'Well_type']]
    decision_df = Tm_df.loc[:, ['Unique_key','Subplot','Final_decision','Error']]
    decision_df['Subplot'] = decision_df['Subplot'].astype(pd.Int64Dtype()).astype(str)
    decision_df['Subplot'] = decision_df['Subplot'].replace('<NA>', 'Original')
    final_curves0 = pd.merge(semifinal_curves,well_type_df, on = 'Unique_key', how = 'inner')
    final_curves1 = pd.merge(final_curves0, Tm_df[['Unique_key','Source_Plate','Ctrl_Tm_z-score']],on = 'Unique_key', how = 'outer')
    final_curves = pd.merge(final_curves1,decision_df, on = ['Unique_key','Subplot'], how = 'outer')
    print("Calculations complete!")

    ################################
    #
    # PART 7: Generating reports
    #
    ################################
    from functools import reduce

    print("Generating reports...")

    failed_df = Tm_df[Tm_df['Final_decision'] == 'Failed']
    all_plates = Tm_df['Assay_Plate'].unique()

    #Count total failed controls per assay plate
    total_failures_series = failed_df[(failed_df['Well_type'] == 'Control')].groupby('Assay_Plate').size()
    total_failures_series = total_failures_series.reindex(all_plates, fill_value = 0)
    total_failures_df = total_failures_series.reset_index(name='Total_failures')
    total_failures_df['Plate_status'] = np.where(total_failures_df['Total_failures'] > args.failed_control_wells, "Failed","OK")

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
    amp_failures_series = failed_df[(failed_df['Well_type'] == 'Control')&(failed_df['Error'].str.contains('Failed: Amplitude too'))].groupby('Assay_Plate').size()
    amp_failures_series = amp_failures_series.reindex(all_plates, fill_value = 0)
    amp_failures_df = amp_failures_series.reset_index(name='Amp_Z-score_failures')

    #List out failed wells
    failed_wells = Tm_df[Tm_df['Well_type'] == 'Control'][Tm_df[Tm_df['Well_type'] == 'Control']['Final_decision'].str.contains('Failed', case=False)].groupby(['Assay_Plate'])['Well'].agg(lambda x: ','.join(x)).reset_index()
    failed_wells= failed_wells.rename(columns={'Well':'Failed ctrl wells'})
    
    #Put it all together
    plate_status_df = Tm_df[['Assay_Plate','Source_Plate']].drop_duplicates()
    list_of_dfs = [plate_status_df,total_failures_df,failed_wells,melting_temp_failures_df,model_failure_df,amp_failures_df]
    plate_report = reduce(lambda  left,right: pd.merge(left,right,on=['Assay_Plate'],how='outer'), list_of_dfs)

    plate_status_only = plate_report.loc[:, ['Assay_Plate', 'Plate_status']]
    Tm_df = pd.merge(Tm_df, plate_status_only, on = 'Assay_Plate') #merge plate status back in for visualization purposes
    Tm_df['Platename'] = np.where(Tm_df['Plate_status'] =='Failed', Tm_df['Source_Plate']+" ("+Tm_df['Plate_status']+")", Tm_df['Source_Plate']) #Add "Failed" to plate name to make it obvious in visulizations

    tmp = Tm_df[Tm_df['Well_type'] == 'Control'].groupby(['Well'],as_index=False)['Error'].apply(lambda x: x[x.str.contains('Failed')].count())
    well_error_count_df = tmp[tmp['Error'] >= 3]
    well_error_count_df = well_error_count_df.rename(columns={'Error':'Wells that repeatedly failed'})

    Tm_df['Unique_key_subplot'] = Tm_df['Unique_key']+"_"+Tm_df['Subplot'].astype(str)

    print("Writing output files")
    Tm_df.to_csv(output_dir_string+"/Final_results.txt",sep="\t",index=False)
    final_curves.to_csv(output_dir_string+"/Final_curves.txt",sep="\t",index=False)
    plate_report.to_csv(output_dir_string+"/Plate_report.txt",sep="\t",index=False)
    well_error_count_df.to_csv(output_dir_string+"/Potential_problems.txt",sep="\t",index=False)
    print("Analysis complete!")
