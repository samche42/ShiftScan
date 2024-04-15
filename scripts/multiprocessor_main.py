import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--concat_file", help="Full file path to concatenated input file")
parser.add_argument("-m", "--metadata", help="Full file path to metadata file")
parser.add_argument("-p", "--processors", help="No. of processors you want to use")
parser.add_argument("-c", "--control_cols", help="No. of processors you want to use")

args = parser.parse_args()

control_string = args.control_cols
control_list = [int(control) for control in control_string.split(",")]

if __name__ == '__main__':

    ################################
    #
    # PART 1: Reading in and parsing the the experimental data 
    #
    ################################
    print("Reading in raw data")

    import pandas as pd
    import numpy as np

    raw_input_df = pd.read_csv(args.concat_file,sep='\t',header=0)
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
    print("Raw data successfully loaded")

    ################################
    #
    # PART 2: Reading in and parsing the the mapping data 
    #
    ################################
    print("Reading in metadata")

    map_raw_df = pd.read_csv(args.metadata,sep='\t',header=0)
    map_raw_df = map_raw_df.rename(columns={'ASSAY_PLATE': 'Assay_Plate','SOURCE_PLATE':'Source_Plate','WELL':'Well','COMPOUND':'Compound','FRACTION':'Fraction'}) 
    map_raw_df['Row'] = map_raw_df['Well'].str[0] #Take first character from Well string and make it Row value 
    map_raw_df['Column'] = map_raw_df['Well'].str[1:] #Take everything beyond first character as Column (includes padded zeros, as temporarily kept as string)
    map_raw_df['Column'] = map_raw_df['Column'].str.replace(r'^(0+)', '',regex=True) #Strip padding zeros
    map_raw_df['Unique_key'] = map_raw_df['Assay_Plate']+'_'+map_raw_df['Row']+map_raw_df['Column'] #Add in new column with generated unique key
    map_raw_df = map_raw_df.astype({'Column':'int'}) #Convert Column column to integer type
    map_raw_df.Fraction = map_raw_df.Fraction.fillna('NA')
    print("Metadata successfully loaded")

    #Now let's subset and do an inner merge
    tmp = map_raw_df[['Assay_Plate','Source_Plate']].drop_duplicates()
    semifinal_df2 = pd.merge(semifinal_df,tmp,on= 'Assay_Plate', how = 'inner')
    final_df = pd.merge(semifinal_df2,map_raw_df[['Unique_key','Compound','Fraction']].drop_duplicates(), on = 'Unique_key', how = 'left') #Left keeps all wells (i.e controls with no compounds)

    #Assign well types
    is_control = final_df['Column'].isin(control_list)
    is_blank = final_df['Compound'].isnull()
    final_df['Well_type'] = np.select([is_control, is_blank],['Control', 'Blank'],default='Experimental')
    print("Control columns assigned")

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
    print("Data normalized")

    ################################
    #
    # PART 4: Processing wells
    #
    ################################
    import multiprocessing
    from multiprocessing import cpu_count
    from functools import partial
    from DSF_functions import slice_list, clean_curve, split_curves,add_curve_data, boltzmann_sigmoid, initial_params, Model_data, process_well, process_wrapper

    #Create empty dataframe for melting temp summary output
    Tm_df = pd.DataFrame(columns = ['Assay_Plate','Source_Plate','Well','Unique_key','Subplot','Well_type','Compound','Fraction','Smooth_Tm','Boltzmann_Tm','Boltzmann_RSS','Amplitude','Curve_height','Error','Warning'])

    all_original_curve_rows = []
    all_tm_rows = []
    all_subcurves = []

    #Process each well
    unique_keys = list(final_df.Unique_key.unique()) #Create list of all unique keys for subsetting
    print("Multiprocessing of wells started")

    with multiprocessing.Pool(processes=int(args.processors)) as pool:
        partial_process_wrapper = partial(process_wrapper, final_df=final_df)
        results = pool.map(partial_process_wrapper, unique_keys)
        for original_curve_rows, sub_curve_rows, Tm_rows in results:
            all_original_curve_rows.extend(original_curve_rows)
            all_subcurves.extend(sub_curve_rows)
            all_tm_rows.extend(Tm_rows)

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
    print("Bad controls identified and removed from downstream analyses")

    ################################
    #
    # PART 6: Calulating results with bad wells removed
    #
    ################################
    import math
    print("Flagging bad experimental wells and calculating z-scores of good wells..")

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
    successful_samples_per_plate = Tm_df[Tm_df['Final_decision'] == 'Pass'].groupby(['Assay_Plate']).size().reset_index(name='count')
    std_dev_df = pd.merge(sum_squares_df,successful_samples_per_plate, on = 'Assay_Plate')
    std_dev_df['Variance'] = std_dev_df['Sum of squares']/(std_dev_df['count']-1)
    std_dev_df['Std_dev_for_plate'] = std_dev_df['Variance'].apply(lambda x: np.sqrt(x))
    Tm_df = pd.merge(Tm_df, std_dev_df, on = 'Assay_Plate')
    Tm_df['No. Std Dev'] = (Tm_df['Final_Tm']-Tm_df['Avg_ctrl_melting_temp'])/Tm_df['Std_dev_for_plate']
    Tm_df['No. Std Dev'] = np.where(Tm_df['Well_type'] == 'Control', np.NaN, Tm_df['No. Std Dev']) #Nullify control std devs
    Tm_df['Row'] = Tm_df['Well'].str[0] #Generate Row column (needed to build visualizations)
    Tm_df['Column'] = (Tm_df['Well'].str[1:]).astype(int) #Generate Column column (needed to build visualizations)
    Tm_df['Relative_amplitude'] = Tm_df['Amplitude']/Tm_df['Ctrl_avg_amplitude']

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
    print("Generating reports...")

    plate_report_tmp = pd.merge(ctrl_error_df, plate_status_df, on = 'Assay_Plate')
    plate_report_tmp = plate_report_tmp.drop('Total_x',axis = 1)
    plate_report = plate_report_tmp.rename(columns={'Total_y':'Total_failures', 'Well':'Wells'})
    plate_report = pd.merge(plate_report, Tm_df[['Assay_Plate','Source_Plate']], on = 'Assay_Plate', how = 'inner')
    plate_report = plate_report[['Source_Plate','Assay_Plate','Plate_status','Total_failures','Wells','Modelling failures','Tm_Z-score_failures','Amp_Z-score_failures']]
    plate_report.drop_duplicates(keep='first', inplace=True)

    tmp = Tm_df[Tm_df['Well_type'] == 'Control'].groupby(['Well'],as_index=False)['Error'].apply(lambda x: x[x.str.contains('Failed')].count())
    well_error_count_df = tmp[tmp['Error'] >= 3]
    well_error_count_df = well_error_count_df.rename(columns={'Error':'No. failed'})

    print("Writing out output files")
    Tm_df.to_csv("Final_results.txt",sep="\t",index=False)
    final_curves.to_csv("Final_curves.txt",sep="\t",index=False)
    plate_report.to_csv("Plate_report.txt",sep="\t",index=False)
    well_error_count_df.to_csv("Potential_problems.txt",sep="\t",index=False)
    print("DSF analysis complete!")
