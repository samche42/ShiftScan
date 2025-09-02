#!/usr/bin/env python3

import dash, math, argparse
from dash import Dash, html, dcc, Input, Output, State, dash_table, callback_context
import pandas as pd
import plotly.express as px
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate
import numpy as np

# Initialize the app
app = Dash(__name__)
app.config.suppress_callback_exceptions = True

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_dir", help="Full file path to output files from previous step")
parser.add_argument("-z", "--port", help="Optional: Port to route to", default = 8050)
parser.add_argument("--mode",choices=['full', 'only_tm', 'dose_response'],default='full',help="Specifies the visualization mode. 'full' is the default.")

args = parser.parse_args()

#Initialize app
app = dash.Dash(__name__)
app.config.suppress_callback_exceptions = True

#Load data
if args.mode == 'full':
    print("Reading in data. This can take a few minutes if processing many plates")
    df_results = pd.read_csv(args.input_dir+'/Final_results.txt', sep='\t', header = 0)
    df_results['Error'].fillna('', inplace=True)
    df_results['Fraction'].fillna(0, inplace=True)
    df_results['group'] = df_results['Assay_Plate'].apply(lambda x: df_results['Assay_Plate'].unique().tolist().index(x) // 16) #Assign group number to every 16 plates
    columns_to_round = ['Well_zscore','Relative_amplitude','Diff from ctrl avg', 'Final_Tm', 'Std. devs from ctrl mean', 'Avg_ctrl_melting_temp','Min_ctrl_zscore_for_plate','Max_ctrl_zscore_for_plate']
    df_results = df_results.assign(**{col: df_results[col].apply(lambda x: round(x, 2)) for col in columns_to_round})
    groups_df = df_results.loc[:, ['Assay_Plate', 'group']].drop_duplicates()
    df_curves = pd.read_csv(args.input_dir+'/Final_curves.txt', sep='\t', header = 0)
    df_curves = pd.merge(df_curves,groups_df, how = 'outer', on='Assay_Plate')
    df_curves['Error'].fillna('', inplace=True)
    plate_report_df = pd.read_csv(args.input_dir+'/Plate_report.txt', sep='\t', header = 0)
    pipette_df = pd.read_csv(args.input_dir+'/Potential_problems.txt', sep='\t', header = 0)
    #Generate default hit table
    original_df  = df_results[(df_results['Well_type'] != 'Control')&(df_results['Plate_status'] == 'OK')&(df_results['Error'] == '')]
    small_df = original_df.loc[:, ['Source_Plate','Well','Subplot','Compound','Fraction','Final_Tm','Well_zscore', 'Std. devs from ctrl mean','Diff from ctrl avg','Relative_amplitude','Max_ctrl_zscore_for_plate','Min_ctrl_zscore_for_plate','Unique_key','Unique_key_subplot','group']]
    small_df.loc[:,"Hit"] = np.where((small_df['Well_zscore'] > small_df['Max_ctrl_zscore_for_plate'])|
                                     (small_df['Well_zscore'] < small_df['Min_ctrl_zscore_for_plate']),
                                    "Hit","Not a hit")
    filtered_df = small_df[((small_df.Hit == 'Hit'))&((small_df.Relative_amplitude <= 1.5)&(small_df.Relative_amplitude >= 0.5))]
    filtered_df.loc[:,'Influence'] = np.where(filtered_df['Diff from ctrl avg'] > 0, "Stabilizer","Destabilizer")
    filtered_df_count = filtered_df.shape[0]
    filtered_stabilizer_count = filtered_df['Influence'].value_counts().get('Stabilizer', 0)
    filtered_destabilizer_count = filtered_df['Influence'].value_counts().get('Destabilizer', 0)
    filtered_df = filtered_df.sort_values('Well_zscore',ascending=False)
    new_cols = [col for col in filtered_df.columns if (col != 'Unique_key_subplot')&(col != 'Unique_key')] + ['Unique_key'] + ['Unique_key_subplot'] #Stick the unique key column at the end
    filtered_df = filtered_df[new_cols]
    print("Data read in. Let's goooo! ")

elif args.mode == 'only_tm':
    print("Reading in data. This can take a few minutes if processing many plates")
    df_results = pd.read_csv(args.input_dir+'/Only_Tm_values.txt', sep='\t', header = 0)
    df_results['Fraction'].fillna(0, inplace=True)
    df_curves = pd.read_csv(args.input_dir+'/Only_Tm_curves.txt', sep='\t', header = 0)
    columns_to_round = ['Smooth_Tm','Amplitude','Curve_height']
    df_results = df_results.assign(**{col: df_results[col].apply(lambda x: round(x, 2)) for col in columns_to_round})
    df_results['Row'] = df_results['Well'].str[0] #Generate Row column (needed to build visualizations)
    df_results['Column'] = (df_results['Well'].str[1:]).astype(int) #Generate Column column (needed to build visualizations)
    df_results['group'] = df_results['Assay_Plate'].apply(lambda x: df_results['Assay_Plate'].unique().tolist().index(x) // 16) #Assign group number to every 16 plates
    groups_df = df_results.loc[:, ['Assay_Plate', 'group']].drop_duplicates()
    print("Data read in. Let's goooo! ")

elif args.mode == 'dose_response':
    print("Reading in data. This can take a few minutes if processing many plates")
    df_results = pd.read_csv(args.input_dir+'/DR_final_results.txt', sep='\t', header = 0)
    df_results['Fraction'].fillna(0, inplace=True)
    df_curves = pd.read_csv(args.input_dir+'/DR_final_curves.txt', sep='\t', header = 0)
    columns_to_round = ['Smooth_Tm','Amplitude','Curve_height']
    df_results = df_results.assign(**{col: df_results[col].apply(lambda x: round(x, 2)) for col in columns_to_round})
    df_results['Row'] = df_results['Well'].str[0] #Generate Row column (needed to build visualizations)
    df_results['Column'] = (df_results['Well'].str[1:]).astype(int) #Generate Column column (needed to build visualizations)
    df_results['group'] = df_results['Assay_Plate'].apply(lambda x: df_results['Assay_Plate'].unique().tolist().index(x) // 16) #Assign group number to every 16 plates
    groups_df = df_results.loc[:, ['Assay_Plate', 'group']].drop_duplicates()
    df_curves = pd.merge(df_curves,groups_df, how = 'outer', on='Assay_Plate')
    df_curves['Error'].fillna('', inplace=True)
    print("Data read in. Let's goooo! ")

#Load layout and associated callbacks
if args.mode == 'full':
    from shiftscan_viewer_functions import generate_default_hit_table, generate_distribution_graph,generate_selected_curve, generate_first_derivative_curve, generate_all_curves,generate_barplot,generate_scatterplot,generate_all_ctrls_graph,generate_Tm_heatmap,generate_zscore_heatmap,update_heatmap,full_viz_layout,full_callbacks
    app.layout = full_viz_layout(plate_report_df,pipette_df, df_results, df_curves, filtered_df, filtered_df_count, filtered_stabilizer_count, filtered_destabilizer_count)
    full_callbacks(app, plate_report_df,pipette_df, df_results, df_curves, filtered_df, filtered_df_count, filtered_stabilizer_count, filtered_destabilizer_count, args.mode)

elif args.mode == 'only_tm':
    from shiftscan_viewer_functions import generate_selected_curve, generate_first_derivative_curve, generate_Tm_heatmap,generate_zscore_heatmap,update_heatmap,only_tm_layout,tm_only_callbacks
    app.layout = only_tm_layout()
    tm_only_callbacks(app, df_results, df_curves, args.mode)

elif args.mode == 'dose_response':
	from shiftscan_viewer_functions import generate_selected_curve, generate_first_derivative_curve, generate_Tm_heatmap,generate_zscore_heatmap,update_heatmap, DR_layout, DR_callbacks
	app.layout = DR_layout(df_results)
	DR_callbacks(app, df_results, df_curves, args.mode)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=args.port, debug=False)