#!/usr/bin/env python3

import dash, math, argparse
from dash import Dash, html, dcc, Input, Output, State, dash_table, callback_context
import pandas as pd
import plotly.express as px
import plotly.colors as pcolors
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from DSF_functions import boltzmann_sigmoid, find_inflection

#color1 = '#B4EDD2'
color1 = '#B4EDD2'
color1_faded = '#e0fff1'
color2 = '#A0CFD3'
color3 = '#8D94BA'
color4 = '#9A7AA0'
color4_faded = '#d1b6d6' 
color5 = '#87677B'
main_header_bkgrd_color = '#e3e3e3'
graph_gray = '#cacacf'

base_style = {'border': '1px black','color': 'black','font-size': '14px','font-weight': 'bold','font-family': 'Arial','padding': '6px'}
base_selected_style = {'border': '3px solid black','color': 'black','font-size': '14px','font-weight': 'bold','font-family': 'Arial','padding':'6px'}
base_disabled = {'border': '1px black','color': 'dark grey','font-size': '14px','font-weight': 'bold','font-family': 'Arial','padding': '6px','background': 'grey'}

tab1_style = {**base_style,'background': color1}
tab1_selected = {**base_selected_style,'background': color1}

tab2_style = {**base_style,'background': color2}
tab2_selected = {**base_selected_style,'background': color2}

tab3_style = {**base_style,'background': color3}
tab3_selected = {**base_selected_style,'background': color3}

tab4_style = {**base_style,'background': color4}
tab4_selected = {**base_selected_style,'background': color4}

tab5_style = {**base_style,'background': color5}
tab5_selected = {**base_selected_style,'background': color5}

data_table_style_data_conditional = [{'if': {'state': 'active'},'backgroundColor': color4_faded,'border': '1px black',},{'if': {'state': 'selected'},'backgroundColor': color4_faded,'border': '1px black',},]
plate_report_style_data_conditional = [{'if': {'state': 'active'},'backgroundColor': color1_faded,'border': '1px black'},{'if': {'state': 'selected'},'backgroundColor': color1_faded,'border': '1px black'},]

def generate_default_hit_table(df_results):
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
    return filtered_df, filtered_df_count, filtered_stabilizer_count, filtered_destabilizer_count

def generate_distribution_graph(y_axis_choice, df_results):
    sub_df = df_results[(df_results['Well_type']=='Control')&(df_results['Final_decision']!='Removed')]
    tm_max = np.nanmax(sub_df['Smooth_Tm'])
    tm_min = np.nanmin(sub_df['Smooth_Tm'])
    figure=px.violin(sub_df, x='Assay_Plate', y=y_axis_choice,color_discrete_sequence =[color2], height = 200,\
                        labels={'Assay_Plate': 'Assay Plate', 'Smooth_Tm': 'Melting temp'},
                        hover_data={'Source_Plate':True})\
                        .update_xaxes(matches=None, tickfont=dict(size=6), tickangle=90)\
                        .update_yaxes(title = 'Tm', range=[tm_min-5,tm_max+5])\
                        .update_layout(margin=dict(l=20, r=20, t=20, b=20))
    return figure

def generate_selected_curve(selected_unique_key, df_curves, mode):
    selected_curve_data = df_curves[df_curves['Unique_key'] == selected_unique_key]
    plate = selected_unique_key.rsplit('_', 1)[0]
    well = selected_unique_key.rsplit('_', 1)[-1]
    original_selected = selected_curve_data[selected_curve_data['Subplot'] == 'Original']

    if mode == 'full':
        source_plate = selected_curve_data['Source_Plate'].unique()[0]
        # Filter subplots for those that passed
        subplots_selected = selected_curve_data[(selected_curve_data['Subplot'] != 'Original') & (selected_curve_data['Final_decision'] == 'Pass')]
        # Define plot title and hover data for the full mode
        plot_title = f'Plate {source_plate}, Well {well}<br><sup>Assay Plate: {plate}</sup>'
        hover_data_config = {'Final_decision': False, 'Assay_Plate': False, 'Temps': False, 'Smooth Fluorescence': False, 'Error': True}

    elif (mode == 'only_tm') | (mode == 'dose_response') :
        # Get all subplots
        subplots_selected = selected_curve_data[selected_curve_data['Subplot'] != 'Original']
        # Define plot title and hover data for the only_tm mode
        plot_title = f'Original curve, Well {well}<br><sup>Assay Plate: {plate}</sup>'
        hover_data_config = None # Use default hover data

    else:
        raise ValueError("Invalid mode specified. Choose 'full' or 'only_tm'.")

    fig = px.scatter(
        subplots_selected,
        x='Temps',
        y='Smooth Fluorescence',
        color='Subplot',
        color_discrete_sequence=[color3, color4, color2, color1],
        labels={'Temps': 'Temperature', 'Smooth Fluorescence': 'Normalized fluorescence'},
        hover_name='Well',
        hover_data=hover_data_config,
        height=500,
        title=plot_title
    )

    # Add the 'original' curve as a dotted line
    fig.add_scatter(
        x=original_selected['Temps'],
        y=original_selected['Smooth Fluorescence'],
        line={'color': graph_gray, 'dash': 'dot'},
        name='Original data'
    )

    fig.update_layout(title={'font': {'size': 15}})
    
    # Reverse plots
    fig.data = fig.data[::-1]
    
    return fig

def generate_first_derivative_curve(selected_unique_key, df_curves, df_results, mode):
    selected_curve_data = df_curves[df_curves['Unique_key'] == selected_unique_key]
    if mode == 'full':
        avg_ctrl_Tm = df_results.loc[df_results['Unique_key'] == selected_unique_key, 'Avg_ctrl_melting_temp'].iloc[0]
        well_Tms = df_results.loc[(df_results['Unique_key'] == selected_unique_key) & (df_results['Final_decision'] == 'Pass'), 'Final_Tm'].unique()
        subplots = selected_curve_data[(selected_curve_data['Subplot'] != 'Original') & (selected_curve_data['Final_decision'] == 'Pass')]
        plot_title = "First derivative <br><sup><span style='color:red'>Red:</span> Avg. ctrl Tm for plate, <span style='color:blue'>Blue:</span> Well Tm</sup>"
    elif (mode == 'only_tm') | (mode == 'dose_response'):
        well_Tms = df_results.loc[df_results['Unique_key'] == selected_unique_key, 'Smooth_Tm'].unique()
        subplots = selected_curve_data[selected_curve_data['Subplot'] != 'Original']
        plot_title = "First derivative <br><sup><span style='color:blue'>Blue:</span> Well Tm</sup>"
    else:
        raise ValueError("Invalid mode specified. Choose 'full' or 'only_tm'.")

    # Original curve derivative
    orig_curve_coords = selected_curve_data[selected_curve_data["Subplot"] == 'Original'].drop_duplicates()
    raw_x = orig_curve_coords["Temps"].values
    raw_y = orig_curve_coords["Smooth Fluorescence"].values
    y_grad = np.gradient(raw_y, raw_x)
    orig_deriv_df = pd.DataFrame({'Temps': raw_x, '1st_deriv': y_grad})
    
    # Subplots derivative
    subplots = subplots.drop_duplicates()
    subplots_x = subplots["Temps"].values
    subplots_y = subplots["Smooth Fluorescence"].values
    
    if len(subplots_y) == 0:
        return None

    subplots_subplots = subplots["Subplot"].values
    subplots_y_grad = np.gradient(subplots_y, subplots_x)
    subplot_deriv_df = pd.DataFrame({'Temps': subplots_x, '1st_deriv': subplots_y_grad, 'Subplot': subplots_subplots})
    #subplot_deriv_df.to_csv(selected_unique_key+"_first_deriv_coord.txt", sep="\t", index=False)

    # Draw the subplots
    fig = px.scatter(subplot_deriv_df, x='Temps', y='1st_deriv', color='Subplot',
                     color_discrete_sequence=[color3, color4, color2, color1],
                     labels={'Temps': 'Temperature', '1st_deriv': 'dF/dT'})

    # Add the original curve as a dotted line
    fig.add_scatter(x=orig_deriv_df['Temps'], y=orig_deriv_df['1st_deriv'],
                    line={'color': graph_gray, 'dash': 'dot'}, name='Original data')

    # Add vertical lines for well Tms
    for value in well_Tms:
        fig.add_vline(x=value, line={'color': 'blue', 'width': 1, 'dash': 'dash'}, name='Tm(s) for selected well')

    # Add control Tm line if in 'full' mode
    if mode == 'full':
        fig.add_vline(x=avg_ctrl_Tm, line={'color': 'red', 'width': 2, 'dash': 'dot'}, name='Avg. ctrl Tm for plate')

    fig.update_layout(title=plot_title, xaxis_title='Temperature', yaxis_title='dF/dT')
    
    # Reverse plots
    fig.data = fig.data[::-1]
    
    return fig

def generate_all_curves(unique_key_list, df_curves):
    hit_curve_data_df_raw = df_curves[df_curves['Unique_key'].isin(unique_key_list)]
    hit_curve_data_df = hit_curve_data_df_raw[~hit_curve_data_df_raw['Final_decision'].str.contains('Failed|Removed', na=False)]
    if len(unique_key_list) < 6: #Preventing squashed images
        plot_num = len(unique_key_list)
        all_curves_figure = px.scatter(hit_curve_data_df, x='Temps', y='Smooth Fluorescence', color='Subplot', color_discrete_sequence=[graph_gray,color3,color4,color2,color1],
                            hover_name='Well', hover_data={'Source_Plate':True, 'Final_decision':False, 'Assay_Plate':False, 'Temps':False, 'Smooth Fluorescence':False, 'Error':False, 'Ctrl_Tm_z-score':False}, #Hover data (Tooltip in Spotfire)
                            category_orders={'Unique_key': unique_key_list},
                            facet_col='Unique_key',facet_col_spacing=0.08)

    else:
        num_rows = math.ceil(len(unique_key_list)/10)
        height = 25*num_rows
        facet_row_max_spacing = 3/height
        all_curves_figure = px.scatter(hit_curve_data_df, x='Temps', y='Smooth Fluorescence', color='Subplot', color_discrete_sequence=[graph_gray,color3,color4,color2,color1],
                            hover_name='Well', hover_data={'Source_Plate':True, 'Final_decision':False, 'Assay_Plate':False, 'Temps':False, 'Smooth Fluorescence':False, 'Error':False, 'Ctrl_Tm_z-score':False}, #Hover data (Tooltip in Spotfire)
                            category_orders={'Unique_key': unique_key_list}, #Maintain row order for facte plots!
                            facet_col='Unique_key', facet_col_wrap=5, facet_col_spacing=0.08,facet_row_spacing = facet_row_max_spacing,#Facet plots by plate and only allow 5 columns. Column spacing had to be adjusted to allow for individual y-axes
                            render_mode = 'auto', height = height) #Height of plots is equal to half the number of plates (coz 2 columns) with each plot 300px high. Width will have to be adjusted
    for annotation in all_curves_figure.layout.annotations:
        text = annotation.text
        unique_key = text.split('=')[1]
        subset = hit_curve_data_df.loc[hit_curve_data_df['Unique_key'] == unique_key, ['Source_Plate', 'Well']].drop_duplicates()
        source = subset['Source_Plate'].iloc[0]
        well = subset['Well'].iloc[0]
        annotation.text = (f"Plate: {source}, Well: {well}")

    return all_curves_figure

def generate_barplot(df,x_axis, color_by):
    df_count = df.groupby(['Platename', 'Fraction']).size().reset_index(name='Count')
    if color_by == 'Fraction':
        if df_count['Fraction'].dtype == 'float64':
            df_count[color_by] = df_count[color_by].astype(int).astype(str)
        else:
            df_count[color_by] = df_count[color_by].astype(str)
    figure = px.bar(df_count, x=x_axis, y='Count', color=color_by, barmode='stack',
                    labels={'Platename': 'Plate'})
    if color_by == 'Fraction':
        figure.update_xaxes(tickfont=dict(size=8), tickangle=90)
    else:
        figure.update_xaxes(tickfont=dict(size=12))
    figure.update_layout(margin=dict(l=20, r=20, t=5, b=5), height=350) #Reduce enormous default margins
    return figure

def generate_scatterplot(hits_df):
    category_order = ['Control', 'Hit', 'Not a hit']
    hits_df['Hit'] = pd.Categorical(hits_df['Hit'], categories=category_order, ordered=True)
    figure = px.strip(hits_df, x = 'Platename',y = 'Final_Tm', color = 'Hit', color_discrete_map={'Control':color2,'Not a hit':color3,'Hit':color4},
                      category_orders={'Hit': category_order},
                      labels={'Final_Tm': 'Melting temp', 'Platename': 'Plate'},
                      hover_data={'Well':True, 'Well_zscore':True, 'Relative_amplitude':True,'Avg_ctrl_melting_temp':True,'Diff from ctrl avg':True, 'Std. devs from ctrl mean':True,'Max_ctrl_zscore_for_plate':True,'Min_ctrl_zscore_for_plate':True})
    figure.update_xaxes(tickfont=dict(size=8), tickangle=90)
    figure.update_layout(margin=dict(l=20, r=20, t=5, b=5), height=350) #Reduce enormous default margins
    figure.update_traces(marker_size=10)
    return figure

def generate_all_ctrls_graph(df):
    no_plates = df['Assay_Plate'].nunique()
    num_rows = math.ceil(no_plates/2)
    height = 25*num_rows
    facet_row_max_spacing = 6/height
    figure=px.scatter(df, x='Temps', y='Smooth Fluorescence', color='Final_decision', color_discrete_map={'Pass':color2,'Failed':color4},\
                                hover_name='Well', hover_data={'Source_Plate':True, 'Final_decision':False, 'Assay_Plate':False, 'Temps':False, 'Smooth Fluorescence':False, 'Error':True, 'Ctrl_Tm_z-score':True},\
                                facet_col='Assay_Plate', facet_col_wrap=2, facet_col_spacing=0.06, facet_row_spacing=facet_row_max_spacing, render_mode = 'auto')
    figure.update_yaxes(matches=None)
    figure.update_xaxes(matches=None)
    figure.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
    figure.for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True))
    figure['layout']['height'] = 250*num_rows
    return figure

def generate_Tm_heatmap(df, facet_row_max_spacing, mode):
    if mode == 'full':
        figure=px.scatter(df, x='Column', y='Row', color='Final_Tm', hover_name='Compound', 
                                    color_continuous_scale=[color1,color2,color3,color4], color_continuous_midpoint=50,
                                    hover_data={'Platename':False, 'Fraction':True, 'Final_Tm':True,'Column':False, 'Row':False, 'Well':True, 'Error':True,'Unique_key':True, 'Well_zscore':True},
                                    facet_col='Platename', facet_col_wrap=1, facet_row_spacing=facet_row_max_spacing, 
                                    render_mode = 'auto')
    elif mode == 'only_tm':
        figure=px.scatter(df, x='Column', y='Row', color='Smooth_Tm', hover_name='Compound', 
                                    color_continuous_scale=[color1,color2,color3,color4], color_continuous_midpoint=50,
                                    hover_data={'Source_Plate':False, 'Fraction':True, 'Smooth_Tm':True,'Column':False, 'Row':False, 'Well':True, 'Error':True,'Unique_key':True},
                                    facet_col='Source_Plate', facet_col_wrap=1, facet_row_spacing=facet_row_max_spacing, 
                                    render_mode = 'auto')
    if mode == 'dose_response':
        figure=px.scatter(df, x='Column', y='Row', color='Final_Tm', hover_name='Compound', 
                                    color_continuous_scale=[color1,color2,color3,color4], color_continuous_midpoint=50,
                                    hover_data={'Source_Plate':False, 'Concentration':True, 'Replicate':True, 'Fraction':True, 'Smooth_Tm':True,'Column':False, 'Row':False, 'Well':True, 'Error':True,'Unique_key':True},
                                    facet_col='Source_Plate', facet_col_wrap=1, facet_row_spacing=facet_row_max_spacing, 
                                    render_mode = 'auto')
    return figure

def generate_zscore_heatmap(df, facet_row_max_spacing):
    figure=px.scatter(df, x='Column', y='Row', color='Well_zscore', hover_name='Compound', 
                                color_continuous_scale=px.colors.diverging.Spectral_r, color_continuous_midpoint=0,
                                hover_data={'Platename':False, 'Fraction':True, 'Final_Tm':True,'Column':False, 'Row':False, 'Well':True, 'Error':True,'Unique_key':True, 'Well_zscore':True}, 
                                facet_col='Platename', facet_col_wrap=1, facet_row_spacing=facet_row_max_spacing, 
                                render_mode = 'auto')
    return figure

def update_heatmap(figure,facet_row_max_spacing,num_rows,size_option):
    figure.update_traces(marker={'symbol': 'square', 'size':size_option, 'line':{'color':'grey','width':1}})
    figure.update_yaxes(autorange='reversed',showgrid=False, zeroline=False) #Get y-axis to go from A to Z, and remove all gridlines
    figure.update_xaxes(showgrid=False, zeroline=False) #Remove all gridlines
    figure.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)'}) #Make plot backgrounds transparent
    figure.update_layout(coloraxis_colorbar = {'len':0.5, 'xanchor':'center', 'yanchor':'bottom','y':1+facet_row_max_spacing,'orientation':'h','ypad':0.8,'thicknessmode':'pixels', 'thickness':20})
    figure.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True)) #Put y-axis labels on all facet plots
    figure.for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True)) #Put x-axis labels on all facet plots
    tab3_height = 500*num_rows
    figure['layout'].update(height = tab3_height,margin=dict(l=10, r=10, t=10, b=10), width = 700)
    return figure

def generate_indiv_Tm_graph(df):
    df['Replicate'] = df['Replicate'].astype(str)
    sub_df = df.dropna(subset=['Final_Tm'])
    fig = go.Figure() #Create empty fig
    color_map = {'1': color1, '2': color2, '3': color3, '4': color5}

    for replicate in sorted(sub_df['Replicate'].unique()):
        replicate_df = sub_df[sub_df['Replicate'] == replicate]
        fig.add_trace(go.Scatter(
            x=replicate_df['Concentration'],
            y=replicate_df['Final_Tm'],
            mode='markers',
            marker=dict(
                color=color_map.get(replicate, 'black'),
                size=8),
            name=f'Replicate {replicate}',
            customdata=replicate_df[['Compound', 'Well', 'zscore', 'Relative_amplitude']],
            hovertemplate=(
                "<b>Compound</b>: %{customdata[0]}<br>" +
                "<b>Concentration</b>: %{x}<br>" +
                "<b>Melting Temp</b>: %{y:.2f}°C<br>" +
                "<b>Well</b>: %{customdata[1]}<br>" +
                "<b>Z-Score</b>: %{customdata[2]:.2f}<br>" +
                "<b>Relative Amplitude</b>: %{customdata[3]:.2f}" +
                "<extra></extra>"
            )
        ))

    # Update the layout of the figure
    fig.update_layout(
        xaxis_title="Concentration",
        yaxis_title="Melting Temperature (°C)",
        margin=dict(l=10, r=10, t=10, b=10),
        height=350,
        legend_title="Replicates",
        legend=dict(orientation="h",yanchor="top",y=-0.2,xanchor="center",x=0.5))
    fig.update_xaxes(autorange="reversed")

    return fig

def generate_delta_Tm_graphs(df,log_it_option):
    df['Replicate'] = df['Replicate'].astype(str)
    sub_df = df.dropna(subset=['Final_Tm']).reset_index()
    control_series = sub_df[sub_df['Concentration'] == 0]['Final_Tm']

    if not control_series.empty:
        ctrl_mean = control_series.values[0]
    else:
        ctrl_mean = 0 
        print("Warning: No control with Concentration=0 found. Delta Tm may be inaccurate.")
    
    sub_df['Delta_tm'] = sub_df['Final_Tm'] - ctrl_mean
    agg_df = sub_df.groupby('Concentration')['Delta_tm'].agg(['mean', 'min', 'max']).reset_index()

    if log_it_option == 'Conc':

        conc_fig = go.Figure() #Create empty fig

        # Add the ribbon
        conc_fig.add_trace(go.Scatter(
            x=agg_df['Concentration'],
            y=agg_df['min'],
            mode='lines',
            line={'width':0},
            showlegend=False))

        conc_fig.add_trace(go.Scatter(
            x=agg_df['Concentration'],
            y=agg_df['max'],
            mode='lines',
            line={'width':0},
            fill='tonexty',
            fillcolor=color4_faded,
            showlegend=False))

        # Plot the deltas
        conc_fig.add_trace(go.Scatter(
            x=agg_df['Concentration'],
            y=agg_df['mean'],
            line={'color':color4},
            showlegend=False))
        
        conc_fig.update_layout(xaxis_title="Conc",yaxis_title="Δ Tm from ctrl (°C)", margin=dict(l=10, r=10, t=10, b=10), height = 300)
        conc_fig.update_xaxes(autorange="reversed")

        return conc_fig

    else:
        #Finding log and bypassing negative values (np.where still sees them)
        agg_df['log(Conc)'] = np.nan
        positive_mask = agg_df['Concentration'] > 0
        agg_df.loc[positive_mask, 'log(Conc)'] = np.log(agg_df.loc[positive_mask, 'Concentration']) #Log masked rows only
        log_fig = go.Figure() #Create empty fig
        
        #Curve fitting and makeing pretty rep curve
        log_conc_x = agg_df['log(Conc)'].dropna().tolist()
        mean_delta_tm_y = agg_df[agg_df['log(Conc)'].notna()]['mean'].tolist()
        popt, _ = curve_fit(boltzmann_sigmoid, log_conc_x, mean_delta_tm_y)
        A_opt, B_opt, C_opt, D_opt = popt
        model_x_plot = np.linspace(min(log_conc_x), max(log_conc_x), 100)
        model_y_plot = boltzmann_sigmoid(model_x_plot, A_opt, B_opt, C_opt, D_opt)
        model_y = boltzmann_sigmoid(log_conc_x, A_opt, B_opt, C_opt, D_opt)
        closest_x_value = min(model_x_plot, key=lambda x: abs(x - C_opt))
        closest_index = model_x_plot.tolist().index(closest_x_value)
        inflection_y = model_y_plot.tolist()[closest_index]

        #Goodness of fit stuff
        ssr = np.sum((mean_delta_tm_y - model_y) ** 2) #Sum of Squared Residuals
        sst = np.sum((mean_delta_tm_y - np.mean(mean_delta_tm_y)) ** 2) #Total Sum of Squares
        r_squared = 1 - (ssr/sst)


        # Plot model
        log_fig.add_trace(go.Scatter(
            x=model_x_plot,
            y=model_y_plot,
            line={'color':color2},
            showlegend=True,
            name = 'Model'))

        # Plot data
        log_fig.add_trace(go.Scatter(
            x=agg_df['log(Conc)'],
            y=agg_df['mean'],
            line={'color':color4},
            showlegend=True,
            name = 'Actual data'))

        # Plot inflection
        log_fig.add_trace(go.Scatter(
            x=[closest_x_value],
            y=[inflection_y],
            mode='markers+text',
            marker={'color':'red', 'size':10},
            text=[f"Infl: {closest_x_value:.4f}<br>SSR: {ssr:.4f}"],
            textposition="middle right",
            showlegend=False))

        log_fig.update_layout(xaxis_title="log(Conc)",
                                yaxis_title="Δ Tm from ctrl (°C)", 
                                margin=dict(l=10, r=10, t=10, b=10), 
                                height = 350,
                                legend=dict(orientation="h",yanchor="top",y=-0.2,xanchor="center",x=0.5))
        log_fig.update_xaxes(autorange="reversed")
        return log_fig

def get_curves_for_avg(df_curves,select_unique_keys_df):
    multi_curve_df = pd.merge(df_curves, select_unique_keys_df, on = 'Unique_key',how='inner')
    #Need to remove points where not all subplots have the same x-axis values (i.e. only one replicate may be presented and creates a messy viz)
    #First, Count how many occurnaces of Smooth Fluoro you have per 'Subplot', 'Concentration', 'Temps' combination
    multi_curve_df['Temp_count'] = multi_curve_df.groupby(['Subplot', 'Concentration', 'Temps'])['Smooth Fluorescence'].transform('count')
    #Next, count how many occurances there are in the controls 
    orig_counts = multi_curve_df[multi_curve_df['Subplot'] == 'Original'].groupby(['Concentration', 'Temps'])['Smooth Fluorescence'].count().rename('Orig_count').reset_index()
    #If the number of temp counts is less than the occurnaces in the controls, boot the row
    multi_curve_df = pd.merge(multi_curve_df, orig_counts, on=['Concentration', 'Temps'], how='left')
    multi_curve_df_avg = multi_curve_df[multi_curve_df['Temp_count'] >= multi_curve_df['Orig_count']].groupby(['Subplot', 'Concentration', 'Temps'])['Smooth Fluorescence'].mean().reset_index()
    return multi_curve_df_avg

def generate_avg_curves(avg_curve_df):
    avg_curve_df = avg_curve_df.rename(columns={'Smooth Fluorescence': 'Avg. Smooth Fluorescence'})
    original_df = avg_curve_df[avg_curve_df['Subplot'] == 'Original']
    subplots_df = avg_curve_df[avg_curve_df['Subplot'] != 'Original'].sort_values('Temps')

    fig = go.Figure()

    #Get colors for concs
    concentrations = original_df['Concentration'].unique()
    continuous_scale = pcolors.make_colorscale([color1,color2,color3,color4,color5])
    colors_required = len(concentrations)
    colors = pcolors.sample_colorscale(continuous_scale,colors_required,low=0.0,high=1.0)
    color_map = {conc: colors[i % len(colors)] for i, conc in enumerate(concentrations)}

    #Trace of 'Original' subplots as dotted lines
    for conc in concentrations:
        df_filtered = original_df[original_df['Concentration'] == conc]
        fig.add_trace(go.Scatter(
            x=df_filtered['Temps'],
            y=df_filtered['Avg. Smooth Fluorescence'],
            mode='lines',
            line_dash = 'dot',
            line=dict(color=color_map[conc], width=2),
            name=f'{conc}' # Legend entry for the line
        ))

    #Add in subplt dot plots
    for conc in concentrations:
        df_filtered = subplots_df[subplots_df['Concentration'] == conc]
        fig.add_trace(go.Scatter(
            x=df_filtered['Temps'],
            y=df_filtered['Avg. Smooth Fluorescence'],
            mode='lines',
            line=dict(color=color_map[conc], width=4),
            name=f'{conc}', # Legend entry for the dots
            showlegend=False 
        ))
    
    fig.update_layout(xaxis_title="Temp (°C)", yaxis_title="Normalized_Fluorescence",margin=dict(l=10, r=10, t=10, b=10), height = 300)
    return fig

def full_viz_layout(plate_report_df,pipette_df, df_results, df_curves, filtered_df, filtered_df_count, filtered_stabilizer_count, filtered_destabilizer_count):
    full_layout = html.Div([
        html.H1('ShiftScan Viewer', style = {'background': main_header_bkgrd_color,'font-family': 'Arial','padding':'20px'}),
        dcc.Tabs([
                #Result overview 
                dcc.Tab([
                    html.Hr(style = {'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}), #Break
                    #Plate report table
                    html.Div([
                        #Header
                        html.H3('Plate report:',style = {'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'0px'}),
                        #Table
                        dash_table.DataTable(
                            data = plate_report_df.to_dict('records'),
                            columns = [{'name': i, 'id': i} for i in plate_report_df.columns],
                            id = 'plate_report_table',
                            sort_action='native',
                            export_format='xlsx',
                            style_data_conditional=plate_report_style_data_conditional, 
                            style_as_list_view=True, style_cell={'fontSize':12, 'font-family':'Arial'}, style_header = {'backgroundColor': color1,'fontWeight': 'bold'})
                        ], style = {'display':'inline-block','vertical-align':'top','width':'75%', 'height': '80%','margin-top':'30px','overflow-y':'scroll', 'margin-right': '50px'}),
                    #Pipetting problems
                    html.Div([
                        #Header
                        html.H3('Possible pipetting issues:',style = {'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'0px'}),
                        #Table
                        dash_table.DataTable(
                                data = pipette_df.to_dict('records'),
                                columns = [{'name': i, 'id': i} for i in pipette_df.columns],
                                id = 'pipette_issues_table',
                                export_format='xlsx',
                                style_data_conditional=plate_report_style_data_conditional, 
                                style_as_list_view=True, style_cell={'fontSize':12, 'font-family':'Arial'}, style_header = {'backgroundColor': color1,'fontWeight': 'bold'})
                        ], style = {'display':'inline-block','vertical-align':'top','width':'20%', 'height': '200px','margin-top':'30px','overflow-y':'scroll'})
                    ],
                label='Plate overview', style=tab1_style, selected_style=tab1_selected, value='tab-1'), #Plate overview tab style

                #CONTROL TAB
                dcc.Tab([
                    html.Hr(style = {'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}), #Break
                    #All ctrl Tm div
                    html.Div([
                            #Header
                            html.H4('Distribution of all control Tms per plate:',style = {'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'0px'}),
                            #Boxplot
                            html.Div([dcc.Loading(dcc.Graph(id = 'all_control_boxes', figure = generate_distribution_graph('Smooth_Tm', df_results),
                                style = {'display':'inline-block','vertical-align':'top','width':'100%','height':200})
                                ,type = 'dot',color = color2)],
                            id='curve_all_boxplots')
                            ],
                        style = {'display':'inline-block','vertical-align':'top','width':'50%','margin-top':'20px','overflow-y':'scroll'}),
                    #Clean ctrl tm div
                    html.Div([
                            #Header
                            html.H4('Distribution of "clean" control Tms per plate:',style = {'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'0px'}),
                            #Boxplot
                            html.Div([dcc.Loading(dcc.Graph(id = 'clean_control_boxes', figure = generate_distribution_graph('Final_Tm',df_results),
                                style = {'display':'inline-block','vertical-align':'top','width':'100%','height':200})
                                ,type = 'dot',color = color2)],
                            id='curve_clean_boxplots')
                            ],
                        style = {'display':'inline-block','vertical-align':'top','width':'50%','margin-top':'20px','overflow-y':'scroll'}), 
                    html.Hr(style = {'margin-bottom':'0px','width': '100%'}), #Break 
                    dcc.Dropdown(id = 'Choice of controls',
                                options=[{'label': 'Plates: 1 - 16', 'value': 0}],
                                value=0,
                                style = {'display':'inline-block','vertical-align':'top','font-family': 'Arial', 'margin-bottom':'0px','margin-top':'5px','fontSize':12, 'width': '50%'}),          
                    #Curves
                    html.H4('Control curves:',style = {'font-family': 'Arial','margin-bottom':'5px', 'margin-top':'5px'}),
                    html.Div([dcc.Loading(dcc.Graph(id = 'control_curves',
                        style = {'display':'inline-block','vertical-align':'top', 'width': '100%','overflow-y':'scroll'})
                    ,type = 'dot',color = color2)],id='curve_graphs')],
                label='Control overview', style=tab2_style, selected_style=tab2_selected, value ='tab-2'),

                #MELTING TEMP/Z-SCORE PLATE OVERVIEW
                dcc.Tab([
                    html.Hr(style = {'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}), #Break
                    html.Div([
                            html.H3('Color plates by:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px','display':'inline-block','vertical-align':'top','width':'50%'}),
                            html.H3('Plate range:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px','display':'inline-block','vertical-align':'top','width':'50%'},id = 'plate_range_header'),
                            dbc.RadioItems(id = 'color_filter', 
                                            options = [{'label':'Melting temp', 'value':'Final_Tm'},
                                                        {'label':'Z-score', 'value':'Well_zscore'}], 
                                                        value = 'Final_Tm',
                                                        style = {'display':'inline-block','vertical-align':'top','font-family': 'Arial', 'margin-bottom':'5px','margin-top':'5px','fontSize':12, 'width': '50%'},
                                                        inline = True),
                            dcc.Dropdown(id = 'Choice of plates',
                                            options=[{'label': 'Plates: 1 - 50', 'value': 0}],
                                                    value=0,
                                                    style = {'display':'inline-block','vertical-align':'top','font-family': 'Arial', 'margin-bottom':'5px','margin-top':'5px','fontSize':12, 'width': '50%'}),
                            dcc.Loading(dcc.Graph(id = 'Tm_plates'),type = 'dot', color = color2),
                                ], style = {'display':'inline-block','vertical-align':'top', 'width': '55%','overflow-y':'scroll', 'height':'600px'}),
                   #Child2: Selected curve (Right of screen)
                    html.Div([
                        html.Div([],id='first_deriv_popup',style = {'margin-bottom':'0px'}),
                        html.Div([],id='popup',style = {'margin-top':'5px'})
                        ],id='graphs_block',style = {'display':'inline-block','vertical-align':'top', 'width': '45%', 'height':'200px'}),
                    ],label='Melting temp overview', style=tab3_style, selected_style=tab3_selected, value ='tab-3'), 

                #HIT TABLE TAB
                dcc.Tab([
                    html.Hr(style = {'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}), #Break
                    #Child A: Data table 
                    html.Div([
                            #Data table
                            html.Div([
                            html.Label("Total hits: "+str(filtered_df_count), style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'font-size': '20px', 'color':'dark green'}, id = 'hits_header'),
                            html.Hr(style = {'margin-bottom':'0px','width': '2%', 'margin-top':'2px', 'border':'0px white'}), #Break
                            html.Label("Stabilizers: "+str(filtered_stabilizer_count), style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'font-size': '16px'}, id = 'stabilizer_header'),
                            html.Hr(style = {'margin-bottom':'0px','width': '2%', 'margin-top':'2px', 'border':'0px white'}), #Break
                            html.Label("Destabilizers: "+str(filtered_destabilizer_count), style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'font-size': '16px'}, id = 'destabilizer_header'),
                            html.Hr(style = {'margin-bottom':'0px','width': '2%', 'margin-top':'2px', 'border':'15px white'}), #Break
                            dash_table.DataTable(
                                data = filtered_df.to_dict('records'),
                                columns = [{'name': i, 'id': i} for i in filtered_df.columns],
                                id = 'results_table_datatable',
                                hidden_columns=['Max_ctrl_zscore_for_plate','Hit','Min_ctrl_zscore_for_plate','group'], 
                                css=[{'selector': '.show-hide', 'rule': 'display: none'},{'selector':'.export','rule': 'margin:5px'}],
                                row_deletable=True,
                                sort_action='custom',
                                sort_mode='single',
                                sort_by=[],
                                export_format='xlsx',
                                style_data_conditional=data_table_style_data_conditional, 
                                style_table = {'height':'250px','overflow-y':'scroll'},
                                style_as_list_view=True, style_cell={'fontSize':12, 'font-family':'Arial'}, style_header = {'backgroundColor': color4,'fontWeight': 'bold'}),#Styling of table
                            ], id = 'results_table'),
                            html.Hr(style = {'margin-bottom':'0px','width': '2%', 'margin-top':'2px', 'border':'15px white'}), #Break
                            #Generate all hit graphs button
                            html.Div([html.Button('Generate all curves', id='generate', n_clicks=None, style = {'width':'200px','margin-top':'10px', 'margin-right': '20px', 'backgroundColor': color1})], style = {'display':'inline-block','vertical-align':'top'}),
                            #Generate summary graphs button
                            html.Div([html.Button('Generate summary graphs', id='summary_graphs', n_clicks=None, style = {'width':'200px','margin-top':'10px', 'margin-right': '20px', 'backgroundColor': color2})],style = {'display':'inline-block','vertical-align':'top'}),
                            #Clear all graphs button
                            html.Div([html.Button('Clear', id='clear', n_clicks=None, style = {'width':'200px','margin-top':'10px', 'margin-right': '80px', 'backgroundColor': color3})],style = {'display':'inline-block','vertical-align':'top'})
                            ], 
                        id='left_panel',
                        style = {'display':'inline-block','vertical-align':'top', 'width': '55%','overflow-y':'scroll','overflow-x':'scroll', 'height':'50%', 'padding':'10px','margin-top':'5px'}), #Styling of table container
                    
                    #Child B: Cut off boxes and pop up graph
                    html.Div([
                        #Child1
                        html.Div([
                                    html.H5('Hit defined by:',style = {'font-family': 'Arial','margin-bottom':'5px', 'margin-top':'0px', 'margin-right':'2px'}),
                                    dcc.RadioItems(id = 'hit_type', 
                                                        options = [
                                                                {'label':'Z-score distribution', 'value':'zscore'},
                                                                {'label':'Std dev from ctrl mean', 'value':'stddev'},
                                                                {'label':'Degrees', 'value':'degrees'},
                                                                ], value = 'zscore', inline = True,
                                                                style = {'font-family': 'Arial', 'margin-bottom':'15px','margin-top':'7px','fontSize':12}),
                                    html.Hr(style = {'margin-bottom':'5px','width': '100%'}),
                                ]),
                        #Child 2: z-score options
                        html.Div([
                                    #Child 1.1: Main z-score header
                                    html.H5('% beyond ctrl z-score range:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'margin-right':'2px'}, id = 'zscore_header'),#
                                    #Child 1.2: Upper z-score cutoff header and box
                                    html.Div([
                                            html.H6('Cutoff:',style = {'font-family': 'Arial','margin-bottom':'5px','margin-top':'7px'}, id = 'zscore_cutoff_header'),#
                                            dcc.Input(id='relative_zscore_cutoff', name = 'Cutoff', type='number', value = 0,style = {'font-family': 'Arial','width':'50px','margin-top':'5px'})],#
                                        style = {'display':'inline-block','vertical-align':'top', 'margin-top':'0px', 'margin-bottom':'5px'}),#Child 1.2 style
                                ],style = {'display':'inline-block','vertical-align':'top', 'width':'33%', 'margin-top':'0px', 'margin-bottom':'5px'}), #Child 2 style
                        #Child 3: Std dev options
                        html.Div([
                                    #Child 1.1: Main stddev header
                                    html.H5('No. of std devs from ctrl mean:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'margin-right':'2px','color':'grey'}, id = 'stddev_header'),#
                                    #Child 1.2: Stddev cutoff header and box
                                    html.Div([
                                            html.H6('Cutoff:',style = {'font-family': 'Arial','margin-bottom':'5px','margin-top':'7px','color':'grey'}, id = 'stddev_cutoff_header'),#
                                            dcc.Input(id='relative_stddev_cutoff', name = 'Cutoff', type='number', value = 3, disabled=True,style = {'font-family': 'Arial','width':'50px','margin-top':'5px'})],#
                                        style = {'display':'inline-block','vertical-align':'top', 'margin-top':'0px', 'margin-bottom':'5px'}),#Child 1.2 style
                                ],style = {'display':'inline-block','vertical-align':'top', 'width':'33%', 'margin-top':'0px', 'margin-bottom':'5px'}), #Child 3 style
                        #Child 4: Degrees options
                        html.Div([
                                    #Child 1.1: Main degree header
                                    html.H5('Degrees from avg. ctrl Tm:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'margin-right':'2px','color':'grey'}, id = 'degree_header'),#
                                    #Child 1.2: Upper z-score cutoff header and box
                                    html.Div([
                                            html.H6('Cutoff:',style = {'font-family': 'Arial','margin-bottom':'5px','margin-top':'7px','color':'grey'}, id = 'degree_cutoff_header'),#
                                            dcc.Input(id='relative_degree_cutoff', name = 'Cutoff', type='number', value = 2, disabled=True,style = {'font-family': 'Arial','width':'50px','margin-top':'5px'})],#
                                        style = {'display':'inline-block','vertical-align':'top', 'margin-top':'0px', 'margin-bottom':'5px'}),#Child 1.2 style
                                ],style = {'display':'inline-block','vertical-align':'top', 'width':'33%', 'margin-top':'0px', 'margin-bottom':'5px'}), #Child 4 style
                        html.Hr(style = {'margin-bottom':'5px','width': '100%'}),
                        #Child 3
                        html.Div([
                                #Child 2.1: Main amplitude header
                                html.H5('Relative amplitude cutoffs:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px'}),
                                #html.H6('*The relative amplitude of the hit must be between these two values',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px'}),
                                #Child 2.2.: Upper amplitude cutoff header and box
                                html.Div([
                                        html.H6('Upper cutoff:',style = {'font-family': 'Arial','margin-bottom':'5px','margin-top':'7px'}),#
                                        dcc.Input(id='upper_limit_amp', name = 'Upper cutoff', type='number', value = 1.5,style = {'font-family': 'Arial','width':'50px','margin-top':'5px'})],#
                                    style = {'display':'inline-block','vertical-align':'top', 'margin-top':'0px', 'margin-bottom':'5px'}),#Child 2.2 style
                                #Child 2.3: Lower amplitude cutoff header and box
                                html.Div([
                                        html.H6('Lower cutoff:',style = {'font-family': 'Arial', 'margin-bottom':'5px','margin-top':'7px'}),#
                                        dcc.Input(id='lower_limit_amp', name = 'Lower cutoff', type='number', value = 0.5, style = {'font-family': 'Arial', 'width':'50px','margin-top':'5px'})],#
                                    style = {'display':'inline-block','vertical-align':'top', 'margin-top':'0px', 'margin-bottom':'5px','margin-left': '30px'}), #Child 2.3 style
                                ],style = {'display':'inline-block','vertical-align':'top', 'width':'40%', 'margin-top':'0px', 'margin-bottom':'5px'}), #Child 2 style
                        #Child 4
                        html.Div([
                                html.H5('Interaction type:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px'}),
                                dbc.RadioItems(id = 'hit_filter', 
                                                options = [
                                                        {'label':'All hits', 'value':'all'},
                                                        {'label':'Only stabilizers', 'value':'stabilizers'},
                                                        {'label':'Only destabilizers', 'value':'destabilizers'},
                                                        ], value = 'all',
                                                        style = {'font-family': 'Arial', 'margin-bottom':'5px','margin-top':'7px','fontSize':12})
                                ],
                                style = {'display':'inline-block','vertical-align':'top', 'width':'30%', 'margin-top':'0px', 'margin-bottom':'5px'}),
                        #Submit button
                        #Child 5: Submit button
                        html.Div([html.Button('Update', id='submit', n_clicks=0)]),
                        #Child 6: Horizontal break
                        html.Hr(style = {'margin-bottom':'0px','width': '100%'}),
                        #Child 7: POPUP GRAPH
                        html.Div([
                            html.Div([],id = 'popup3',style = {'display':'inline-block','vertical-align':'top', 'width':'50%'}),
                            html.Div([],id = 'popup4',style = {'display':'inline-block','vertical-align':'top', 'width':'50%'}),
                            html.Div([html.Button('Clear', id='clear_popup', n_clicks=0)],id='clear_pop_up_div',hidden=True)
                            ],id='popup3_div',style = {'width': '100%', 'margin-top':'0px', 'padding':'5px'})
                            ],
                        id='control_boxes', #Child B ID
                        style = {'display':'inline-block','vertical-align':'top', 'width': '40%', 'overflow-x':'scroll', 'margin-top': '20px'}), #Child B style
                    
                    #Child C: All graphs for hits
                    dcc.Store(id='current_table_data_store'), # Use dcc.Store to pass data
                    html.Div(id='hits_dropdown_container', children=[
                        html.Div([
                            html.Label('Select hit range:'),
                            dcc.Dropdown(id='hits_dropdown', options=[], value=0, clearable=False,
                                         style={'width': '150px', 'margin-top': '5px','fontSize':12})
                        ], style={'font-family': 'Arial', 'margin-bottom': '5px'})
                    ], style={'display': 'none'}), # Initially hide the dropdown container
                    html.Div(id='all_graphs_div', style={'width': '95%'})
                    ],label='Hit list', style=tab4_style, selected_style=tab4_selected, value ='tab-4'), #Tab div
                ],id ='tabs',value = 'tab-0',vertical =False, style = {'width':'100%'}) #All tabs style and other details
    ],id = 'page')
    return full_layout

def only_tm_layout():
    only_tm_layout = html.Div([
        html.H1('ShiftScan Viewer', style = {'background': main_header_bkgrd_color,'font-family': 'Arial','padding':'20px'}),
        dcc.Tabs([
                #Result overview 
                dcc.Tab([],
                label='Plate overview', style=tab1_style, selected_style=tab1_selected, disabled_style = base_disabled, value='tab-1', disabled=True), #Plate overview tab style

                #CURVE TAB
                dcc.Tab([],
                label='Control overview', style=tab2_style, selected_style=tab2_selected, disabled_style = base_disabled, value ='tab-2', disabled=True),

                #MELTING TEMP/Z-SCORE PLATE OVERVIEW
                dcc.Tab([
                    html.Hr(style = {'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}), #Break
                    html.Div([
                            html.H3('Color plates by:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px','display':'inline-block','vertical-align':'top','width':'50%'}),
                            html.H3('Plate range:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px','display':'inline-block','vertical-align':'top','width':'50%'},id = 'plate_range_header'),
                            dcc.Dropdown(id = 'Choice of plates',
                                            options=[{'label': 'Plates: 1 - 50', 'value': 0}],
                                                    value=0,
                                                    style = {'display':'inline-block','vertical-align':'top','font-family': 'Arial', 'margin-bottom':'5px','margin-top':'5px','fontSize':12, 'width': '50%'}),
                            dcc.Loading(dcc.Graph(id = 'Tm_plates'),type = 'dot', color = color2),
                                ], style = {'display':'inline-block','vertical-align':'top', 'width': '50%','overflow-y':'scroll', 'height':'600px'}),
                   #Child2: Selected curve (Right of screen)
                    html.Div([
                        html.Div([],id='first_deriv_popup',style = {}),
                        html.Div([],id='popup',style = {})
                        ],id='graphs_block',style = {'display':'inline-block','vertical-align':'top', 'width': '50%'}),
                    ],label='Melting temp overview', style=tab3_style, selected_style=tab3_selected, value ='tab-3'), 

                #HIT TABLE TAB
                dcc.Tab([],label='Hit list', style=tab4_style, selected_style=tab4_selected, disabled_style = base_disabled, value ='tab-4', disabled=True), #Tab div
                ],id ='tabs',value = 'tab-3',vertical =False, style = {'width':'100%'}) #All tabs style and other details
    ],id = 'page')
    return only_tm_layout

def DR_layout(df_results):
    dr_layout = html.Div([
        html.H1('ShiftScan Viewer', style = {'background': main_header_bkgrd_color,'font-family': 'Arial','padding':'20px'}),
        dcc.Tabs([
                #Result overview 
                dcc.Tab([],
                label='Plate overview', style=tab1_style, selected_style=tab1_selected, disabled_style = base_disabled, value='tab-1', disabled=True), #Plate overview tab style

                #CONTROL TAB
                dcc.Tab([
                    html.Hr(style = {'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}), #Break
                    #All ctrl Tm div
                    html.Div([
                            #Header
                            html.H4('Distribution of all control Tms per plate:',style = {'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'0px'}),
                            #Boxplot
                            html.Div([dcc.Loading(dcc.Graph(id = 'all_control_boxes', figure = generate_distribution_graph('Smooth_Tm', df_results),
                                style = {'display':'inline-block','vertical-align':'top','width':'100%','height':200})
                                ,type = 'dot',color = color2)],
                            id='curve_all_boxplots')
                            ],
                        style = {'display':'inline-block','vertical-align':'top','width':'50%','margin-top':'20px','overflow-y':'scroll'}),
                    #Clean ctrl tm div
                    html.Div([
                            #Header
                            html.H4('Distribution of "clean" control Tms per plate:',style = {'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'0px'}),
                            #Boxplot
                            html.Div([dcc.Loading(dcc.Graph(id = 'clean_control_boxes', figure = generate_distribution_graph('Final_Tm',df_results),
                                style = {'display':'inline-block','vertical-align':'top','width':'100%','height':200})
                                ,type = 'dot',color = color2)],
                            id='curve_clean_boxplots')
                            ],
                        style = {'display':'inline-block','vertical-align':'top','width':'50%','margin-top':'20px','overflow-y':'scroll'}), 
                    html.Hr(style = {'margin-bottom':'0px','width': '100%'}), #Break 
                    dcc.Dropdown(id = 'Choice of controls',
                                options=[{'label': 'Plates: 1 - 16', 'value': 0}],
                                value=0,
                                style = {'display':'inline-block','vertical-align':'top','font-family': 'Arial', 'margin-bottom':'0px','margin-top':'5px','fontSize':12, 'width': '50%'}),          
                    #Curves
                    html.H4('Control curves:',style = {'font-family': 'Arial','margin-bottom':'5px', 'margin-top':'5px'}),
                    html.Div([dcc.Loading(dcc.Graph(id = 'control_curves',
                        style = {'display':'inline-block','vertical-align':'top', 'width': '100%','overflow-y':'scroll'})
                    ,type = 'dot',color = color2)],id='curve_graphs')],
                label='Control overview', style=tab2_style, selected_style=tab2_selected, value ='tab-2'),

                #MELTING TEMP/Z-SCORE PLATE OVERVIEW
                dcc.Tab([
                    html.Hr(style = {'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}), #Break
                    html.Div([
                            html.H3('Plate range:',style = {'font-family': 'Arial','margin-bottom':'5px', 'margin-top':'5px','display':'inline-block','vertical-align':'top','width':'20%'},id = 'plate_range_header'),
                            dcc.Dropdown(id = 'Choice of plates',
                                            options=[{'label': 'Plates: 1 - 50', 'value': 0}],
                                                    value=0,
                                                    style = {'display':'inline-block','vertical-align':'top','font-family': 'Arial', 'margin-bottom':'5px','margin-top':'5px','fontSize':12, 'width': '50%'}),
                            dcc.Loading(dcc.Graph(id = 'Tm_plates'),type = 'dot', color = color2),
                                ], style = {'display':'inline-block','vertical-align':'top', 'width': '50%','overflow-y':'scroll', 'height':'600px'}),
                   #Child2: Selected curve (Right of screen)
                    html.Div([
                        html.Div([],id='first_deriv_popup',style = {}),
                        html.Div([],id='popup',style = {})
                        ],id='graphs_block',style = {'display':'inline-block','vertical-align':'top', 'width': '50%'}),
                    ],label='Melting temp overview', style=tab3_style, selected_style=tab3_selected, value ='tab-3'), 

                #DOSE RESPONSE VIZ
                dcc.Tab([
                    html.Hr(style={'margin-bottom': '0px', 'width': '100%', 'margin-top': '10px'}), #Break
                    html.Div([
                        html.H4('Choose compound:', style={'font-family': 'Arial', 'margin-bottom': '0px', 'margin-top': '5px'}, id='Compound_choice_header'),
                        dcc.Dropdown(id='Compound_choice',
                                     options=[{'label': 'Placeholder', 'value': 0}],
                                     value=0,
                                     style={'font-family': 'Arial', 'margin-bottom': '0px', 'margin-top': '5px', 'margin-left': '5px', 'fontSize': 12, 'width':'50%'}),
                            ],style={'display': 'flex', 'align-items': 'center'}),
                    #Std dev block
                    html.Div([
                        html.Div([
                            html.H4('Replicate std. deviation cutoff:', style={'font-family': 'Arial', 'margin-bottom': '5px', 'margin-top': '5px', 'display': 'inline-block', 'vertical-align': 'top', 'width': '70%'}, id='Std_dev_choice_header'),
                            dcc.Input(id='std_dev_cutoff_input', name = 'Std. dev cutoff', type='number', value = 3, disabled=False,style = {'font-family': 'Arial','width':'50px','margin-top':'5px'}),
                            ], id = 'std_dev_block',style={'flex': 1, 'padding': '0 5px', 'min-width': 0}),
                        
                        #Amp upper block
                        html.Div([
                            html.H4('Rel. amplitude upper cutoff:', style={'font-family': 'Arial', 'margin-bottom': '5px', 'margin-top': '5px', 'display': 'inline-block', 'vertical-align': 'top', 'width': '70%'}, id='Amp_upper_choice_header'),
                            dcc.Input(id='amp_upper_cutoff_input', name = 'Rel. amplitude upper cutoff', type='number', value = 2, disabled=False,style = {'font-family': 'Arial','width':'50px','margin-top':'5px'}),
                            ], id = 'amp_upper_block',style={'flex': 1, 'padding': '0 5px', 'min-width': 0}),
                        
                        #Amp lower block
                        html.Div([
                            html.H4('Rel. amplitude lower cutoff:', style={'font-family': 'Arial', 'margin-bottom': '5px', 'margin-top': '5px', 'display': 'inline-block', 'vertical-align': 'top', 'width': '70%'}, id='Amp_lower_choice_header'),
                            dcc.Input(id='amp_lower_cutoff_input', name = 'Rel. amplitude lower cutoff', type='number', value = 0.2, disabled=False,style = {'font-family': 'Arial','width':'50px','margin-top':'5px'}),
                            ], id = 'amp_lower_block',style={'flex': 1, 'padding': '0 5px', 'min-width': 0}),
                    ],style={'display': 'flex', 'flex-direction': 'row','width': '100%', 'margin-top': '10px', 'margin-bottom': '15px'}, id = 'cutoff_block'),
                    
                    html.Hr(style={'margin-bottom': '0px', 'width': '100%', 'margin-top': '5px'}), #Break

                    # Parent container for the graphs with flex styling
                    html.Div([
                        # All points
                        html.Div([
                            html.H4('Replicate melting temps', style={'font-family': 'Arial', 'margin-bottom': '5px', 'margin-top': '5px', 'width': '100%'}),
                            dcc.Loading(dcc.Graph(id='DR_indiv_Tm_graph', style={'width': '100%', 'height': '100%'} ), type='dot', color=color2)
                        ], id='DR_plots_indiv_tm', style={'flex': 1, 'padding': '0 10px', 'min-width': 0}),

                        # Average Tm per conc with ribbon
                        html.Div([
                                html.H4('Delta Tm vs. Concentration', style={'font-family': 'Arial', 'margin-bottom': '5px', 'margin-top': '5px'}),
                                dcc.Loading(dcc.Graph(id='DR_avg_Tm_graph', style={'width': '100%', 'height': '100%'} ), type='dot', color=color2),
                                html.H6('X-axis option', style={'font-family': 'Arial', 'margin-bottom': '5px', 'margin-top': '2px','margin-left':'50px'}),
                                dbc.RadioItems(id = 'conc_choice', 
                                options = [{'label':'Conc', 'value':'Conc'},
                                            {'label':'log(Conc)', 'value':'log'}], 
                                            value = 'Conc',
                                            style = {'font-family': 'Arial', 'margin-bottom':'5px','margin-top':'5px','margin-left':'50px','fontSize':12},
                                            inline = True)
                        ], id='DR_plots_raw_tm', style={'flex': 1, 'padding': '0 10px', 'min-width': 0}),

                        # Average curves
                        html.Div([
                            html.H4('Avg. curve per conc', style={'font-family': 'Arial', 'margin-bottom': '5px', 'margin-top': '5px', 'width': '100%'}),
                            dcc.Loading(dcc.Graph(id='DR_delta_Tm_graph', style={'width': '100%', 'height': '100%'} ), type='dot', color=color2)
                        ], id='DR_plots_delta_tm', style={'flex': 1, 'padding': '0 10px', 'min-width': 0})
                    ], style={'display': 'flex', 'flex-direction': 'row','margin-top': '15px'}), #flex style for equal spacing

                ], label='Dose response', style=tab4_style, selected_style=tab4_selected, value='tab-4'),  # Tab div
                ],id ='tabs',value = 'tab-4',vertical =False, style = {'width':'100%'}) #All tabs style and other details
    ],id = 'page')
    return dr_layout

def full_callbacks(app, plate_report_df,pipette_df, df_results, df_curves, filtered_df, filtered_df_count, filtered_stabilizer_count, filtered_destabilizer_count, mode):
    color1 = '#B4EDD2'
    color1_faded = '#e0fff1'
    color2 = '#A0CFD3'
    color3 = '#8D94BA'
    color4 = '#9A7AA0'
    color4_faded = '#d1b6d6' 
    color5 = '#87677B'
    main_header_bkgrd_color = '#e3e3e3'
    graph_gray = '#cacacf'

    base_style = {'border': '1px black','color': 'black','font-size': '14px','font-weight': 'bold','font-family': 'Arial','padding': '6px'}
    base_selected_style = {'border': '3px solid black','color': 'black','font-size': '14px','font-weight': 'bold','font-family': 'Arial','padding':'6px'}

    tab1_style = {**base_style,'background': color1}
    tab1_selected = {**base_selected_style,'background': color1}

    tab2_style = {**base_style,'background': color2}
    tab2_selected = {**base_selected_style,'background': color2}

    tab3_style = {**base_style,'background': color3}
    tab3_selected = {**base_selected_style,'background': color3}

    tab4_style = {**base_style,'background': color4}
    tab4_selected = {**base_selected_style,'background': color4}

    tab5_style = {**base_style,'background': color5}
    tab5_selected = {**base_selected_style,'background': color5}

    data_table_style_data_conditional = [{'if': {'state': 'active'},'backgroundColor': color4_faded,'border': '1px black',},{'if': {'state': 'selected'},'backgroundColor': color4_faded,'border': '1px black',},]
    plate_report_style_data_conditional = [{'if': {'state': 'active'},'backgroundColor': color1_faded,'border': '1px black'},{'if': {'state': 'selected'},'backgroundColor': color1_faded,'border': '1px black'},]

    ####################
    #TAB 1 visualizations
    ####################

    #Highlight entire row in plate report when selected
    @app.callback(
        Output('plate_report_table', 'style_data_conditional'),
        [Input('plate_report_table', 'active_cell')],
        State('plate_report_table', 'data'))
    def update_plate_report_row_color(active,table_data):
        style_plate = plate_report_style_data_conditional.copy()
        if active:
            style_plate.append({'if': {'row_index': active['row']}, 'backgroundColor': color1_faded,'border': '1px black'},)
            selected_row_data = table_data[active['row']]
            return style_plate

    #Highlight entire row in plate report when selected
    @app.callback(
        Output('pipette_issues_table', 'style_data_conditional'),
        [Input('pipette_issues_table', 'active_cell')],
        State('pipette_issues_table', 'data'))
    def update_pipette_issues_row_color(active,table_data):
        style_plate = plate_report_style_data_conditional.copy()
        if active:
            style_plate.append({'if': {'row_index': active['row']}, 'backgroundColor': color1_faded,'border': '1px black'},)
            selected_row_data = table_data[active['row']]
            return style_plate

    ####################
    #TAB 2 visualizations
    ####################

    @app.callback(
        Output('Choice of controls', 'options'),
        Input('tabs', 'value'))
    def enable_options(tabValue):
        if tabValue == 'tab-2':
            no_plates = df_results['Assay_Plate'].nunique()
            options = []
            for i in range(0, no_plates, 16):
                label = f"Plates: {i + 1} - {min(i + 16, no_plates)}"
                value = i // 16
                options.append({'label': label, 'value': value})
            return options
        else:
            raise PreventUpdate

    #Generate ctrl curves
    @app.callback(
        Output('control_curves','figure'),
        [Input('tabs','value'),
        Input('Choice of controls', 'value')])
    def generate_curves(tabValue, dropdown_choice):
        if tabValue == 'tab-2':
            subdata = df_curves[(df_curves['Well_type']=='Control')&(df_curves['Subplot']!='Original')&(df_curves['Final_decision']!='Removed')]
            list_of_subdataframes = [d for _, d in subdata.groupby('group')] #Separate plates into groups of 16 plates
            chosen_df = list_of_subdataframes[dropdown_choice]
            figure = generate_all_ctrls_graph(chosen_df)
            return figure
        else:
            raise PreventUpdate

    ####################
    #TAB 3 visualizations
    ####################

    @app.callback(
        Output('Choice of plates', 'options'),
        Input('tabs', 'value'))
    def enable_options(tabValue):
        if tabValue == 'tab-3':
            no_plates = df_results['Assay_Plate'].nunique()
            options = []
            for i in range(0, no_plates, 16):
                label = f"Plates: {i + 1} - {min(i + 16, no_plates)}"
                value = i // 16
                options.append({'label': label, 'value': value})
            return options
        else:
            raise PreventUpdate

    #Generate heatmaps of plates
    @app.callback(
        Output('Tm_plates','figure'),
        [Input('color_filter','value'),
        Input('tabs','value'),
        Input('Choice of plates', 'value')])
    def generate_heatmaps(color_choice,tabValue, dropdown_choice):
        if tabValue == 'tab-3':
            no_plates = df_results['Assay_Plate'].nunique()
            no_wells = df_results['Well'].nunique()
            if no_wells <= 96:
                size = 50
            else:
                size = 25
            subdata = df_results[(df_results['Final_decision']!='Removed')]
            if no_plates <= 16:
                num_rows = math.ceil(no_plates)
                height = 50*num_rows
                facet_row_max_spacing = 6/height
                if color_choice == "Final_Tm":
                    figure = generate_Tm_heatmap(subdata, facet_row_max_spacing, mode)
                elif color_choice == "Well_zscore":
                    figure= generate_zscore_heatmap(subdata, facet_row_max_spacing)
                updated_figure = update_heatmap(figure, facet_row_max_spacing,num_rows,size)
                return updated_figure
            else:
                list_of_subdataframes = [d for _, d in subdata.groupby('group')] #Separate plates into groups of 16 plates
                chosen_df = list_of_subdataframes[dropdown_choice]
                sub_no_plates = chosen_df['Assay_Plate'].nunique()
                sub_num_rows = math.ceil(sub_no_plates)
                sub_height = 50*sub_num_rows
                facet_row_max_spacing = 6/sub_height
                if color_choice == "Final_Tm":
                    figure = generate_Tm_heatmap(chosen_df, facet_row_max_spacing, mode)
                elif color_choice == "Well_zscore":
                    figure= generate_zscore_heatmap(chosen_df, facet_row_max_spacing)
                updated_figure = update_heatmap(figure, facet_row_max_spacing,sub_num_rows,size)
                return updated_figure
        else:
            raise PreventUpdate

    #Render curve graph for selected data in Melting temp tab
    @app.callback(
        Output('popup', 'children'),
        Output('first_deriv_popup','children'),
        Input('Tm_plates', 'clickData'))
    def update_popup(clickData):
        if clickData:
            selected_unique_key = clickData['points'][0]['customdata'][5]
            curve_figure = generate_selected_curve(selected_unique_key, df_curves, mode)
            curve_figure.update_layout(margin=dict(l=10, r=10, t=40, b=10), height=300)
            first_deriv_curve = generate_first_derivative_curve(selected_unique_key, df_curves, df_results, mode)
            orig_graph_object = dcc.Graph(id = 'selected_curve',figure=curve_figure)
            if first_deriv_curve == None:
                return orig_graph_object, None
            else:
                first_deriv_curve.update_layout(margin=dict(l=10, r=10, t=40, b=10), height=300)
                first_deriv_graph_object = dcc.Graph(id = 'first_deriv_curve',figure=first_deriv_curve)
                return orig_graph_object, first_deriv_graph_object
        else:
            return None, None

    ####################
    #TAB 4 visualizations
    ####################

    #Custom Sorting
    @app.callback(
        Output('results_table_datatable', 'data'),
        Input('results_table_datatable', 'sort_by'),
        State('results_table_datatable', 'data')
    )
    def sort_table(sort_by, rows):
        if not sort_by or rows is None:
            return rows

        df = pd.DataFrame(rows)
        # Apply all sorts in reverse order so that the first sort takes precedence
        for sort in reversed(sort_by):
            df.sort_values(
                by=sort['column_id'],
                ascending=sort['direction'] == 'asc',
                inplace=True,
                kind='mergesort'  # stable sort preserves prior order
            )
        return df.to_dict('records') 

    # En/Dis-able options based on hit ID choice
    @app.callback(
        Output('relative_zscore_cutoff', 'disabled'),
        Output('relative_stddev_cutoff', 'disabled'),
        Output('relative_degree_cutoff', 'disabled'),
        Output('zscore_header', 'style'),
        Output('stddev_header', 'style'),
        Output('degree_header', 'style'),
        Output('zscore_cutoff_header', 'style'),
        Output('stddev_cutoff_header', 'style'),
        Output('degree_cutoff_header', 'style'),
        Input('hit_type', 'value'),
        prevent_initial_call=True)
    def hit_choice_option_enabled(hit_choice_value):
        header_disabled_style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'margin-right':'2px','color':'grey'}
        header_enabled_style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'margin-right':'2px','color':'black'}
        cutoff_enabled_style = {'font-family': 'Arial','margin-bottom':'5px','margin-top':'7px'}
        cutoff_disabled_style = {'font-family': 'Arial','margin-bottom':'5px','margin-top':'7px', 'color': 'grey'}
        if hit_choice_value == 'zscore':
            return False, True, True, header_enabled_style, header_disabled_style, header_disabled_style, cutoff_enabled_style, cutoff_disabled_style, cutoff_disabled_style
        elif hit_choice_value == 'stddev':
            return True, False, True, header_disabled_style, header_enabled_style, header_disabled_style, cutoff_disabled_style, cutoff_enabled_style, cutoff_disabled_style
        elif hit_choice_value == 'degrees':
            return True, True, False, header_disabled_style, header_disabled_style, header_enabled_style, cutoff_disabled_style, cutoff_disabled_style, cutoff_enabled_style
        else:
            return True, True, True, header_disabled_style, cutoff_disabled_style, cutoff_disabled_style, cutoff_disabled_style

    #Update hit table based on input box values
    @app.callback(
        Output('results_table', 'children'),
        Input('submit', 'n_clicks'),
        Input('results_table_datatable', 'data_previous'),
        State('results_table_datatable', 'data'),
        Input('hit_type', 'value'),
        State('relative_zscore_cutoff', 'value'),
        State('relative_stddev_cutoff', 'value'),
        State('relative_degree_cutoff', 'value'),
        State('lower_limit_amp', 'value'),
        State('upper_limit_amp', 'value'),
        State('hit_filter','value'),
        State('results_table_datatable', 'derived_viewport_data'),
        prevent_initial_call=True)
    def update_hit_table(n_clicks, previous, current, hit_choice_value, zscore_cutoff, stddev_cutoff, degree_cutoff, lower_value_amp, upper_value_amp, hit_option, view):
        if ((n_clicks < 1)&(previous is None)) | ((n_clicks >= 1)&(previous is None)) :
            original_df  = df_results[(df_results['Well_type'] != 'Control')&(df_results['Plate_status'] == 'OK')&(df_results['Error'] == '')]
            #small_df = original_df.loc[:, ['Source_Plate','Well','Subplot','Compound','Fraction','Final_Tm','Well_zscore', 'Relative_amplitude','Max_ctrl_zscore_for_plate','Min_ctrl_zscore_for_plate','Diff from ctrl avg','Std. devs from ctrl mean','Unique_key','Unique_key_subplot']]
            small_df = original_df.loc[:, ['Source_Plate','Well','Subplot','Compound','Fraction','Final_Tm','Well_zscore', 'Std. devs from ctrl mean','Diff from ctrl avg','Relative_amplitude','Max_ctrl_zscore_for_plate','Min_ctrl_zscore_for_plate','Unique_key','Unique_key_subplot','group']]
            if hit_choice_value == 'zscore':
                cutoff = zscore_cutoff/100
                small_df.loc[:,"Hit"] = np.where((small_df['Well_zscore'] > small_df['Max_ctrl_zscore_for_plate'] + ((small_df['Max_ctrl_zscore_for_plate'].abs())*cutoff))|
                                                 (small_df['Well_zscore'] < small_df['Min_ctrl_zscore_for_plate'] - ((small_df['Min_ctrl_zscore_for_plate'].abs())*cutoff)),
                                                "Hit","Not a hit")
            elif hit_choice_value == 'stddev':
                small_df.loc[:,"Hit"] = np.where((small_df['Std. devs from ctrl mean'] > stddev_cutoff)|
                                                 (small_df['Std. devs from ctrl mean'] < -stddev_cutoff),
                                                "Hit","Not a hit")
            elif hit_choice_value == 'degrees':
                small_df.loc[:,"Hit"] = np.where((small_df['Diff from ctrl avg'] > degree_cutoff)|
                                                 (small_df['Diff from ctrl avg'] < -degree_cutoff),
                                                "Hit","Not a hit")
            else:
                pass
            filtered_df = small_df[((small_df.Hit == 'Hit'))&((small_df.Relative_amplitude <= upper_value_amp)&(small_df.Relative_amplitude >= lower_value_amp))]
            filtered_df.loc[:,'Influence'] = np.where(filtered_df['Diff from ctrl avg'] > 0, "Stabilizer","Destabilizer")
            if hit_option == 'stabilizers':
                filtered_df = filtered_df[filtered_df['Influence'] == 'Stabilizer']
            elif hit_option == 'destabilizers':
                filtered_df = filtered_df[filtered_df['Influence'] == 'Destabilizer']
            else:
                pass
            filtered_df_count = filtered_df.shape[0]
            filtered_stabilizer_count = filtered_df['Influence'].value_counts().get('Stabilizer', 0)
            filtered_destabilizer_count = filtered_df['Influence'].value_counts().get('Destabilizer', 0)
            filtered_df = filtered_df.sort_values('Well_zscore',ascending=False)
            new_cols = [col for col in filtered_df.columns if (col != 'Unique_key_subplot')&(col != 'Unique_key')] + ['Unique_key'] + ['Unique_key_subplot'] #Stick the unique key column at the end
            filtered_df = filtered_df[new_cols]
            data_table = dash_table.DataTable(filtered_df.to_dict('records'), [{'name': i, 'id': i} for i in filtered_df.columns], 
                    id = 'results_table_datatable',
                    hidden_columns=['Max_ctrl_zscore_for_plate','Hit','Min_ctrl_zscore_for_plate','group'], 
                    style_as_list_view=True, 
                    row_deletable=True, 
                    sort_action='custom',
                    sort_mode='single',
                    sort_by=[],
                    export_format='xlsx',
                    style_table = {'height':'250px','overflow-y':'scroll'},
                    style_cell={'fontSize':12, 'font-family':'Arial'},
                    css=[{"selector": ".show-hide", "rule": "display: none"}], 
                    style_header = {'backgroundColor': color4,'fontWeight': 'bold'})
            header = html.Label("Total hits: "+str(filtered_df_count), style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'font-size': '20px', 'color':'dark green'}, id = 'hits_header')
            break1 = html.Hr(style = {'margin-bottom':'0px','width': '2%', 'margin-top':'2px', 'border':'0px white'})
            stabilizers = html.Label("Stabilizers: "+str(filtered_stabilizer_count), style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'font-size': '16px'}, id = 'stabilizer_header')
            break2 = html.Hr(style = {'margin-bottom':'0px','width': '2%', 'margin-top':'2px', 'border':'0px white'})
            destabilizers = html.Label("Destabilizers: "+str(filtered_destabilizer_count), style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'font-size': '16px'}, id = 'destabilizer_header')
            return header,break1,stabilizers,break2,destabilizers,data_table
        else: #User change to datatable (e.g. delete row)
            current_df = pd.DataFrame(current)
            current_df.loc[:,'Influence'] = np.where(current_df['Diff from ctrl avg'] > 0, "Stabilizer","Destabilizer")
            filtered_df_count = current_df.shape[0]
            filtered_stabilizer_count = current_df['Influence'].value_counts().get('Stabilizer', 0)
            filtered_destabilizer_count = current_df['Influence'].value_counts().get('Destabilizer', 0)
            data_table = dash_table.DataTable(current_df.to_dict('records'), [{'name': i, 'id': i} for i in current_df.columns], 
                    id = 'results_table_datatable',
                    hidden_columns=['Max_ctrl_zscore_for_plate','Hit','Min_ctrl_zscore_for_plate','group'], 
                    style_as_list_view=True, 
                    row_deletable=True, 
                    sort_action='custom',
                    sort_mode='single',
                    sort_by=[],
                    export_format='xlsx',
                    style_table = {'height':'250px','overflow-y':'scroll'},
                    style_cell={'fontSize':12, 'font-family':'Arial'},
                    css=[{"selector": ".show-hide", "rule": "display: none"}],  
                    style_header = {'backgroundColor': color4,'fontWeight': 'bold'})
            header = html.Label("Total hits: "+str(filtered_df_count), style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'font-size': '20px', 'color':'dark green'}, id = 'hits_header')
            break1 = html.Hr(style = {'margin-bottom':'0px','width': '2%', 'margin-top':'2px', 'border':'0px white'})
            stabilizers = html.Label("Stabilizers: "+str(filtered_stabilizer_count), style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'font-size': '16px'}, id = 'stabilizer_header')
            break2 = html.Hr(style = {'margin-bottom':'0px','width': '2%', 'margin-top':'2px', 'border':'0px white'})
            destabilizers = html.Label("Destabilizers: "+str(filtered_destabilizer_count), style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'font-size': '16px'}, id = 'destabilizer_header')
            return header,break1,stabilizers,break2,destabilizers,data_table

    #Highlight entire row in results Datatable when selected and produce pop up graph
    @app.callback(
        [Output('results_table_datatable', 'style_data_conditional'),
        Output('popup3', 'children'),
        Output('popup4', 'children'),
        Output('clear_pop_up_div','hidden')],
        [Input('results_table_datatable', 'active_cell'),
        Input('clear_popup','n_clicks')],
        State('results_table_datatable', 'derived_viewport_data'))
    def update_selected_row_color(active,n_clicks, table_data):
        style = data_table_style_data_conditional.copy()
        ctx = callback_context
        if not ctx.triggered:
            raise PreventUpdate
        trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
        if trigger_id == 'clear_popup':
            if n_clicks > 0:
                style.append({'if': {'row_index': active['row']}, 'backgroundColor': color4_faded,'border': '1px black'},)
                return style, None, None, True
        elif trigger_id =='results_table_datatable':
            if active:
                style.append({'if': {'row_index': active['row']}, 'backgroundColor': color4_faded,'border': '1px black'},)
                current_data_table = pd.DataFrame(table_data)
                selected_row_data = table_data[active['row']]
                selected_unique_key = selected_row_data['Unique_key']
                curve_figure = generate_selected_curve(selected_unique_key, df_curves, mode)
                curve_figure.update_layout(title={'yref': 'paper','y' : 1,'yanchor' : 'bottom','font':{'size':12}}, title_pad_b = 25, margin=dict(l=20, r=20, t=50, b=20), height = 200)
                curve_figure.update_layout(font = {'size':8}, showlegend=False)
                graph_object = dcc.Graph(id = 'selected_curve',figure=curve_figure, style = {'fontSize':12, 'font-family':'Arial'})
                first_deriv_curve = generate_first_derivative_curve(selected_unique_key, df_curves, df_results, mode)
                first_deriv_curve.update_layout(title={'yref': 'paper','y' : 1,'yanchor' : 'bottom','font':{'size':12}}, title_pad_b = 25, margin=dict(l=20, r=20, t=50, b=20), height = 200)
                first_deriv_curve.update_layout(font = {'size':8}, showlegend=False)
                graph_object2 = dcc.Graph(id = 'first_deriv_curve',figure=first_deriv_curve, style = {'fontSize':12, 'font-family':'Arial'})
                return style, graph_object, graph_object2, False
            else:
                return style, None, None, True
        raise PreventUpdate


    # Callback to store the table data (assuming it's generated elsewhere)
    @app.callback(
        Output('current_table_data_store', 'data'),
        [Input('generate', 'n_clicks'),
         Input('submit', 'n_clicks')],
        [State('results_table_datatable', 'derived_virtual_data')] # Assuming this exists
    )
    def store_table_data(generate_clicks, submit_clicks, current_table_data):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if current_table_data is not None:
            return current_table_data
        raise PreventUpdate

    # Callback to manage the visibility and options of the hits_dropdown
    @app.callback(
        [Output('hits_dropdown_container', 'style'),
         Output('hits_dropdown', 'options'),
         Output('hits_dropdown', 'value')],
        [Input('current_table_data_store', 'data'),
         Input('submit', 'n_clicks'),
         Input('clear', 'n_clicks'),
         Input('summary_graphs', 'n_clicks'),
         Input('generate', 'n_clicks')], # Also trigger on generate to reset value
        [State('hits_dropdown', 'value')] # Keep the current value if possible
    )
    def manage_hits_dropdown(current_table_data, update_clicks, clear_clicks, summary_clicks, generate_clicks, current_dropdown_value):
        ctx = dash.callback_context
        if not ctx.triggered_id or not current_table_data:
            return {'display': 'none'}, [], 0

        triggered_id = ctx.triggered_id

        # If 'submit', 'clear', or 'summary_graphs' were most recently clicked,
        # or if there's no data, hide and clear the dropdown.
        if triggered_id in ['submit', 'clear', 'summary_graphs']:
            return {'display': 'none'}, [], 0

        key_list = [d['Unique_key'] for d in current_table_data]
        max_hits = 50

        if len(key_list) > max_hits:
            num_pages = math.ceil(len(key_list) / max_hits)
            dropdown_options = []
            for i in range(num_pages):
                start_item = i * max_hits + 1
                end_item = min((i + 1) * max_hits, len(key_list))
                dropdown_options.append({'label': f'Hits {start_item} - {end_item}', 'value': i})

            # Set the value to 0 if generating, otherwise try to keep the current value
            ctx = dash.callback_context
            trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
            if trigger_id == 'generate':
                chosen_value = 0
            else:
                # Ensure the current_dropdown_value is still a valid option
                if current_dropdown_value in [opt['value'] for opt in dropdown_options]:
                    chosen_value = current_dropdown_value
                else:
                    chosen_value = 0 # Default to first page if current value is invalid

            return {'display': 'block'}, dropdown_options, chosen_value
        else:
            return {'display': 'none'}, [], 0

    # Generating or clearing graphs for all hits
    @app.callback(
        Output('all_graphs_div', 'children'),
        [Input('generate', 'n_clicks'),
         Input('clear', 'n_clicks'),
         Input('summary_graphs', 'n_clicks'),
         Input('submit', 'n_clicks'),
         Input('hits_dropdown', 'value')],  # This is now reliably present when visible
        [State('current_table_data_store', 'data'),
        State('results_table_datatable', 'derived_virtual_data')] # Use the stored data
    )
    def update_or_clear_graphs(generate_clicks, clear_clicks, summary_clicks, update_clicks, hit_choice, current_table_data,hit_data):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]

        max_hits = 50

        if trigger_id == 'generate' or trigger_id == 'hits_dropdown':
            if current_table_data is None:
                raise PreventUpdate

            key_list = [d['Unique_key'] for d in hit_data]

            if not key_list:
                return None  # Return nothing if there's no data

            if len(key_list) > max_hits:
                current_page = hit_choice if hit_choice is not None else 0 # Use the value from the dropdown
                start_index = current_page * max_hits
                end_index = start_index + max_hits
                paginated_key_list = key_list[start_index:end_index]

                new_figure = generate_all_curves(paginated_key_list, df_curves)
                new_figure.update_yaxes(matches=None)
                new_figure.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
                new_figure.update_layout(height=200 * len(paginated_key_list) / 5)

                graph_object = dcc.Graph(id='all_graphs_object', figure=new_figure, style={'font-family': 'Arial'})

                # The dropdown is now managed by a separate callback, so we just return the graph
                return graph_object

            elif len(key_list) < 6:
                new_figure = generate_all_curves(key_list, df_curves)
                new_figure.update_yaxes(matches=None)
                new_figure.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
                new_figure.update_layout(height=250)
                graph_object = dcc.Graph(id='all_graphs_object', figure=new_figure, style={'font-family': 'Arial'})
                return graph_object

            elif len(key_list) < 10:
                new_figure = generate_all_curves(key_list, df_curves)
                new_figure.update_yaxes(matches=None)
                new_figure.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
                new_figure.update_layout(height=400)
                graph_object = dcc.Graph(id='all_graphs_object', figure=new_figure, style={'font-family': 'Arial'})
                return graph_object

            else:
                new_figure = generate_all_curves(key_list, df_curves)
                new_figure.update_yaxes(matches=None)
                new_figure.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
                new_figure.update_layout(height=200 * len(key_list) / 5)
                graph_object = dcc.Graph(id='all_graphs_object', figure=new_figure, style={'font-family': 'Arial'})
                return graph_object

        elif trigger_id == 'summary_graphs':
            if (summary_clicks is not None) and (hit_data is not None):
                hits_key_list = [d['Unique_key_subplot'] for d in hit_data]
                df_results.loc[:, 'Hit'] = np.where(df_results['Unique_key_subplot'].isin(hits_key_list), 'Hit', 'Not a hit')
                df_results.loc[:, 'Hit'] = np.where(df_results['Well_type'] == 'Control', "Control", df_results['Hit'])
                df_hits_raw = df_results[df_results['Hit'] == 'Hit']
                df_hits_sub = df_hits_raw.loc[:, ['Platename','Fraction','Well','Hit']]
                df_hits = df_hits_sub.drop_duplicates()
                Hit_distribution_figure = generate_scatterplot(df_results)
                Hits_per_plate_figure = generate_barplot(df_hits,'Platename', 'Fraction')
                Hits_per_fraction_figure = generate_barplot(df_hits, 'Fraction', 'Platename')
            else:
                Hit_distribution_figure = None
                Hits_per_plate_figure = None
                Hits_per_fraction_figure = None

            graph_Div = html.Div([
                html.Hr(style={'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}),
                html.Div([
                    html.H3('Distribution of hits:',style={'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'20px'}),
                    dcc.Loading(dcc.Graph(id='Hit_distribution_graph',figure=Hit_distribution_figure),type='dot', color=color5)
                ], id='hit_plots_div1', style={'display':'inline-block','vertical-align':'top', 'width': '100%'}),
                html.Hr(style={'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}),
                html.Div([
                    html.H3('Hits per plate colored by fraction:',style={'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'20px'}),
                    dcc.Loading(dcc.Graph(id='Hits_per_plate',figure=Hits_per_plate_figure),type='dot', color=color5)
                ], id='hit_plots_div2', style={'display':'inline-block','vertical-align':'top', 'width': '50%'}),
                html.Div([
                    html.H3('Hits per fraction colored by plate:',style={'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'20px'}),
                    dcc.Loading(dcc.Graph(id='Hits_per_fraction',figure=Hits_per_fraction_figure),type='dot', color=color5)
                ], id='hit_plots_div3', style={'display':'inline-block','vertical-align':'top', 'width': '50%'})
            ])
            return graph_Div

        elif trigger_id == 'clear':
            if (clear_clicks is not None):
                return None

        elif trigger_id == 'submit':
            if update_clicks is not None:
                return None

        raise PreventUpdate

def tm_only_callbacks(app, df_results, df_curves, mode):  
    @app.callback(
        Output('Choice of plates', 'options'),
        Input('tabs', 'value'))
    def enable_options(tabValue):
        if tabValue == 'tab-3':
            no_plates = df_results['Assay_Plate'].nunique()
            options = []
            for i in range(0, no_plates, 16):
                label = f"Plates: {i + 1} - {min(i + 16, no_plates)}"
                value = i // 16
                options.append({'label': label, 'value': value})
            return options
        else:
            raise PreventUpdate

    #Generate heatmaps of plates
    @app.callback(
        Output('Tm_plates','figure'),
        [Input('tabs','value')])
    def generate_heatmaps(tabValue):
        if tabValue == 'tab-3':
            no_plates = df_results['Assay_Plate'].nunique()
            no_wells = df_results['Well'].nunique()
            if no_wells <= 96:
                size = 50
            else:
                size = 25
            if no_plates <= 16:
                num_rows = math.ceil(no_plates)
                height = 50*num_rows
                facet_row_max_spacing = 6/height
                figure = generate_Tm_heatmap(df_results, facet_row_max_spacing, mode)
                updated_figure = update_heatmap(figure, facet_row_max_spacing,num_rows,size)
                return updated_figure
            else:
                list_of_subdataframes = [d for _, d in subdata.groupby('group')] #Separate plates into groups of 16 plates
                chosen_df = list_of_subdataframes[dropdown_choice]
                sub_no_plates = chosen_df['Assay_Plate'].nunique()
                sub_num_rows = math.ceil(sub_no_plates)
                sub_height = 50*sub_num_rows
                facet_row_max_spacing = 6/sub_height
                figure = generate_Tm_heatmap(chosen_df, facet_row_max_spacing, mode)
                updated_figure = update_heatmap(figure, facet_row_max_spacing,sub_num_rows,size)
                return updated_figure
        else:
            raise PreventUpdate

    #Render curve graph for selected data in Melting temp tab
    @app.callback(
        Output('popup', 'children'),
        Output('first_deriv_popup','children'),
        Input('Tm_plates', 'clickData'))
    def update_popup(clickData):
        if clickData:
            selected_unique_key = clickData['points'][0]['customdata'][5]
            curve_figure = generate_selected_curve(selected_unique_key, df_curves, mode)
            first_deriv_curve = generate_first_derivative_curve(selected_unique_key, df_curves, df_results, mode)
            orig_graph_object = dcc.Graph(id = 'selected_curve',figure=curve_figure)
            if first_deriv_curve == None:
                return orig_graph_object, None
            else:
                first_deriv_graph_object = dcc.Graph(id = 'first_deriv_curve',figure=first_deriv_curve)
                return orig_graph_object, first_deriv_graph_object
        else:
            return None, None

def DR_callbacks(app, df_results, df_curves, mode): 

    ####################
    #TAB 2 visualizations
    ####################

    @app.callback(
        Output('Choice of controls', 'options'),
        Input('tabs', 'value'))
    def enable_options(tabValue):
        if tabValue == 'tab-2':
            no_plates = df_results['Assay_Plate'].nunique()
            options = []
            for i in range(0, no_plates, 16):
                label = f"Plates: {i + 1} - {min(i + 16, no_plates)}"
                value = i // 16
                options.append({'label': label, 'value': value})
            return options
        else:
            raise PreventUpdate

    #Generate ctrl curves
    @app.callback(
        Output('control_curves','figure'),
        [Input('tabs','value'),
        Input('Choice of controls', 'value')])
    def generate_curves(tabValue, dropdown_choice):
        if tabValue == 'tab-2':
            subdata = df_curves[(df_curves['Well_type']=='Control')&(df_curves['Subplot']!='Original')&(df_curves['Final_decision']!='Removed')]
            list_of_subdataframes = [d for _, d in subdata.groupby('group')] #Separate plates into groups of 16 plates
            chosen_df = list_of_subdataframes[dropdown_choice]
            figure = generate_all_ctrls_graph(chosen_df)
            return figure
        else:
            raise PreventUpdate

    ####################
    #TAB 3 visualizations
    ####################

    @app.callback(
        Output('Choice of plates', 'options'),
        Input('tabs', 'value'))
    def enable_options(tabValue):
        if tabValue == 'tab-3':
            no_plates = df_results['Assay_Plate'].nunique()
            options = []
            for i in range(0, no_plates, 16):
                label = f"Plates: {i + 1} - {min(i + 16, no_plates)}"
                value = i // 16
                options.append({'label': label, 'value': value})
            return options
        else:
            raise PreventUpdate

    #Generate heatmaps of plates
    @app.callback(
        Output('Tm_plates','figure'),
        [Input('tabs','value'),
        Input('Choice of plates', 'value')])
    def generate_heatmaps(tabValue, dropdown_choice):
        if tabValue == 'tab-3':
            no_plates = df_results['Assay_Plate'].nunique()
            no_wells = df_results['Well'].nunique()
            if no_wells <= 96:
                size = 50
            else:
                size = 25
            if no_plates <= 16:
                num_rows = math.ceil(no_plates)
                height = 50*num_rows
                facet_row_max_spacing = 6/height
                figure = generate_Tm_heatmap(df_results, facet_row_max_spacing, mode)
                updated_figure = update_heatmap(figure, facet_row_max_spacing,num_rows,size)
                return updated_figure
            else:
                subdata = df_results[(df_results['Final_decision']!='Removed')]
                list_of_subdataframes = [d for _, d in subdata.groupby('group')] #Separate plates into groups of 16 plates
                chosen_df = list_of_subdataframes[dropdown_choice]
                sub_no_plates = chosen_df['Assay_Plate'].nunique()
                sub_num_rows = math.ceil(sub_no_plates)
                sub_height = 50*sub_num_rows
                facet_row_max_spacing = 6/sub_height
                figure = generate_Tm_heatmap(chosen_df, facet_row_max_spacing, mode)
                updated_figure = update_heatmap(figure, facet_row_max_spacing,sub_num_rows,size)
                return updated_figure
        else:
            raise PreventUpdate

    #Render curve graph for selected data in Melting temp tab
    @app.callback(
        Output('popup', 'children'),
        Output('first_deriv_popup','children'),
        Input('Tm_plates', 'clickData'))
    def update_popup(clickData):
        if clickData:
            selected_unique_key = clickData['points'][0]['customdata'][7]
            curve_figure = generate_selected_curve(selected_unique_key, df_curves, mode)
            first_deriv_curve = generate_first_derivative_curve(selected_unique_key, df_curves, df_results, mode)
            orig_graph_object = dcc.Graph(id = 'selected_curve',figure=curve_figure)
            if first_deriv_curve == None:
                return orig_graph_object, None
            else:
                first_deriv_graph_object = dcc.Graph(id = 'first_deriv_curve',figure=first_deriv_curve)
                return orig_graph_object, first_deriv_graph_object
        else:
            return None, None

    ####################
    #TAB 4 visualizations
    ####################

    @app.callback(
        [Output('Compound_choice', 'options'),
         Output('Compound_choice', 'value')], # Chaining
        Input('tabs', 'value')
    )
    def enable_options(tabValue):
        if tabValue == 'tab-4':
            raw_compounds = df_results['Compound'].unique()
            compounds = raw_compounds[raw_compounds != 'Control']
            
            # Prevent error if there are no plates
            if len(compounds) == 0:
                return [], None

            # Create the options for the dropdown
            options = [{'label': i, 'value': i} for i in compounds]
            
            # Set the dropdown's value to the FIRST plate in the list
            initial_value = compounds[0]
            
            return options, initial_value 
        else:
            raise PreventUpdate


    @app.callback(
        [Output('DR_indiv_Tm_graph', 'figure'),
         Output('DR_avg_Tm_graph', 'figure'),
         Output('DR_delta_Tm_graph', 'figure')],
        [Input('tabs', 'value'),
         Input('Compound_choice', 'value'),
         Input('std_dev_cutoff_input', 'value'),
         Input('amp_upper_cutoff_input', 'value'),
         Input('amp_lower_cutoff_input', 'value'),
         Input('conc_choice','value')]
    )
    def generate_avg_Tm_plots(tabValue, compound_choice, std_dev_cutoff, upper_amp_cutoff, lower_amp_cutoff, log_it_option):
        if not compound_choice or tabValue != 'tab-4' or std_dev_cutoff is None or upper_amp_cutoff is None or lower_amp_cutoff is None:
            raise PreventUpdate

        chosen_cmpd_df = df_results[(df_results['Compound'] == compound_choice)] 
        chosen_cmpd_plate = chosen_cmpd_df['Assay_Plate'].unique()[0]
        plate_controls_for_chosen_compound = df_results[(df_results['Assay_Plate'] == chosen_cmpd_plate)&(df_results['Well_type'] == 'Control')&(df_results['Final_decision'] == 'Pass')].copy()
        compound_choice_df = pd.concat([chosen_cmpd_df, plate_controls_for_chosen_compound], axis=0, ignore_index=True)
        compound_choice_df['zscore'] = compound_choice_df.groupby(['Concentration'])['Final_Tm'].transform(stats.zscore).abs()
        filtered_df = compound_choice_df[(compound_choice_df['Relative_amplitude'] < upper_amp_cutoff)&(compound_choice_df['Relative_amplitude'] > lower_amp_cutoff)&(compound_choice_df['zscore'] < std_dev_cutoff)].copy()
        getting_curves = filtered_df.loc[:,['Unique_key','Concentration']]
        multi_curve_df_avg = get_curves_for_avg(df_curves, getting_curves) #Data crunching for avg curves
        
        #Generate three figs
        indiv_figure = generate_indiv_Tm_graph(filtered_df)
        delta_figure = generate_delta_Tm_graphs(filtered_df,log_it_option)
        multi_line_avg_figure = generate_avg_curves(multi_curve_df_avg)
        
        return indiv_figure, delta_figure, multi_line_avg_figure









