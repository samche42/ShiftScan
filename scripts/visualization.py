#!/usr/bin/env python3

import dash
from dash import Dash, html, dcc, Input, Output, State, dash_table, callback_context
import pandas as pd
import plotly.express as px
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate
import numpy as np
import argparse
import math

# Initialize the app
app = Dash(__name__)
app.config.suppress_callback_exceptions = True

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_dir", help="Full file path to output files from previous step")

args = parser.parse_args()

####################
#
# STYLES
#
####################


yellow_tab_style = {
    'background': '#FFEF79',
    'border': '1px solid #FFDA0A',
    'color': 'black',
    'font-size': '14px',
    'font-weight': 600,
    'font-family': 'Arial',
    'border-radius': '4px',
    'padding':'6px'
}

orange_tab_style = {
    'background': '#F9C159',
    'border': '1px solid orange',
    'color': 'black',
    'font-size': '14px',
    'font-weight': 600,
    'font-family': 'Arial',
    'border-radius': '4px',
    'padding':'6px',
}

red_tab_style = {
    'background': '#FE667E',
    'border': '1px solid red',
    'color': 'black',
    'font-size': '14px',
    'font-weight': 600,
    'font-family': 'Arial',
    'border-radius': '4px',
    'padding':'6px',
}

purple_tab_style = {
    'background': '#D397FF',
    'border': '1px solid purple',
    'color': 'black',
    'font-size': '14px',
    'font-weight': 600,
    'font-family': 'Arial',
    'border-radius': '4px',
    'padding':'6px',
}

blue_tab_style = {
    'background': '#A7D1FE',
    'color': 'black',
    'border': '1px solid navy',
    'font-size': '14px',
    'font-weight': 600,
    'font-family': 'Arial',
    'border-radius': '4px',
    'padding':'6px',
}

green_tab_style = {
    'background': '#a6ed6f',
    'color': 'black',
    'border': '1px solid green',
    'font-size': '14px',
    'font-weight': 600,
    'font-family': 'Arial',
    'border-radius': '4px',
    'padding':'6px',
}

tab_selected_style = {
    'background': '#D3D3D3',
    'color': 'black',
    'border': '1px solid black',
    'font-size': '14px',
    'font-weight': 600,
    'font-family': 'Arial',
    'border-radius': '4px',
    'padding':'6px',
}

plate_report_style_data_conditional = [
    {
        'if': {'state': 'active'},
        'backgroundColor': '#FFFBB2',
        'border': '1px solid #FFFBB2'
    },
    {
        'if': {'state': 'selected'},
        'backgroundColor': '#FFFBB2',
        'border': '1px solid #FFFBB2'
    },
]


pipette_style_data_conditional = [
    {
        'if': {'state': 'active'},
        'backgroundColor': '#FFFBB2',
        'border': '1px solid #FFFBB2'
    },
    {
        'if': {'state': 'selected'},
        'backgroundColor': '#FFFBB2',
        'border': '1px solid #FFFBB2'
    },
]

data_table_style_data_conditional = [
    {
        'if': {'state': 'active'},
        'backgroundColor': '#cef786',
        'border': '1px solid green',
    },
    {
        'if': {'state': 'selected'},
        'backgroundColor': '#cef786',
        'border': '1px solid green',
    },
]


#Fiddle with input data to correct format
df_results = pd.read_csv(args.input_dir+'/Final_results.txt', sep='\t', header = 0)
df_results['Error'].fillna('', inplace=True)
df_results['Well_zscore']=df_results['Well_zscore'].apply(lambda x:round(x,2))
df_results['Relative_amplitude']=df_results['Relative_amplitude'].apply(lambda x:round(x,2))
df_curves = pd.read_csv(args.input_dir+'/Final_curves.txt', sep='\t', header = 0)
df_curves['Error'].fillna('', inplace=True)
df_plate_report = pd.read_csv(args.input_dir+'/Plate_report.txt', sep='\t', header = 0)
df_pipette = pd.read_csv(args.input_dir+'/Potential_problems.txt', sep='\t', header = 0)

####################
#
# GENERATING DEFAULT DATA TABLE (TAB5)
#
####################

default_datatable = df_results[(df_results['Well_type'] != 'Control')&(df_results['Plate_status'] == 'OK')&(df_results['Error'] == '')]
default_datatable = default_datatable.rename(columns = {'Well_zscore':'zscore'})
default_datatable['zscore'] = default_datatable['zscore'].abs() #Convert values to absolute values for % calculations
default_datatable['Max_ctrl_zscore_for_plate'] = default_datatable['Max_ctrl_zscore_for_plate'].abs()
default_datatable = default_datatable[((default_datatable.zscore >= -default_datatable.Max_ctrl_zscore_for_plate))&((default_datatable.Relative_amplitude <= 1.5)&(default_datatable.Relative_amplitude >= 0.5))]
default_datatable.sort_values('zscore',ascending=False, inplace= True)
default_datatable = default_datatable.rename(columns = {'zscore':'Well_zscore'})
default_datatable['Max_ctrl_zscore_for_plate']=default_datatable['Max_ctrl_zscore_for_plate'].apply(lambda x:round(x,2))
default_datatable_count = default_datatable['Unique_key'].nunique()

####################
#
# WEBPAGE LAYOUT
#
####################

app.layout = html.Div([html.H1('DSF Analysis Visualizations', style = {'background': '#e3e3e3','font-family': 'Arial','border-radius': '4px','border': f'1px solid #70706d','padding':'10px'}),
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
                        df_plate_report.to_dict('records'), 
                        [{'name': i, 'id': i} for i in df_plate_report.columns],
                        id = 'plate_report_table',
                        sort_action='native',
                        export_format='xlsx',
                        style_data_conditional=plate_report_style_data_conditional, 
                        style_as_list_view=True, style_cell={'fontSize':12, 'font-family':'Arial'}, style_header = {'backgroundColor': '#FFEF79','fontWeight': 'bold'})
                    ], style = {'display':'inline-block','vertical-align':'top','width':'60%', 'height': '80%','margin-top':'30px','overflow-y':'scroll', 'margin-right': '50px'}),
                #Pipetting problems
                html.Div([
                    #Header
                    html.H3('Possible pipetting issues:',style = {'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'0px'}),
                    #Table
                    dash_table.DataTable(
                        df_pipette.to_dict('records'), 
                        [{'name': i, 'id': i} for i in df_pipette.columns],
                        id = 'pipette_issues_table',
                        export_format='xlsx',
                        style_data_conditional=pipette_style_data_conditional, 
                        style_as_list_view=True, style_cell={'fontSize':12, 'font-family':'Arial'}, style_header = {'backgroundColor': '#FFEF79','fontWeight': 'bold'})
                    ], style = {'display':'inline-block','vertical-align':'top','width':'20%', 'height': '80%','margin-top':'30px','overflow-y':'scroll'})
                ],
            label='Plate overview', style=orange_tab_style, selected_style=tab_selected_style, value='tab-1'), #Plate overview tab style
            #CURVE TAB
            dcc.Tab([
                html.Hr(style = {'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}), #Break
                #All ctrl Tm div
                html.Div([
                        #Header
                        html.H3('Distribution of all control Tms per plate:',style = {'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'0px'}),
                        #Boxplot
                        html.Div([dcc.Graph(id = 'all_control_boxes',
                        figure=px.violin(df_results[(df_results['Well_type']=='Control')], x='Assay_Plate', y='Smooth_Tm',height = 200,
                        labels={'Assay_Plate': 'Assay Plate', 'Smooth_Tm': 'Melting temp'}) 
                        .update_xaxes(matches=None)
                        .update_yaxes(range=[35,70])
                        .update_layout(margin=dict(l=20, r=20, t=20, b=20)), #Free x-axis
                        style = {'display':'inline-block','vertical-align':'top','width':'100%','height':200})],
                        id='curve_all_boxplots')
                        ],
                    style = {'display':'inline-block','vertical-align':'top','width':'50%','margin-top':'20px','overflow-y':'scroll'}),
                #Clean ctrl tm div
                html.Div([
                        #Header
                        html.H3('Distribution of "clean" control Tms per plate:',style = {'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'0px'}),
                        #Boxplot
                        html.Div([dcc.Graph(id = 'clean_control_boxes',
                        figure=px.violin(df_results[(df_results['Well_type']=='Control')], x='Assay_Plate', y='Final_Tm',height = 200,
                        labels={'Assay_Plate': 'Assay Plate', 'Final_Tm': 'Melting temp'}) 
                        .update_xaxes(matches=None)
                        .update_yaxes(range=[35,70])
                        .update_layout(margin=dict(l=20, r=20, t=20, b=20)), #Free x-axis
                        style = {'display':'inline-block','vertical-align':'top','width':'100%','height':200})],
                        id='curve_clean_boxplots')
                        ],
                    style = {'display':'inline-block','vertical-align':'top','width':'50%','margin-top':'20px','overflow-y':'scroll'}), 
                html.Hr(style = {'margin-bottom':'0px','width': '100%'}), #Break           
                #Curves
                html.H3('Control curves:',style = {'font-family': 'Arial','margin-bottom':'10px', 'margin-top':'20px'}),
                html.Div([dcc.Graph(id = 'control_curves',
                    figure=px.scatter(df_curves[(df_curves['Well_type']=='Control')&(df_curves['Subplot']!='Original')], x='Temps', y='Smooth Fluorescence', color='Final_decision', color_discrete_sequence=['#FABF73','#F06D4E','#3489eb'],
                        hover_name='Well', hover_data={'Final_decision':False, 'Assay_Plate':False, 'Temps':False, 'Smooth Fluorescence':False, 'Error':True, 'Ctrl_Tm_z-score':True}, #Hover data (Tooltip in Spotfire)
                        facet_col='Assay_Plate', facet_col_wrap=2, facet_col_spacing=0.06, #Facet plots by plate and only allow 2 columns. Column spacing had to be adjusted to allow for individual y-axes
                        render_mode = 'auto', height = 250*(len(df_curves['Assay_Plate'].unique()))) #Height of plots is equal to half the number of plates (coz 2 columns) with each plot 500px high. Width will have to be adjusted
                    .update_yaxes(matches=None) #Free y-axis
                    .update_xaxes(matches=None) #Free x-axis
                    .for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True)) #Put y-axis labels on all facet plots
                    .for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True)), #Put x-axis labels on all facet plots
                    style = {'display':'inline-block','vertical-align':'top', 'width': '100%','overflow-y':'scroll'})
                ],id='curve_graphs')],
            label='Control overview', style=red_tab_style, selected_style=tab_selected_style, value ='tab-2'),

            #MELTING TEMP PLATE OVERVIEW
            dcc.Tab([
                html.Hr(style = {'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}), #Break
                #Child 1: Plate maps of melting temps (Left of screen)
                html.Div([dcc.Graph(id = 'Tm_plates',
                    figure=px.scatter(df_results, x='Column', y='Row', color='Final_Tm', hover_name='Compound', 
                        color_continuous_scale=px.colors.sequential.Oranges, 
                        hover_data={'Final_Tm':True,'Column':False, 'Row':False, 'Well':True, 'Error':True,'Unique_key':True, 'Well_zscore':True}, 
                        facet_col='Spotfire_Platename', facet_col_wrap=1, facet_row_spacing=0.03, 
                        render_mode = 'auto', height = 500*(len(df_results['Spotfire_Platename'].unique())))
                    .update_traces(marker={'symbol': 'square', 'size':25, 'line':{'color':'grey','width':1}})
                    .update_yaxes(autorange='reversed',showgrid=False, zeroline=False) #Get y-axis to go from A to Z, and remove all gridlines
                    .update_xaxes(showgrid=False, zeroline=False) #Remove all gridlines
                    .update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)'}) #Make plot backgrounds transparent
                    .update_layout(coloraxis_colorbar = {'len':0.5, 'xanchor':'center', 'orientation':'h','ypad':5,'thicknessmode':'pixels', 'thickness':20})
                    .for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True)) #Put y-axis labels on all facet plots
                    .for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True)) #Put x-axis labels on all facet plots
                    )],id='Tm_graphs',style = {'display':'inline-block','vertical-align':'top', 'width': '50%','overflow-y':'scroll', 'height':500}), #Plate graph style
                #Child2: Selected curve (Right of screen)
                html.Div([],id='popup',style = {'display':'inline-block','vertical-align':'top', 'width': '50%'})
                ],label='Melting temp overview', style=purple_tab_style, selected_style=tab_selected_style, value ='tab-3'), #Melting temp tab style and other details
            
            #Z-SCORE TAB
            dcc.Tab([
                html.Hr(style = {'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}), #Break
                #Child 1: Plate maps of z-scores (Left of screen)
                html.Div([dcc.Graph(id = 'Z_plates',
                    figure=px.scatter(df_results, x='Column', y='Row', color='Well_zscore', hover_name='Compound', 
                        color_continuous_scale=px.colors.diverging.Spectral_r, color_continuous_midpoint=0,
                        hover_data={'Final_Tm':True,'Column':False, 'Row':False, 'Well':True, 'Error':True,'Unique_key':True, 'Well_zscore':True}, 
                        facet_col='Spotfire_Platename', facet_col_wrap=1, facet_row_spacing=0.03, 
                        render_mode = 'auto', height = 500*(len(df_results['Spotfire_Platename'].unique())))
                    .update_traces(marker={'symbol': 'square', 'size':25, 'line':{'color':'grey','width':1}})
                    .update_yaxes(autorange='reversed',showgrid=False, zeroline=False) #Get y-axis to go from A to Z, and remove all gridlines
                    .update_xaxes(showgrid=False, zeroline=False) #Remove all gridlines
                    .update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)'}) #Make plot backgrounds transparent
                    .update_layout(coloraxis_colorbar = {'len':0.5, 'xanchor':'center', 'orientation':'h','ypad':5,'thicknessmode':'pixels', 'thickness':20})
                    .for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True)) #Put y-axis labels on all facet plots
                    .for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True)) #Put x-axis labels on all facet plots
                    )],id='Z_graphs',style = {'display':'inline-block','vertical-align':'top', 'width': '50%','overflow-y':'scroll', 'height':500}), #Plate graph style
                #Child2: Selected curve (Right of screen)
                html.Div([],id='popup2',style = {'display':'inline-block','vertical-align':'top', 'width': '50%'})
                ],label='Z-score overview', style=blue_tab_style, selected_style=tab_selected_style, value ='tab-4'), #Melting temp tab style and other details
            
            #HIT TABLE TAB
            dcc.Tab([
                html.Hr(style = {'margin-bottom':'0px','width': '100%', 'margin-top':'10px'}), #Break
                #Child A: Data table 
                html.Div([
                        #Data table
                        html.Div([
                        html.Label("Total hits: "+str(default_datatable_count), style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'font-size': '20px', 'color':'dark green'}, id = 'hits_header'),
                        dash_table.DataTable(
                            default_datatable.to_dict('records'), 
                            [{'name': i, 'id': i} for i in default_datatable.loc[:, ['Source_Plate','Well','Subplot','Compound','Fraction','Final_Tm','Well_zscore','Relative_amplitude','Max_ctrl_zscore_for_plate','Unique_key']]],
                            id = 'results_table_datatable',
                            hidden_columns=['Unique_key','Max_ctrl_zscore_for_plate'],
                            css=[{'selector': '.show-hide', 'rule': 'display: none'},{'selector':'.export','rule': 'margin:5px'}],
                            row_deletable=True,
                            sort_action='native',
                            export_format='xlsx',
                            style_data_conditional=data_table_style_data_conditional, 
                            style_table = {'height':'250px','overflow-y':'scroll'},
                            style_as_list_view=True, style_cell={'fontSize':12, 'font-family':'Arial'}, style_header = {'backgroundColor': '#a6ed6f','fontWeight': 'bold'}),#Styling of table
                        ], id = 'results_table'),
                        #Generate all hit graphs button
                        html.Div([html.Button('Generate all graphs', id='generate', n_clicks=None, style = {'margin-top':'20px', 'margin-right': '20px'})], style = {'display':'inline-block','vertical-align':'top'}),
                        #Clear all graphs button
                        html.Div([html.Button('Clear all graphs', id='clear', n_clicks=0, style = {'margin-top':'20px'})],style = {'display':'inline-block','vertical-align':'top'}),
                        ], 
                    id='left_panel',
                    style = {'display':'inline-block','vertical-align':'top', 'width': '50%','overflow-y':'scroll','overflow-x':'scroll', 'height':'50%', 'padding':'10px','margin-top':'5px'}), #Styling of table container
                
                #Child B: Cut off boxes and pop up graph
                html.Div([
                    #Child 1
                    html.Div([
                                #Child 1.1: Main z-score header
                                html.H4('% beyond ctrl z-score range:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px'}),#
                                #Child 1.2: Upper z-score cutoff header and box
                                html.Div([
                                        html.H5('Cutoff:',style = {'font-family': 'Arial','margin-bottom':'5px','margin-top':'7px'}),#
                                        dcc.Input(id='relative_zscore_cutoff', name = 'Cutoff', type='number', value = 0,style = {'font-family': 'Arial','width':'50px','margin-top':'5px'})],#
                                    style = {'display':'inline-block','vertical-align':'top', 'margin-top':'0px', 'margin-bottom':'5px'}),#Child 1.2 style
                            ],style = {'display':'inline-block','vertical-align':'top', 'width':'49%', 'margin-top':'0px', 'margin-bottom':'15px'}), #Child 1 style
                    #Child 2
                    html.Div([
                            #Child 2.1: Main amplitude header
                            html.H4('Relative amplitude cutoffs:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px'}),
                            #Child 2.2.: Upper amplitude cutoff header and box
                            html.Div([
                                    html.H5('Upper cutoff:',style = {'font-family': 'Arial','margin-bottom':'5px','margin-top':'7px'}),#
                                    dcc.Input(id='upper_limit_amp', name = 'Upper cutoff', type='number', value = 1.5,style = {'font-family': 'Arial','width':'50px','margin-top':'5px'})],#
                                style = {'display':'inline-block','vertical-align':'top', 'margin-top':'0px', 'margin-bottom':'5px'}),#Child 2.2 style
                            #Child 2.3: Lower amplitude cutoff header and box
                            html.Div([
                                    html.H5('Lower cutoff:',style = {'font-family': 'Arial', 'margin-bottom':'5px','margin-top':'7px'}),#
                                    dcc.Input(id='lower_limit_amp', name = 'Lower cutoff', type='number', value = 0.5, style = {'font-family': 'Arial', 'width':'50px','margin-top':'5px'})],#
                                style = {'display':'inline-block','vertical-align':'top', 'margin-top':'0px', 'margin-bottom':'5px','margin-left': '30px'}), #Child 2.3 style
                            ],style = {'display':'inline-block','vertical-align':'top', 'width':'49%', 'margin-top':'0px', 'margin-bottom':'5px'}),#Child 2 style
                        #Submit button
                    #Child 3: Submit button
                    html.Div([html.Button('Update', id='submit', n_clicks=0)]),
                    #Child 4: Horizontal break
                    html.Hr(style = {'margin-bottom':'0px','width': '100%'}),
                    #Child 5: POPUP GRAPH
                    html.Div([
                        html.Div([],id = 'popup3',style = {'display':'inline-block','vertical-align':'top', 'width':'50%'}),
                        html.Div([],id = 'popup4',style = {'display':'inline-block','vertical-align':'top', 'width':'50%'}),
                        html.Div([html.Button('Clear', id='clear_popup', n_clicks=0)],id='clear_pop_up_div',hidden=True)
                        ],id='popup3_div',style = {'width': '100%', 'margin-top':'0px', 'padding':'5px'})
                        ],
                    id='control_boxes', #Child B ID
                    style = {'display':'inline-block','vertical-align':'top', 'width': '45%', 'overflow-x':'scroll', 'margin-top': '20px'}), #Child B style
                
                #Child C: All graphs for hits
                html.Div([],id = 'all_graphs_div', style = {'height':'500px','overflow-y':'scroll', 'width': '45%'})

                ],label='Hit list', style=green_tab_style, selected_style=tab_selected_style, value ='tab-5') #Tab div
            ],id ='tabs',value = 'tab-0',vertical =False, style = {'width':'100%'}) #All tabs style and other details
])

####################
#
# FUNCTIONS
#
####################
def generate_selected_curve(selected_unique_key):
    selected_curve_data = df_curves[df_curves['Unique_key'] == selected_unique_key]
    plate = selected_unique_key.split('_')[0]+'_'+selected_unique_key.split('_')[1]+'_'+selected_unique_key.split('_')[2]
    well = selected_unique_key.split('_')[3]
    source_plate = selected_curve_data['Source_Plate'].unique()[0]
    original_selected = selected_curve_data[selected_curve_data['Subplot'] == 'Original']
    subplots_selected = selected_curve_data[selected_curve_data['Subplot'] != 'Original']
    selected_curve_figure=px.scatter(subplots_selected, x='Temps', y='Smooth Fluorescence', color='Subplot', color_discrete_sequence=["#64960c","#abd95b","#4d7804"],
                 labels={'Temps': 'Temperature', 'Smooth Fluorescence': 'Normalized fluorescence'},
                hover_name='Well', hover_data={'Final_decision':False, 'Assay_Plate':False, 'Temps':False, 'Smooth Fluorescence':False, 'Error':True}, #Hover data (Tooltip in Spotfire)
                render_mode = 'auto', height = 500, title = 'Selected: Plate '+source_plate+', Well '+well+"<br><sup> Assay Plate: "+plate+"</sup>")
    selected_curve_figure.add_scatter(x=original_selected['Temps'], y=original_selected['Smooth Fluorescence'],line={'color':'grey','dash':'dot'}, name = 'Original data')
    selected_curve_figure.update_layout(title = {'font':{'size':15}})
    return selected_curve_figure

def generate_first_derivative_curve(selected_unique_key):
    selected_curve_data = df_curves[df_curves['Unique_key'] == selected_unique_key]
    avg_ctrl_Tm = df_results[df_results['Unique_key'] == selected_unique_key]['Avg_ctrl_melting_temp'].unique()[0]
    well_Tm = df_results[df_results['Unique_key'] == selected_unique_key]['Final_Tm'].unique()[0]
    plate = selected_unique_key.split('_')[0]+'_'+selected_unique_key.split('_')[1]+'_'+selected_unique_key.split('_')[2]
    well = selected_unique_key.split('_')[3]
    source_plate = selected_curve_data['Source_Plate'].unique()[0]
    raw_x = selected_curve_data["Temps"].values.tolist()
    raw_y = selected_curve_data["Smooth Fluorescence"].values.tolist()
    subplot = selected_curve_data["Subplot"].values.tolist()
    y_grad = list(np.gradient(raw_y, raw_x))
    deriv_df = pd.DataFrame({'Temps': raw_x,'Smooth': raw_y,'1st_deriv': y_grad,'Subplot':subplot})
    #Drawing the plot
    # 1. Seprate out data (original to be drawn as line, others as points)
    original_1st_deriv = deriv_df[deriv_df['Subplot'] == 'Original']
    subplots_1st_deriv = deriv_df[deriv_df['Subplot'] != 'Original']
    # 2. Draw the 'original' curve as a line
    fig = px.scatter(subplots_1st_deriv, x='Temps', y='1st_deriv', color='Subplot', 
                  color_discrete_sequence=["#64960c","#abd95b","#4d7804"],
                  labels={'Temps': 'Temperature', '1st_deriv': 'DF/DT'})
    # 3. Add in spliced subplots
    fig.add_scatter(x=original_1st_deriv['Temps'], y=original_1st_deriv['1st_deriv'],line={'color':'grey','dash':'dot'}, name = 'Original data')
    # 4. Add Avg ctrl Tm line
    fig.add_vline(x=avg_ctrl_Tm, line=dict(color='red', width=1, dash='dash'), name = 'Avg. ctrl Tm for plate')
    # 5. Add well final Tm
    fig.add_vline(x=well_Tm, line=dict(color='blue', width=1, dash='dash'), name = 'Tm for selected well')
    # 6. Reorder figure data so the line is at the back
    fig.data = (fig.data[1],fig.data[0])
    # 7. Add labels
    fig.update_layout(title = "First derivative <br><sup> <span style='color:red'>Red:</span> Avg. ctrl Tm for plate, <span style='color:blue'>Blue:</span> Well Tm", xaxis_title='Temperature',yaxis_title='DF/DT')
    return fig

def generate_all_curves(unique_key_list):
    hit_curve_data_df = df_curves[df_curves['Unique_key'].isin(unique_key_list)]
    num_rows = math.ceil(len(unique_key_list)/4)
    height = 50*num_rows
    facet_row_max_spacing = 6/height
    all_curves_figure = px.scatter(hit_curve_data_df, x='Temps', y='Smooth Fluorescence', color='Subplot', color_discrete_sequence=['#FABF73','#F06D4E'],
                        hover_name='Well', hover_data={'Final_decision':False, 'Assay_Plate':False, 'Temps':False, 'Smooth Fluorescence':False, 'Error':True, 'Ctrl_Tm_z-score':True}, #Hover data (Tooltip in Spotfire)
                        facet_col='Unique_key', facet_col_wrap=2, facet_col_spacing=0.08,facet_row_spacing = facet_row_max_spacing,#Facet plots by plate and only allow 2 columns. Column spacing had to be adjusted to allow for individual y-axes
                        render_mode = 'auto', height = height) #Height of plots is equal to half the number of plates (coz 2 columns) with each plot 300px high. Width will have to be adjusted
    return all_curves_figure


####################
#
# CALLBACKS
#
####################

#Render curve graph for selected data in Melting temp tab
@app.callback(
    Output('popup', 'children'),
    [Input('Tm_plates', 'clickData')]
)
def update_popup(clickData):
    if clickData:
        selected_unique_key = clickData['points'][0]['customdata'][3]
        curve_figure = generate_selected_curve(selected_unique_key)
        graph_object = dcc.Graph(id = 'selected_curve',figure=curve_figure)
        return graph_object
    else:
        return None

#Render curve graph for selected data in Z score tab
@app.callback(
    Output('popup2', 'children'),
    [Input('Z_plates', 'clickData')]
)
def update_popup2(clickData):
    if clickData:
        selected_unique_key = clickData['points'][0]['customdata'][3]
        curve_figure = generate_selected_curve(selected_unique_key)
        graph_object = dcc.Graph(id = 'selected_curve',figure=curve_figure)
        return graph_object
    else:
        return None


#Highlight entire row in plate report when selected
@app.callback(
    Output('plate_report_table', 'style_data_conditional'),
    [Input('plate_report_table', 'active_cell')],
    State('plate_report_table', 'data')
)
def update_plate_report_row_color(active,table_data):
    style_plate = plate_report_style_data_conditional.copy()
    if active:
        style_plate.append({'if': {'row_index': active['row']}, 'backgroundColor': '#FFFBB2','border': '1px solid #FFFBB2'},)
        seleted_row_data = table_data[active['row']]
        return style_plate

#Highlight entire row in results Datatable when selected and produce pop up graph
@app.callback(
    [Output('results_table_datatable', 'style_data_conditional'),
    Output('popup3', 'children'),
    Output('popup4', 'children'),
    Output('clear_pop_up_div','hidden')],
    [Input('results_table_datatable', 'active_cell'),
    Input('clear_popup','n_clicks')],
    State('results_table_datatable', 'data')
)
def update_selected_row_color(active,n_clicks, table_data):
    style = data_table_style_data_conditional.copy()
    ctx = callback_context
    if not ctx.triggered:
        raise PreventUpdate
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if trigger_id == 'clear_popup':
        if n_clicks > 0:
            style.append({'if': {'row_index': active['row']}, 'backgroundColor': '#cef786','border': '1px solid green'},)
            return style, None, None, True
    elif trigger_id =='results_table_datatable':
        if active:
            style.append({'if': {'row_index': active['row']}, 'backgroundColor': '#cef786','border': '1px solid green'},)
            seleted_row_data = table_data[active['row']]
            selected_unique_key = seleted_row_data['Unique_key']
            curve_figure = generate_selected_curve(selected_unique_key)
            curve_figure.update_layout(title={'yref': 'paper','y' : 1,'yanchor' : 'bottom','font':{'size':12}}, title_pad_b = 25, margin=dict(l=20, r=20, t=50, b=20), height = 200)
            curve_figure.update_layout(font = {'size':8}, showlegend=False)
            graph_object = dcc.Graph(id = 'selected_curve',figure=curve_figure, style = {'fontSize':12, 'font-family':'Arial'})
            first_deriv_curve = generate_first_derivative_curve(selected_unique_key)
            first_deriv_curve.update_layout(title={'yref': 'paper','y' : 1,'yanchor' : 'bottom','font':{'size':12}}, title_pad_b = 25, margin=dict(l=20, r=20, t=50, b=20), height = 200)
            first_deriv_curve.update_layout(font = {'size':8}, showlegend=False)
            graph_object2 = dcc.Graph(id = 'first_deriv_curve',figure=first_deriv_curve, style = {'fontSize':12, 'font-family':'Arial'})
            return style, graph_object, graph_object2, False
        else:
            return style, None, None, True
    raise PreventUpdate

#Generating or clearing for all hits
@app.callback(
    Output('all_graphs_div', 'children'),
    [Input('generate', 'n_clicks'),
     Input('clear', 'n_clicks')],
    [State('results_table_datatable', 'data')])
def update_or_clear_graphs(generate_clicks, clear_clicks, current_table_data):
    ctx = callback_context
    if not ctx.triggered:
        raise PreventUpdate
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if trigger_id == 'generate':
        if generate_clicks > 0:
            key_list = [d['Unique_key'] for d in current_table_data]
            if len(key_list) > 50:
                cropped_list = key_list[:50]
                new_figure = generate_all_curves(cropped_list)
                new_figure.update_yaxes(matches=None)
                new_figure.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
                new_figure.update_layout(height=200 * len(cropped_list)/2, title = "WARNING: Too many figures to plot. <br><sup> The first 50 have been generated. <br>  ", title_font_color="red") #Each plot is 400px high (i.e. 200 is half of 400)
                graph_object = dcc.Graph(id='all_graphs_object', figure=new_figure, style={'font-family': 'Arial'})
                return graph_object
            else:
                new_figure = generate_all_curves(key_list)
                new_figure.update_yaxes(matches=None)
                new_figure.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
                new_figure.update_layout(height=200 * len(key_list)/2) #Each plot is 400px high (i.e. 200 is half of 400)
                graph_object = dcc.Graph(id='all_graphs_object', figure=new_figure, style={'font-family': 'Arial'})
                return graph_object
    elif trigger_id == 'clear':
        if clear_clicks > 0:
            return None
    raise PreventUpdate

#Update hit table based on input box values
@app.callback(
    Output('results_table', 'children'),
    Input('submit', 'n_clicks'),
    State('relative_zscore_cutoff', 'value'),
    State('lower_limit_amp', 'value'),
    State('upper_limit_amp', 'value'),
    prevent_initial_call=True
)
def update_hit_table(n_clicks, zscore_cutoff, lower_value_amp, upper_value_amp):
    if n_clicks is None:
        raise PreventUpdate
    else:
        original_df  = df_results[(df_results['Well_type'] != 'Control')&(df_results['Plate_status'] == 'OK')&(df_results['Error'] == '')]
        small_df = original_df.loc[:, ['Source_Plate','Well','Subplot','Compound','Fraction','Final_Tm','Well_zscore','Unique_key', 'Relative_amplitude','Max_ctrl_zscore_for_plate']]
        small_df = small_df.rename(columns = {'Well_zscore':'zscore'})
        cutoff = zscore_cutoff/100
        small_df["Hit"] = np.where(small_df['zscore'].abs() > 
                        (small_df['Max_ctrl_zscore_for_plate'].abs()+
                         (small_df['Max_ctrl_zscore_for_plate'].abs()*cutoff))
                        ,"Hit","Not a hit")
        filtered_df = small_df[((small_df.Hit == 'Hit'))&((small_df.Relative_amplitude <= upper_value_amp)&(small_df.Relative_amplitude >= lower_value_amp))]
        filtered_df.sort_values('zscore',ascending=False, inplace= True)
        filtered_df = filtered_df.rename(columns = {'zscore':'Well_zscore'})
        filtered_df_count = filtered_df['Unique_key'].nunique()
        header = html.Label("Total hits: "+str(filtered_df_count), style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px', 'font-size': '20px', 'color':'dark green'}, id = 'hits_header')
        data_table = dash_table.DataTable(filtered_df.to_dict('records'), [{'name': i, 'id': i} for i in filtered_df.columns], 
                id = 'results_table_datatable',
                hidden_columns=['Unique_key','Max_ctrl_zscore_for_plate','Hit'], 
                style_as_list_view=True, 
                row_deletable=True, 
                sort_action='native',
                export_format='xlsx',
                style_table = {'height':'250px','overflow-y':'scroll'},
                style_cell={'fontSize':12, 'font-family':'Arial'}, 
                style_header = {'backgroundColor': '#a6ed6f','fontWeight': 'bold'})
    return header, data_table


#Callback to fix figures auto adjusting to tiny
#Taken from: https://github.com/plotly/dash-core-components/issues/922
@app.callback(
    [Output('control_curves','figure'),
    Output('Tm_plates','figure'),
    Output('Z_plates','figure')],
    Input('tabs','value'),
    State('control_curves','figure'),
    State('Tm_plates','figure'),
    State('Z_plates','figure')
)

def controlGraphSizes(tabValue,figure1,figure2, figure3):
    tab2_height = 250*(len(df_curves['Assay_Plate'].unique()))
    tab3_height = 500*(len(df_curves['Assay_Plate'].unique()))
    tab4_height = 500*(len(df_curves['Assay_Plate'].unique()))
    if tabValue == 'tab-1':
        return dash.no_update
    elif tabValue == 'tab-2':
        figure1['layout']['height'] = tab2_height
        return [figure1,dash.no_update,dash.no_update]
    elif tabValue == 'tab-3':
        figure2['layout'].update(plot_bgcolor = '#FFFFFF',paper_bgcolor = '#FFFFFF', height = tab3_height, width = '50vw', margin=dict(l=30, r=30, t=50, b=0))
        return [dash.no_update,figure2,dash.no_update]
    elif tabValue == 'tab-4':
        figure3['layout'].update(plot_bgcolor = '#FFFFFF',paper_bgcolor = '#FFFFFF', height = tab4_height, width = '50vw', margin=dict(l=30, r=30, t=50, b=0))
        return [dash.no_update,dash.no_update,figure3]
    elif tabValue == 'tab-5':
        return dash.no_update
    else:
        raise PreventUpdate

# Run the app
if __name__ == '__main__':
    app.run(host ='0.0.0.0',debug=True)

