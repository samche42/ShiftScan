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

args = parser.parse_args()

####################
#
# STYLES
#
####################

color1 = '#B4EDD2'
color1_faded = '#e0fff1'
color2 = '#A0CFD3'
color3 = '#8D94BA'
color4 = '#9A7AA0'
color4_faded = '#d1b6d6' 
color5 = '#87677B'
main_header_bkgrd_color = '#e3e3e3'
graph_gray = '#cacacf'

base_style = {
    'border': '1px black',
    'color': 'black',
    'font-size': '14px',
    'font-weight': 'bold',
    'font-family': 'Arial',
    'padding': '6px'}

base_selected_style = {
    'border': '3px solid black',
    'color': 'black',
    'font-size': '14px',
    'font-weight': 'bold',
    'font-family': 'Arial',
    'padding':'6px'}

base_disabled = {
    'border': '1px black',
    'color': 'dark grey',
    'font-size': '14px',
    'font-weight': 'bold',
    'font-family': 'Arial',
    'padding': '6px',
    'background': 'grey'}

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

####################
#
# Read in files and format
#
####################

print("Reading in data. This can take a few minutes if processing many plates")
df_results = pd.read_csv(args.input_dir+'/Only_Tm_values.txt', sep='\t', header = 0)
df_curves = pd.read_csv(args.input_dir+'/Only_Tm_curves.txt', sep='\t', header = 0)

columns_to_round = ['Smooth_Tm','Amplitude','Curve_height']
df_results = df_results.assign(**{col: df_results[col].apply(lambda x: round(x, 2)) for col in columns_to_round})
df_results['Row'] = df_results['Well'].str[0] #Generate Row column (needed to build visualizations)
df_results['Column'] = (df_results['Well'].str[1:]).astype(int) #Generate Column column (needed to build visualizations)

df_results['group'] = df_results['Assay_Plate'].apply(lambda x: df_results['Assay_Plate'].unique().tolist().index(x) // 16) #Assign group number to every 16 plates
groups_df = df_results.loc[:, ['Assay_Plate', 'group']].drop_duplicates()

####################
#
# FUNCTIONS
#
####################

def generate_selected_curve(selected_unique_key):
    selected_curve_data = df_curves[df_curves['Unique_key'] == selected_unique_key]
    plate = selected_unique_key.rsplit('_', 1)[0]
    well = selected_unique_key.rsplit('_', 1)[-1]
    original_selected = selected_curve_data[selected_curve_data['Subplot'] == 'Original']
    subplots_selected = selected_curve_data[(selected_curve_data['Subplot'] != 'Original')]
    selected_curve_figure=px.scatter(subplots_selected, x='Temps', y='Smooth Fluorescence', color='Subplot', 
                color_discrete_sequence=[color3,color4,color2,color1],
                labels={'Temps': 'Temperature', 'Smooth Fluorescence': 'Normalized fluorescence'},
                hover_name='Well',
                render_mode = 'auto', height = 500, title = ', Well '+well+"<br><sup> Assay Plate: "+plate+"</sup>")
    selected_curve_figure.add_scatter(x=original_selected['Temps'], y=original_selected['Smooth Fluorescence'],line={'color':graph_gray,'dash':'dot'}, name = 'Original data')
    selected_curve_figure.update_layout(title = {'font':{'size':15}})
    selected_curve_figure.data = selected_curve_figure.data[::-1]
    return selected_curve_figure

def generate_Tm_heatmap(df, facet_row_max_spacing):
    figure=px.scatter(df, x='Column', y='Row', color='Smooth_Tm', hover_name='Compound', 
                                color_continuous_scale=[color1,color2,color3,color4], color_continuous_midpoint=50,
                                facet_col='Assay_Plate', facet_col_wrap=1, facet_row_spacing=facet_row_max_spacing, 
                                hover_data={'Unique_key':True},
                                render_mode = 'auto')
    return figure

def update_heatmap(figure,facet_row_max_spacing,num_rows):
    figure.update_traces(marker={'symbol': 'square', 'size':25, 'line':{'color':'grey','width':1}})
    figure.update_yaxes(autorange='reversed',showgrid=False, zeroline=False) #Get y-axis to go from A to Z, and remove all gridlines
    figure.update_xaxes(showgrid=False, zeroline=False) #Remove all gridlines
    figure.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)'}) #Make plot backgrounds transparent
    figure.update_layout(coloraxis_colorbar = {'len':0.5, 'xanchor':'center', 'yanchor':'bottom','y':1+facet_row_max_spacing,'orientation':'h','ypad':0.8,'thicknessmode':'pixels', 'thickness':20})
    figure.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True)) #Put y-axis labels on all facet plots
    figure.for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True)) #Put x-axis labels on all facet plots
    tab3_height = 500*num_rows
    figure['layout'].update(height = tab3_height,margin=dict(l=10, r=10, t=10, b=10), width = 700)
    return figure

####################
#
# WEBPAGE LAYOUT
#
####################

app.layout = html.Div([
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
                html.Div([],id='popup',style = {'display':'inline-block','vertical-align':'top', 'width': '50%'})
                ],label='Melting temp overview', style=tab3_style, selected_style=tab3_selected, value ='tab-3'), 

            #HIT TABLE TAB
            dcc.Tab([],label='Hit list', style=tab4_style, selected_style=tab4_selected, disabled_style = base_disabled, value ='tab-4', disabled=True), #Tab div
            ],id ='tabs',value = 'tab-3',vertical =False, style = {'width':'100%'}) #All tabs style and other details
],id = 'page')

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
    [Input('tabs','value')])
def generate_heatmaps(tabValue):
    if tabValue == 'tab-3':
        no_plates = df_results['Assay_Plate'].nunique()
        if no_plates <= 16:
            num_rows = math.ceil(no_plates)
            height = 50*num_rows
            facet_row_max_spacing = 6/height
            figure = generate_Tm_heatmap(df_results, facet_row_max_spacing)
            updated_figure = update_heatmap(figure, facet_row_max_spacing,num_rows)
            return updated_figure
        else:
            list_of_subdataframes = [d for _, d in subdata.groupby('group')] #Separate plates into groups of 16 plates
            chosen_df = list_of_subdataframes[dropdown_choice]
            sub_no_plates = chosen_df['Assay_Plate'].nunique()
            sub_num_rows = math.ceil(sub_no_plates)
            sub_height = 50*sub_num_rows
            facet_row_max_spacing = 6/sub_height
            figure = generate_Tm_heatmap(chosen_df, facet_row_max_spacing)
            updated_figure = update_heatmap(figure, facet_row_max_spacing,sub_num_rows)
            return updated_figure
    else:
        raise PreventUpdate

#Render curve graph for selected data in Melting temp tab
@app.callback(
    Output('popup', 'children'),
    Input('Tm_plates', 'clickData'))
def update_popup(clickData):
    if clickData:
        selected_unique_key = clickData['points'][0]['customdata'][0]
        curve_figure = generate_selected_curve(selected_unique_key)
        graph_object = dcc.Graph(id = 'selected_curve',figure=curve_figure)
        return graph_object
    else:
        return None

# Run the app
if __name__ == '__main__':
    app.run(host ='0.0.0.0',port = args.port, debug=False)