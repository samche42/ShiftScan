import dash
from dash import Dash, html, dcc, Input, Output, State, dash_table
import pandas as pd
import plotly.express as px
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

# Initialize the app
app = Dash(__name__)


#Fiddle with input data to correct format
df_results = pd.read_csv('Final_results.txt', sep='\t', header = 0)
df_results['Error'].fillna('', inplace=True)
df_results['No. Std Dev']=df_results['No. Std Dev'].apply(lambda x:round(x,2))
df_curves = pd.read_csv('Final_curves.txt', sep='\t', header = 0)
df_curves['Error'].fillna('', inplace=True)


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

data_table_style_data_conditional = [
    {
        'if': {'state': 'active'},
        'backgroundColor': '#f3e3ff',
        'border': '1px solid purple',
    },
    {
        'if': {'state': 'selected'},
        'backgroundColor': '#f3e3ff',
        'border': '1px solid purple',
    },
]

####################
#
# GENERATING DEFAULT DATA TABLE (TAB5)
#
####################

default_datatable = df_results[(df_results['Well_type'] != 'Control')&(df_results['Plate_status'] == 'OK')&(df_results['Error'] == '')]
default_datatable = default_datatable.rename(columns = {'No. Std Dev':'zscore'})
default_datatable = default_datatable[(default_datatable.zscore <= -3)|(default_datatable.zscore >= 3)]
default_datatable.sort_values('zscore',ascending=False, inplace= True)
default_datatable = default_datatable.rename(columns = {'zscore':'No. Std Dev'})

####################
#
# WEBPAGE LAYOUT
#
####################

app.layout = html.Div([html.H1('DSF Analysis Visualizations', style = {'background': '#e3e3e3','font-family': 'Arial','border-radius': '4px','border': f'1px solid #70706d','padding':'20px'}),
    dcc.Tabs([
            #Home tab. Currently empty
            dcc.Tab([],label='Home', style=yellow_tab_style, selected_style=tab_selected_style, value='tab-1'),

            #CURVE TAB
            dcc.Tab([
                html.Div([dcc.Graph(id = 'control_curves',
                    figure=px.scatter(df_curves[(df_curves['Well_type']=='Control')&(df_curves['Subplot']!='Original')], x='Temps', y='Smooth Fluorescence', color='Final_decision', color_discrete_sequence=['#5CBF12','#F06D4E'],
                        hover_name='Well', hover_data={'Final_decision':False, 'Assay_Plate':False, 'Temps':False, 'Smooth Fluorescence':False, 'Error':True}, #Hover data (Tooltip in Spotfire)
                        facet_col='Assay_Plate', facet_col_wrap=2, facet_col_spacing=0.06, #Facet plots by plate and only allow 2 columns. Column spacing had to be adjusted to allow for individual y-axes
                        render_mode = 'auto', height = 250*(len(df_curves['Assay_Plate'].unique()))) #Height of plots is equal to half the number of plates (coz 2 columns) with each plot 500px high. Width will have to be adjusted
                    .update_yaxes(matches=None) #Free y-axis
                    .update_xaxes(matches=None) #Free x-axis
                    .for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True)) #Put y-axis labels on all facet plots
                    .for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True)), #Put x-axis labels on all facet plots
                    style = {'display':'inline-block','vertical-align':'top', 'width': '100%','overflow-y':'scroll'})
                ],id='curve_graphs')],
            label='Control curve overview', style=orange_tab_style, selected_style=tab_selected_style, value ='tab-2'),

            #MELTING TEMP PLATE OVERVIEW
            dcc.Tab(
                [
                #Child 1: Plate maps of melting temps (Left of screen)
                html.Div([dcc.Graph(id = 'Tm_plates',
                    figure=px.scatter(df_results, x='Column', y='Row', color='Final_Tm', hover_name='Compound', 
                        color_continuous_scale=px.colors.sequential.Oranges, 
                        hover_data={'Final_Tm':True,'Column':False, 'Row':False, 'Well':True, 'Error':True,'Unique_key':True, 'No. Std Dev':True}, 
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
                ],label='Melting temp overview', style=red_tab_style, selected_style=tab_selected_style, value ='tab-3'), #Melting temp tab style and other details
            
            #Z-SCORE TAB
            dcc.Tab(
                [
                #Child 1: Plate maps of z-scores (Left of screen)
                html.Div([dcc.Graph(id = 'Z_plates',
                    figure=px.scatter(df_results, x='Column', y='Row', color='No. Std Dev', hover_name='Compound', 
                        color_continuous_scale=px.colors.diverging.Spectral_r, color_continuous_midpoint=0,
                        hover_data={'Final_Tm':True,'Column':False, 'Row':False, 'Well':True, 'Error':True,'Unique_key':True, 'No. Std Dev':True}, 
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
                ],label='Z-score overview', style=purple_tab_style, selected_style=tab_selected_style, value ='tab-4'), #Melting temp tab style and other details
            
            #HIT TABLE TAB
            dcc.Tab(
                [
                #Child 1: Data table 
                html.Div([dash_table.DataTable(
                        default_datatable.to_dict('records'), 
                        [{'name': i, 'id': i} for i in default_datatable.loc[:, ['Source_Plate','Well','Subplot','Compound','Fraction','Final_Tm','No. Std Dev','Unique_key']]],
                        id = 'results_table_datatable',
                        hidden_columns=['Unique_key'],
                        css=[{'selector': '.show-hide', 'rule': 'display: none'},{'selector':'.export','rule': 'margin:5px'}],
                        row_deletable=True,
                        sort_action='native',
                        export_format='xlsx',
                        style_data_conditional=data_table_style_data_conditional, 
                        style_as_list_view=True, style_cell={'fontSize':12, 'font-family':'Arial'}, style_header = {'backgroundColor': '#eed9ff','fontWeight': 'bold'})], #Styling of table
                    id='results_table',
                    style = {'display':'inline-block','vertical-align':'top', 'width': '50%','overflow-y':'scroll','overflow-x':'scroll', 'height':400, 'padding':'10px','margin-top':'5px'}), #Styling of table container
                
                #Child 2: Cut off boxes and pop up graph
                html.Div([
                        #Header
                        html.H4('Define hit z-score cutoffs here:',style = {'font-family': 'Arial','margin-bottom':'1px', 'margin-top':'0px'}),#
                        #Hit limit inputs
                        html.Div([ #Upper cutoff layout component (Holds header and input box)
                                html.H5('Upper cutoff:',style = {'font-family': 'Arial','margin-bottom':'5px','margin-top':'7px'}),#
                                dcc.Input(id='upper_limit', name = 'Upper cutoff', type='number', value = 3,style = {'font-family': 'Arial','width':'30px','margin-top':'5px'})],#
                            style = {'display':'inline-block','vertical-align':'top', 'width':'30%', 'margin-top':'0px', 'margin-bottom':'15px'}),
                        html.Div([ #lower cutoff layout component (Holds header and input box)
                                html.H5('Lower cutoff:',style = {'font-family': 'Arial', 'margin-bottom':'5px','margin-top':'7px'}),#
                                dcc.Input(id='lower_limit', name = 'Lower cutoff', type='number', value = -3, style = {'font-family': 'Arial', 'width':'30px','margin-top':'5px'})],#
                            style = {'display':'inline-block','vertical-align':'top', 'width':'30%', 'margin-top':'0px', 'margin-bottom':'15px'}),
                        #Submit button
                        html.Div([html.Button('Update', id='submit', n_clicks=0)]),
                        #Horizontal break followed by submit button
                        html.Hr(style = {'margin-bottom':'0px'}),
                        #POPUP GRAPH
                        html.Div([],id='popup3',style = {'width': '100%', 'margin-top':'0px'})
                        ],
                    id='control_boxes',
                    style = {'display':'inline-block','vertical-align':'top', 'width': '40%','overflow-x':'scroll', 'height':400,  'padding':'20px'}),
                ],label='Hit list', style=blue_tab_style, selected_style=tab_selected_style, value ='tab-5') #Melting temp tab style and other details
            ],id ='tabs',vertical =False, style = {'width':'100%'}) #All tabs style and other details
])

#FUNCTIONS
def generate_selected_curve(selected_unique_key):
    selected_curve_data = df_curves[df_curves['Unique_key'] == selected_unique_key]
    plate = selected_unique_key.split('_')[0]+'_'+selected_unique_key.split('_')[1]+'_'+selected_unique_key.split('_')[2]
    well = selected_unique_key.split('_')[3]
    selected_curve_figure=px.scatter(selected_curve_data, x='Temps', y='Smooth Fluorescence', color='Subplot', color_discrete_sequence=px.colors.qualitative.Vivid,
                hover_name='Well', hover_data={'Final_decision':False, 'Assay_Plate':False, 'Temps':False, 'Smooth Fluorescence':False, 'Error':True}, #Hover data (Tooltip in Spotfire)
                render_mode = 'auto', height = 500, title = 'Selected: Plate '+plate+', Well '+well)
    return selected_curve_figure


#CALLBACKS

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

#Update hit table based on input box values
@app.callback(
    Output('results_table', 'children'),
    Input('submit', 'n_clicks'),
    State('lower_limit', 'value'),
    State('upper_limit', 'value'),
    prevent_initial_call=True
)
def update_hit_table(n_clicks, lower_value, upper_value):
    if n_clicks is None:
        raise PreventUpdate
    else:
        original_df  = df_results[(df_results['Well_type'] != 'Control')&(df_results['Plate_status'] == 'OK')&(df_results['Error'] == '')]
        small_df = original_df.loc[:, ['Source_Plate','Well','Subplot','Compound','Fraction','Final_Tm','No. Std Dev','Unique_key']]
        small_df = small_df.rename(columns = {'No. Std Dev':'zscore'})
        filtered_df = small_df[(small_df.zscore <= float(lower_value))|(small_df.zscore >= float(upper_value))]
        filtered_df.sort_values('zscore',ascending=False, inplace= True)
        filtered_df = filtered_df.rename(columns = {'zscore':'No. Std Dev'})
    return dash_table.DataTable(filtered_df.to_dict('records'), [{'name': i, 'id': i} for i in filtered_df.columns], 
                id = 'results_table_datatable',
                hidden_columns=['Unique_key'], 
                style_as_list_view=True, 
                row_deletable=True, 
                sort_action='native',
                export_format='xlsx',
                style_cell={'fontSize':12, 'font-family':'Arial'}, 
                style_header = {'backgroundColor': '#eed9ff','fontWeight': 'bold'})

#Highlight entire row in Datatable when selected and produce pop up graph
@app.callback(
    [Output('results_table_datatable', 'style_data_conditional'),
    Output('popup3', 'children')],
    [Input('results_table_datatable', 'active_cell')],
    State('results_table_datatable', 'data')
)
def update_selected_row_color(active,table_data):
    style = data_table_style_data_conditional.copy()
    if active:
        style.append({'if': {'row_index': active['row']}, 'backgroundColor': '#f3e3ff','border': '1px solid purple'},)
        seleted_row_data = table_data[active['row']]
        selected_unique_key = seleted_row_data['Unique_key']
        curve_figure = generate_selected_curve(selected_unique_key)
        curve_figure.update_layout(title={'yref': 'paper','y' : 1,'yanchor' : 'bottom'}, title_pad_b = 10, margin=dict(l=20, r=20, t=40, b=20), height = 300)
        graph_object = dcc.Graph(id = 'selected_curve',figure=curve_figure, style = {'fontSize':10, 'font-family':'Arial'})
        return style, graph_object
    else:
        return style, None


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
    app.run(debug=True)

