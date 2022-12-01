# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:37:24 2022

@author: Hasan
"""
from dash import Dash, dcc, html, Input, Output
import plotly.express as px
import pandas as pd
import math as m
import plotly.graph_objects as go

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
df_url='https://docs.google.com/spreadsheets/d/e/2PACX-1vShcA7rp1YWrsac3P_QcPPpjzvOFjanT9egjRevOX4OezGGsLPuRJXSrDT5meQzvMeJN3XbQVpn2xIt/pub?output=xlsx'
app = Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

df = pd.read_excel(df_url,sheet_name='InputSep')
dfMonth = pd.read_excel(df_url,sheet_name='Sheet3')
dfVar = pd.read_excel(df_url,sheet_name='VarDist')
M = pd.read_excel(df_url,sheet_name='Months')['Month'].tolist()
MI = pd.read_excel(df_url,sheet_name='Months')['MonthIndex'].tolist()
start, stop, MonthInit,h, duration = 1,1273, 15, 300, 80
DN = ['N-Terminal Domain', 'RBD', 'Post-RBD','Subunit-2']
DR = [[1,325],[326,525],[526,725],[726,1273]]
colors = {
    'background': 'white123',
    'text': 'gold',
    'text2': 'black',
}
mark1 = {DR[i][1]: {'label': str(DR[i][1]), 'style': {'color': colors['text2']}} for i in range(len(DN))}
mark2 = {m.ceil( (DR[i][0]+DR[i][1])/2 ): {'label': DN[i], 'style': {'color': colors['text2']}} for i in range(len(DN))}
mark1.update(mark2)

fig = px.scatter(df, x="Sites", y="MonthIndex",
                  size="delemus_score", hover_name="Residues",
                  hover_data=['Residues','Sites','delemus_score','SAPs','Variants','MonthIndex','Month','Domain'],
                  size_max=15,width=1500, height=900,
                  marginal_x="histogram",
                  color = 'Domain', 
                  # animation_frame="MonthIndex", animation_group="Residues",#pd.Series(['N-Terminal Domain','RBD','Post-RBD','Subunit-2','VoC/VoI Sites'])
                  # range_y=[1-2,34+2],range_x=[1-2.5,1273+2.5],
                  )

figVar = go.Figure()
for col in dfVar.columns:
    if col != 'Month':
        figVar.add_trace(go.Scatter(x=dfVar['Month'], y=dfVar[col],
                        line_shape='spline', name=col))

app.layout = html.Div(style={'backgroundColor': colors['background']}, children=[
    html.H1(
        children='Dynamic Leading Mutations in Spike Glycoprotein of SARS-CoV-2',
        style={
            'font_family': 'Fantasy',
            'font_size': '40px',
            'textAlign': 'center',
            'color': colors['text']
    }),
    html.H1(children='(DeLeMus)',
        style={
            'textAlign': 'center',
            'color': colors['text'],
            'font_size': '40px',
            'font_family': 'cursive',
    }),
    html.Br(),
    dcc.Slider(
        dfMonth['MonthIndex'].min(),
        dfMonth['MonthIndex'].max(),
        step = None,
        value = MonthInit,
        id = 'month--slider',
        marks = {str(MI[i]): {'label':str(M[i]) , 'style': {"transform": "rotate(15deg)",'color': colors['text2']}} for i in range(len(dfMonth))}
    ),
    
    html.Div(style={'backgroundColor': colors['background'],'width': '78%', 'display': 'inline-block', 'padding': '0 20'},children=[
        dcc.Graph(
            id = 'monthly-leading-sites',
            hoverData={'points': [{'customdata':['D614','G(11027),N(4),S(1)','α,β,γ',MonthInit]}]}
        )
    ]),
    
    html.Div([
        dcc.Graph(id='SAP-Dist')
    ], style={'width': '22%','display': 'inline-block'}),
    
    html.Div(dcc.RangeSlider(
        id='range-slider',
        min=start, max=stop, step=1,
        marks=mark1,
        tooltip={"placement": "bottom", "always_visible": True},
        value=[start, 1234])
    , style={'width': '90%', 'padding': '0px 20px 20px 20px'}
    ),
    
    html.H1(
        children='Leading Mutations Map in Spike Glycoprotein of SARS-CoV-2',
        style={
            'font_family': 'cursive',
            'font_size': '10px',
            'textAlign': 'center',
            'color': colors['text']
    }),
    html.Div([
        dcc.Graph(
            id = 'All-leading-sites',
            figure = fig)
    ]),

    html.H1(
        children='Variant Distributions of SARS-CoV-2',
        style={
            'font_family': 'cursive',
            'font_size': '10px',
            'textAlign': 'center',
            'color': colors['text']
    }),
    
    html.Div([
        dcc.Graph(
            id = 'Var-Dist',
            figure = figVar)
    ])
])

@app.callback(
    Output('monthly-leading-sites', 'figure'),
    Input('month--slider', 'value'),
    Input('range-slider', 'value'),
    )

def update_figure(selected_Month,slider_range):
    indexMonth = dfMonth[(dfMonth['MonthIndex'] == selected_Month)].index.values
    Title = dfMonth.at[indexMonth[0],'Month']
    filtered_df = df.loc[df['MonthIndex'] == selected_Month]
    low, high = slider_range
    fig = px.scatter(filtered_df,
                x = 'Sites',
                y = 'Month',
                color='Domain', #size='delemus_score', 
                hover_name = 'Residues',
                size_max=15,
                hover_data=['Residues','SAPs','Variants','MonthIndex'],
                range_x=[low, high], height=h,
                title=Title)
    fig.update_layout(transition_duration=duration) #, hovermode='closest'
    fig.update_yaxes(title='')
    # fig.add_vline(x=614, line_width=1.3, line_dash="dash",line_color="black",row = 1)
    return fig
    
@app.callback(
    Output('SAP-Dist', 'figure'),
    Input('monthly-leading-sites', 'hoverData'),
    )
    
def update_SAP_Dist(hoverData):
    Residues = hoverData['points'][0]['customdata'][0]
    SAPs = hoverData['points'][0]['customdata'][1]
    Variants = hoverData['points'][0]['customdata'][2]
    MonthIndex = hoverData['points'][0]['customdata'][3]
    index = df[(df['Residues'] == Residues)&
               (df['SAPs'] == SAPs)&
               (df['Variants'] == Variants)&
               (df['MonthIndex'] == MonthIndex)].index.values
    SAP0 = df.at[index[0],'SAPs']
    SAP =SAP0.replace(")","").replace("(","").split(",")
    if "(" in SAP0:
        AA = [mut[0] for mut in SAP]
        Freq_log10 = [m.log10(int(mut[1:])+1) for mut in SAP] 
        title = '<b>{}</b>'.format(Residues)
        dfSAP = pd.DataFrame({'AA':AA,'#Mutations (LogScale)':Freq_log10})
        fig = px.bar(dfSAP,x = 'AA',y = '#Mutations (LogScale)', height=h)
        fig.update_layout(transition_duration=duration, hovermode='closest')
        fig.add_annotation(x=0,y=1,xanchor='right',yanchor='bottom',
                        xref='paper',yref='paper',showarrow = False,align='right',
                        text=title)
        return fig
    else:
        AA = ['G','N','S']
        Freq_log10 = [m.log10(mut) for mut in [1102,4,1]]
        title = '<b>{}</b>'.format('D614')
        dfSAP = pd.DataFrame({'AA':AA,'#Mutations (LogScale)':Freq_log10})
        fig = px.bar(dfSAP,x = 'AA',y = '#Mutations (LogScale)', height=h)
        fig.update_layout(transition_duration=duration, hovermode='closest')
        fig.add_annotation(x=0,y=1,xanchor='right',yanchor='bottom',
                        xref='paper',yref='paper',showarrow = False,align='right',
                        text=title)
        return fig

if __name__ == '__main__':
    app.run_server(debug=True)