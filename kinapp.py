import dash
import dash_table
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State

import dash_cytoscape as cyto

import pickle
from scipy import sparse
import pandas as pd
import re

import diffusion
from kinapp_helper import InputValidator

cyto.load_extra_layouts()
# network as a coo matrix
network = pickle.load(open('results/kinase_matrix.pkl', 'rb'))

def make_nodes(network):
    """
    Constructs cytoscape elements (node) dict from COO matrix
    """
    elements = []
    node_i, node_j = network.nonzero()
    track = []
    for index, label_i in enumerate(node_i):
        label_j = node_j[index] 
        if label_i not in track:
            elements.append({'data': {'id': label_i, 'label': label_i}})
            track.append(label_i)
        if label_j not in track:
            elements.append({'data': {'id': label_j, 'label': label_j}})
            track.append(label_j)
        elements.append({'data': {'source': label_i, 'target': label_j}})
        #if index > 500:
        #    break
    print(len(elements))
    return elements


app = dash.Dash(__name__)

app.layout = html.Div(
    [
#        cyto.Cytoscape(
#            id='kin-map',
#            layout={'name': 'circle'},
#            #layout={'name': 'cose-bilkent', 'numIter': 50, 'padding': 30,
#            #        'nodeDimensionsIncludeLabels': 'true'},
#            #layout={'name': 'cola', 'numIter': 10},
#            #layout={'name': 'klay', 'numIter': 10},
#            #layout={'name': 'dagre', 'numIter': 10},
#            style={'width': '100%', 'height': '1200px'},
#            elements=make_nodes(network)
#        ),
        
        html.Div(
            dcc.Textarea(
                id='kinase-list',
                placeholder='Enter kinase IDs in HUGO format (ex: CDC7, AURAB)',
                style={'margin':'auto', 'width': '75%', 'height': '20%'}
                #value='ABC'
            ),
            style={'display':'flex', 'justifyContent':'center', 'padding':'500px 0px 0px 0px'}
        ),

        html.Div(
            html.Button('SUBMIT', id='submit-button', n_clicks=None),
            style={'display':'flex', 'justifyContent':'center'}
        ),
        
        html.Div(id='response', children="some text"),

        html.Div(id='output', style={'display':'none'})        
    ]
)


def parse_input(text):
    re.findall('\w+', text)


@app.callback([
    Output('output', 'children'),
    Output('output', 'style'),
    Output('response', 'children')],
    [Input('submit-button', 'n_clicks')],
    [State('kinase-list', 'value')]
)
def print_output(n_clicks, value):
    if n_clicks is not None and value is not None:
        #validate input
        labels = re.findall('\w+', value)
        valid_labels = InputValidator().validate(labels)
        if len(valid_labels) == 0:
            msg = """Enter some valid kinase IDs to diffuse. 
                    The app expects HUGO (hgnc) symbols. See https://www.genenames.org/"""
            return None, None, msg
        
        dif = diffusion.Diffusion(network, 0)
        res = dif.diffuse()
        print(res)
        df = pd.DataFrame({'a':[1,2,3], 'b':[1,2,3]})    
        
        msg = f"Diffusing from nodes # {', '.join([str(x) for x in valid_labels])}"

        return dash_table.DataTable(
                id='hi',
                columns=[{"name": i, "id": i} for i in df.columns],
                data = df.to_dict('records'),
                export_columns='all', export_format='csv'
            ), {'display':'flex', 'justifyContent':'center'}, msg
    
    return None, None, None

if __name__ == '__main__':
    app.run_server(debug=True)
