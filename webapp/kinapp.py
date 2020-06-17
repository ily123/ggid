import dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State

import dash_cytoscape as cyto
import pickle
from scipy import sparse

cyto.load_extra_layouts()
# network as a coo matrix
network = pickle.load(open('../../ggid/results/kinase_matrix.pkl', 'rb'))

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
        cyto.Cytoscape(
            id='kin-map',
            layout={'name': 'circle'},
            #layout={'name': 'cose-bilkent', 'numIter': 50, 'padding': 30,
            #        'nodeDimensionsIncludeLabels': 'true'},
            #layout={'name': 'cola', 'numIter': 10},
            #layout={'name': 'klay', 'numIter': 10},
            #layout={'name': 'dagre', 'numIter': 10},
            style={'width': '100%', 'height': '1200px'},
            elements=make_nodes(network)
        ),
        
        html.Div(
            dcc.Textarea(
                id='kinase-list',
                placeholder='Enter kinase IDs (ex: TP53, AURAB)',
                style={'margin':'auto', 'width': '75%', 'height': '20%'}
                #value='ABC'
            ),
            style={'display':'flex'}
        ),

        html.Div(
            html.Button('SUBMIT', id='submit-button', n_clicks=0),
            style={'display':'flex', 'justifyContent':'center'}
        ),

        html.Div(id='dummy', style={'display':'none'})        
    ]
)

@app.callback(
    Output('dummy', 'children'),
    [Input('submit-button', 'n_clicks')],
    [State('kinase-list', 'value')]
)
def print_output(n_clicks, value):
    if n_clicks > 0:
        print(value)
    print(1) 

if __name__ == '__main__':
    app.run_server(debug=True)
