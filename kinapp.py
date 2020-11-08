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
import color_gradient

pd.options.display.float_format = '${:,.2f}'.format

cyto.load_extra_layouts()

image_url = 'url(https://image.freepik.com/free-photo/elegant-black-handmade-technique-aquarelle_23-2148300751.jpg)'

# network as a coo matrix
network = pickle.load(open('results/kinase_matrix.pkl', 'rb'))

def make_nodes(network):
    """
    Constructs cytoscape elements (node) dict from COO matrix
    """
    elements = []
    node_i, node_j = network.adj_matrix.nonzero()
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
        if index > 50:
            break
    print(len(elements))
    return elements

def make_nodes2(sim_matrix, labels):
    """
    Constructs cytoscape elements (node) dict from COO matrix

    {'data':{'id': label_i, 'label': label_i}}
    {'data':{'source': label_i, 'target': label_j}}
    """

    elements = []
    node_added = []
    for label in labels:
        edges = sim_matrix.get_edges_for_protein(label)
        for edge in edges:
            edge_notation = {'data': {'source': label, 'target': edge}}
            elements.append(edge_notation)
        node_added += [node for node in [label] + edges if node not in node_added]

    for node in node_added:
        islabel=0
        if node in labels:
            islabel=1
            print(node)
        elements.append({'data': {'id': node, 'label': node, 'label_flag':islabel}})

    return elements

def create_cytoscape_div(sim_matrix, labels, diffusion_result):
    """Constructs updated cytoscape div."""
    elements = make_nodes2(sim_matrix, labels)
    cyto_div = cyto.Cytoscape(id='kin-map',
                              layout = {'name': 'circle'},
                              #layout = {'name': 'klay'}, # Good-ish
                              #layout = {'name': 'spread'}, # BAD
                              #layout = {'name': 'cose'}, # this one is okay-ish
                              #layout = {'name': 'breadthfirst', 'roots': '#CDK1, #AURKB, #AURKA'}, # this one is good
                              #layout={'name': 'cose-bilkent', 'numIter': 50, 'padding': 30,
                              #        'nodeDimensionsIncludeLabels': 'true'},
                              style={'width':'75%', 'height': '500px', 'display': 'flex'},
                              stylesheet=get_cyto_stylesheet(diffusion_result),
                              elements = elements)
    return cyto_div

def get_cyto_stylesheet(diffusion_result):
    stylesheet=[
        {
            'selector': 'node',
            'style': {
                'label': 'data(label)'
            }
        },
        {
            'selector':'[label_flag=1]',
            'style':{'shape':'diamond', 'background-color':'red', 'line-color':'black', 'border-color':'black', 'border-width':'2'}
        },
        {
            'selector':'[label_flag=0]',
            'style':{'border-width':'1', 'border-color':'black'}
        }
        #{
        #    'selector':'[label_flag=0]',
        #    'style':{'background-color':'MediumPurple'}#'rgba(135, 206, 250, 0.5)'}#'gray'}#08306b'}
        #},
        #{
        #    'selector':'[label_flag=0]',
        #    'style':{'shape':'triangle'}#'rgba(135, 206, 250, 0.5)'}#'gray'}#08306b'}
        #}
    ]
    print(diffusion_result.zscore)
    colors = get_colors(diffusion_result)
    color_selectors = make_selectors_colors(colors)
    print(color_selectors)
    return stylesheet + color_selectors

def make_selectors_colors(colors):

    selectors = []
    for protein, color in colors.items():
        selector = {'selector':'[label="%s"]' % protein,
                    'style':{'background-color':color}}
        selectors.append(selector)
    return selectors

def get_colors(diffusion_result):
    xcolor_gradient = color_gradient.ColorGradientGenerator()
    xcolor_gradient.create_color_map2(base_color=[128,128,128])
    #xcolor_gradient.create_color_map(palette='Greens')
    colors = xcolor_gradient.map_colors(diffusion_result.zscore)
    color_dict = dict(zip(diffusion_result.protein, colors))
    print(color_dict)
    return color_dict

app = dash.Dash(__name__)
app.title = 'GGid'

app.layout = html.Div(
        #style={'background-image': 'url(https://upload.wikimedia.org/wikipedia/commons/2/22/North_Star_-_invitation_background.png)'},
    children=[
        html.Div(
            dcc.Textarea(
                id='kinase-list',
                placeholder='Enter kinase IDs in HUGO format (ex: CDC7, AURAB)',
                style={'margin':'auto', 'width': '75%', 'height': '20%'},
                value='CDK1'
            ),
            style={'display':'flex', 'justifyContent':'center', 'padding':'100px 0px 0px 0px'}
        ),
        html.Div(
            html.Button('SUBMIT', id='submit-button', n_clicks=None),
            style={'display':'flex', 'justifyContent':'center'}
        ),
        html.Div(id='kin-map-container',
            children=cyto.Cytoscape(
                id='kin-map',
                #layout={'name': 'circle'},
                layout={'name': 'cose-bilkent', 'numIter': 50, 'padding': 30,
                        'nodeDimensionsIncludeLabels': 'true'},
                #layout={'name': 'cola', 'numIter': 10},
                #layout={'name': 'klay', 'numIter': 10},
                #layout={'name': 'dagre', 'numIter': 10},
                style={'width': '75%', 'height': '500px', 'display':'flex'},
                elements=make_nodes(network)
            ),
            style={'justifyContent':'center', 'display':'flex'}
        ),
        html.Div(id='response', children="some text"),
        html.Div(id='output', style={'display':'none'})
    ]
    #style={'background-image':image_url, 'background-size':'cover', 'height':'100%', 'position':'fixed', 'width':'100%', 'top':'0px', 'left':'0px'}
)


def parse_input(text):
    re.findall('\w+', text)


#from dash_table.Format import Format, Group, Scheme, Symbol
#frmt = Format(scheme=Scheme.fixed, precision=2, group=Group.yes,
#                groups=3,
#                group_delimiter='.',
#                decimal_delimiter=',')


import dash_table.FormatTemplate as FormatTemplate
from dash_table.Format import Format, Scheme, Sign, Symbol

@app.callback([
    Output('output', 'children'),
    Output('output', 'style'),
    Output('response', 'children'),
    Output('kin-map-container', 'children')],
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
            return None, None, msg, None

        dif = diffusion.Diffusion(network, valid_labels)
        res = dif.diffuse()

        diffusion_result = res.get_as_pandas_df_without_labels()
        new_net = create_cytoscape_div(network, valid_labels, diffusion_result)
        print(new_net)
        msg = f"Diffusing from nodes # {', '.join([str(x) for x in valid_labels])}"
        print(diffusion_result.head())

        return ([dash_table.DataTable(
                id='hi',
                columns=[{"name": i, "id": i} for i in diffusion_result.columns],
                data = diffusion_result.to_dict('records'),
                export_columns='all',
                export_format='csv',
                page_size=20)], {'display':'flex', 'justifyContent':'center'}, msg, [new_net])

    return None, None, None, None

if __name__ == '__main__':
    app.run_server(debug=True)
