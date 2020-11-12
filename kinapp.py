import pickle
import re
import time

import dash
import dash_core_components as dcc
import dash_cytoscape as cyto
import dash_html_components as html
import dash_table
import dash_table.FormatTemplate as FormatTemplate
import pandas as pd
import plotly.graph_objects as go
from dash.dependencies import Input, Output, State
from dash_table.Format import Format, Scheme, Sign, Symbol
from scipy import sparse

import color_gradient
import cross_validation
import diffusion
from kinapp_helper import InputValidator

pd.options.display.float_format = "${:,.2f}".format

cyto.load_extra_layouts()

image_url = "url(https://image.freepik.com/free-photo/elegant-black-handmade-technique-aquarelle_23-2148300751.jpg)"

# network as a coo matrix
network = pickle.load(open("network/kinase_matrix.pkl", "rb"))


def get_diffusion_result(labeled_kinases):
    """Returns formatted diffusion results."""
    experiment = diffusion.Diffusion(network, labeled_kinases)
    result = experiment.diffuse()
    # make z-score table
    zscore_table = result.get_as_pandas_df_without_labels()
    # make updated graph
    graph_nodes, node_styling = create_cytoscape_div(
        network, labeled_kinases, zscore_table
    )
    return zscore_table, graph_nodes, node_styling


def convert_to_dash_table(zscore_table):
    """Converts pandas zscore table to Dash data table."""

    dash_table_ = dash_table.DataTable(
        id="hi",
        columns=[{"name": i, "id": i} for i in zscore_table.columns],
        data=zscore_table.to_dict("records"),
        export_columns="all",
        export_format="csv",
        page_size=20,
    )
    return dash_table_


def draw_roc_curve(tpr, fpr):
    """Creates plotly scatter for the ROC curve."""
    fig = go.Figure()
    trace = go.Scatter(x=tpr, y=fpr)
    fig.add_trace(trace)
    fig.update_layout(width=500, height=500)
    return fig


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
            elements.append({"data": {"id": label_i, "label": label_i}})
            track.append(label_i)
        if label_j not in track:
            elements.append({"data": {"id": label_j, "label": label_j}})
            track.append(label_j)
        elements.append({"data": {"source": label_i, "target": label_j}})
    # print(len(elements))
    return elements


def make_nodes2(sim_matrix, labels, proteins=[]):
    """
    Constructs cytoscape elements (node) dict from COO matrix

    {'data':{'id': label_i, 'label': label_i}}
    {'data':{'source': label_i, 'target': label_j}}
    """

    elements = []
    node_added = []
    for label in labels + proteins:
        edges = sim_matrix.get_edges_for_protein(label)
        for edge in edges:
            if edge not in proteins + labels:
                continue
            edge_notation = {"data": {"source": label, "target": edge}}
            elements.append(edge_notation)
        # node_added += [node for node in [label] + edges if node not in node_added]
    node_added = labels + proteins
    for node in node_added:
        islabel = 0
        if node in labels:
            islabel = 1
            # print(node)
        elements.append({"data": {"id": node, "label": node, "label_flag": islabel}})

    return elements


def create_cytoscape_div(sim_matrix, labels, diffusion_result):
    """Constructs updated cytoscape div."""
    top_hits = diffusion_result[diffusion_result.zscore >= 2]
    top_hits = list(top_hits.protein.values)
    elements = make_nodes2(sim_matrix, labels, top_hits)
    new_elements = elements
    new_style = get_cyto_stylesheet(diffusion_result)
    return new_elements, new_style


def get_cyto_stylesheet(diffusion_result):
    stylesheet = [
        {"selector": "node", "style": {"label": "data(label)"}},
        {
            "selector": "[label_flag=1]",
            "style": {
                "shape": "diamond",
                "background-color": "red",
                "line-color": "black",
                "border-color": "black",
                "border-width": "2",
            },
        },
        {
            "selector": "[label_flag=0]",
            "style": {"border-width": "1", "border-color": "black"},
        }
        # {
        #    'selector':'[label_flag=0]',
        #    'style':{'background-color':'MediumPurple'}#'rgba(135, 206, 250, 0.5)'}#'gray'}#08306b'}
        # },
        # {
        #    'selector':'[label_flag=0]',
        #    'style':{'shape':'triangle'}#'rgba(135, 206, 250, 0.5)'}#'gray'}#08306b'}
        # }
    ]
    # print(diffusion_result.zscore)
    colors = get_colors(diffusion_result)
    color_selectors = make_selectors_colors(colors)
    # print(color_selectors)
    return stylesheet + color_selectors


def make_selectors_colors(colors):

    selectors = []
    for protein, color in colors.items():
        selector = {
            "selector": '[label="%s"]' % protein,
            "style": {"background-color": color},
        }
        selectors.append(selector)
    return selectors


def get_colors(diffusion_result):
    xcolor_gradient = color_gradient.ColorGradientGenerator()
    xcolor_gradient.create_color_map2(base_color=[128, 128, 128])
    # xcolor_gradient.create_color_map(palette='Greens')
    colors = xcolor_gradient.map_colors(diffusion_result.zscore)
    color_dict = dict(zip(diffusion_result.protein, colors))
    # print(color_dict)
    return color_dict


app = dash.Dash(
    __name__, external_stylesheets=["https://codepen.io/chriddyp/pen/bWLwgP.css"]
)
app.title = "GGid"

app.layout = html.Div(
    # style={'background-image': 'url(https://upload.wikimedia.org/wikipedia/commons/2/22/North_Star_-_invitation_background.png)'},
    children=[
        html.Div(
            dcc.Textarea(
                id="kinase-list",
                placeholder="Enter kinase IDs in HUGO format (ex: CDC7, AURAB)",
                style={"margin": "auto", "width": "75%", "height": "20%"},
                value="CDK1",
            ),
            style={
                "display": "flex",
                "justifyContent": "center",
                "padding": "100px 0px 0px 0px",
            },
        ),
        html.Div(
            html.Button("SUBMIT", id="submit-button", n_clicks=None),
            style={"display": "flex", "justifyContent": "center"},
        ),
        html.Div(
            dcc.Checklist(
                id="loo-switch",
                options=[{"label": "conduct LOO validation", "value": "loo_true"}],
                value=[],
            ),
            style={"display": "flex", "justifyContent": "center"},
        ),
        html.Div(
            id="kin-map-container",
            children=[
                cyto.Cytoscape(
                    id="kin-map",
                    layout={"name": "circle"},
                    style={"width": "75%", "height": "500px", "display": "flex"},
                    elements=make_nodes(network),
                ),
                dcc.RadioItems(
                    id="network-layout",
                    options=[
                        {"label": "Circle", "value": "circle"},
                        {"label": "Klay", "value": "klay"},
                        {"label": "Dagre", "value": "dagre"},
                        {"label": "Cola", "value": "cola"},
                        {"label": "Cose", "value": "cose"},
                        {"label": "Cose-Bilkenet", "value": "cose-bilkent"},
                        {"label": "Euler", "value": "euler"},
                        {"label": "Spread", "value": "spread"},
                        {"label": "Concentric", "value": "concentric"},
                    ],
                    value="circle",
                ),
            ],
            style={"justifyContent": "center", "display": "flex"},
        ),
        html.Div(id="status-line"),
        html.Div(id="results-container", children=[]),
        # line below is a hidden utility switch for chaining callbacks
        # if user presents valid kinase list for diffusion, the switch
        # value is set to 'valid' and diffusion callback fires
        dcc.RadioItems(
            id="diffusion-switch", value="invalid", style={"display": "none"}
        ),
    ]
)


def parse_input(text):
    re.findall("\w+", text)


def create_results_div():
    """Wraps diffusion results as children for the results-container."""
    pass


@app.callback(
    [
        Output("status-line", "children"),
        Output("diffusion-switch", "value"),
    ],
    [Input("submit-button", "n_clicks")],
    [State("kinase-list", "value")],
    prevent_initial_call=True,
)
def check_validity(n_clicks, protein_list):
    """Validates protein list input."""
    message = "Status: "
    diffusion_switch = "invalid"
    proteins = re.findall("\w+", protein_list)

    if len(proteins) == 0:
        message = "You haven't entered anything in the textbox"
    else:
        valid_kinases = InputValidator().validate(proteins)
        invalid_kinases = [kinase for kinase in proteins if kinase not in valid_kinases]
        if len(valid_kinases) == 0:
            message = """
                None of the strings you entered were recognized as valid kinase ids.
                The app expects HUGO (hgnc) symbols for human kinases.
                See list here https://www.genenames.org/
            """
        if len(invalid_kinases) > 0:
            message += (
                "These items were not recognized as human kinases {strings}".format(
                    strings=invalid_kinases
                )
            )
        if len(valid_kinases) > 0:
            message += "Diffusing valid kinases {kinases}".format(kinases=valid_kinases)
            diffusion_switch = "valid"
    return message, diffusion_switch


@app.callback(
    [
        Output("results-container", "children"),
        Output("kin-map", "elements"),
        Output("kin-map", "stylesheet"),  # 'stylesheet' is different from 'style'
    ],
    [Input("diffusion-switch", "value")],
    [State("kinase-list", "value")],
    prevent_initial_call=True,
)
def diffuse(diffusion_switch, protein_list):
    proteins = re.findall("\w+", protein_list)
    valid_kinases = InputValidator().validate(proteins)

    loo_switch = "off"

    if loo_switch == "off":
        print("hello")
        zscore_table, graph_nodes, node_styling = get_diffusion_result(valid_kinases)
    else:
        pass
        # do LOO and average results

    # create results container
    # add table + figure, or just table to
    return [convert_to_dash_table(zscore_table)], graph_nodes, node_styling

    ## print(diffusion_result.head())
    # fig = None
    # if len(valid_labels) > 1 and "loo_true" in loo_switch:
    #    print("loo true")
    #    print(network)
    #    loo_experiment = cross_validation.LOOValitation(network, valid_labels)
    #    res = loo_experiment.run_validation()
    #    tpr, fpr, auc = loo_experiment.get_roc()
    #    fig = draw_roc_curve(tpr, fpr)
    #    msg_ = """Leave-one-out cross-validation result: {loo_auc}""".format(
    #        loo_auc=auc
    #    )
    #    msg += msg_
    # return (
    #    [
    #        dash_table.DataTable(
    #            id="hi",
    #            columns=[{"name": i, "id": i} for i in diffusion_result.columns],
    #            data=diffusion_result.to_dict("records"),
    #            export_columns="all",
    #            export_format="csv",
    #            page_size=20,
    #        )
    #    ],
    #    {"display": "flex", "justifyContent": "center"},
    #    msg,
    #    elements,
    #    stylesheet,
    #    fig,
    # )


@app.callback(
    Output("kin-map", "layout"),
    Input("network-layout", "value"),
)
def change_graph_layout(layout_type):
    """Changes layout of the network graph."""
    print("network layout radio item executed", time.time())
    layout = {"name": layout_type}
    return layout


if __name__ == "__main__":
    app.run_server(debug=True)
