import pickle
import re
import time

import dash
import dash_bootstrap_components as dbc
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


def get_cross_validation_result(labeled_kinases):
    """Does LOO validation and returns container with formatted results."""
    loo_experiment = cross_validation.LOOValitation(network, labeled_kinases)
    zscore_table = loo_experiment.run_validation()
    # get ROC figure
    tpr, fpr, auc = loo_experiment.get_roc()
    fig = draw_roc_curve(tpr, fpr)
    # define results div
    result_div = html.Div(
        children=[
            html.Details(
                children=[
                    html.Summary(
                        "Cross-validation AUC is %1.2f (expand for details)" % auc
                    ),
                    html.Ul(
                        children=[
                            html.Li(
                                """
                                AUC ~= 0.50: information does not diffuse between the
                                given kinases.
                                This is because the kinases are not closely connected in
                                the network.
                                Either the network is not modeling their
                                relationship properly, or these kinases are not
                                very similar (i.e. do not share low-level GO terms).
                                You should diffuse each kinase on its own
                                to find its functional neighbors.
                                """
                            ),
                            html.Li(
                                """
                                0.60 < AUC < 0.80: there is some signal! The kinases in
                                the set diffuse some information to each other,
                                meaning they are somewhat connected in the network.

                                Alternatively, there could
                                be two independent kinase clusters of in your set.
                                Kinases within the same cluster are predictive of each
                                other, but not of the kinases in the other cluster.
                                Each cluster may be responsible for some unique
                                sub-function relevant to the overall process.
                                Identify the clusters by looking at the graph visual
                                above, and diffuse each cluster on its own to find
                                its functional neighbors.
                                """
                            ),
                            html.Li(
                                """
                                AUC > 0.80: kinases are strongly predictive of each
                                other. They are well connected in the network, and thier
                                neighbors may be involved in the same type of process.
                                For the list of their closest neigbors see table below.
                                """
                            ),
                        ]
                    ),
                    dcc.Graph(id="loo-roc", figure=fig),
                ],
            ),
            convert_to_dash_table(zscore_table),
        ]
    )
    # make updated graph
    graph_nodes, node_styling = create_cytoscape_div(
        network, labeled_kinases, zscore_table
    )
    return result_div, graph_nodes, node_styling


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
    trace = go.Scatter(x=fpr, y=tpr)
    fig.add_trace(trace)
    fig.update_layout(width=500, height=500)
    fig.layout.title = "Leave-one-out cross-validation ROC"
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
    top_hits = diffusion_result[diffusion_result.zscore >= 1]
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
    mask = diffusion_result.zscore >= 1
    diffusion_result = diffusion_result.loc[mask, :]
    xcolor_gradient = color_gradient.ColorGradientGenerator()
    xcolor_gradient.create_color_map2(base_color=[128, 128, 128])
    # xcolor_gradient.create_color_map(palette='Greens')
    colors = xcolor_gradient.map_colors(diffusion_result.zscore)
    color_dict = dict(zip(diffusion_result.protein, colors))
    # print(color_dict)
    return color_dict


def get_main_tab():
    """Returns main tab content."""

    return [
        html.Div(
            dbc.Textarea(
                id="kinase-list",
                placeholder="Enter kinase IDs in HUGO format (ex: CDC7, AURAB)",
                style={"margin": "auto", "width": "75%", "height": "20%"},
                value="CDK1",
            ),
            style={
                "display": "flex",
                "justifyContent": "center",
                "padding": "50px 0px 0px 0px",
            },
        ),
        html.Div(
            dbc.Button("SUBMIT", id="submit-button", color="dark", n_clicks=None),
            style={"display": "flex", "justifyContent": "center", "padding-top": "5px"},
        ),
        html.Div(
            dbc.Checklist(
                id="loo-switch",
                options=[
                    {
                        "label": "cross-validate (for 2+ kinases)",
                        "value": "on",
                    }
                ],
                value=[],  # when you click checkmark its value is added to this list
            ),
            style={"display": "flex", "justifyContent": "center"},
        ),
        dbc.FormGroup(
            [
                dbc.Label("Graph layout"),
                dbc.RadioItems(
                    id="network-layout-toggle",
                    options=[
                        {"label": "Circle", "value": "circle"},
                        {"label": "Concentric", "value": "concentric"},
                        {"label": "Spread", "value": "spread"},
                        {"label": "Force-directed", "value": "cola"},
                    ],
                    value="circle",
                    inline=True,
                ),
            ]
        ),
        html.Div(
            id="kin-map-container",
            children=[
                # dbc package broke zoom and panning, so I had to use
                # raw pixels to set size of network area, and
                # specify pan and zoom manually -- not a fan!
                cyto.Cytoscape(
                    id="kin-map",
                    layout={"name": "circle", "fit": True},
                    style={
                        "height": "600px",
                        "width": "1000px",
                        "border": "1px solid grey",
                    },
                    zoom=0.1,
                    elements=make_nodes(network),
                    pan={"x": 500, "y": 300},
                ),
            ],
            # style={"justifyContent": "center", "display": "flex"},
        ),
        html.Div(id="status-line"),  # style={"width": "100%", "margin": "auto"}),
        html.Div(id="results-container", children=[]),
        # line below is a hidden utility switch for chaining callbacks
        # if user presents valid kinase list for diffusion, the switch
        # value is set to 'valid' and diffusion callback fires
        dcc.RadioItems(
            id="diffusion-switch", value="invalid", style={"display": "none"}
        ),
    ]


def get_theory_tab():
    """Returns content of theory tab."""
    return "This is the theory tab."


def get_example_tab():
    """Returns contents of the example tab."""
    return "This is the example tab."


def parse_input(text):
    re.findall("\w+", text)


def create_results_div():
    """Wraps diffusion results as children for the results-container."""
    pass


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "GGid"

tab_main = dbc.Tab(label="Diffusion Tool", children=[])
tab_theory = dbc.Tab(label="How it works", children=[])
tab_example = dbc.Tab(label="Example", children=[])

tab_main.children = get_main_tab()
tab_theory.children = get_theory_tab()
tab_example.children = get_example_tab()

app.layout = dbc.Container(
    id="main-container",
    children=[
        dbc.CardHeader(
            html.H1([dbc.Badge("GGID: Human Kinome", color="dark")]),
            style={"background-color": "#343a40"},
        ),
        html.Br(),
        dbc.Container(
            id="content",
            children=[
                html.P(
                    "Diffuse information over Gene Ontology-derived network of human kinases."
                ),
                dbc.Tabs(children=[tab_main, tab_theory, tab_example]),
            ],
        ),
        html.Br(),
        html.Br(),
        dbc.CardFooter(
            id="footer",
            children=[
                html.P(
                    "Ilya Novikov 2020 // github: BLAH BLAH",
                    style={
                        "text-align": "center",
                        "color": "white",
                        "padding-top": "1rem",
                    },
                )
            ],
            style={"background-color": "#343a40"},
        ),
    ],
)


@app.callback(
    [
        Output("status-line", "children"),
        Output("diffusion-switch", "value"),
    ],
    [Input("submit-button", "n_clicks")],
    [State("kinase-list", "value")],
    prevent_initial_call=True,
)
def validate_inputs(n_clicks, protein_list):
    """Validates protein list submitted by user."""
    message = []  # html.H4("Experiment Result:")]
    # diffusion switch, by default, is set to return no-update signal (ie: do nothing)
    diffusion_switch = dash.no_update
    proteins = re.findall("\w+", protein_list)

    if len(proteins) == 0:
        message.append(
            dbc.Alert("You haven't entered anything in the textbox", color="danger")
        )
    else:
        valid_kinases = InputValidator().validate(proteins)
        invalid_kinases = [kinase for kinase in proteins if kinase not in valid_kinases]
        if len(valid_kinases) == 0:
            message.append(
                dbc.Alert(
                    """
                    None of the strings you entered were recognized as valid kinase ids.
                    The app expects HUGO (hgnc) symbols for human kinases.
                    See list here https://www.genenames.org/
                    """,
                    color="danger",
                )
            )
        if len(valid_kinases) > 0:
            message.append(
                dbc.Alert(
                    "Diffusing kinases {kinases}".format(
                        kinases=", ".join(valid_kinases)
                    ),
                    color="success",
                )
            )
            diffusion_switch = "valid"
        if len(invalid_kinases) > 0:
            message.append(
                dbc.Alert(
                    "These items were not recognized as human kinases: {kins}".format(
                        kins=", ".join(invalid_kinases)
                    ),
                    color="dark",
                ),
            )
    print(message)
    # wrap text in Markdown element before returning;
    # block diffusion switch from updating (and triggering diffusion)
    # if value of switch hasn't changed to "valid"
    return message, diffusion_switch


@app.callback(
    [
        Output("results-container", "children"),
        Output("kin-map", "elements"),
        Output("kin-map", "stylesheet"),  # 'stylesheet' is different from 'style'
    ],
    [Input("diffusion-switch", "value")],
    [State("kinase-list", "value"), State("loo-switch", "value")],
    prevent_initial_call=True,
)
def diffuse(diffusion_switch, protein_list, loo_switch):
    """Conducts diffusion experiment with provided kinases."""
    proteins = re.findall("\w+", protein_list)
    valid_kinases = InputValidator().validate(proteins)
    print("state of the diffusion switch:", diffusion_switch)
    print("diffusion started.")
    if "on" in loo_switch:
        # averaged post-diffusion results produced via LOO validation
        result_div, graph_nodes, node_styling = get_cross_validation_result(
            valid_kinases
        )
    else:
        # plain diffusion, with no LOO validation and averaging
        zscore_table, graph_nodes, node_styling = get_diffusion_result(valid_kinases)
        result_div = [convert_to_dash_table(zscore_table)]

    return result_div, graph_nodes, node_styling


@app.callback(
    Output("kin-map", "layout"),
    Input("network-layout-toggle", "value"),
)
def change_graph_layout(layout_type):
    """Changes layout of the network graph."""
    print("network layout radio item executed", time.time())
    layout = {"name": layout_type}
    return layout


if __name__ == "__main__":
    app.run_server(debug=True)
