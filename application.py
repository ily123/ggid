import pickle
import re

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_cytoscape as cyto
import dash_html_components as html
import dash_table
import pandas as pd
import plotly.graph_objects as go
from dash.dependencies import Input, Output, State

import app_text
import color_gradient
import cross_validation
import diffusion
from kinapp_helper import InputValidator

# laod data and set some defaults:
network = pickle.load(open("network/kinase_matrix.pkl", "rb"))
pd.options.display.float_format = "{:,.2f}".format
cyto.load_extra_layouts()


# ---
# diffusion and cytoscape logic:


def get_diffusion_result(labeled_kinases, zscore_cutoff):
    """Returns components of the diffusion results."""
    experiment = diffusion.Diffusion(network, labeled_kinases)
    result = experiment.diffuse()
    # make z-score table
    zscore_table = result.get_result_df_with_zscore()
    # make updated graph
    graph_nodes, node_styling = create_cytoscape_div(
        network, labeled_kinases, zscore_table, zscore_cutoff
    )

    return zscore_table, graph_nodes, node_styling


def get_cross_validation_result(labeled_kinases, zscore_cutoff):
    """Does LOO validation and returns container with formatted results."""
    loo_experiment = cross_validation.LOOValitation(network, labeled_kinases)
    zscore_table = loo_experiment.run_validation()
    # get ROC figure
    tpr, fpr, auc = loo_experiment.get_roc()
    roc_fig = draw_roc_curve(tpr, fpr, auc)
    # define results div
    result_div = html.Div(
        children=[
            html.Details(
                children=[
                    html.Summary(
                        "Cross-validation AUC is %1.2f (expand for details)" % auc,
                    ),
                    dcc.Graph(id="loo-roc", figure=roc_fig),
                    html.Ul(
                        children=[
                            html.Li(app_text.auc_explanation["random"]),
                            html.Li(app_text.auc_explanation["moderate"]),
                            html.Li(app_text.auc_explanation["high"]),
                        ]
                    ),
                ],
            ),
            html.Details(
                children=[
                    html.Summary("Full z-score and ranks table (expand for details)"),
                    convert_to_dash_table(zscore_table),
                ]
            ),
        ]
    )
    # make updated graph
    graph_nodes, node_styling = create_cytoscape_div(
        network, labeled_kinases, zscore_table, zscore_cutoff
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


def draw_roc_curve(tpr, fpr, auc):
    """Creates plotly scatter for the ROC curve."""
    fig = go.Figure()
    loo_trace = go.Scatter(
        x=fpr, y=tpr, mode="lines", line_color="red", name="input set (AUC=%1.2f)" % auc
    )
    random_trace = go.Scatter(
        x=[0, 1],
        y=[0, 1],
        mode="lines",
        line_color="gray",
        line_dash="dot",
        line_width=1,
        name="random",
    )
    fig.add_traces([loo_trace, random_trace])
    fig.update_layout(
        width=500,
        height=500,
        legend_x=0.98,
        legend_y=0.02,
        legend_xanchor="right",
        title="Leave-one-out cross-validation of the input proteins",
        title_x=0.5,
        xaxis_title="False Positive Rate",
        yaxis_title="True Positive Rate",
        xaxis_range=[-0.01, 1.01],
        yaxis_range=[-0.01, 1.01],
    )
    return fig


def get_network_elements(network):
    """Constructs cytoscape view of the entire network."""
    # add all proteins and their edges to the view
    elements = []
    protein_name = network.proteins
    nodes_i, nodes_j = network.network.nonzero()
    # iterate over all edges
    for index, node_i_index in enumerate(nodes_i):
        node_j_index = nodes_j[index]
        # get protein names that make up the edge
        node_i_name = protein_name[node_i_index]
        node_j_name = protein_name[node_j_index]
        # construct cytoscape elements
        element_i = {"data": {"id": node_i_name, "label": node_i_name}}
        element_j = {"data": {"id": node_j_name, "label": node_j_name}}
        edge_ij = {"data": {"source": node_i_name, "target": node_j_name}}
        edge_ji = {"data": {"source": node_j_name, "target": node_i_name}}
        # the matrix is symmetrical, so we need to account for redundant edges
        if element_i not in elements:
            elements.append(element_i)
        if element_j not in elements:
            elements.append(element_j)
        if (edge_ij not in elements) and (edge_ji not in elements):
            elements.append(edge_ij)
    return elements


def get_cluster_elements(network, input_proteins, top_hits, diffusion_result):
    """Constructs cytoscape view for the post-diffusion result."""
    # The goal is to make a graph view of the input proteins and
    # their most closely connected post-diffusion hits. (Together,
    # these form a cluster.) We only want to display the cluster
    # proteins and the connections they form to each other.
    elements = []
    cluster = input_proteins + top_hits
    # find connections within cluster
    for protein in cluster:
        edges = network.get_edges_for_protein(protein)
        cluster_edges = [edge for edge in edges if edge in cluster]
        for edge in cluster_edges:
            edge_ij = {"data": {"source": protein, "target": edge}}
            edge_ji = {
                "data": {"source": edge, "target": protein}
            }  # no redundant edges
            if (edge_ij not in elements) and (edge_ji not in elements):
                elements.append(edge_ij)
    # now add nodes to the view, and attach experimental data to the nodes
    zscores = dict(zip(diffusion_result.protein, diffusion_result.zscore))
    ranks = dict(zip(diffusion_result.protein, diffusion_result["rank"]))
    for node in cluster:
        node_notation = {
            "data": {
                "id": node,
                "label": node,
                "input_label_flag": 0 if node not in input_proteins else 1,
                "zscore": None if node not in zscores else zscores[node],
                "rank": None if node not in ranks else ranks[node],
            }
        }
        elements.append(node_notation)
    return elements


def create_cytoscape_div(network, labels, diffusion_result, zscore_cutoff):
    """Constructs updated cytoscape div."""
    # regardless of experiment type, we want to show all label nodes
    # plus all the top-scoring unlabeled nodes
    top_hits = diffusion_result[diffusion_result.zscore >= zscore_cutoff]
    top_hits_nodes = list(top_hits.protein.values)
    graph_elements = get_cluster_elements(
        network, labels, top_hits_nodes, diffusion_result
    )
    # however, the coloring of the nodes needs to be handled differently
    # based on experiment type (if it's LOO we want gradient applied to labels too)
    is_loo = len(network.proteins) == len(diffusion_result)
    if is_loo:
        labels_result = diffusion_result[diffusion_result.initial_state == 1]
        top_hits_and_labels = pd.concat([top_hits, labels_result])
        # to get an even gradient, rerank within the set
        top_hits_and_labels["rank"] = top_hits_and_labels["rank"].rank()
        graph_style = get_cyto_stylesheet(top_hits_and_labels, is_loo)
    else:
        graph_style = get_cyto_stylesheet(top_hits, is_loo)
    return graph_elements, graph_style


def get_cyto_stylesheet(diffusion_result, is_loo):
    """Style nodes according to post-diffusion results."""
    stylesheet = [
        {"selector": "node", "style": {"label": "data(label)"}},
        {
            "selector": "[input_label_flag=1]",
            "style": {
                "shape": "diamond",
                "background-color": "white" if is_loo else "red",
                "line-color": "black",
                "border-color": "black",
                "border-width": "2",
            },
        },
        {
            "selector": "[input_label_flag=0]",
            "style": {"border-width": "1", "border-color": "black"},
        },
    ]
    # color nodes according to their post-diffusion z-score
    # this will overwrite default white of labels in a loo experiment
    protein_colors = get_colors(diffusion_result)
    color_selectors = make_selector_colors(protein_colors)
    return stylesheet + color_selectors


def make_selector_colors(protein_colors):
    """Color proteins according to their post-diffusion rank."""
    selectors = []
    for protein, color in protein_colors.items():
        selector = {
            "selector": '[label="%s"]' % protein,
            "style": {"background-color": color},
        }
        selectors.append(selector)
    return selectors


def get_colors(diffusion_result):
    """Generate color gradient for displayed nodes."""
    xcolor_gradient = color_gradient.ColorGradientGenerator()
    xcolor_gradient.create_color_map2(base_color=[255, 51, 51])
    colors = xcolor_gradient.map_colors(
        1 - diffusion_result["rank"] / len(diffusion_result)
    )
    color_dict = dict(zip(diffusion_result.protein, colors))
    print(diffusion_result)
    print(diffusion_result["rank"])
    print(color_dict)
    return color_dict


def get_graph_layout_options(isdisabled=True):
    """Returns layout options for the radio form.

    Parameters
    ----------
    isdiabled : bool
        If True, turns off the CPU-intensive options
    """
    options = [
        {"label": "Circle", "value": "circle"},
        {"label": "Concentric", "value": "concentric"},
        {
            "label": "Spread",
            "value": "spread",
            "disabled": isdisabled,
        },
        {
            "label": "Force-directed",
            "value": "cola",
            "disabled": isdisabled,
        },
    ]
    return options


def get_main_tab():
    """Returns main tab content."""
    return [
        html.Div(
            dbc.Textarea(
                id="input-kinase-list",
                placeholder="Enter kinase IDs in HUGO format (ex: CDC7, AURKB)",
                value="CDK1, CDC7",
            ),
        ),
        html.Div(
            id="submit-button-div",
            children=dbc.Button(
                "DIFFUSE", id="submit-button", color="dark", n_clicks=None
            ),
        ),
        dbc.Form(
            id="input-toggles-form",
            children=[
                dbc.FormGroup(
                    dbc.Checklist(
                        id="loo-switch",
                        options=[
                            {
                                "label": "cross-validate (for 2+ kinases)",
                                "value": "on",
                            }
                        ],
                        value=[],  # when you checkmark, value is added to this list
                    ),
                    className="mr-3",
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("graph hits with zscore of >="),
                        dbc.Input(id="zscore-cutoff", value=2, type="number"),
                    ],
                ),
            ],
            inline=True,
        ),
        html.Br(),
        html.Div(id="status-line"),
        html.Br(),
        dbc.Form(
            [
                dbc.FormGroup(
                    id="kin-map-layout-selector",
                    children=[
                        dbc.Label("Adjust graph layout"),
                        dbc.RadioItems(
                            id="network-layout-toggle",
                            options=get_graph_layout_options(),
                            value="circle",
                            inline=True,
                        ),
                    ],
                )
            ]
        ),
        html.Hr(),
        html.Div(
            id="kin-map-container",
            children=[
                # DBC package (?) broke zoom and panning, so I had to use
                # raw pixels to set size of network area, and
                # to specify pan and zoom manually -- not a fan!
                # Do not move it to a separate style sheet.
                cyto.Cytoscape(
                    id="kin-map",
                    layout={"name": "circle", "fit": True},
                    style={
                        "height": "600px",
                        "width": "1000px",
                    },
                    zoom=0.1,
                    elements=get_network_elements(network),
                    pan={"x": 500, "y": 300},
                ),
            ],
        ),
        html.Hr(),
        html.Div(
            id="node-info-wrapper",
            children=[
                html.Table(id="node-info", children=construct_empty_table()),
                html.P(
                    """Note: input proteins are excluded from
                          ranking unless cross-validated."""
                ),
            ],
            style={"display": "None"},
        ),
        html.Div(id="results-container", children=[]),
        # Line below is a hidden utility switch for chaining callbacks
        # if user presents valid kinase list for diffusion, the switch
        # value is set to 'valid' and diffusion callback fires.
        # Do not move it to a separate style sheet.
        dcc.RadioItems(
            id="diffusion-switch", value="invalid", style={"display": "none"}
        ),
    ]


def get_theory_tab():
    """Returns content of theory tab."""
    content = dcc.Markdown(app_text.tabs["theory_tab"])
    return content


def get_example_tab():
    """Returns contents of the example tab."""
    content = dcc.Markdown(app_text.tabs["example_tab"])
    return content


def get_footer_elements():
    """Returns children of the footer div."""
    footer_elems = [
        html.P(className="footer-text", children="\u00A9 2020"),
        html.A(
            className="footer-text",
            href="https://www.ilyanovikov.io",
            children="Ilya Novikov",
        ),
        html.A(
            href="https://github.com/ily123/ggid",
            children=html.Img(
                className="github-logo", src="assets/GitHub-Mark-Light-32px.png"
            ),
        ),
        html.A(
            href="https://linkedin.com/in/ilyabnovikov",
            children=html.Img(className="github-logo", src="assets/linkedin-3-32.png"),
        ),
    ]
    return footer_elems


def construct_empty_table():
    """Creates empty scaffold for node info table."""
    empty_table = [
        html.Tr(
            children=[
                html.Th("Protein"),
                html.Th("Input Protein"),
                html.Th("Rank"),
                html.Th("Zscore"),
            ]
        ),
        html.Tr(
            children=[
                html.Td("Click on a node"),
                html.Td(),
                html.Td(),
                html.Td(),
            ]
        ),
    ]
    return empty_table


# ---
# main app layout:


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server
app.title = "GGid"

tab_main = dbc.Tab(label="Diffusion Tool", children=get_main_tab())
tab_theory = dbc.Tab(label="Theory", children=get_theory_tab())
tab_example = dbc.Tab(label="Example", children=get_example_tab())

app.layout = dbc.Container(
    id="main-container",
    children=[
        dbc.CardHeader(
            id="header",
            children=html.H1([dbc.Badge("GGID: Human Kinome", color="dark")]),
        ),
        html.Br(),
        dbc.Container(
            id="content",
            children=[
                html.P(
                    """Find clusters of connected kinases in the
                    GO network of the human kinome."""
                ),
                dbc.Tabs(children=[tab_main, tab_example, tab_theory]),
            ],
        ),
        html.Br(),
        html.Br(),
        dbc.CardFooter(id="footer", children=get_footer_elements()),
    ],
)


# ---
# app callbacks:


@app.callback(
    Output("node-info", "children"),
    Input("kin-map", "tapNodeData"),
    State("diffusion-switch", "value"),
    prevent_initial_call=True,
)
def display_tapped_node_data(data, diffusion_switch):
    """Displays node ranks in table below graph."""
    if diffusion_switch != "valid":
        return None  # unless diffusion experiment has been done, return nothing
    out = construct_empty_table()
    out[1] = html.Tr(
        children=[
            html.Td(data["id"]),
            html.Td("True" if data["input_label_flag"] == 1 else "False"),
            html.Td(data["rank"] if data["rank"] else "None"),
            html.Td("%1.2f" % data["zscore"] if data["zscore"] else "None"),
        ]
    )
    return out


@app.callback(
    [
        Output("status-line", "children"),
        Output("diffusion-switch", "value"),
    ],
    [Input("submit-button", "n_clicks")],
    [State("input-kinase-list", "value")],
    prevent_initial_call=True,
)
def validate_inputs(n_clicks, protein_list):
    """Validates protein list submitted by user."""
    message = []
    # diffusion switch, by default, is set to return no-update signal (ie: do nothing)
    diffusion_switch = dash.no_update
    proteins = re.findall("\w+", protein_list)
    proteins = [protein.upper() for protein in proteins]

    if len(proteins) == 0:
        message.append(
            dbc.Alert("You haven't entered anything in the textbox", color="danger")
        )
        return message, diffusion_switch
    else:
        valid_kinases = InputValidator().validate(proteins)
        valid_kinases = [vk for vk in valid_kinases if vk in network.proteins]
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
            return message, diffusion_switch
        if len(valid_kinases) > 0:
            message.append(
                dbc.Alert(
                    children=[
                        html.P(
                            """
                            Diffusing kinases in the input set: {kinases}.
                            See diffusion results (graph and z-score table) below.
                            """.format(
                                kinases=", ".join(valid_kinases)
                            )
                        ),
                        html.P(
                            """
                            Higher z-score indicates closer connectivity of the protein
                            to the input set. Only proteins with z-score > 2 are
                            show in the graph view.
                            """
                        ),
                    ],
                    color="success",
                )
            )
            diffusion_switch = "valid"
        if len(invalid_kinases) > 0:
            message.append(
                dbc.Alert(
                    "These items were not found in network: {kins}".format(
                        kins=", ".join(invalid_kinases)
                    ),
                    color="dark",
                ),
            )
    # wrap alerts into a collapsible
    message_wrap = html.Details(
        children=[
            html.Summary("Experiment summary (expand for details)"),
            html.Div(children=message),
        ]
    )
    return message_wrap, diffusion_switch


@app.callback(
    [
        Output("results-container", "children"),
        Output("kin-map", "elements"),
        Output("kin-map", "stylesheet"),  # 'stylesheet' is different from 'style'
        Output("network-layout-toggle", "options"),
        Output("node-info-wrapper", "style"),
    ],
    [Input("diffusion-switch", "value")],
    [
        State("input-kinase-list", "value"),
        State("loo-switch", "value"),
        State("zscore-cutoff", "value"),
    ],
    prevent_initial_call=True,
)
def diffuse(diffusion_switch, protein_list, loo_switch, zscore_cutoff):
    """Conducts diffusion experiment with input kinases."""
    proteins = re.findall("\w+", protein_list)
    # have to validate again -- I am not good at passing data b/w callbacks
    valid_kinases = InputValidator().validate(proteins)
    valid_kinases = [vk for vk in valid_kinases if vk in network.proteins]
    if "on" in loo_switch and len(valid_kinases) >= 2:
        # averaged post-diffusion results produced via LOO validation
        result_div, graph_nodes, node_styling = get_cross_validation_result(
            valid_kinases, zscore_cutoff
        )
    else:
        # plain diffusion, with no LOO validation
        zscore_table, graph_nodes, node_styling = get_diffusion_result(
            valid_kinases, zscore_cutoff
        )
        result_div = (
            html.Details(
                children=[
                    html.Summary("Full z-score and ranks table (expand for details)"),
                    convert_to_dash_table(zscore_table),
                ]
            ),
        )
    graph_options = get_graph_layout_options(isdisabled=False)
    node_info_style = {}  # this div is hidden prior to diffusion
    return result_div, graph_nodes, node_styling, graph_options, node_info_style


@app.callback(
    Output("kin-map", "layout"),
    Input("network-layout-toggle", "value"),
)
def change_graph_layout(layout_type):
    """Changes layout of the network graph."""
    layout = {"name": layout_type}
    return layout


if __name__ == "__main__":
    app.run_server(debug=True)
