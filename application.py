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


def get_diffusion_result(labeled_kinases, zscore_cutoff):
    """Returns container with formatted diffusion results."""
    experiment = diffusion.Diffusion(network, labeled_kinases)
    result = experiment.diffuse()
    # make z-score table
    zscore_table = result.get_as_pandas_df_without_labels()
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


def create_cytoscape_div(sim_matrix, labels, diffusion_result, zscore_cutoff):
    """Constructs updated cytoscape div."""
    top_hits = diffusion_result[diffusion_result.zscore >= zscore_cutoff]
    top_hits_nodes = list(top_hits.protein.values)
    elements = make_nodes2(sim_matrix, labels, top_hits_nodes)
    new_elements = elements
    new_style = get_cyto_stylesheet(top_hits)
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
        },
    ]
    colors = get_colors(diffusion_result)
    color_selectors = make_selectors_colors(colors)
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
    """Generate color gradient for displayed nodes."""
    xcolor_gradient = color_gradient.ColorGradientGenerator()
    xcolor_gradient.create_color_map2(base_color=[255, 51, 51])
    print(diffusion_result["rank"])
    colors = xcolor_gradient.map_colors(
        1 - diffusion_result["rank"] / len(diffusion_result)
    )
    color_dict = dict(zip(diffusion_result.protein, colors))
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
            children=[
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
                dbc.Input(id="zscore-cutoff", value=2, type="number"),
            ],
            style={"display": "flex", "justifyContent": "center"},
        ),
        html.Div(
            id="kin-map-container",
            children=[
                dbc.FormGroup(
                    [
                        dbc.Label("Adjust graph layout"),
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
            style={
                "display": "grid",
                "justify-content": "center",
                "align-content": "center",
            },
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

    content = [
        dcc.Markdown(
            """
    ### The problem
    ---
    Proteins don't work in isolation. They form networks and work in tandem.

    This means that once you've identified
    a protein associated with some process (ex: a disease), you can query that
    protein's network to find additional proteins that are involved.
    Once you understand the network and the role of all relevant protein players,
    you can design effective theurapeutics.

    Construction of such networks, and learning from them, is one of the
    subfields in computational biology. The industry standard for collecting and
    visualizing these data is the [STRING](https://string-db.org/) database. There
    you can find all sorts of protein networks build on top of various data (physical
    interactions, co-expression of proteins, co-expression of protein mRNA, etc).

    ### GO term similarity network
    ---

    In this dashboard, I created a novel way of presenting another type of
    information: GO terms. GO terms come from the [Gene Ontology](http://geneontology.org/)
    database. They are tags assigned to proteins in order to describe proteins job within the cell (very similar to how
    IMDB assigns tags to movies). Proteins that do similar jobs have similar GO terms. We can build a network from that!

    To build a network, we compare proteins to each other, in an all-vs-all manner, to
    produce a similarity matrix that show how similar any pair of proteins are. We then
    reduce the similarity matrix to a graph (network) by dropping out connections
    that are weak (ie the two proteins are not similar). What we end up with is a network.

    The circle pic you see loaded in the first tab is the GO term network for 350 human
    kinases, build exactly like I described. (For those who want the nitty gritty details:
    similarity between two proteins was measured using Resnik similarity, the network
    was constructed by dropping out all edges outside of the top 10. The proteins, beforehand,
    were sanitized by removing those with fewer than 10 annotations.)


    ### Diffusion
    ---

    To extract information from this network, we can simply examine it by hand. Given a protein, what
    are its closest neighbors? That works, but gets out of hand really quickly when you
    want to look at more than just 1 protein. (Let's say you have 10 proteins involved
    in some disease, and you want to see all of their neighbors.)

    To help us, there are a bunch of algorithms that
    simplify the process of finding neighors within networks. The one I am using
    is called Information Diffusion. It was originally developed in the 1990s to
    model diffusion of heat through metal, but has since been deployed to protein
    networks as well.

    The "heat" in our case is information (a binary label: 1 if a protein does X, and 0 otherwise).
    We diffuse the label over the network of proteins. Post-diffusion, the proteins
    that have the highest label content are the ones most closely connected to the
    source.

    ### Kinases
    ---

    Ok, that's all great. But this app only includes a network of kinases.
    why did you build a kinase network specifically?

    I chose kinases for this demo app for technical reasons:

    * I wasn't sure the browser can render a full network of 20,000 human proteins
    * Meanwhile, the number of kinases is realtively small (just 500),
    so they are easy to model
    * Plus, kinases are very important in human disease (they drive signaling
    by phosphorylating other proteins), and as a result well-studied and annotated

    You can apply the framework described here to the entire human proteome (or other actors
    like RNAs).


    ### Does this tool actually work for predicting novel associations?
    ---

    This tool is good at presenting networks of existing information (ie cross-validation
    is pretty decent >= 0.80). That being said, prediction of future interactions
    from present data is pretty iffy (retrospective auc ~0.6-0.7). I am not showing any data here,
    just recalling the experimental results from the waning days of my training.
    I think some of it is in my thesis.

    As to why this model doesn't generalize well (i.e. it overtrains)  - I think this is
    because of the binary nature of GO term annotations. If two proteins are found
    to be associated, but their GO terms are completely different before the discovery
    is made, the model isn't going to be aware the relationship. The model predicts two proteins
    that already share GO terms as functionally related. If there is no similarit to
    cease on, the model simily won't create an edge in th network.

    ALL THIS NEEDS EDITING // FIRST DRAFT


    **If the tool can't predict the future, wth did you make this dashboard??**

    I wanted to write a network app. I have multiple ideas for using the
    basic approach in other domains. And now I have an app for that. Just
    have to change the underlying network and edit some text, and I am good to go!
    """
        )
    ]

    return content


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
tab_theory = dbc.Tab(label="Theory", children=[])
tab_example = dbc.Tab(label="Example", children=[])

tab_main.children = get_main_tab()
tab_theory.children = get_theory_tab()
tab_example.children = get_example_tab()

app.layout = dbc.Container(
    id="main-container",
    children=[
        dbc.CardHeader(
            id="header",
            children=html.H1([dbc.Badge("GGID: Human Kinome", color="dark")]),
            style={"width": "100%", "background-color": "#343a40"},
        ),
        html.Br(),
        dbc.Container(
            id="content",
            children=[
                html.P(
                    """Find clusters of connected kinases by diffusing information over
                    GO term network of the human kinome."""
                ),
                dbc.Tabs(children=[tab_main, tab_example, tab_theory]),
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
    [
        State("kinase-list", "value"),
        State("loo-switch", "value"),
        State("zscore-cutoff", "value"),
    ],
    prevent_initial_call=True,
)
def diffuse(diffusion_switch, protein_list, loo_switch, zscore_cutoff):
    """Conducts diffusion experiment with provided kinases."""
    proteins = re.findall("\w+", protein_list)
    valid_kinases = InputValidator().validate(proteins)
    if "on" in loo_switch:
        # averaged post-diffusion results produced via LOO validation
        result_div, graph_nodes, node_styling = get_cross_validation_result(
            valid_kinases, zscore_cutoff
        )
    else:
        # plain diffusion, with no LOO validation and averaging
        zscore_table, graph_nodes, node_styling = get_diffusion_result(
            valid_kinases, zscore_cutoff
        )
        result_div = [convert_to_dash_table(zscore_table)]

    # add z-score explanation
    result_div.insert(
        0,
        dbc.Alert(
            """
            Showing proteins with z-score >= 2 in the graph view.
            See full results in table below.
            Higher z-score indicates closer connectivity
            of the protein to the input set.
            """
        ),
    )
    result_div.insert(1, html.H3(dbc.Badge("Diffusion Result:", color="secondary")))
    return result_div, graph_nodes, node_styling


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
