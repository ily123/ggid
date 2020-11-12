import pickle
import re

import dash
import dash_core_components as dcc
import dash_cytoscape as cyto
import dash_html_components as html
import dash_table
import dash_table.FormatTemplate as FormatTemplate
import pandas as pd
from dash.dependencies import Input, Output, State
from dash_table.Format import Format, Scheme, Sign, Symbol
from scipy import sparse

import diffusion
from kinapp_helper import InputValidator

pd.options.display.float_format = "${:,.2f}".format

cyto.load_extra_layouts()

image_url = "url(https://image.freepik.com/free-photo/elegant-black-handmade-technique-aquarelle_23-2148300751.jpg)"

# network as a coo matrix
network = pickle.load(open("network/kinase_matrix.pkl", "rb"))


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
        # elements.append({"data": {"source": label_i, "target": label_j}})
        print(index)
        # if index > 65:  # bugs start at 65
        #    break
    print(len(elements))
    return elements


app = dash.Dash(__name__)
app.title = "GGid"

app.layout = html.Div(
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
        html.Details(
            children=[
                html.Summary("Config Menu"),
                dcc.Checklist(
                    options={"label": "conduct LOO validation", "value": True},
                    value=True,
                ),
            ],
            id="config-menu",
        ),
        style={"display": "flex", "justifyContent": "center"},
    ),
    html.Div(
        id="kin-map-container",
        children=cyto.Cytoscape(
            id="kin-map",
            # layout={'name': 'circle'},
            layout={
                "name": "cose-bilkent",
                "numIter": 50,
                "padding": 30,
                "nodeDimensionsIncludeLabels": "true",
            },
            # layout={'name': 'cola', 'numIter': 10},
            # layout={'name': 'klay', 'numIter': 10},
            # layout={'name': 'dagre', 'numIter': 10},
            style={"width": "75%", "height": "500px", "display": "flex"},
            elements=make_nodes(network),
        ),
        style={"justifyContent": "center", "display": "flex"},
    ),
    html.Div(id="response", children="some text"),
    html.Div(id="output", style={"display": "none"}),
    # style={'background-image':image_url, 'background-size':'cover', 'height':'100%', 'position':'fixed', 'width':'100%', 'top':'0px', 'left':'0px'}
)

# from dash_table.Format import Format, Group, Scheme, Symbol
# frmt = Format(scheme=Scheme.fixed, precision=2, group=Group.yes,
#                groups=3,
#                group_delimiter='.',
#                decimal_delimiter=',')


@app.callback(
    [
        Output("output", "children"),
        Output("output", "style"),
        Output("response", "children"),
    ],
    [Input("submit-button", "n_clicks")],
    [State("kinase-list", "value")],
)
def print_output(n_clicks, value):
    if n_clicks is not None and value is not None:
        # validate input
        labels = re.findall("\w+", value)
        valid_labels = InputValidator().validate(labels)
        if len(valid_labels) == 0:
            msg = """Enter some valid kinase IDs to diffuse.
                    The app expects HUGO (hgnc) symbols. See https://www.genenames.org/"""
            return None, None, msg

        dif = diffusion.Diffusion(network, valid_labels)
        res = dif.diffuse()

        df = res.get_as_pandas_df_without_labels()
        # df['z_score']=df['z_score'].map("{:,.3f}".format)
        # df['final_label']=df['final_label'].map("{:,.4}".format)

        msg = f"Diffusing from nodes # {', '.join([str(x) for x in valid_labels])}"
        print(df.head())
        return (
            [
                dash_table.DataTable(
                    id="hi",
                    #                columns=[{"name": i, "id": i, 'type':'numeric','format':Format(precision=3)} for i in df.columns],
                    columns=[{"name": i, "id": i} for i in df.columns],
                    data=df.to_dict("records"),
                    export_columns="all",
                    export_format="csv",
                    page_size=20,
                )
                # html.Div(dcc.Graph(id='test_graph', figure={'data':[{'x':[1,2,3], 'y':[1,2,3], 'type':'bar'}]}))
                # html.Div(html.Img(src=app.get_asset_url('post_diffusion_pic.png'))),
                # html.Img(src=app.get_asset_url('post_diffusion_pic.png')),
            ],
            {"display": "flex", "justifyContent": "center"},
            msg,
        )

    return None, None, None


if __name__ == "__main__":
    app.run_server(debug=True)
