# GGID

Protein annotation via Information Diffusion over GO-derived protein Graphs.
The deployed app lives [here](https://ggid.herokuapp.com/).

![Python](https://img.shields.io/badge/python-3.6%7C3.7-blue.svg)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Table of Content
1. [Graphical Abstract](#graphical-abstract)
2. [Overview](#overview)
3. [Installation](#installation)
4. [Deployment](#deployment)
5. [License](#license)

## Graphical Abstract 
<img src="assets/graphical_abstract.png" alt="drawing" width="700"/>

## Overview
A detailed description of the method is available as part of the deployed application (THEORY tab).
However, here is a brief overview:

* We can annotate proteins by propagating labels through protein networks

* In this repo, I created a new type of kinase network, one based purely on GO term annotations

* To propagate labels through this network, I also created a web-app interface
    using the [Dash](https://dash.plotly.com/) framework combined with the [Cytoscape](https://github.com/plotly/dash-cytoscape) library.

## License
