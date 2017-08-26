# -*- coding: utf-8 -*-
import locale
import sys

import click
import dash
from dash.dependencies import Output, Input
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import pandas as pd
import plotly.graph_objs as go


from .gene_vs_gene import GeneVsGene
from .subset import SubsetGroup
from .scatter import Scatter

# Use commas to separate thousands
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')


def run_singlecell_dash(cell_metadata, counts, group_col, javascript=None):
# def run_singlecell_dash(javascript=None):
    app = dash.Dash()

    if javascript is not None:
        app.scripts.append_script({"external_url": javascript})

    # Necessary with a modular layout since we add callbacks without adding a
    # layout
    app.config.supress_callback_exceptions = True

    # Add nice looking css
    app.css.append_css({"external_url":
                            "https://codepen.io/chriddyp/pen/bWLwgP.css"})

    # creating a new MyBlock will register all callbacks
    # gene_vs_gene = GeneVsGene(app, cell_metadata, counts, group_col)
    subset = SubsetGroup(app, cell_metadata[group_col].unique())
    scatter = Scatter(app, cell_metadata, counts, group_col)

    # import pdb; pdb.set_trace()
    # now insert this component into the app's layout
    app.layout = html.Div([html.H1('Single Cell Dashboard App'),
                           subset.layout,
                           scatter.layout
                           ],
                          className='ten columns offset-by-one')
    return app


