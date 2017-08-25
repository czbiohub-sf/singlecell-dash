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


from .gene_vs_gene import GeneVsGeneBlock

# Use commas to separate thousands
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')


def run_singlecell_dash(cell_metadata, counts, group_col):
    app = dash.Dash()

    # Necessary with a modular layout since we add callbacks without adding a
    # layout
    app.config.supress_callback_exceptions = True

    # creating a new MyBlock will register all callbacks
    block = GeneVsGeneBlock(app, cell_metadata, counts, group_col)

    import pdb; pdb.set_trace()
    # now insert this component into the app's layout
    app.layout = html.Div([html.H1('My App'), block.layout],
                          className='ten columns offset-by-one')
    return app


