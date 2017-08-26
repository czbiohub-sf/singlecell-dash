
from dash.dependencies import Output, Input
import dash_html_components as html
import dash_core_components as dcc
import pandas as pd
import plotly.graph_objs as go

from .base import BaseBlock, CONFIG_DICT


class SubsetGroup(BaseBlock):

    ID = 'group_name'

    def __init__(self, app, groups, name='Plate'):
        self.groups = groups
        self.name = name
        super().__init__(app)

    @property
    def layout(self):
        return html.Div([  # plate selector
            html.Label(f'{self.name}:',
                       htmlFor=self.ID,
                       className='one columns offset-by-one'),
            dcc.Dropdown(
                id=self.ID,
                options=[{'label': i, 'value': i} for i in self.groups],
                value=self.groups[0],
                className='four columns',
                clearable=False
            ),
        ],
            className='row'
        )