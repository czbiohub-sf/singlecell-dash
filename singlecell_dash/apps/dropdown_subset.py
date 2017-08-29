"""Dropdown selector for a particular group"""
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


class SubsetBase(BaseBlock):
    """Class to inherit from for using the subsetted group name"""

    SUBSET_ID = SubsetGroup.ID

    def __init__(self, app, cell_metadata, dropdown_col):
        self.cell_metadata = cell_metadata
        self.dropdown_col = dropdown_col

        self.metadata_grouped = self.cell_metadata.groupby(self.dropdown_col)

        super().__init__(app)

    def _get_dropdown_barcodes(self, group_name):
        """Given the name of a group, get the cells/barcodes of that group"""
        # TODO: add logic to deal with the special case of "All" samples
        return self.metadata_grouped.groups[group_name]