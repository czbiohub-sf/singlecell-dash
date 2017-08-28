from dash.dependencies import Output, Input
import dash_html_components as html
import dash_core_components as dcc
import pandas as pd
import plotly.graph_objs as go

from .base import BaseBlock, CONFIG_DICT
from .dropdown_subset import SubsetGroup


class UMIsVsGenesGate(BaseBlock):

    ID = 'umis_vs_genes'
    
    def __init__(self, app, cell_metadata, group_col,
                 n_molecules_col='total_reads',
                 n_genes_col='n_genes'):
        """Show number of molecules vs number of genes and 'gate' the cells
        
        Parameters
        ----------
        app
        cell_metadata
        n_molecules_col : str
            Column name in "cell_metadata" indicating the total number of 
            unique molecules counted per cell
        n_genes_col : str
            Column name in "cell_metadata" indicating the total number of genes 
            detected per cell
        """
        self.cell_metadata = cell_metadata
        self.group_col = group_col
        self.n_molecules_col = n_molecules_col
        self.n_genes_col = n_genes_col
        self.metadata_grouped = self.cell_metadata.groupby(self.group_col)

        super().__init__(app)

    @property
    def layout(self):
        return html.Div([
            html.H4('Reads vs genes', className='row',
                    style={'padding-top': '20px',
                           'padding-bottom': '60px'}),
            dcc.Graph(id=self.ID, config=CONFIG_DICT.copy()),
        ],
            className='six columns')

    def callbacks(self, app):

        @app.callback(
            Output(self.ID, 'figure'),
            [Input(SubsetGroup.ID, 'value')])
        def update_reads_vs_genes(group_name):
            """When a group is selected, update the reads v genes scatter"""

            group_barcodes = self.metadata_grouped.groups[group_name]
            group_metadata_subset = self.cell_metadata.loc[group_barcodes]

            x = group_metadata_subset[self.n_molecules_col]
            y = group_metadata_subset[self.n_genes_col]

            return {
                "data": [go.Scatter(x=x,
                                    y=y,
                                    mode='markers', hoverinfo='text',
                                    text=group_metadata_subset.index,
                                    customdata=group_metadata_subset.index)],
                "layout": go.Layout(xaxis={'title': 'Reads per cell'},
                                    yaxis={'title': "Genes per cell"},
                                    margin={'b': 40, 't': 10, 'r': 0},
                                    hovermode='closest', dragmode='select'),
            }