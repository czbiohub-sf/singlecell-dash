from dash.dependencies import Output, Input
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import pandas as pd
import plotly.graph_objs as go

from .base import BaseBlock, CONFIG_DICT
from .color_by import ColorByGeneExpression, ColorByMetadata
from .dropdown_subset import SubsetGroup
from .umis_vs_genes import UMIsVsGenesGate

AXIS_TYPES = 'Linear', 'Log'


class GeneVsGene(BaseBlock):
    ID = 'gene_vs_gene_plot'

    def __init__(self, app, cell_metadata, counts, group_col,
                 gene_x='Actb', gene_y='Eef1a1'):
        self.cell_metadata = cell_metadata

        self.counts = counts
        self.group_col = group_col
        self.genes = sorted(counts.columns)

        self.gene_x = gene_x
        self.gene_y = gene_y

        self.metadata_grouped = self.cell_metadata.groupby(self.group_col)

        super().__init__(app)

    @property
    def layout(self):
        return html.Div([  # Gene v Gene plots
            html.H4('Gene vs gene comparison', className='row',
                    style={'padding-top': '20px'}),
            html.Div([
                html.Div([  # gene selector 1
                    dcc.Dropdown(
                        id='xaxis-column',
                        options=[{'label': i, 'value': i} for i in
                                 self.genes],
                        value=self.gene_x,
                    ),
                    dcc.RadioItems(
                        id='xaxis-type',
                        options=[{'label': i, 'value': i} for i in AXIS_TYPES],
                        value='Linear',
                        labelStyle={'display': 'inline-block'}
                    )
                ],
                    className='six columns'
                ),
                html.Div([  # gene selector 2
                    dcc.Dropdown(
                        id='yaxis-column',
                        options=[{'label': i, 'value': i} for i in
                                 self.genes],
                        value=self.gene_y
                    ),
                    dcc.RadioItems(
                        id='yaxis-type',
                        options=[{'label': i, 'value': i} for i in AXIS_TYPES],
                        value='Linear',
                        labelStyle={'display': 'inline-block'}
                    )
                ],
                    className='six columns'
                ),

            ],
                className='row'),

            # gene v gene plot
            dcc.Graph(id=self.ID,
                      config=CONFIG_DICT.copy()),
        ],
            className='six columns'
        )

    def callbacks(self, app):
        @app.callback(
            Output(self.ID, 'figure'),
            [Input('xaxis-column', 'value'),
             Input('yaxis-column', 'value'),
             Input('xaxis-type', 'value'),
             Input('yaxis-type', 'value'),
             Input(SubsetGroup.ID, 'value'),
             Input(UMIsVsGenesGate.ID, 'selectedData')
             ])
        def update_gene_vs_gene_scatter(xaxis_col, yaxis_col,
                                        xaxis_type, yaxis_type,
                                        group_name, selectedData):
            """Update the gene vs gene scatter plot"""

            group_barcodes = self.metadata_grouped.groups[group_name]

            group_counts = self.counts.loc[group_barcodes]
            alpha = pd.Series(1.0, index=group_counts.index)
            hovertext = group_counts[[xaxis_col, yaxis_col]].apply(
                lambda x: '{}: {}, {}'.format(x.name, x[0], x[1]), 1)

            if selectedData and selectedData['points']:
                barcodes = {d['customdata'] for d in selectedData['points']}
                alpha.loc[~alpha.index.isin(barcodes)] = 0.1
                hovertext[~hovertext.index.isin(barcodes)] = ''

            return {
                'data': [
                    go.Scatter(x=group_counts[xaxis_col],
                               y=group_counts[yaxis_col],
                               mode='markers',
                               marker={'opacity': alpha},
                               hovertext=hovertext.values)
                ],
                'layout': go.Layout(
                    xaxis={'title': xaxis_col + ' ' + group_name,
                           'type': xaxis_type.lower()},
                    yaxis={'title': yaxis_col,
                           'type': yaxis_type.lower(),
                           'scaleanchor': 'x'},
                    margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
                    hovermode='closest', dragmode='select')
            }
