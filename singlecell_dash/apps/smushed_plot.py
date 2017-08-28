"""Visualize 2d embedding (e.g. TSNE) of cell space"""
from dash.dependencies import Output, Input
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import pandas as pd
import plotly.graph_objs as go

from .base import BaseBlock, CONFIG_DICT
from .dropdown_subset import SubsetGroup
from .umis_vs_genes import UMIsVsGenesGate


class SmushedPlot(BaseBlock):
    ID = 'smushed'
    SELECTED_GENE_ID = 'selected_gene'

    config_dict = CONFIG_DICT.copy()

    def __init__(self, app, cell_metadata, group_col, smushed, counts,
                 top_genes):
        """Visualize two-dimensional embedding of gene expression data"""
        self.cell_metadata = cell_metadata
        self.group_col = group_col
        self.smushed = smushed
        self.counts = counts
        self.genes = sorted(counts.columns)
        self.top_genes = top_genes

        self.metadata_grouped = self.cell_metadata.groupby(self.group_col)

        super().__init__(app)

    @property
    def layout(self):
        return html.Div([
            html.H4('Cell TSNE', className='row',
                    style={'padding-top': '20px'}),

            html.Div([  # gene selector
                html.Label('Gene: ',
                           htmlFor=self.SELECTED_GENE_ID,
                           className='two columns offset-by-four'
                           ),
                dcc.Dropdown(
                    id=self.SELECTED_GENE_ID,
                    options=[{'label': i, 'value': i} for i in
                             self.genes],
                    value='',
                    className='six columns'
                ),
            ],
                className='row',
            ),

            # 2d embedding/TSNE plot
            dcc.Graph(id=self.ID,
                      config=self.config_dict.copy())
        ],
            className='six columns'
        )

    def _top_genes_hovertext(self, barcode):
        """Return an HTML string of the top genes for each cell (barcode)"""
        first_line = '{}<br>Top genes:<br>'.format(barcode)
        top_genes = ['{}. {}'.format(i + 1, gene)
                     for i, gene in enumerate(self.top_genes[barcode])]
        return first_line + '<br>'.join(top_genes)

    def callbacks(self, app):

        @app.callback(
            Output(self.ID, 'figure'),
            [Input(SubsetGroup.ID, 'value'),
             Input(UMIsVsGenesGate.ID, 'selectedData'),
             Input(self.SELECTED_GENE_ID, 'value')])
        def update_cell_tsne(group_name, selectedData=None,
                             selected_gene=None):
            smushed = self.smushed[group_name]
            alpha = pd.Series(1.0, index=smushed.index)

            if selected_gene:
                group_barcodes = self.metadata_grouped.groups[group_name]
                log_gene_data = np.log10(self.counts.loc[group_barcodes][
                                             selected_gene] + 1)
                hovertext = log_gene_data.to_frame().apply(
                    lambda x: '{}: {:.1f}'.format(x.name, x[0]), 1
                )
            else:
                log_gene_data = None
                hovertext = smushed.index.to_series()
                hovertext = hovertext.map(self._top_genes_hovertext)

            if selectedData and selectedData['points']:
                barcodes = {d['customdata'] for d in selectedData['points']}

                alpha.loc[~alpha.index.isin(barcodes)] = 0.1
                hovertext[~hovertext.index.isin(barcodes)] = ''

            return {
                "data": [go.Scatter(x=smushed[0],
                                    y=smushed[1],
                                    marker=dict(color=log_gene_data,
                                                cmin=0.0, cmax=6.0,
                                                colorscale='Viridis',
                                                opacity=alpha,
                                                showscale=bool(selected_gene),
                                                colorbar=dict(thickness=20,
                                                              title='log10 KPM',
                                                              titleside='right')
                                                ),
                                    mode='markers', hoverinfo='text',
                                    text=hovertext.values,
                                    customdata=smushed.index)],
                "layout": go.Layout(
                    xaxis={'title': 'TSNE 1', 'showticklabels': False},
                    yaxis={'title': "TSNE 2", 'showticklabels': False},
                    hovermode='closest', dragmode='select'),
            }
