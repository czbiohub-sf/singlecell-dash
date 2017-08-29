"""Visualize 2d embedding (e.g. TSNE) of cell space"""
from dash.dependencies import Output, Input
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import pandas as pd
import plotly.graph_objs as go

from ..common import diff_exp
from .base import BaseBlock, CONFIG_DICT
from .dropdown_subset import SubsetBase
from .smushed_plot import SmushedPlot
from .umis_vs_genes import UMIsVsGenesGate


class DifferentialExpression(SubsetBase):
    ID = 'diff_expr'
    DIFFERENCE_TYPE = 'difference_type'

    XAXIS_TITLE = 'Log10 Expression'
    YAXIS_TITLE = 'Gene'

    config_dict = CONFIG_DICT.copy()

    def __init__(self, app, cell_metadata, dropdown_col, counts):
        """Visualize two-dimensional embedding of gene expression data"""

        self.counts = counts
        self.genes = sorted(counts.columns)

        self.log_counts = np.log10(self.counts + 1.0)

        super().__init__(app, cell_metadata, dropdown_col)

    @property
    def layout(self):
        return html.Div([  # differential gene plots
                    html.H4('Differential Expression', className='row',
                            style={'padding-top': '20px'}),
                    html.P('Use the lasso or box to select cells on the left '
                           'to see which genes are high- or low-expressed '
                           'compared to the remaining cells. Double-click the '
                           'plot to reset the selection.'),
                    html.Div([
                        dcc.RadioItems(
                            id=self.DIFFERENCE_TYPE,
                            options=[{'label': i, 'value': i} for i in
                                     ['High Expression', 'Low Expression']],
                            value='High Expression',
                            labelStyle={'display': 'inline-block'}
                        )
                    ],
                        className='row'
                    ),

                    # differential expression plot
                    dcc.Graph(id=self.ID,
                              config=self.config_dict.copy()),
                ],
                    className='six columns')

    def callbacks(self, app):

        @app.callback(
            Output(self.DIFFERENCE_TYPE, 'options'),
            [Input(SmushedPlot.ID, 'selectedData')])
        def update_diff_exp_ctrl(selectedDataTSNE=None):
            if selectedDataTSNE is not None:
                return [{'label': i, 'value': i} for i in
                        ['High Expression', 'Low Expression']]
            else:
                return [{'label': i, 'value': i, 'disabled': True} for i in
                        ['High Expression', 'Low Expression']]

        @app.callback(
            Output(self.ID, 'figure'),
            [Input(self.SUBSET_ID, 'value'),
             Input(UMIsVsGenesGate.ID, 'selectedData'),
             Input(SmushedPlot.ID, 'selectedData'),
             Input(self.DIFFERENCE_TYPE, 'value')])
        def update_diff_exp(group_name, selectedDataQC=None,
                            selectedDataTSNE=None, difference_type=None):
            group_barcodes = self._get_dropdown_barcodes(group_name)
            log_counts = self.log_counts.loc[group_barcodes]

            if selectedDataQC and selectedDataQC['points']:
                group_barcodes = [d['customdata'] for d in
                                  selectedDataQC['points']]

            if selectedDataTSNE and selectedDataTSNE['points']:
                selected_barcodes = {d['customdata'] for d in
                                     selectedDataTSNE['points']}
                selected_barcodes = selected_barcodes.intersection(
                    group_barcodes)
            else:
                selected_barcodes = set()

            # if the plate changed without updating the TSNE selection somehow,
            # this will be empty
            if selectedDataTSNE is not None:
                unselected_barcodes = [b for b in group_barcodes if
                                       b not in selected_barcodes]
                selected_barcodes = list(selected_barcodes)
                diff_stats = diff_exp(log_counts, selected_barcodes,
                                      unselected_barcodes)

                # Bonferroni cutoff
                bonferonni_cutoff = diff_stats['p'] < (0.05 / len(diff_stats))
                z_scores = diff_stats[bonferonni_cutoff]['z']

                if difference_type == "High Expression":
                    z_scores = z_scores[z_scores > 0]
                    genes_to_show = list(z_scores.nlargest(5).index)[::-1]
                else:
                    z_scores = z_scores[z_scores < 0]
                    genes_to_show = list(z_scores.nsmallest(5).index)[::-1]

                if not genes_to_show:
                    return {
                        "data": [],
                        "layout": go.Layout(
                            title="No differentially expressed genes!")
                    }

                x1 = np.concatenate(
                    [log_counts.loc[selected_barcodes, g] for g in
                     genes_to_show])
                y1 = np.concatenate(
                    [[g] * len(selected_barcodes) for g in genes_to_show])
                x2 = np.concatenate(
                    [log_counts.loc[unselected_barcodes, g] for g in
                     genes_to_show])
                y2 = np.concatenate(
                    [[g] * len(unselected_barcodes) for g in genes_to_show])

                return {
                    "data": [go.Box(x=x1, y=y1,
                                    name='Selected', orientation='h',
                                    jitter=0.5),
                             go.Box(x=x2, y=y2,
                                    name='Unselected', orientation='h',
                                    jitter=0.5),
                             ],

                    "layout": go.Layout(title="Differential Expression",
                                        xaxis={'title': self.XAXIS_TITLE},
                                        yaxis={'title': self.YAXIS_TITLE},
                                        boxmode='group'
                                        )
                }
            else:
                genes_to_show = log_counts.mean().nlargest(5).index

                x1 = np.concatenate(
                    [log_counts.loc[group_barcodes, g] for g in genes_to_show])
                y1 = np.concatenate(
                    [[g] * len(group_barcodes) for g in genes_to_show])

                return {
                    "data": [go.Box(x=x1, y=y1,
                                    name='Selected', orientation='h',
                                    jitter=0.5)],
                    "layout": go.Layout(title="Top Gene Expression",
                                        xaxis={'title': self.XAXIS_TITLE},
                                        yaxis={'title': self.YAXIS_TITLE},
                                        ),
                }

