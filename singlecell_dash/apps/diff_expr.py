"""Visualize 2d embedding (e.g. TSNE) of cell space"""
from collections import OrderedDict

from dash.dependencies import Output, Input
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from scipy import stats

from ..common import diff_exp
from .base import BaseBlock, CONFIG_DICT
from .dropdown_subset import SubsetBase
from .smushed_plot import SmushedPlot
from .umis_vs_genes import UMIsVsGenesGate


def sparse_diff_exp(matrix, group1, group2, index):
    """Computes differential expression between group 1 and group 2
    for each column in the dataframe counts.

    Parameters
    ----------
    matrix : SparseDataFrame
        Sparse gene expression data
    group1 : list of str
        Ids for the samples in group 1
    group2 : list of str
        IDs of the samples in group 2
    index : list of str
        Names of the genes in the data

    Returns a dataframe of Z-scores and p-values."""

    g1 = matrix[group1, :]
    g2 = matrix[group2, :]

    g1mu = np.asarray(g1.mean(0)).flatten()
    g2mu = np.asarray(g2.mean(0)).flatten()

    mean_diff = g1mu - g2mu
    # E[X^2] - (E[X])^2
    pooled_sd = np.sqrt(((g1.power(2)).mean(0) - g1mu ** 2) / len(group1)
                        + ((g2.power(2)).mean(0) - g2mu ** 2) / len(group2))
    pooled_sd = np.asarray(pooled_sd).flatten()

    z_scores = np.zeros_like(pooled_sd)
    nz = pooled_sd > 0
    z_scores[nz] = np.nan_to_num(mean_diff[nz] / pooled_sd[nz])

    # t-test
    p_vals = np.clip((1 - stats.norm.cdf(np.abs(z_scores))) * 2 * matrix.shape[1], 0, 1)

    df = pd.DataFrame(OrderedDict([('z', z_scores), ('p', p_vals),
                                   ('mean1', g1mu), ('mean2', g2mu)]),
                      index=index)

    return df


class DifferentialExpression(SubsetBase):
    ID = 'diff_expr'
    DIFFERENCE_TYPE = 'difference_type'

    XAXIS_TITLE = 'UMI Counts'
    YAXIS_TITLE = 'Gene'

    config_dict = CONFIG_DICT.copy()

    def __init__(self, app, cell_metadata, dropdown_col, counts):
        """Visualize two-dimensional embedding of gene expression data"""

        self.counts = counts
        self.genes = sorted(counts.columns)

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

            # counts = self.counts[group_barcodes, :]

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
                diff_stats = sparse_diff_exp(self.counts, selected_barcodes,
                                             unselected_barcodes,
                                             self.counts.columns)

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
                    [np.ravel(self.counts[selected_barcodes, g].todense())
                     for g in genes_to_show])
                y1 = np.concatenate(
                    [[g] * len(selected_barcodes) for g in genes_to_show])
                x2 = np.concatenate(
                    [np.ravel(self.counts[unselected_barcodes, g].todense())
                     for g in genes_to_show])
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
                # genes_to_show = counts.mean().nlargest(5).index
                #
                # x1 = np.concatenate(
                #     [counts[group_barcodes, g] for g in genes_to_show])
                # y1 = np.concatenate(
                #     [[g] * len(group_barcodes) for g in genes_to_show])
                # data = [go.Box(x=x1, y=y1,
                #                     name='Selected', orientation='h',
                #                     jitter=0.5)]
                return {
                    "data": [],
                    "layout": go.Layout(
                        title="Select cells to see differential expression",
                                        xaxis={'title': self.XAXIS_TITLE},
                                        yaxis={'title': self.YAXIS_TITLE},
                                        ),
                }

