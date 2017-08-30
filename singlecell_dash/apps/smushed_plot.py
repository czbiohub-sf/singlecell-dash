"""Visualize 2d embedding (e.g. TSNE) of cell space"""
from dash.dependencies import Output, Input
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import pandas as pd
import plotly.graph_objs as go

from .base import BaseBlock, CONFIG_DICT
from .color_by import ColorByGeneExpression, ColorByMetadata
from .dropdown_subset import SubsetBase
from .umis_vs_genes import UMIsVsGenesGate


class SmushedPlot(SubsetBase):
    ID = 'smushed'

    config_dict = CONFIG_DICT.copy()

    def __init__(self, app, cell_metadata, dropdown_col, smushed, counts,
                 top_genes=None):
        """Visualize two-dimensional embedding of gene expression data"""
        self.smushed = smushed
        self.counts = counts
        self.genes = sorted(counts.columns)

        self.top_genes = top_genes

        self.log_counts = np.log10(self.counts + 1.0)
        self.log_counts_max = self.log_counts.max().max()

        self.metadata_columns = cell_metadata.columns.difference(
            [dropdown_col])

        super().__init__(app, cell_metadata, dropdown_col)

    @property
    def layout(self):
        return dcc.Graph(id=self.ID,
                      config=self.config_dict.copy())

    def _top_genes_hovertext(self, barcode):
        """Return an HTML string of the top genes for each cell (barcode)"""
        if self.top_genes is not None:
            first_line = '{}<br>Top genes:<br>'.format(barcode)
            top_genes = ['{}. {}'.format(i + 1, gene)
                         for i, gene in enumerate(self.top_genes[barcode])]
            return first_line + '<br>'.join(top_genes)
        else:
            return barcode

    def _scatter(self, df, color=None, data_type=None, opacity=None,
                 showscale=False, colorbar_title=None, text=None,
                 customdata=None, name=None, **kwargs):
        cmin = None
        cmax = None

        if data_type is not None and 'expression'.startswith(data_type):
            cmin = 0.0
            cmax = self.log_counts_max

        kwargs['marker'] = dict(color=color, cmin=cmin, cmax=cmax,
                                colorscale='Viridis', opacity=opacity,
                                showscale=showscale,
                                colorbar=dict(thickness=20,
                                              title=colorbar_title,
                                              titleside='right'))

        scatter = go.Scatter(x=df[0], y=df[1], mode='markers',
                             hoverinfo='text', text=text,
                             customdata=customdata, name=name, **kwargs)
        return scatter

    @staticmethod
    def _values_hovertext(series):
        ids = series.index
        values = series.values
        strings = ['{}: {:.1f}'.format(i, v) for i, v in zip(ids, values)]
        return pd.Series(strings, index=ids)

    def callbacks(self, app):

        @app.callback(
            Output(self.ID, 'figure'),
            [Input(self.SUBSET_ID, 'value'),
             Input(UMIsVsGenesGate.ID, 'selectedData'),
             Input(ColorByGeneExpression.ID, 'value'),
             Input(ColorByMetadata.ID, 'value')])
        def update_cell_tsne(group_name, selectedDataQC=None,
                             selected_gene=None, selected_metadata=None):
            smushed = self.smushed[group_name]
            alpha = pd.Series(1.0, index=smushed.index)

            if selected_gene and not selected_metadata:
                group_barcodes = self._get_dropdown_barcodes(group_name)
                log_gene_data = self.log_counts.loc[group_barcodes,
                                                    selected_gene]
                hovertext = self._values_hovertext(log_gene_data)
            else:
                log_gene_data = None
                hovertext = smushed.index.to_series()
                hovertext = hovertext.map(self._top_genes_hovertext)

            if selectedDataQC and selectedDataQC['points']:
                barcodes = {d['customdata'] for d in selectedDataQC['points']}

                alpha.loc[~alpha.index.isin(barcodes)] = 0.1
                hovertext[~hovertext.index.isin(barcodes)] = ''

            if selected_metadata:
                groupby = self.cell_metadata[selected_metadata]
                try:
                    # Numeric data: use a color scale
                    groupby = pd.to_numeric(groupby)
                    hovertext = self._values_hovertext(groupby)

                    scatter = self._scatter(smushed, color=groupby,
                                            opacity=alpha, showscale=True,
                                            colorbar_title=selected_metadata,
                                            text=hovertext.values,
                                            customdata=smushed.index)
                    scatters = [scatter]
                except ValueError:
                    scatters = []
                    # Categorical data: plot every group
                    for name, df in smushed.groupby(groupby):
                        scatter = self._scatter(df, name=name,
                                                opacity=alpha.loc[df.index],
                                                text=hovertext.loc[df.index],
                                                customdata=df.index)
                        scatters.append(scatter)
            else:
                scatter = self._scatter(smushed, color=log_gene_data,
                                        showscale=bool(selected_gene),
                                        colorbar_title='log10 KPM',
                                        data_type='expression',
                                        text=hovertext.values,
                                        customdata=smushed.index,
                                        opacity=alpha)
                scatters = [scatter]

            return {
                "data": scatters,
                "layout": go.Layout(
                    xaxis={'title': 'TSNE 1', 'showticklabels': False},
                    yaxis={'title': "TSNE 2", 'showticklabels': False},
                    hovermode='closest', dragmode='lasso'),
            }
