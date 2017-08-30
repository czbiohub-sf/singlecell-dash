"""Color the scatter plots by gene expression or metadata columns"""

from dash.dependencies import Output, Input
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import pandas as pd

from .dropdown_subset import SubsetBase


class ColorByGeneExpression(SubsetBase):

    GENE_ID = 'selected_gene'

    def __init__(self, app, counts, cell_metadata, dropdown_col):
        self.genes = sorted(counts.columns)

        super().__init__(app, cell_metadata, dropdown_col)

    @property
    def layout(self):
        return html.Div([  # gene selector
                html.Label('Color by gene expression: ',
                           htmlFor=self.GENE_ID,
                           className='four columns offset-by-two'
                           ),
                dcc.Dropdown(
                    id=self.GENE_ID,
                    options=[{'label': i, 'value': i} for i in
                             self.genes],
                    value='',
                    className='six columns'
                ),
            ],
                className='row',
            )

    def _single_gene_data(self, selected_gene, barcodes=None, log=True):
        """Get gene expression values for a single gene"""
        if barcodes:
            gene_data = pd.Series(self.counts[barcodes, selected_gene],
                                  index=barcodes)
        else:
            gene_data = pd.Series(self.counts[:, selected_gene],
                                  index=self.cell_metadata.index)
        if log:
            gene_data = np.log10(gene_data + 1)
        return gene_data


class ColorByMetadata(SubsetBase):

    METADATA_ID = 'selected_metadata'

    def __init__(self, app, cell_metadata, dropdown_col):

        super().__init__(app, cell_metadata, dropdown_col)

        self.metadata_columns = \
            self.cell_metadata.columns.difference([self.dropdown_col])

    @property
    def layout(self):
        return html.Div([  # gene selector
                html.Label('Color by cell metadata: ',
                           htmlFor=self.METADATA_ID,
                           className='four columns offset-by-two'
                           ),
                dcc.Dropdown(
                    id=self.METADATA_ID,
                    options=[{'label': i, 'value': i} for i in
                             self.metadata_columns],
                    value='',
                    className='six columns'
                ),
            ],
                className='row',
            )