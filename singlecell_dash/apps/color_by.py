"""Color the scatter plots by gene expression or metadata columns"""

from dash.dependencies import Output, Input
import dash_html_components as html
import dash_core_components as dcc

from .dropdown_subset import SubsetBase


class ColorByGeneExpression(SubsetBase):

    ID = 'selected_gene'

    def __init__(self, app, counts):
        self.genes = sorted(counts.columns)

        super().__init__(app)

    @property
    def layout(self):
        return html.Div([  # gene selector
                html.Label('Color by gene expression: ',
                           htmlFor=self.ID,
                           className='four columns offset-by-two'
                           ),
                dcc.Dropdown(
                    id=self.ID,
                    options=[{'label': i, 'value': i} for i in
                             self.genes],
                    value='',
                    className='six columns'
                ),
            ],
                className='row',
            )


class ColorByMetadata(SubsetBase):

    ID = 'selected_metadata'

    def __init__(self, app, cell_metadata):

        super().__init__(app)

        self.metadata_columns = \
            self.cell_metadata.columns.difference([self.dropdown_col])

    @property
    def layout(self):
        return html.Div([  # gene selector
                html.Label('Color by cell metadata: ',
                           htmlFor=self.ID,
                           className='four columns offset-by-two'
                           ),
                dcc.Dropdown(
                    id=self.ID,
                    options=[{'label': i, 'value': i} for i in
                             self.metadata_columns],
                    value='',
                    className='six columns'
                ),
            ],
                className='row',
            )