# -*- coding: utf-8 -*-
import locale
import sys

import click
import dash
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import pandas as pd
import plotly.graph_objs as go

from common import diff_exp, Plates


# Use commas to separate thousands
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')


@click.command()
@click.option('--data-folder', default='data', help='Location of data files')
@click.option('--metadata', default='data/plate_metadata.csv',
              help='Full path of metadata file describing each plate')
@click.option('--verbose', help="Print the filenames as they are being read",
              is_flag=True)
@click.option('--port', help="Changes the port where the app is being "
                             "hosted from", default='8050')
@click.option('--host', help="Changes the host from which the app is being "
                             "hosted (i.e. global or local host). Default is "
                             "None (localhost). Change to '0.0.0.0' for "
                             "global host", default=None)
@click.option('--debug', help="Run the Dash server in debug mode",
              is_flag=True)
@click.version_option(version='v0.1.0')
def cli(data_folder, metadata, verbose, port, host, debug):
    """Run a dashboard showing sequencing QC of single-cell RNA-seq plates"""
    plates = Plates(data_folder, metadata, verbose)

    app = dash.Dash()
    app.css.append_css({"external_url":
                            "https://codepen.io/chriddyp/pen/bWLwgP.css"})


    # these aren't all working at the moment, may be a bug in Dash
    # note: always pass a copy of this dict, as it seems like it gets modified
    config_dict = {'modeBarButtonsToRemove': ['sendDataToCloud',
                                              'pan2d',
                                              'zoomIn2d',
                                              'zoomOut2d',
                                              'autoScale2d',
                                              'resetScale2d',
                                              'hoverCompareCartesian',
                                              'hoverClosestCartesian',
                                              'toggleSpikelines',
                                              ],
                   'displaylogo'           : False
                   }

    app.layout = html.Div([
        html.H1('MACA Dashboard', className='row'),

        html.Div([
            html.P('''This is a dashboard for preliminary evaluation of sequenced plates.
        The first half gives a view of all the plates at once so you can catch
        plates with significantly fewer reads or genes per cell, or with an odd
        proportion of ERCCs or mapped reads. You can also compare the average
        gene expression between plates in one TSNE. If
        a plate clusters with those of the wrong tissue, then it might be mislabelled.
        '''),

            html.P('''The second half lets you drill down on an individual plate. You can see
        its metadata, including tissue, mouse ID, age, and quality stats for the
        alignment. The plots allow you to gate out bad cells, analyze gene expression,
        and look for differential expression between groups of your choosing.
        ''')
        ],
                className='row'
        ),

        html.H2('QC of all plates', className='row',
                style={'padding-top': '20px'}),

        html.Div([  # all plates section
            html.Div([  # mean genes vs mean reads
                dcc.Graph(
                        id='plate_summary1',
                        figure={
                            "data"  : [go.Scatter(
                                    x=plates.plate_summaries[
                                        plates.MEAN_READS_PER_CELL],
                                    y=plates.plate_summaries[
                                        plates.MEDIAN_GENES_PER_CELL],
                                    mode='markers',
                                    hoverinfo='text',
                                    text=plates.plate_summaries.index)],
                            "layout": go.Layout(
                                    title="Plate reads and genes",
                                    xaxis={
                                        'title': 'Mean reads per cell'},
                                    yaxis={
                                        'title': "Median genes per cell"},
                                    hovermode='closest'
                            ),
                        },
                        config=config_dict.copy()
                )
            ],
                    className='six columns'
            ),

            html.Div([  # percent ERCC
                dcc.Graph(
                        id='plate_summary2',
                        figure={"data": [go.Scatter(
                                x=plates.plate_summaries[
                                    plates.PERCENT_ERCC],
                                y=plates.plate_summaries[
                                    plates.PERCENT_MAPPED_READS],
                                mode='markers', hoverinfo='text',
                                text=plates.plate_summaries.index)],
                            "layout"  : go.Layout(
                                    title="Plate ERCC and mapped reads",
                                    xaxis={
                                        'title': 'Percent ERCC'},
                                    yaxis={
                                        'title': "Percent mapped reads"},
                                    hovermode='closest'
                            ),
                        },
                        config=config_dict.copy()
                )
            ],
                    className='six columns'
            ),
        ],
                className='row'
        ),

        html.H2('TSNE of all plates', className='row',
                style={'padding-top': '20px'}),

        html.Div([  # TSNE of all plates
            html.Div([  # metadata feature selector
                html.Div([
                    html.Label('Color by:',
                               htmlFor='plate_metadata_feature',
                               className='two columns offset-by-one'),
                    dcc.Dropdown(
                            id='plate_metadata_feature',
                            options=[{'label': i, 'value': i} for i in
                                     plates.plate_metadata_features],
                            value='tissue_subtissue',
                            className='six columns'
                    ),
                ],
                        className='row',
                ),
                # all-plates TSNE
                dcc.Graph(
                        id='tsne_all_plates', figure={
                            "data"  : [go.Scatter(x=plates.bulk_smushed[0],
                                                  y=plates.bulk_smushed[1],
                                                  mode='markers',
                                                  hoverinfo='text',
                                                  text=plates.plate_summaries.index)],
                            "layout": go.Layout(
                                    title="Plate TSNE (in silico 'bulk' gene expression)",
                                    xaxis={'title'         : 'TSNE 1',
                                           'showticklabels': False},
                                    yaxis={'title'         : "TSNE 2",
                                           'showticklabels': False},
                                    hovermode='closest'
                            ),
                        },
                        config=config_dict.copy()
                )
            ],
                    className='eight columns'
            ),
        ],
                className='row'
        ),

        # plate statistics
        html.H2('Single plate statistics', className='row',
                style={'padding-top': '20px'}),

        html.Div([  # plate selector
            html.Label('Plate:',
                       htmlFor='plate_name',
                       className='one columns offset-by-one'),
            dcc.Dropdown(
                    id='plate_name',
                    options=[{'label': i, 'value': i} for i in
                             plates.plate_summaries.index],
                    value=plates.plate_summaries.index[0],
                    className='four columns'
            ),
        ],
                className='row'
        ),

        html.Div(id='plate_stats_table', className='row',
                 style={'column-width': '400px'}),

        # plate QC and other plots
        html.H2('Single plate visualizations', className='row',
                style={'padding-top': '20px'}),

        html.Div([
            html.P('''In the graphs below, each point represents a single cell.
        In the 'Reads vs Genes' graph, use the box selector to 'gate' cells
        based on number of distinct genes and number of reads.
        This will propagate to the other graphs.
        If you lasso a group of cells in the TSNE, then differential expression for
        those cells relative to the rest of the cells on the plate will appear to
        the right.
        ''')
        ],
                className='row'
        ),

        html.Div([
            # single plate reads vs genes
            html.Div([
                html.H4('Reads vs genes', className='row',
                        style={'padding-top': '20px',
                               'padding-bottom': '60px'}),

                dcc.Graph(id='single_plate_reads_vs_genes',
                          config=config_dict.copy()),
            ],
                    className='six columns'
            ),

            html.Div([  # Gene v Gene plots
                html.H4('Gene vs gene comparison', className='row',
                        style={'padding-top': '20px'}),

                html.Div([
                    html.Div([  # gene selector 1
                        dcc.Dropdown(
                                id='xaxis-column',
                                options=[{'label': i, 'value': i} for i in
                                         plates.genes.columns],
                                value='Actb'
                        ),
                        dcc.RadioItems(
                                id='xaxis-type',
                                options=[{'label': i, 'value': i} for i in
                                         ['Linear', 'Log']],
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
                                         plates.genes.columns],
                                value='Malat1'
                        ),
                        dcc.RadioItems(
                                id='yaxis-type',
                                options=[{'label': i, 'value': i} for i in
                                         ['Linear', 'Log']],
                                value='Linear',
                                labelStyle={'display': 'inline-block'}
                        )
                    ],
                            className='six columns'
                    ),
                ],
                        className='row'
                ),

                # gene v gene plot
                dcc.Graph(id='gene_gene_plot',
                          config=config_dict.copy()),
            ],
                    className='six columns'
            ),
        ],
                className='row'
        ),

        html.Div([
            # TSNE
            html.Div([
                html.H4('Cell TSNE', className='row',
                        style={'padding-top': '20px'}),

                html.Div([  # gene selector
                    html.Label('Gene: ',
                               htmlFor='selected_gene',
                               className='two columns offset-by-four'
                               ),
                    dcc.Dropdown(
                            id='selected_gene',
                            options=[{'label': i, 'value': i} for i in
                                     plates.genes.columns],
                            value='',
                            className='six columns'
                    ),
                ],
                        className='row',
                ),

                # TSNE plot
                dcc.Graph(id='single_plate_tsne',
                          config=config_dict.copy())
            ],
                    className='six columns'
            ),

            html.Div([  # differential gene plots
                html.H4('Differential Expression', className='row',
                        style = {'padding-top': '20px'}),

                html.Div([
                    dcc.RadioItems(
                            id='expression-type',
                            options=[{'label': i, 'value': i} for i in
                                     ['High Expression', 'Low Expression']],
                            value='High Expression',
                            labelStyle={'display': 'inline-block'}
                    )
                ],
                        className='row'
                ),

                # differential expression plot
                dcc.Graph(id='diff_exp',
                          config=config_dict.copy()),
            ],
                    className='six columns'
            ),

        ],
            className='row'
        )
    ],
            className='ten columns offset-by-one'
    )

    @app.callback(
            dash.dependencies.Output('single_plate_reads_vs_genes', 'figure'),
            [dash.dependencies.Input('plate_name', 'value')])
    def update_single_plate_reads_vs_genes(plate_name):
        """Callback when a plate is selected, update the reads v genes scatter plot
        """

        plate_barcodes = plates.cell_metadata.groupby('WELL_MAPPING').groups[
            plate_name]
        cell_metadata_subset = plates.cell_metadata.loc[plate_barcodes]

        return {
            "data"  : [go.Scatter(x=cell_metadata_subset['total_reads'],
                                  y=cell_metadata_subset['n_genes'],
                                  mode='markers', hoverinfo='text',
                                  text=cell_metadata_subset.index,
                                  customdata=cell_metadata_subset.index)],
            "layout": go.Layout(xaxis={'title': 'Reads per cell'},
                                yaxis={'title': "Genes per cell"},
                                margin={'b': 40, 't': 10, 'r': 0},
                                hovermode='closest', dragmode='select'),
        }


    @app.callback(
            dash.dependencies.Output('single_plate_tsne', 'figure'),
            [dash.dependencies.Input('plate_name', 'value'),
             dash.dependencies.Input('single_plate_reads_vs_genes',
                                     'selectedData'),
             dash.dependencies.Input('selected_gene', 'value')])
    def update_cell_tsne(plate_name, selectedData=None, selected_gene=None):
        smushed = plates.cell_smushed[plate_name]
        alpha = pd.Series(1.0, index=smushed.index)

        if selected_gene:
            plate_barcodes = plates.cell_metadata.groupby('WELL_MAPPING').groups[
                plate_name]
            log_gene_data = plates.counts_per_million.loc[plate_barcodes][
                selected_gene].map(lambda x: np.log10(x + 1.0))
            hovertext = log_gene_data.to_frame().apply(
                    lambda x: '{}: {:.1f}'.format(x.name, x[0]), 1
            )
        else:
            log_gene_data = None
            hovertext = smushed.index.to_series()
            hovertext = hovertext.map(
                lambda x: '{}<br>Top genes:<br>'.format(x)
                          + '<br>'.join(['{}. {}'.format(i+1, gene) for i, gene
                                         in enumerate(plates.top_genes[x])]))

        if selectedData and selectedData['points']:
            barcodes = [d['customdata'] for d in selectedData['points']]

            alpha.loc[~alpha.index.isin(barcodes)] = 0.1
            hovertext[~hovertext.index.isin(barcodes)] = ''

        return {
            "data"  : [go.Scatter(x=smushed[0],
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
            "layout": go.Layout(xaxis={'title': 'TSNE 1', 'showticklabels': False},
                                yaxis={'title': "TSNE 2", 'showticklabels': False},
                                hovermode='closest', dragmode='select'),
        }


    @app.callback(
            dash.dependencies.Output('expression-type', 'options'),
            [dash.dependencies.Input('single_plate_tsne', 'selectedData')])
    def update_diff_exp_ctrl(selectedDataTSNE=None):
        if selectedDataTSNE is not None:
            return [{'label': i, 'value': i} for i in
                    ['High Expression', 'Low Expression']]
        else:
            return [{'label': i, 'value': i, 'disabled': True} for i in
                    ['High Expression', 'Low Expression']]


    @app.callback(
            dash.dependencies.Output('diff_exp', 'figure'),
            [dash.dependencies.Input('plate_name', 'value'),
             dash.dependencies.Input('single_plate_reads_vs_genes',
                                     'selectedData'),
             dash.dependencies.Input('single_plate_tsne', 'selectedData'),
             dash.dependencies.Input('expression-type', 'value')])
    def update_diff_exp(plate_name, selectedDataQC=None,
                        selectedDataTSNE=None, expression_type=None):
        all_barcodes = plates.cell_metadata.groupby('WELL_MAPPING').groups[
            plate_name]

        if selectedDataQC and selectedDataQC['points']:
            all_barcodes = [d['customdata'] for d in selectedDataQC['points']]

        log_counts = np.log10(plates.counts_per_million.loc[all_barcodes] + 1.0)

        if selectedDataTSNE and selectedDataTSNE['points']:
            selected_barcodes = {d['customdata'] for d in
                                 selectedDataTSNE['points']}.intersection(
                all_barcodes)
        else:
            selected_barcodes = set()

        # if the plate changed without updating the TSNE selection somehow,
        # this will be empty
        if selected_barcodes:
            unselected_barcodes = [b for b in all_barcodes if
                                   b not in selected_barcodes]
            selected_barcodes = list(selected_barcodes)
            diff_stats = diff_exp(log_counts, selected_barcodes,
                                  unselected_barcodes)

            # Bonferroni cutoff
            z_scores = diff_stats[diff_stats['p'] < 0.05 / len(diff_stats)]['z']

            if expression_type == "High Expression":
                z_scores = z_scores[z_scores > 0]
                genes_to_show = list(z_scores.nlargest(5).index)[::-1]
            else:
                z_scores = z_scores[z_scores < 0]
                genes_to_show = list(z_scores.nsmallest(5).index)[::-1]

            if not genes_to_show:
                return {
                    "data"  : [],
                    "layout": go.Layout(
                            title="No differentially expressed genes")
                }

            x1 = np.concatenate(
                    [log_counts.loc[selected_barcodes, g] for g in genes_to_show])
            y1 = np.concatenate(
                    [[g] * len(selected_barcodes) for g in genes_to_show])
            x2 = np.concatenate([log_counts.loc[unselected_barcodes, g] for g in
                                 genes_to_show])
            y2 = np.concatenate(
                    [[g] * len(unselected_barcodes) for g in genes_to_show])

            return {
                "data"  : [go.Box(x=x1, y=y1,
                                  name='Selected', orientation='h',
                                  jitter=0.5),
                           go.Box(x=x2, y=y2,
                                  name='Unselected', orientation='h',
                                  jitter=0.5),
                           ],

                "layout": go.Layout(title="Differential Expression",
                                    xaxis={'title': 'Log expression'},
                                    yaxis={'title': "Gene"},
                                    boxmode='group'
                                    )
            }
        else:
            genes_to_show = log_counts.mean().nlargest(5).index

            x1 = np.concatenate(
                    [log_counts.loc[all_barcodes, g] for g in genes_to_show])
            y1 = np.concatenate([[g] * len(all_barcodes) for g in genes_to_show])

            return {
                "data"  : [go.Box(x=x1, y=y1,
                                  name='Selected', orientation='h',
                                  jitter=0.5)
                           ],

                "layout": go.Layout(title="Top Gene Expression",
                                    xaxis={'title': 'Log Expression'},
                                    yaxis={'title': "Gene"},
                                    ),
            }


    @app.callback(
            dash.dependencies.Output('gene_gene_plot', 'figure'),
            [dash.dependencies.Input('plate_name', 'value'),
             dash.dependencies.Input('xaxis-column', 'value'),
             dash.dependencies.Input('yaxis-column', 'value'),
             dash.dependencies.Input('xaxis-type', 'value'),
             dash.dependencies.Input('yaxis-type', 'value'),
             dash.dependencies.Input('single_plate_reads_vs_genes',
                                     'selectedData')])
    def update_gene_vs_gene_scatter(plate_name, xaxis_column_name,
                                    yaxis_column_name,
                                    xaxis_type, yaxis_type, selectedData=None):
        """Update the gene vs gene scatter plot"""
        plate_barcodes = plates.cell_metadata.groupby('WELL_MAPPING').groups[
            plate_name]
        genes_subset = plates.genes.loc[plate_barcodes]

        if selectedData and selectedData['points']:
            barcodes = [d['customdata'] for d in selectedData['points']]
            genes_subset = genes_subset.loc[barcodes, :]

        return {
            'data'  : [go.Scatter(
                    x=genes_subset[xaxis_column_name],
                    y=genes_subset[yaxis_column_name],
                    text=genes_subset.index,
                    mode='markers',
            )],
            'layout': go.Layout(
                    xaxis={
                        'title'    : xaxis_column_name,
                        'type'     : ('log', 'linear')[xaxis_type == 'Linear'],
                    },
                    yaxis={
                        'title'      : yaxis_column_name,
                        'type'       : ('log', 'linear')[yaxis_type == 'Linear'],
                        'scaleanchor': 'x'
                    },
                    margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
                    hovermode='closest')
        }


    @app.callback(
            dash.dependencies.Output('tsne_all_plates', 'figure'),
            [dash.dependencies.Input('plate_metadata_feature', 'value')]
    )
    def color_tsne_all_plates(plate_metadata_feature):
        """Update the all-plate TSNE plot when a new feature is selected"""
        feature = plates.plate_metadata[plate_metadata_feature]

        scatters = []

        for name, df in plates.bulk_smushed.groupby(feature):
            scatter = go.Scatter(x=df[0],
                                 y=df[1],
                                 mode='markers',
                                 name=name,
                                 hoverinfo='text',
                                 text=df.index)
            scatters.append(scatter)

        return {
            "data"  : scatters,
            "layout": go.Layout(
                    title="Plate TSNE (in silico 'bulk' gene expression)",
                    xaxis={'title'         : 'TSNE 1',
                           'showticklabels': False},
                    yaxis={'title'         : "TSNE 2",
                           'showticklabels': False},
                    hovermode='closest',
                    showlegend=True)
        }

    def maybe_format(item):
        """Pretty-format a string, integer, float, or percent
    
        Parameters
        ----------
        item : pandas.Series
            A single-item series containing a .name attribute and a value in the
            first (0th) index
        """
        value = item[0]
        if pd.isnull(value):
            return 'N/A'
        elif isinstance(value, str):
            return value
        elif 'percent' in item.name.lower():
            return '{:.2g}%'.format(value)
        elif (isinstance(value, float)  # this must go before ints!
              or np.issubdtype(value, np.number)):
            if value >= 1e3:
                return locale.format("%d", int(value), grouping=True)
            else:
                return locale.format("%.3g", value, grouping=True)
        elif (isinstance(value, int)
              or np.issubdtype(value, np.integer)):
            return locale.format("%d", value, grouping=True)
        else:
            raise TypeError


    @app.callback(
        dash.dependencies.Output('plate_stats_table', 'children'),
        [dash.dependencies.Input('plate_name', 'value')])
    def update_plate_stats_table(plate_name):
        summary = plates.plate_summaries.loc[plate_name]
        metadata = plates.plate_metadata.loc[plate_name]

        metadata = metadata.append(summary)

        row_list = [html.Tr([html.Th(name), html.Td(item)])
                    for name, item in
                    metadata.to_frame().apply(maybe_format, 1).items()]

        return html.Table(row_list, className='six columns')

    # this is where the magic happens
    app.run_server(host=host, debug=debug, port=port)


if __name__ == '__main__':
    cli()
