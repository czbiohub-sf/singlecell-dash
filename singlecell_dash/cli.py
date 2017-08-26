import click

from singlecell_dash.common import TenX_Runs
from singlecell_dash.app import run_singlecell_dash


@click.command()
@click.option('--data-folder', default='data', help='Location of data files')
@click.option('--metadata', default='data/plate_metadata.csv',
              help='Full path of metadata file describing each plate')
@click.option('--genes-to-drop', default='Rn45s',
              help='Gene symbols to remove from the counts matrix for '
                   'calculating TSNE. If more than one, must be '
                   'comma-separated ')
@click.option('--verbose', help="Print the filenames as they are being read",
              is_flag=True)
@click.option('--port', help="Changes the port where the app is being "
                             "hosted from", default=8050, type=int)
@click.option('--host', help="Changes the host from which the app is being "
                             "hosted (i.e. global or local host). Default is "
                             "None (localhost). Change to '0.0.0.0' for "
                             "global host", default=None)
@click.option('--javascript',
              help="Location of an arbitrary javacsript file you want to add"
                   " to the Dash app", default=None)
@click.option('--debug', help="Run the Dash server in debug mode",
              is_flag=True)
@click.version_option(version='v0.1.0')
def cli(data_folder, metadata, genes_to_drop, verbose, port, host, javascript,
        debug):
    """Run a dashboard showing sequencing QC of single-cell RNA-seq plates"""
    # plates = Plates(data_folder, metadata, genes_to_drop=genes_to_drop,
    #                 verbose=verbose)

    tenx_runs = TenX_Runs(data_folder, genes_to_drop=genes_to_drop,
                          verbose=verbose)

    app = run_singlecell_dash(tenx_runs.cell_metadata,
                              tenx_runs.counts_per_million,
                              tenx_runs.SAMPLE_MAPPING)


    # app = run_singlecell_dash()
    # this is where the magic happens
    app.run_server(host=host, debug=debug, port=port)


if __name__ == '__main__':
    cli()
