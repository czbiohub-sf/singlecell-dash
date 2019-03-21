# ARCHIVED

This repository was a work-in-progress during the best of times and it is no longer under develop. It may still contain some kernel of useful code, but it is not guaranteed to work at all.

# singlecell-dash

Dashboard for visualizing per-plate sequencing quality metrics and basic
analysis for a single-cell RNA-seq project.


## Installation

To install:

```
git clone git@github.com:czbiohub/singlecell-dash
cd singlecell-dash
conda env create --file environment.yml
```

This will create an environment called `singlecell-dash`, which you can then activate
using `source activate singlecell-dash`. This has plotly, dash, etc components installed.

## Usage

To run the dash app, do: `python app.py`, which should show the following output:

```
$ python app.py
 * Running on http://127.0.0.1:8050/ (Press CTRL+C to quit)
 * Restarting with stat
 * Debugger is active!
 * Debugger PIN: 118-904-628
```


NB: If you CTRL+Z or CTRL+X instead of CTRL+C, the 8050 port will still be in
use :( -- so you'll have to close your terminal window/tab and start anew.


To get help, check out the options with `--help`:

```
$ python app.py --help
Usage: app.py [OPTIONS]

  Run a dashboard showing sequencing QC of single-cell RNA-seq plates

Options:
  --data-folder TEXT  Location of data files
  --metadata TEXT     Full path of metadata file describing each plate
  --verbose           Print the filenames as they are being read
  --port INTEGER      Changes the port where the app is being hosted from
  --host TEXT         Changes the host from which the app is being hosted
                      (i.e. global or local host). Default is None
                      (localhost). Change to '0.0.0.0' for global host
  --javascript TEXT   Location of an arbitrary javacsript file you want to add
                      to the Dash app
  --debug             Run the Dash server in debug mode
  --version           Show the version and exit.
  --help              Show this message and exit.
```
