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