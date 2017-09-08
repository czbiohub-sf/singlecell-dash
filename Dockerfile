FROM jupyter/datascience-notebook

MAINTAINER Joshua Batson <joshua.batson@gmail.com>

USER root

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
 && apt-get install -y graphviz-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

USER jovyan

RUN conda config --system --add channels bioconda

RUN conda install --yes \
    'cycler' \
    'decorator' \
    'freetype' \
    'icu' \
    'ipython_genutils' \
    'jsonschema' \
    'jupyter_core' \
    'libpng' \
    'matplotlib' \
    'mkl' \
    'nbformat' \
    'networkx' \
    'numpy' \
    'openssl' \
    'pandas' \
    'pip' \
    'plotly' \
    'pyparsing' \
    'pyqt' \
    'python=3.6' \
    'python-dateutil' \
    'pytz' \
    'qt' \
    'readline' \
    'requests' \
    'seaborn' \
    'setuptools' \
    'sip' \
    'six' \
    'sqlite' \
    'scikit-learn' \
    'tk' \
    'traitlets' \
    'wheel' \
    'xz' \
    'zlib'

RUN /bin/bash -c 'source activate root'

RUN yes | pip install \
    'click' \
    'colorlover' \
    'dash>=0.17.7' \
    'dash_core_components>=0.11' \
    'dash-html-components>=0.6.2' \
    'dash-renderer>=0.7.3' \
    'fastcluster' \
    'flask' \
    'flask-compress' \
    'flask-seasurf' \
    'ipython-genutils' \
    'itsdangerous' \
    'hermione' \
    'jinja2' \
    'jupyter-core' \
    'markupsafe' \
    'werkzeug'

RUN conda install --yes graphviz

RUN pip install --upgrade git+https://github.com/pygraphviz/pygraphviz.git \
    --install-option="--library-path=/opt/conda/pkgs/graphviz-2.38.0-4/lib" \
    --install-option="--include-path=/opt/conda/pkgs/graphviz-2.38.0-4/include"

RUN /bin/bash  -c 'source deactivate'

RUN mkdir -p /home/$NB_USER/shared_scratch \
 && touch /home/$NB_USER/shared_scratch/.keep \
 && chown -R $NB_USER:users /home/$NB_USER/shared_scratch \
 && mkdir -p /home/$NB_USER/scratch \
 && touch /home/$NB_USER/scratch/.keep \
 && chown -R $NB_USER:users /home/$NB_USER/scratch

RUN mkdir -p /home/$NB_USER/singlecell_dash \
 && chown -R $NB_USER:users /home/$NB_USER/singlecell_dash

COPY tissue_analysis.py /home/$NB_USER/
COPY singlecell_dash /home/$NB_USER/singlecell_dash
COPY ["Tissue Exploration Notebook eXtravaganza.ipynb", "/home/$NB_USER/"]
