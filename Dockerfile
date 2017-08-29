FROM jupyter/minimal-notebook

ADD environment.txt /home/jovyan/environment.txt

RUN conda install -y --file /home/jovyan/environment.txt

RUN conda install -c anaconda graphviz
