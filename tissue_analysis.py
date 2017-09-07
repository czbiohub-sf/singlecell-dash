
import argparse
import os
from collections import Counter, OrderedDict

import fastcluster
import matplotlib.colors
import matplotlib.figure

import networkx
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy
import scipy.stats as stats
import seaborn as sns

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Circle

from networkx.drawing.nx_agraph import graphviz_layout
from sklearn.neighbors import NearestNeighbors

import singlecell_dash.common as common


TISSUES = {'Tongue', 'Liver', 'Bladder', 'Kidney', 'Spleen', 'Marrow',
           'Lung', 'Muscle', 'Heart', 'Thymus', 'Mammary'}

AGE_3M_SAMPLES = {'10X_P4_0', '10X_P4_1', '10X_P4_2', '10X_P4_3',
                  '10X_P4_4', '10X_P4_5', '10X_P4_6', '10X_P4_7',
                  '10X_P6_0', '10X_P6_1', '10X_P6_2', '10X_P6_3',
                  '10X_P6_4', '10X_P6_5', '10X_P6_6', '10X_P6_7',
                  '10X_P7_0', '10X_P7_1', '10X_P7_2', '10X_P7_3',
                  '10X_P7_4', '10X_P7_5', '10X_P7_6', '10X_P7_7',
                  '10X_P7_8', '10X_P7_9', '10X_P7_10', '10X_P7_11',
                  '10X_P7_12', '10X_P7_13', '10X_P7_14', '10X_P7_15'}


def diff_exp_clusters(Z, expression_df, clusters):
    cluster_sizes = dict(Counter(clusters).most_common())
    n_clusters = len(cluster_sizes)

    cluster_sum_umi = expression_df.groupby(clusters).sum()
    cluster_ssq_umi = expression_df.groupby(clusters).apply(
            lambda df: (df ** 2).sum())

    root, rd = scipy.cluster.hierarchy.to_tree(Z, rd=True)

    def de(lbl_1, lbl_2, group1, group2):
        print(f'Comparing {group1} to {group2}')

        group1_n_cells = sum(cluster_sizes[c] for c in group1)
        group2_n_cells = sum(cluster_sizes[c] for c in group2)

        group1_mean = cluster_sum_umi[group1, :].sum(axis=0) / group1_n_cells
        group2_mean = cluster_sum_umi[group2, :].sum(axis=0) / group2_n_cells

        mean_diff = group1_mean - group2_mean

        group1_var = (cluster_ssq_umi[group1, :].sum(axis=0)
                      / group1_n_cells - group1_mean ** 2)
        group2_var = (cluster_ssq_umi[group2, :].sum(axis=0)
                      / group2_n_cells - group2_mean ** 2)

        pooled_sd = np.sqrt(group1_var / group1_n_cells
                            + group2_var / group2_n_cells)

        z_scores = np.zeros_like(pooled_sd)
        nz = pooled_sd > 0
        z_scores[nz] = np.nan_to_num(mean_diff[nz] / pooled_sd[nz])

        # t-test
        p_vals = np.clip((1 - stats.norm.cdf(np.abs(z_scores)))
                         * 2 * z_scores.shape[0], 0, 1)

        df = pd.DataFrame(OrderedDict([('z', z_scores),
                                       ('p', p_vals),
                                       ('group1', group1_mean),
                                       ('group2', group2_mean)]),
                          index=cluster_expression_df.index)

        df = df[df['p'] < 0.001]
        df['diff'] = df['group1'] - df['group2']

        df.sort_values('diff', ascending=False, inplace=True)

        name = f'differential_gene_expression_{lbl_1}_v_{lbl_2}'

        df.to_csv(file_format.format(name, 'csv'))


    for i in range(0, 2 * n_clusters - 1):
        if i >= n_clusters:
            left_child = rd[i].get_left()
            left_clusters = (left_child.pre_order(lambda x: x.id))

            right_child = rd[i].get_right()
            right_clusters = (right_child.pre_order(lambda x: x.id))

            # don't calculate if it's redundant with a 1-vs-all comp
            if i == 2 * n_clusters - 2 and (len(left_clusters) == 1
                                            or len(right_clusters) == 1):
                continue

            de(left_child.id, right_child.id, left_clusters, right_clusters)

        if i < 2 * n_clusters - 2:
            below = rd[i].pre_order(lambda x: x.id)
            above = [j for j in range(len(cluster_sizes)) if j not in below]

            # don't calculate redundant comparison
            if len(above) == 1:
                continue

            de(i, 'all', below, above)

    group_list = [(i, rd[i].pre_order(lambda x: x.id))
                  for i in range(0, 2 * n_clusters - 1)]
    group_list[-1] = ('total', group_list[-1][1])

    return group_list



def batch_plots(tenx, coords):
    """Generate heatmaps to diagnose batch effects across the clustering"""

    cell_metadata = tenx.cell_metadata.join(coords)

    fig,ax = plt.subplots(1, 2, figsize=(18,6))

    sns.heatmap(np.log10(pd.crosstab(cell_metadata['cluster'],
                                     cell_metadata['Mouse'])+1),
                annot=True, fmt="d", ax=ax[0])
    sns.heatmap(np.log10(pd.crosstab(cell_metadata['cluster'],
                                     cell_metadata['Sex'])+1),
                annot=True, fmt="d", ax=ax[1])

    plt.show()
