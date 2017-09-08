
import os
from collections import Counter, OrderedDict

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy
import scipy.stats as stats
import seaborn as sns

sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

import matplotlib.pyplot as plt

import singlecell_dash.common as common
import singlecell_dash.network_util as nutil

import hermione as hm


__all__ = ['TISSUES', 'load_tissue', 'cluster_tissue',
           'diff_exp_clusters', 'batch_plots',
           'subset_exp', 'plot_gene_joyplot']

TISSUES = {'Tongue', 'Liver', 'Bladder', 'Kidney', 'Spleen', 'Marrow',
           'Lung', 'Muscle', 'Heart', 'Thymus', 'Mammary'}


def load_tissue(data_folder, tissue):
    """Load all the 3-month data for the chosen tissue"""

    # load in data with appropriate subsets
    tenx = common.TenX_Runs(data_folder,
                            filters={'Age': [3], 'Tissue': [tissue]},
                            verbose=True)

    exp_df = pd.DataFrame(tenx.genes.matrix.todense(),
                          index=tenx.genes.rows,
                          columns=tenx.genes.columns)

    knn_cache_file = os.path.join(data_folder, 'knncache',
                                  f'{tissue}-k-500.npy')

    knn_cache = nutil.KNNCache(knn_cache_file)

    return tenx,exp_df,knn_cache


def cluster_tissue(exp_df:pd.DataFrame,
                   knn_cache:nutil.KNNCache, k:int):
    knn_graph = knn_cache.get_knn_graph(k=25)  # change k here

    coords = nutil.network_layout(exp_df.index, knn_graph)
    coords['cluster'], Z = nutil.label_propagation(exp_df, knn_graph)

    return coords, Z


def diff_exp_clusters(Z, expression_df, clusters, verbose=False):
    cluster_sizes = dict(Counter(clusters).most_common())
    n_clusters = len(cluster_sizes)

    cluster_sum_umi = expression_df.groupby(clusters).sum().values
    cluster_ssq_umi = expression_df.groupby(clusters).apply(
            lambda df: (df ** 2).sum()).values

    root, rd = scipy.cluster.hierarchy.to_tree(Z, rd=True)

    def de(group1, group2):
        if verbose:
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
                          index=expression_df.columns)

        df = df[df['p'] < 0.001]
        df['diff'] = df['group1'] - df['group2']

        df.sort_values('diff', ascending=False, inplace=True)

        return df


    de_dict = dict()

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

            de_dict[left_child.id, right_child.id] = de(left_clusters,
                                                        right_clusters)

        if i < 2 * n_clusters - 2:
            below = rd[i].pre_order(lambda x: x.id)
            above = [j for j in range(len(cluster_sizes)) if j not in below]

            # don't calculate redundant comparison
            if len(above) == 1:
                continue

            de_dict[i, 'all'] = de(below, above)

    return de_dict


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


def subset_exp(tenx, exp_df, filters, knncache):
    samples = tenx.cell_metadata.query(
            ' & '.join(f'{k} in {filters[k]}' for k in filters)
    ).index

    exp_subset = exp_df.loc[samples]

    if knncache is not None:
        ix = np.where(tenx.cell_metadata.index.isin(samples))[0]
        knn_subset = knncache.subset_cache(ix)
        return exp_subset, knn_subset
    else:
        return exp_subset


def plot_gene_joyplot(exp_df:pd.DataFrame, coords:pd.DataFrame,
                      gene:str, Z:np.ndarray):
    hm.joyplot(np.log2(exp_df + 1).join(coords), gene, 'cluster',
               scipy.cluster.hierarchy.leaves_list(Z))
    plt.show()
