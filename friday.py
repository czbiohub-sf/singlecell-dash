
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

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Circle

from networkx.drawing.nx_agraph import graphviz_layout
from sklearn.neighbors import NearestNeighbors

import singlecell_dash.common as common


class HandlerCircle(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = Circle(xy=center, radius=width / 4.0, alpha=0.4)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]


def plot_labelprop_mpl(coords, communities, file_name=None, title=''):
    u_community = np.unique(communities)

    cmap = matplotlib.cm.tab20
    cmap.set_over('black')

    ix = np.random.permutation(np.arange(coords.shape[0], dtype=int))
    x = coords[ix, 0]
    y = coords[ix, 1]

    fig = matplotlib.figure.Figure(figsize=(12, 12))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    ax.scatter(x, y, s=60, alpha=0.8, linewidth=0,
               color=cmap(communities[ix]))
    ax.tick_params(left='off', labelleft='off', bottom='off', labelbottom='off')

    ax.set_title(title)

    lbl_rects = [(Circle((0, 0), 1, color=cmap(c)), c) for c in u_community]

    fig.legend(*zip(*lbl_rects), **{'handler_map': {Circle: HandlerCircle()},
                                    'loc': 7, 'fontsize': 'large'})

    if file_name:
        FigureCanvasAgg(fig).print_figure(file_name)


def label_propagation(G, verbose=False):
    node_labels = {node: i for i, node in enumerate(G.nodes())}

    n_changes = 1
    while n_changes:
        n_changes = 0
        for node in G.nodes():
            neighbor_labels = Counter([node_labels[n] for n in G.neighbors(node)])
            pop_label = neighbor_labels.most_common(1)[0][0]
            if node_labels[node] != pop_label:
                node_labels[node] = pop_label
                n_changes += 1
        if verbose:
            print("Round with ", n_changes, " to labels.")

    label_renames = {label: i for i, (label, c)
                     in enumerate(Counter(node_labels.values()).most_common())}

    for node in node_labels:
        node_labels[node] = label_renames[node_labels[node]]

    if verbose:
        print("Most common labels, in the form label, frequency")
        print(Counter(node_labels.values()).most_common())

    return node_labels


def network_layout(matrix, k=30):
    nbrs = NearestNeighbors(k, algorithm='brute',
                            metric='cosine', n_jobs=-1).fit(matrix)
    G = networkx.from_scipy_sparse_matrix(nbrs.kneighbors_graph(matrix))

    node_labels = label_propagation(G, verbose=True)
    communities_labelprop = np.array([node_labels[i] for i in range(matrix.shape[0])])

    pos = graphviz_layout(G, prog="sfdp")
    coords = np.array([pos[i] for i in range(len(pos))])
    print(coords.shape)

    return coords, communities_labelprop


def expression(matrix, group):
    g = matrix[group,:].tocsc()
    mu = np.asarray(g.mean(0)).flatten()
    std = np.sqrt(np.asarray((g.power(2)).mean(0)).flatten() - mu ** 2)
    percent_nz = 100*np.asarray((g > 0).mean(0)).flatten()

    return mu, std, percent_nz


def cluster_expression(tenx, clusters, skip=1):
    df = pd.DataFrame(index=tenx.genes.columns)

    for c in np.unique(clusters):
        mu, std, percent_nz = expression(tenx.genes[::skip,:], clusters == c)
        df[f'Cluster {c} mean UMI'] = mu
        df[f'Cluster {c} std UMI'] = std
        df[f'Cluster {c} % present'] = percent_nz

    return df


def load_tissue(tissue, data_folder, channels_to_use = None):

    genes_to_drop = 'Malat1|Rn45s|Rpl10|Rpl10a|Rpl10l|Rpl11|Rpl12|Rpl13|Rpl13a|Rpl14|Rpl15|Rpl17|Rpl18|Rpl18a|Rpl19|Rpl21|Rpl22|Rpl22l1|Rpl23|Rpl23a|Rpl24|Rpl26|Rpl27|Rpl27a|Rpl28|Rpl29|Rpl3|Rpl30|Rpl31|Rpl31-ps12|Rpl32|Rpl34|Rpl34-ps1|Rpl35|Rpl35a|Rpl36|Rpl36a|Rpl36al|Rpl37|Rpl37a|Rpl38|Rpl39|Rpl39l|Rpl3l|Rpl4|Rpl41|Rpl5|Rpl6|Rpl7|Rpl7a|Rpl7l1|Rpl8|Rpl9|Rplp0|Rplp1|Rplp2|Rplp2-ps1|Rps10|Rps11|Rps12|Rps13|Rps14|Rps15|Rps15a|Rps15a-ps4|Rps15a-ps6|Rps16|Rps17|Rps18|Rps19|Rps19-ps3|Rps19bp1|Rps2|Rps20|Rps21|Rps23|Rps24|Rps25|Rps26|Rps27|Rps27a|Rps27l|Rps28|Rps29|Rps3|Rps3a|Rps4x|Rps4y2|Rps5|Rps6|Rps6ka1|Rps6ka2|Rps6ka3|Rps6ka4|Rps6ka5|Rps6ka6|Rps6kb1|Rps6kb2|Rps6kc1|Rps6kl1|Rps7|Rps8|Rps9|Rpsa'.split(
        '|')

    tenx = common.TenX_Runs(data_folder, tissue=tissue, verbose=True, genes_to_drop=genes_to_drop,
                            channels_to_use=channels_to_use)

    return tenx


def cluster(tenx, skip, file_format, k, tissue):

    coords, communities_labelprop = network_layout(tenx.genes.matrix[::skip], k=k)

    coords_df = pd.DataFrame({'0': coords[:, 0], '1': coords[:, 1], 'cluster': communities_labelprop},
                             index=tenx.genes.rows[::skip])

    file_name = file_format.format('smushed', 'csv')
    coords_df.to_csv(file_name)

    file_name = file_format.format('embedding','png')
    plot_labelprop_mpl(coords, communities_labelprop, title=f'{tissue}: Graph layout of clusters',
                       file_name=file_name)

    return communities_labelprop


def diff_exp_clusters(cluster_expression_df, cluster_sizes, file_format):
    n_clusters = len(cluster_sizes)

    cluster_sum_umi = np.vstack(
            [cluster_sizes[c] * cluster_expression_df[f'Cluster {c} mean UMI'].values
             for c in range(n_clusters)]
    )

    cluster_ssq_umi = np.vstack(
            [cluster_sizes[c] * (cluster_expression_df[f'Cluster {c} std UMI'].values ** 2
                     + cluster_expression_df[f'Cluster {c} mean UMI'].values ** 2)
             for c in range(n_clusters)]
    )

    Z = fastcluster.linkage(cluster_sum_umi, method='average', metric='cosine')

    fig = matplotlib.figure.Figure(figsize=(12, 12))

    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    scipy.cluster.hierarchy.dendrogram(Z, ax=ax,
                                       color_threshold=0, above_threshold_color='grey')

    ax.set_title('Hierarchical structure of cell-type clusters')
    ax.set_xlabel('Cluster Label')
    ax.tick_params(labelleft='off')

    FigureCanvasAgg(fig).print_figure(file_format.format('dendrogram', 'png'))


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


if __name__ == '__main__':
    all_tissues = ['Tongue', 'Liver', 'Bladder', 'Kidney', 'Spleen', 'Marrow',
                   'Lung', 'Muscle', 'Heart', 'Thymus', 'Mammary']

    channels_to_use = ['10X_P4_0', '10X_P4_1', '10X_P4_2', '10X_P4_3',
                       '10X_P4_4', '10X_P4_5', '10X_P4_6', '10X_P4_7',
                       '10X_P6_0', '10X_P6_1', '10X_P6_2', '10X_P6_3',
                       '10X_P6_4', '10X_P6_5', '10X_P6_6', '10X_P6_7',
                       '10X_P7_0', '10X_P7_1', '10X_P7_2', '10X_P7_3',
                       '10X_P7_4', '10X_P7_5', '10X_P7_6', '10X_P7_7',
                       '10X_P7_8', '10X_P7_9', '10X_P7_10', '10X_P7_11',
                       '10X_P7_12', '10X_P7_13', '10X_P7_14', '10X_P7_15']

    parser = argparse.ArgumentParser()

    parser.add_argument('--data_folder')
    parser.add_argument('--n_samples', type=int, default=-1)
    parser.add_argument('--tissues', nargs='*', default=None)
    parser.add_argument('--k', type=int, default=25)

    args = parser.parse_args()

    if args.tissues is None:
        args.tissues = all_tissues

    for tissue in args.tissues:
        print(f'Processing {tissue}...')
        tenx = load_tissue(tissue, args.data_folder,
                           channels_to_use=channels_to_use)


        if not os.path.exists(os.path.join(f'{args.data_folder}',
                                           '10x_data', 'tissues', tissue)):
            os.mkdir(os.path.join(f'{args.data_folder}',
                                  '10x_data', 'tissues', tissue))

        if args.n_samples < 1:
            skip = 1
            file_format = os.path.join(args.data_folder, '10x_data', 'tissues',
                                       tissue, '{}.{}')
        else:
            skip = tenx.genes.matrix.shape[0] // args.n_samples
            skip = max(skip, 1)

            file_format = os.path.join(
                    args.data_folder, '10x_data', 'tissues', tissue,
                    f'{{}}-{tissue}-{args.n_samples}-{args.k}.{{}}'
            )

        clusters = cluster(tenx, skip, file_format=file_format,
                           k=args.k, tissue=tissue)

        print('Computing cluster expression...')

        # genewise mean expression and percent non-zero for each cluster
        cluster_expression_df = cluster_expression(tenx, clusters, skip)

        # drop zeros
        cluster_expression_df = cluster_expression_df.loc[
            cluster_expression_df.max(axis=1) != 0
        ]

        # round for readability and output to csv
        cluster_expression_df = np.round(cluster_expression_df, 2)
        cluster_expression_df.to_csv(file_format.format('cluster-expresion', 'csv'))

        cluster_sizes = dict(Counter(clusters).most_common())

        print('Computing differential expression...')

        group_list = diff_exp_clusters(cluster_expression_df, cluster_sizes,
                                       file_format)

        with open(file_format.format('summary', 'txt'), 'w') as OUT:
            print('Clustered {} cells into {} clusters'.format(
                    sum(cluster_sizes.values()), len(cluster_sizes)),
                  file=OUT)
            print('\t'.join(('group', 'n_cells', 'member_clusters')), file=OUT)
            for i,gl in group_list:
                print('{}\t{}\t{}'.format(i, sum(cluster_sizes[j] for j in gl),
                                          ', '.join(map(str, gl))), file=OUT)
