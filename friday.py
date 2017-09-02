
import itertools
from collections import OrderedDict, Counter
import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt # necessary?
# import seaborn as sns

import scipy.spatial.distance as sdist
import scipy.stats as stats
import scipy.sparse as sparse
import scipy.cluster.hierarchy

from sklearn.neighbors import NearestNeighbors, kneighbors_graph

import fastcluster

import networkx
from networkx.drawing.nx_agraph import graphviz_layout

import plotly
import plotly.graph_objs as go

import singlecell_dash.common as common
import sparse_dataframe
import matplotlib.figure
import matplotlib.colors
from matplotlib.patches import Circle
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.legend_handler import HandlerPatch


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



def plot_labelprop_mpl_tissues(coords, tissue_list, file_name=None, title=''):
    u_tissue = set(tissue_list)

    tissue_d = {t:i for i,t in enumerate(u_tissue)}

    cmap = matplotlib.cm.tab20
    cmap.set_over('black')

    ix = np.random.permutation(np.arange(coords.shape[0], dtype=int))
    x = coords[ix, 0]
    y = coords[ix, 1]

    fig = matplotlib.figure.Figure(figsize=(12, 12))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    ax.scatter(x, y, s=60, alpha=0.8, linewidth=0,
               color=[cmap(tissue_d[tissue_list[i]]) for i in ix])
    ax.tick_params(left='off', labelleft='off', bottom='off', labelbottom='off')

    ax.set_title(title)

    lbl_rects = [(Circle((0, 0), 1, color=cmap(tissue_d[t])), t) for t in u_tissue]

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
        if (verbose):
            print("Round with ", n_changes, " to labels.")

    labels = list(set(node_labels.values()))
    label_renames = {label: i for i, (label, c)
                     in enumerate(Counter(node_labels.values()).most_common())}

    for node in node_labels:
        node_labels[node] = label_renames[node_labels[node]]

    if (verbose):
        print("Most common labels, in the form label, frequency")
        print(Counter(node_labels.values()).most_common())

    return node_labels


def network_layout(matrix, k=30):
    nbrs = NearestNeighbors(k, algorithm='brute', metric='cosine').fit(matrix)
    G = networkx.from_scipy_sparse_matrix(nbrs.kneighbors_graph(matrix))

    node_labels = label_propagation(G, verbose=True)
    communities_labelprop = np.array([node_labels[i] for i in range(matrix.shape[0])])

    pos = graphviz_layout(G, prog="sfdp")
    coords = np.array([pos[i] for i in range(len(pos))])
    print(coords.shape)

    return coords, communities_labelprop


def plot_labelprop(coords, communities, filename=None, title=''):
    scatters = []

    n_communities = len(np.unique(communities))

    scatter = go.Scatter(x=coords[0],
                         y=coords[1],
                         mode='markers',
                         name=communities,
                         hoverinfo='text',
                         text=communities,
                         marker=dict(color=communities,
                                     cmin=0.0, cmax=n_communities - 1,
                                     colorscale='Viridis'
                                     ))

    fig = {
        'data': [scatter],
        "layout": go.Layout(title=f'Graph layout of cell type clusters: {title}',
                            xaxis={'showticklabels': False},
                            yaxis={'showticklabels': False},
                            hovermode='closest', dragmode='select',
                            show_legend=True)
    }

    if filename is None:
        plotly.offline.iplot(fig)
    else:
        plotly.offline.plot(fig, filename=filename, image='png')


def tenx_diff_exp(matrix, columns, communities):
    diff_expr_dfs = []

    for c in np.unique(communities):
        group1 = (communities == c)
        group2 = (communities != c)
        diff_expr_df = sparse_diff_exp(matrix,
                                       group1, group2,
                                       columns).sort_values('p')
        diff_expr_df['community'] = c
        diff_expr_dfs.append(diff_expr_df)

    diff_expr = pd.concat(diff_expr_dfs)
    print(diff_expr.shape)
    print(diff_expr.head())
    return diff_expr


def tenx_diff_exp_all(tenx_data, communities):
    diff_expr_dfs = []

    for c1, c2 in itertools.combinations(np.unique(communities), 2):
        group1 = (communities == c1)
        group2 = (communities == c2)
        diff_expr_df = sparse_diff_exp(tenx_data.genes.matrix,
                                       group1, group2, tenx_data.genes.columns).sort_values('p')
        diff_expr_df['community1'] = c1
        diff_expr_df['community2'] = c2
        diff_expr_dfs.append(diff_expr_df)

    diff_expr = pd.concat(diff_expr_dfs)
    print(diff_expr.shape)
    print(diff_expr.head())
    return diff_expr


def diff_exp(counts, group1, group2):
    """Computes differential expression between group 1 and group 2
    for each column in the dataframe counts.
    Returns a dataframe of Z-scores and p-values."""

    mean_diff = counts.loc[group1].mean() - counts.loc[group2].mean()
    pooled_sd = np.sqrt(counts.loc[group1].var() / len(group1)
                        + counts.loc[group2].var() / len(group2))
    z_scores = mean_diff / pooled_sd
    z_scores = z_scores.fillna(0)

    # t-test
    p_vals = (1 - stats.norm.cdf(np.abs(z_scores))) * 2

    df = pd.DataFrame({'z': z_scores})
    df['p'] = p_vals

    return df

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

def sparse_diff_exp(matrix, group1, group2, index):
    """Computes differential expression between group 1 and group 2
    for each column in the dataframe counts.

    Returns a dataframe of Z-scores and p-values."""

    g1 = matrix[group1, :]
    g2 = matrix[group2, :]

    g1mu = np.asarray(g1.mean(0)).flatten()
    g2mu = np.asarray(g2.mean(0)).flatten()

    mean_diff = g1mu - g2mu
    # E[X^2] - (E[X])^2
    pooled_sd = np.sqrt(((g1.power(2)).mean(0) - g1mu ** 2) / len(group1)
                        + ((g2.power(2)).mean(0) - g2mu ** 2) / len(group2))
    pooled_sd = np.asarray(pooled_sd).flatten()

    z_scores = np.zeros_like(pooled_sd)
    nz = pooled_sd > 0
    z_scores[nz] = np.nan_to_num(mean_diff[nz] / pooled_sd[nz])

    # t-test
    p_vals = np.clip((1 - stats.norm.cdf(np.abs(z_scores))) * 2 * matrix.shape[1], 0, 1)

    df = pd.DataFrame(OrderedDict([('z', z_scores), ('p', p_vals),
                                   ('mean1', g1mu), ('mean2', g2mu)]),
                      index=index)

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
    cluster_sizes_dict = dict(cluster_sizes)

    cluster_sum_umi = np.vstack([size*cluster_expression_df[f'Cluster {c} mean UMI'].values
                                 for c, size in cluster_sizes])

    cluster_ssq_umi = np.vstack([size*(
                                     cluster_expression_df[f'Cluster {c} std UMI'].values**2 +
                                     cluster_expression_df[f'Cluster {c} mean UMI'].values**2)
                                 for c, size in cluster_sizes])

    Z = fastcluster.linkage(cluster_sum_umi, method='average', metric='cosine')

    scipy.cluster.hierarchy.dendrogram(Z)
    plt.savefig(file_format.format('dendrogram','png'))

    root, rd = scipy.cluster.hierarchy.to_tree(Z, rd=True)

    def de(group1, group2):
        print(f'Comparing {group1} to {group2}')

        left_n_cells = sum(cluster_sizes_dict[c] for c in group1)
        right_n_cells = sum(cluster_sizes_dict[c] for c in group2)

        left_mean = cluster_sum_umi[group1, :].sum(axis=0) / left_n_cells
        right_mean = cluster_sum_umi[group2, :].sum(axis=0) / right_n_cells

        mean_diff = left_mean - right_mean

        # here is where things went wrong

        left_var = cluster_ssq_umi[group1, :].sum(axis=0) / left_n_cells - left_mean ** 2
        right_var = cluster_ssq_umi[group2, :].sum(axis=0) / right_n_cells - right_mean ** 2

        pooled_sd = np.sqrt(left_var / left_n_cells + right_var / right_n_cells)

        z_scores = np.zeros_like(pooled_sd)
        nz = pooled_sd > 0
        z_scores[nz] = np.nan_to_num(mean_diff[nz] / pooled_sd[nz])

        # t-test
        p_vals = np.clip((1 - stats.norm.cdf(np.abs(z_scores))) * 2 * z_scores.shape[0], 0, 1)

        df = pd.DataFrame(OrderedDict([('z', z_scores), ('p', p_vals),
                                       ('left_mean', left_mean), ('right_mean', right_mean)]),
                          index=cluster_expression_df.index)

        df = df[df['p'] < 0.001]
        df['diff'] = df['left_mean'] - df['right_mean']

        df.sort_values('diff', ascending=False, inplace=True)
        # df.drop('scale', axis=1, inplace=True)

        name = 'dge_' + '-'.join(map(str, group1)) + '_vs_' + '-'.join(map(str, group2))
        df.to_csv(file_format.format(name, 'csv'))
        

    for i in range(len(cluster_sizes), 2*len(cluster_sizes)-1):
        left_clusters = (rd[i].get_left().pre_order(lambda x: x.id))

        right_clusters = (rd[i].get_right().pre_order(lambda x: x.id))
        de(left_clusters, right_clusters)

    for i in range(0, 2 * len(cluster_sizes) - 2):

        below = rd[i].pre_order(lambda x: x.id)
        above = [j for j in range(len(cluster_sizes)) if j not in below]

        de(below, above)


def all_tissues(tissues, data_folder, samples, ks, channels_to_use=None):
    genes_to_drop = 'Malat1|Rn45s|Rpl10|Rpl10a|Rpl10l|Rpl11|Rpl12|Rpl13|Rpl13a|Rpl14|Rpl15|Rpl17|Rpl18|Rpl18a|Rpl19|Rpl21|Rpl22|Rpl22l1|Rpl23|Rpl23a|Rpl24|Rpl26|Rpl27|Rpl27a|Rpl28|Rpl29|Rpl3|Rpl30|Rpl31|Rpl31-ps12|Rpl32|Rpl34|Rpl34-ps1|Rpl35|Rpl35a|Rpl36|Rpl36a|Rpl36al|Rpl37|Rpl37a|Rpl38|Rpl39|Rpl39l|Rpl3l|Rpl4|Rpl41|Rpl5|Rpl6|Rpl7|Rpl7a|Rpl7l1|Rpl8|Rpl9|Rplp0|Rplp1|Rplp2|Rplp2-ps1|Rps10|Rps11|Rps12|Rps13|Rps14|Rps15|Rps15a|Rps15a-ps4|Rps15a-ps6|Rps16|Rps17|Rps18|Rps19|Rps19-ps3|Rps19bp1|Rps2|Rps20|Rps21|Rps23|Rps24|Rps25|Rps26|Rps27|Rps27a|Rps27l|Rps28|Rps29|Rps3|Rps3a|Rps4x|Rps4y2|Rps5|Rps6|Rps6ka1|Rps6ka2|Rps6ka3|Rps6ka4|Rps6ka5|Rps6ka6|Rps6kb1|Rps6kb2|Rps6kc1|Rps6kl1|Rps7|Rps8|Rps9|Rpsa'.split('|')

    file_prefix = data_folder + '/10x_data/tissues/'

    all_genes = sparse_dataframe.SparseDataFrame()
    tissue_list = []

    for tissue in tissues:
        tenx = common.TenX_Runs(data_folder, tissue=tissue, verbose=True,
                                genes_to_drop=genes_to_drop,
                                channels_to_use=channels_to_use)


        skip = tenx.genes.matrix.shape[0] // samples
        skip = max(skip, 1)
        print(tenx.genes.matrix.shape, skip)

        if all_genes.matrix is None:
            all_genes.columns = tenx.genes.columns[:]
            all_genes.rows = tenx.genes.rows[::skip]
            all_genes.matrix = tenx.genes.matrix[::skip]
        else:
            all_genes.rows.extend(tenx.genes.rows[::skip])
            all_genes.matrix = sparse.vstack((all_genes.matrix, tenx.genes.matrix[::skip]),
                                             format='csc')

        tissue_list.extend(tissue for i in range(0, tenx.genes.matrix.shape[0], skip))

    print(all_genes.matrix.shape)

    print(len(all_genes.rows))

    for k in ks:
        file_suffix = f'-all-{samples}-{k}'

        coords, communities_labelprop = network_layout(all_genes.matrix, k=k)

        coords_df = pd.DataFrame({'0': coords[:, 0],
                                  '1': coords[:, 1],
                                  'cluster': communities_labelprop},
                                 index=all_genes.rows)

        coords_df.to_csv(file_prefix + 'smushed' + file_suffix + '.csv')

        plot_labelprop_mpl_tissues(coords, tissue_list, title='Graph layout of tissues',
                           file_name=file_prefix + 'tissue-embedding' + file_suffix + '.png')

        plot_labelprop_mpl(coords, communities_labelprop,
                           title='Graph layout of clusters',
                           file_name=file_prefix + 'embedding' + file_suffix + '.png')

        diff_expr_df = tenx_diff_exp(all_genes.matrix,
                                     all_genes.columns,
                                     communities_labelprop)

        de = diff_expr_df[diff_expr_df['p'] < 0.001]
        de = de[de['z'] > 0]
        de['scale'] = np.abs(de['mean1'] - de['mean2'])

        bigten = de.groupby('community')['mean1'].nlargest(50)
        bigten = pd.DataFrame(bigten)
        bigten.columns = ['Mean UMIs']
        bigten.to_csv(file_prefix + 'diff exp' + file_suffix + '.csv')



if __name__ == '__main__':
    tissues = ['Tongue', 'Liver', 'Bladder', 'Kidney', 'Spleen', 'Marrow', 'Lung',
       'Muscle', 'Heart'] #, 'Thymus', 'Mammary']

    channels_to_use = ['10X_P4_0',
                       '10X_P4_1',
                       '10X_P4_2',
                       '10X_P4_3',
                       '10X_P4_4',
                       '10X_P4_5',
                       '10X_P4_6',
                       '10X_P4_7',
                       '10X_P6_0',
                       '10X_P6_1',
                       '10X_P6_2',
                       '10X_P6_3',
                       '10X_P6_4',
                       '10X_P6_5',
                       '10X_P6_6',
                       '10X_P6_7',
                       '10X_P7_0',
                       '10X_P7_1',
                       '10X_P7_2',
                       '10X_P7_3',
                       '10X_P7_4',
                       '10X_P7_5',
                       '10X_P7_6',
                       '10X_P7_7',
                       '10X_P7_8',
                       '10X_P7_9',
                       '10X_P7_10',
                       '10X_P7_11',
                       '10X_P7_12',
                       '10X_P7_13',
                       '10X_P7_14',
                       '10X_P7_15']

    import sys
    for tissue in tissues[6:7]:
        tissue = 'Mammary'
        print(f'Processing {tissue}...')
        samples = int(sys.argv[1])
        k = 25
        data_folder = '/data1/maca'
        tenx = load_tissue(tissue, data_folder, channels_to_use=channels_to_use)

        skip = tenx.genes.matrix.shape[0] // samples
        skip = max(skip, 1)

        if samples < 1:
            skip = 1

        file_suffix = f'-{tissue}-{samples}-{k}'
        if not os.path.exists(f'{data_folder}/10x_data/tissues/{tissue}'):
            os.mkdir(f'{data_folder}/10x_data/tissues/{tissue}')
        file_format = data_folder + '/10x_data/tissues/' + tissue + '/{}' + file_suffix + '.{}'

        if samples < 1:
            file_format = data_folder + '/10x_data/tissues/' + tissue + '/{}' + '.{}'

        clusters = cluster(tenx, skip, file_format=file_format, k=k, tissue=tissue)

        print('Computing cluster expression.')
        cluster_expression_df = cluster_expression(tenx, clusters, skip)
        cluster_expression_df = cluster_expression_df.loc[cluster_expression_df.max(axis=1) != 0]
        cluster_expression_df = np.round(cluster_expression_df, 2)

        #cluster_expression_df.to_csv(data_folder + f'/10x_data/tissues/{tissue}-cluster-expression.csv')

        cluster_sizes = Counter(clusters).most_common()
        pd

        print('Computing differential expression.')

        diff_exp_clusters(cluster_expression_df, cluster_sizes, file_format)
        #diff_exp_clusters(tenx, skip, clusters, file_format)

    #all_tissues(tissues, '/data1/maca', 500, ks=(25,50), channels_to_use=channels_to_use)


