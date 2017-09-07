
import os
from collections import Counter

import fastcluster
import matplotlib
import matplotlib.pyplot as plt
import networkx
import numpy as np
import pandas as pd
import scipy
import scipy.cluster.hierarchy
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Circle
from networkx.drawing.nx_agraph import graphviz_layout
from sklearn.metrics.pairwise import cosine_distances


class KNNCache(object):
    """Cache for creating kNNs from a pre-computed bigger kNN array

    Create a cache from tenx data
    > knn_cache = KNNCache(tenx.genes.matrix, k=200)
    Save for future use
    > np.save(knn_cache_filename, knn_cache.knn_array)
    Load from file
    > knn_cache = KNNCache(knn_cache_filename)
    Get KNN graphs of various k < max_k
    > G_25 = knn_cache.get_knn_graph(25)
    > G_50 = knn_cache.get_knn_graph(50)
    """

    def __init__(self, data_or_filename, max_k:int=200):

        if isinstance(data_or_filename, scipy.sparse.csr_matrix):
            self.knn_array = cosine_distances(
                    data_or_filename).argsort(axis=1)[:,:max_k]
            self.max_k = max_k
        elif isinstance(data_or_filename, np.ndarray):
            self.knn_array = data_or_filename
            self.max_k = self.knn_array.shape[1]
        elif isinstance(data_or_filename, str):
            if os.path.exists(data_or_filename):
                self.knn_array = np.load(data_or_filename)
                self.max_k = self.knn_array.shape[1]
            else:
                raise ValueError(f'File {data_or_filename} not found')
        elif data_or_filename is not None:
            raise ValueError('Invalid value given for data_or_filename')

    def __repr__(self):
        n_samples = self.knn_array.shape[0]
        return f'Cached kNN array of {n_samples} samples, max k = {self.max_k}'

    def get_knn_graph(self, k):
        if k > self.max_k:
            raise ValueError(f'k={k} is too large for this cache')

        return networkx.from_dict_of_lists({
            i:self.knn_array[i,:k] for i in range(self.knn_array.shape[0])
        })

    def subset_cache(self, index):
        index = sorted(index)
        index_d = {j:i for i,j in enumerate(index)}

        # subset and translate rows
        knn_subset = []
        for i in index:
            knn_subset.append([index_d[j] for j in self.knn_array[i,:]
                               if j in index_d])

        max_k = min(len(ks) for ks in knn_subset)

        knn_subset = np.vstack([ks[:max_k] for ks in knn_subset])

        return KNNCache(knn_subset)


def label_propagation(exp_df:pd.DataFrame,
                      G:networkx.Graph, verbose:bool=False):
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

    communities_labelprop = np.array([node_labels[i] for i in range(G.order())])

    cluster_sum_umi = exp_df.groupby(communities_labelprop).sum()

    Z = fastcluster.linkage(cluster_sum_umi, method='average', metric='cosine')

    return communities_labelprop, Z


def network_layout(index, G:networkx.Graph):
    pos = graphviz_layout(G, prog="sfdp")
    coords = np.array([pos[i] for i in range(len(pos))])
    print(coords.shape)

    coords_df = pd.DataFrame(coords, columns=['0', '1'], index=index)

    return coords_df


class HandlerCircle(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = Circle(xy=center, radius=width / 4.0, alpha=0.4)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]


def plot_labelprop(coords_df:pd.DataFrame, Z,
                   color_by='cluster', file_name=None):
    unique_clusters = np.unique(coords_df[color_by])

    cmap = matplotlib.cm.tab20
    cmap.set_over('black')

    # no particular reason to this if the coords aren't sorted by cluster
    ix = np.random.permutation(np.arange(coords_df.shape[0], dtype=int))
    coords2 = coords_df.iloc[ix]

    fig,ax = plt.subplots(1, 2, figsize=(18,6), gridspec_kw={'wspace': 0.05})

    ax[0].scatter(coords2['0'], coords2['1'], s=60, alpha=0.8, linewidth=0,
                  color=cmap(coords2[color_by]))
    ax[0].tick_params(left='off', labelleft='off', bottom='off', labelbottom='off')

    ax[0].set_title('Network Layout')

    ddata = scipy.cluster.hierarchy.dendrogram(Z, ax=ax[1], color_threshold=0,
                                               above_threshold_color='grey')

    for i, (ic, dc) in enumerate(zip(ddata['icoord'][:-1],
                                     ddata['dcoord'][:-1])):
        x = 0.5 * sum(ic[1:3])
        y = dc[1]

        ax[1].plot(x, y, 'o', c='grey')
        ax[1].annotate(i + len(unique_clusters), (x, y), xytext=(0, -5),
                       textcoords='offset points', fontsize='large', va='top',
                       ha='center')

    ax[1].set_title('Hierarchical structure of cell clusters')
    ax[1].tick_params(labelleft='off', left='off')
    for s in ('top', 'bottom', 'left', 'right'):
        ax[1].spines[s].set_visible(False)

    lbl_rects = [(Circle((0, 0), 1, color=cmap(c)), c)
                 for c in unique_clusters]

    fig.legend(*zip(*lbl_rects), **{'handler_map': {Circle: HandlerCircle()},
                                    'loc': 7, 'fontsize': 'large'})

    if file_name:
        plt.savefig(file_name)

    plt.show()
