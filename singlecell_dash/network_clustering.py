
from collections import Counter, OrderedDict

import matplotlib.colors
import matplotlib.figure

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Circle

import networkx
from networkx.drawing.nx_agraph import graphviz_layout

from sklearn.neighbors import NearestNeighbors



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
                            metric='cosine').fit(matrix)
    G = networkx.from_scipy_sparse_matrix(nbrs.kneighbors_graph(matrix))

    node_labels = label_propagation(G, verbose=True)
    communities_labelprop = np.array([node_labels[i] for i in range(matrix.shape[0])])

    pos = graphviz_layout(G, prog="sfdp")
    coords = np.array([pos[i] for i in range(len(pos))])
    print(coords.shape)

    return coords, communities_labelprop


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
