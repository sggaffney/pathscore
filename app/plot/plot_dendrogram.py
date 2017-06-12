
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

def plot_dendrogram(m, out_svg='test_dendrogram.svg'):
    """
    Dendrogram plot, using distance matrix, calculating 'average' linkage.
    Args:
        m (np.array): Square distance matrix.
        out_svg (str): Path for dendrogram
    Returns:
        leaf_indices (list of ints): reordering after clustering
    """
    n_pways = 50
    leaf_sep = 21  # pixels, given 96 ppi in css
    dpi = 96
    px_width = 250
    fig_width = float(px_width) / dpi
    fig_height = float(leaf_sep) * n_pways / dpi

    # Hierarchical clustering
    c = squareform(m)
    Z = linkage(c, method='average')

    hf, hax = plt.subplots(1, figsize=(fig_width, fig_height))
    d = dendrogram(
            Z,
            leaf_font_size=8,  # font size for the x axis labels. default None.
            orientation='left',  # default top

            distance_sort=False,
            ax=hax,

            color_threshold=None,
            get_leaves=True,
            labels=None,
            count_sort='descendent',
            show_leaf_counts=True,
            no_plot=False,
            no_labels=False,
            leaf_rotation=None,
            leaf_label_func=None,
            link_color_func=None,

            above_threshold_color='b',
            p=30,
            truncate_mode=None,
            show_contracted=False,
    )

    hax.set_position([0, 0, 0.98, 1])  # allow margin to show vertical lines
    hax.set_axis_off()
    # hf.savefig('test_dendrogram.png', dpi=dpi)
    hf.savefig(out_svg, dpi=dpi)
    leaf_indices = [int(i) for i in d['ivl']][::-1]
    return leaf_indices
