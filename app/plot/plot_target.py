"""angles are clockwise rotation from North, unless _e specified (for anticlockwise from East."""
import numpy as np

# import matplotlib as mpl
# mpl.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

from .utils import get_target_bbox

N_CUTOFF = 1000000  # max effective pathway size before downscaling.


def plot_target(n_actual=None, n_effective=None, n_max=None, pc_dict=None, 
                n_pway_genes=None, exclusive_genes=list(), n_all_patients=None,
                ax=None, cmap_name='hot_r'):
    """Plot wedges, given radius of inner ring and gene coverages."""
    if ax is None:
        hfig, ax = plt.subplots(figsize=(5.5, 5.5),
                                subplot_kw=dict(xlim=(0, 1), ylim=(0, 1),
                                                position=(0, 0, 1, 1)))
    else:
        hfig = ax.get_figure()
    seg_width = 0.01 * hfig.get_figwidth()
    rot_scale = float(20)
    rot_max = float(55)
    border_width = float(2)  # pt. for wedge border.
    fontsize = 8

    n_effective = min(n_effective, n_max)
    r_actual, r_effective = get_target_radii(n_actual, n_effective, n_max,
                                             n_cutoff=N_CUTOFF)

    gene_tuples = sorted(list(pc_dict.items()), key=(lambda v: v[1]),
                         reverse=True)
    gene_tuples.sort(key=(lambda v: v[0] in exclusive_genes), reverse=True)
    gene_list = [i[0] for i in gene_tuples]
    pc_list = np.array([i[1] for i in gene_tuples])
    n_mutated = len(pc_list)  # 1 'gap' for each mutated gene

    delta = 1.5  # default angle separation between segments
    k = 2  # arbitrary scalar greater than 1. delta inversely proportional to k.
    delta = min(delta, 360 / (k * n_mutated))  # ensure gaps don't sum>360.
    angle_per_pc = (360 - n_mutated * float(delta)) / pc_list.sum()
    pway_pc = float(n_mutated) / n_pway_genes * 100

    # RADIUS
    offset = 0.0025  # to hide inner wedge border
    r_max = r_effective + seg_width - offset
    r_text = r_max + offset + 0.001 * hfig.get_figwidth()
    
    # THETA
    angles = pc_list * angle_per_pc
    # in [0, 360) clockwise from N
    theta_starts = np.insert((angles + delta).cumsum(), 0, np.array([0]))[:-1]
    theta_ends = angles + theta_starts
    assert isinstance(theta_starts, np.ndarray)
    assert isinstance(theta_ends, np.ndarray)
    theta_text = 0.5 * (theta_starts + theta_ends)
    # in alt coords: radians and anti-clockwise from E
    theta_starts_e = -theta_ends + 90
    theta_ends_e = -theta_starts + 90
    theta_text_rad = theta_text * np.pi / 180
    theta_text_e = -theta_text + 90
    theta_text_erad = theta_text_e * np.pi / 180

    # Text positions. center on (0.5, 0.5)
    textpos, coverage_pos = dict(), dict()
    textpos['X'] = r_text * np.cos(theta_text_erad) + 0.5
    textpos['Y'] = r_text * np.sin(theta_text_erad) + 0.5
    rot_test = rot_scale * np.tan(-np.pi / 2 - theta_text_rad)
    textpos['rot'] = np.maximum(np.minimum(rot_max, rot_test), -rot_max)
    textpos['h'], textpos['v'] = zip(*[get_text_hv(a) for a in theta_text])
    
    r_coverage = r_max - float(seg_width) / 2
    coverage_pos['X'] = r_coverage * np.cos(theta_text_erad) + 0.5
    coverage_pos['Y'] = r_coverage * np.sin(theta_text_erad) + 0.5
    
    # PLOT
    patches = []
    tab_cols = []
    for i in range(len(pc_list)):
        ec = 'r' if gene_list[i] in exclusive_genes else 'b'
        tab_cols.append(ec)
        wedge = mpatches.Wedge((0.5, 0.5), r_max, theta_starts_e[i], 
                               theta_ends_e[i], width=seg_width,
                               joinstyle='miter')
        patches.append(wedge)
    cmap = plt.cm.get_cmap(cmap_name)
    wcollection = PatchCollection(patches, cmap=cmap, alpha=1, clim=(0, 100),
                                  lw=border_width, edgecolor=tab_cols,
                                  match_original=True)
    wcollection.set_array(np.array(pc_list))
    ax.add_collection(wcollection)
    ax.axis('off')
    obj_list = list()  # will hold objects with bounding boxes to inspect
    obj_list.append(wcollection)

    # text: COVERAGES
    for ind, cov in enumerate(pc_list):
        col = 'w' if cov > 60 else 'k'
        plt.text(coverage_pos['X'][ind], coverage_pos['Y'][ind], str(int(round(cov))),
                 ha='center', va='center', fontsize=fontsize, color=col)
    # text: GENE NAMES
    for ind, gene in enumerate(gene_list):
        t = plt.text(textpos['X'][ind], textpos['Y'][ind], gene, 
                     ha=textpos['h'][ind], va=textpos['v'][ind],
                     fontsize=fontsize, rotation=textpos['rot'][ind])
        obj_list.append(t)

    # RED PIE
    plt.pie([pway_pc, 100 - pway_pc], colors=['red', '#FA6864'], startangle=90, 
            counterclock=False, radius=r_effective, center=(0.5, 0.5), 
            wedgeprops=dict(lw=0))

    # BLACK DOT
    c = mpatches.Circle((0.5, 0.5), radius=r_actual, color='k')
    ax.add_artist(c)
    plt.text(0.5, 0.5, str(n_pway_genes), color='w', ha='center', va='center',
             fontsize=fontsize, fontweight='bold')

    # BLACK OUTLINE TO RED REGION
    c = mpatches.Circle((0.5, 0.5), radius=r_effective, edgecolor='k', 
                        facecolor='none', lw=0.5, capstyle='projecting')
    ax.add_artist(c)
    ax.set_position([0, 0, 1, 1])
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    bbox = get_target_bbox(obj_list, ax=ax, hfig=hfig)
    ax.set_xlim([bbox.xmin, bbox.xmax])
    ax.set_ylim([bbox.ymin, bbox.ymax])
    hfig.set_figwidth(hfig.get_figwidth() * bbox.width)
    hfig.set_figheight(hfig.get_figheight() * bbox.height)

    return ax
    

def get_target_radii(n_actual, n_effective, n_max,
                     n_cutoff=None):
    """Convert pathway sizes size to radii, using smaller of two cutoffs.

    Cutoffs are n_max (largest observed effective size or genome size)
    and n_cutoff (largest absolute effective size before rescaling happens). Radii are
    scaled down to stop growing at smaller cutoff, keeping areas in proportion.

    Assumes target is plotted in 1 unit by 1 unit area, and n_effective > n_actual."""
    n_limit = min(n_max, n_cutoff)
    if n_effective > n_limit:
        r_effective = 0.5
        # area_scale_factor = float(n_limit) / n_effective  # downscaling; < 1.
        r_scale_factor = np.sqrt(float(n_actual) / n_effective)
        r_actual = r_effective * r_scale_factor
    else:
        r_effective = 0.5 * np.sqrt(float(n_effective) / n_limit)
        r_actual = 0.5 * np.sqrt(float(n_actual) / n_limit)
    return r_actual, r_effective


def get_text_hv(theta_text):
    if 0 <= theta_text < 85:
        h, v = 'left', 'bottom'
    elif 85 <= theta_text < 95:
        h, v = 'left', 'center'
    elif 95 <= theta_text < 180:
        h, v = 'left', 'top'
    elif 180 <= theta_text < 265:
        h, v = 'right', 'top'
    elif 265 <= theta_text < 275:
        h, v = 'right', 'center'
    elif 275 <= theta_text < 360:
        h, v = 'right', 'bottom'
    else:
        h, v = 'left', 'center'
    return h, v


def norm_angle(theta):
    if theta < 0:
        theta += 360
    if theta > 360:
        theta = 360 - theta
    return theta
