import numpy as np
from matplotlib.transforms import Bbox


def get_target_bbox(obj_list, ax=None, hfig=None):
    """Get bounding box for iterable of plot objects.

    SEE:
    # BBOX: https://stackoverflow.com/questions/24581194/matplotlib-text-bounding-box-dimensions
    # https://stackoverflow.com/questions/22667224/matplotlib-get-text-bounding-box-independent-of-backend
    """
    transf = ax.transData.inverted()
    renderer = hfig.canvas.get_renderer()
    bb_list = []
    for obj in obj_list:
        bb_list.append(obj.get_tightbbox(renderer=renderer).transformed(transf))
    extents = np.vstack([bb.extents for bb in bb_list])
    # print(extents)

    x0 = max(extents[:, 0]) if ax.xaxis_inverted() else min(extents[:, 0])
    y0 = max(extents[:, 1]) if ax.yaxis_inverted() else min(extents[:, 1])
    x1 = min(extents[:, 2]) if ax.xaxis_inverted() else max(extents[:, 2])
    y1 = min(extents[:, 3]) if ax.yaxis_inverted() else max(extents[:, 3])

    bbox = Bbox(np.array([[x0, y0], [x1, y1]]))
    # print(bbox)
    return bbox
