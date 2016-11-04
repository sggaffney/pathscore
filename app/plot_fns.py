import os
from collections import OrderedDict

import numpy as np
from bokeh.resources import Resources
from bokeh.embed import components
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models.tools import HoverTool

import naming_rules
import misc
from get_effective_pathways import load_pathway_list_from_file


def get_scatter_dict(upload_obj):
    """Create dictionary required for scatter view."""
    detail_path = naming_rules.get_detailed_path(upload_obj)
    all_pathways = load_pathway_list_from_file(detail_path)
    data_pways, data_pvals, data_effect, data_d = [], [], [], []
    for p in all_pathways:
        pval = float(p.p_value)
        if pval >= 0.05:
            continue
        data_pvals.append(pval)
        data_effect.append(float(p.n_effective) / p.n_actual)
        data_pways.append(p)
        data_d.append(p.D)
    x = np.log2(np.array(data_effect))  # effect size
    y = -np.log10(np.array(data_pvals))  # p-value
    d_vals = np.array(data_d)  # alternatively: np.log2...
    # adjust zero pvalues. e-15.9 seems to be minimum.
    max_y = max([np.ceil(max([i for i in x if i != np.inf])),
                 np.float64(17)])
    y[y == np.inf] = max_y
    pnames = [misc.strip_contributors(p.nice_name) for p in data_pways]
    xyvalues = ColumnDataSource({'effect': x,
                                 'pvals': y,
                                 'D': d_vals,
                                 'pname': pnames})
    tools = "resize,crosshair,pan,wheel_zoom,box_zoom,reset,tap," \
            "box_select,hover"  # poly_select,lasso_select, previewsave
    plot_config = dict(plot_height=400, plot_width=600, logo=None,
                       tools=tools, title_text_font_size='14pt',
                       toolbar_location='right')

    plot1 = figure(title='Effect size vs p-value',
                   x_axis_label="log2 fold change",
                   y_axis_label="-log10 p-value",
                   **plot_config)
    plot1.xaxis.axis_label_text_font_size = "12pt"
    plot1.yaxis.axis_label_text_font_size = "12pt"
    plot1.scatter('effect', 'pvals', source=xyvalues, size=10, color="red",
                  alpha=0.1, marker="circle", line_color="firebrick",
                  line_alpha=0.5)

    hover = plot1.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([("name", "@pname")])

    resources = Resources(mode="cdn")
    script, div = components(plot1, resources)
    js_name = naming_rules.get_js_name(upload_obj)
    # IDS
    all_ids = [p.path_id for p in all_pathways if float(p.p_value) < 0.05]
    js_ids = [p.path_id for p in all_pathways if p.gene_set]
    # INDICES IN JS_OBJECT OF PLOT POINTS. plot->js
    js_inds = []
    for i in all_ids:
        try:
            js_inds.append(js_ids.index(i))
        except ValueError:
            js_inds.append(-1)
    # INDICES IN PLOT OF JS_OBJECT ITEMS (A SUBSET)
    plot_inds = [all_ids.index(i) for i in js_ids]
    proj_dir = naming_rules.get_project_folder(upload_obj)
    if os.path.exists(os.path.join(proj_dir, 'matrix_svg_cnv')):
        has_cnv = True
    else:
        has_cnv = False

    return {
        'js_name': js_name,
        'js_inds': js_inds,
        'plot_inds': plot_inds,
        'has_cnv': has_cnv,
        'bokeh_script': script,
        'bokeh_div': div,
        'resources': resources
        }


def get_tree_data(upload_obj):
    # load ordered dictionary of path_ids : pathway_display_name
    names_path = naming_rules.get_tree_score_paths(upload_obj)[1]
    # tree_path = naming_rules.get_apache_path(tree_path)
    names_odict = OrderedDict()  # ordered dictionary of path_id: name
    with open(names_path, 'rU') as f:
        for line in f:
            vals = line.strip('\n').split('\t')
            if len(vals) != 2:
                continue
            names_odict[vals[0]] = vals[1]
    return names_odict