import os
from collections import OrderedDict

import matplotlib as mpl
import pandas as pd
import numpy as np
import bokeh
import bokeh.plotting as bpl  # figure, ColumnDataSource, show, output_notebook
from bokeh.models.tools import HoverTool
from bokeh.resources import Resources, CDN
from bokeh.embed import components
from bokeh.models.renderers import GlyphRenderer
from bokeh.models.markers import Circle

import naming_rules
import misc
from get_effective_pathways import load_pathway_list_from_file

tools = ("hover,tap,resize,previewsave,pan,wheel_zoom,"
         "box_zoom,box_select,reset,crosshair")  # poly_select,lasso_select
plot_config = dict(plot_height=400, plot_width=600, logo=None, tools=tools,
                   toolbar_location='right',
                   min_border=0, outline_line_width=0)
scatter_config = dict(name='scattered', line_alpha=0.9, alpha=0.7)


resources = CDN  # Resources(mode="cdn")


def get_bokeh_components(plot_obj):
    """Get script and div for interactive plot, from Bokeh plot object."""
    script, div = components(plot_obj, resources)
    return script, div


def get_span(location=0.5, dim='width', **style_kws):
    """Plot span element.
    Args:
        location: axis location value (e.g. float or categorical)
        dim (str): 'width' or 'height'
    """
    style = dict(line_width=2, line_color="blue", line_dash=[6, 6])
    if style_kws:
        style.update(style_kws)
    return bokeh.models.Span(location=location, dimension=dim, **style)


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
    xyvalues = bpl.ColumnDataSource({'effect': x,
                                     'pvals': y,
                                     'D': d_vals,
                                     'pname': pnames})
    # tools = "resize,crosshair,pan,wheel_zoom,box_zoom,reset,tap," \
    #         "box_select,hover"  # poly_select,lasso_select, previewsave
    # plot_config = dict(plot_height=400, plot_width=600, logo=None,
    #                    tools=tools, title_text_font_size='14pt',
    #                    toolbar_location='right')
    plot1 = bpl.figure(title='Effect size vs p-value',
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
    script, div = get_bokeh_components(plot1)  # BOKEH OBJECTS

    js_name = naming_rules.get_js_name(upload_obj)
    # IDS. all_ids and js_ids are pathway identifiers.
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


class MDSPlotter(object):
    """Build interactive MDS plot."""
    def __init__(self, upload_obj):
        """Initialise with upload object and matrix dataframe.

        2d MDS points will be determined by metric/algorithm choice.

        Args:
            upload_obj (models.UserFile): project object
        """
        self.proj_dir = naming_rules.get_project_folder(upload_obj)
        self.proj_name = upload_obj.get_local_filename()
        self.df, self.use_paths = None, None
        # LOAD MATRIX DF
        df_path = naming_rules.get_mds_df_path(upload_obj)
        self.df = pd.read_csv(df_path, sep='\t', index_col=0)
        self.paths_loaded = False
        self.detail_path = naming_rules.get_detailed_path(upload_obj)

    def load_pathways(self):
        # BUILD USED PATHWAY LIST

        all_paths = load_pathway_list_from_file(self.detail_path)
        self.js_ids = [int(p.path_id) for p in all_paths if p.gene_set]
        path_tuples = [(int(p.path_id), p) for p in all_paths]
        self.use_ids = list(self.df.index.values)
        path_tuples = [i for i in path_tuples if i[0] in set(self.use_ids)]
        path_tuples.sort(key=lambda p: self._find_ix(p[0]))
        assert ([p[0] for p in path_tuples] == self.use_ids)
        self.use_paths = [i[1] for i in path_tuples]

        # GATHER PATHWAY ATTRIBUTE LISTS
        self.nice_names = [p.nice_name for p in self.use_paths]
        self.n_patients = self.df.sum(axis=1).values
        self.gene_sets = [MDSPlotter.line_split(list(p.gene_coverage),
                                                n=4, max_lines=1) for p in
                          self.use_paths]
        self.effect_ratios = [i.n_effective / float(i.n_actual) for i in
                              self.use_paths]
        self.n_genes_mut = [i.n_genes_mutated for i in self.use_paths]
        n_genes_frac = MDSPlotter.rescale(np.array(self.n_genes_mut), new_min=0,
                                          new_max=1, use_log=True)
        self.n_genes_cols = [bokeh.colors.RGB(r, g, b) for (r, g, b) in
                             255 * mpl.cm.viridis(n_genes_frac)[:, :3]]
        self.paths_loaded = True


    def _find_ix(self, id):
        return self.use_ids.index(id)

    def _load_mds_pts(self, metric=None, mds_alg=None):
        # LOAD SAVED POINTS (NUMPY ARRAY)
        # initial dir: /Users/sgg/Dropbox/Townsend/mds_pathways
        path_str = '{proj_dir}/npy_pts/{metric}_{alg}.npy'
        pts_path = path_str.format(proj_dir=self.proj_dir, metric=metric,
                                   alg=mds_alg)
        pts = np.load(pts_path)
        x, y = pts[:, 0], pts[:, 1]
        return x, y

    def build_mds_plot(self, metric='jaccard', mds_alg='MDS',
                       size_from_cov=False):
        """Build bokeh MDS plot. Return Bokeh 'Figure' object.

        Args:
            metric (str): distance metric for mds calculation
            mds_alg (str): MDS or NMDS
            size_from_cov (bool): sizing based on patient_cov or effect ratio."""
        if not self.paths_loaded:
            self.load_pathways()

        x, y = self._load_mds_pts(metric=metric, mds_alg=mds_alg)

        xyvalues = bpl.ColumnDataSource(
            {'x': x, 'y': y, 'nice_name': self.nice_names,
             'n_patients': self.n_patients,
             'mut_genes': self.n_genes_mut,
             'gene_set': self.gene_sets,
             'cov_size': MDSPlotter.rescale(
                 self.n_patients, new_min=5, new_max=30),
             'effect_ratio': self.effect_ratios,
             'effect_ratio_size': MDSPlotter.rescale(
                 np.array(self.effect_ratios), new_min=5, new_max=40),
             # alt size: rescale_area(np.log(effect_ratios), new_min=10),
             'n_genes_col': self.n_genes_cols
             })
        title = '{} - {} ({})'.format(self.proj_name, metric, mds_alg)
        plot1 = bpl.figure(title=title, **plot_config)
        size_attr = 'cov_size' if size_from_cov else 'effect_ratio_size'
        plot1.circle('x', 'y', size=size_attr, color='n_genes_col',
                     source=xyvalues, **scatter_config)
        # prettify_axes(plot1)
        # embolden_nonselected(plot1)

        # N_TICKS
        plot1.xaxis[0].ticker.desired_num_ticks = 0
        plot1.yaxis[0].ticker.desired_num_ticks = 0
        # HOVER TOOL
        hover = plot1.select(dict(type=HoverTool))
        from collections import OrderedDict
        hover.tooltips = OrderedDict([("name", "@nice_name"),
                                      ("n_patients", "@n_patients"),
                                      ("mut_genes", "@mut_genes"),
                                      ("effect_ratio", "@effect_ratio"),
                                      ("gene_set", "@gene_set")
                                      ])
        plot1.xaxis[0].major_tick_line_width = 0
        return plot1, title

    @staticmethod
    def rescale(vec, new_min=None, new_max=None, use_log=False):
        if use_log:
            vec = np.log(vec)
        else:
            vec = np.array(vec)
        old_min = min(vec)
        old_max = max(vec)
        old_range = old_max - old_min
        new_range = new_max - new_min
        vec_norm = (vec - old_min) / float(old_range)
        vec_new = vec_norm * new_range + new_min
        return vec_new

    @staticmethod
    def rescale_area(vec, new_min=None, new_max=None, use_log=False):
        if use_log:
            vec = np.log(vec)
        else:
            vec = np.array(vec, 'float64')
        old_min = min(vec)
        old_max = max(vec)
        old_ratio = old_max / old_min

        if new_min is not None and new_max is not None:
            raise Exception(
                'Provide only one of new_min and new_max to allow rescaling.')
        if new_min is not None:
            new_max = new_min * np.sqrt(old_ratio)
        if new_max is not None:
            new_min = new_max / np.sqrt(old_ratio)
        old_range = old_max - old_min
        new_range = new_max - new_min
        vec_norm = (vec - old_min) / float(old_range)
        vec_new = vec_norm * new_range + new_min
        return vec_new

    @staticmethod
    def line_split(vals, n=3, max_lines=3):
        """Create multi-line comma-separated string with n items per line."""
        vals = list(vals)
        out_lines = []
        for group in range(len(vals)/n + 1):
            start = group * n
            line_str = ', '.join(vals[start:start + n])
            out_lines.append(line_str)
        out_str = '\n'.join(out_lines[:max_lines])
        if len(out_lines) > max_lines:
            out_str += '...'
        return out_str


def get_mds_dict(upload_obj, metric='jaccard', mds_alg='NMDS'):
    """Create dictionary required for scatter view."""

    mds = MDSPlotter(upload_obj)
    plot1, title = mds.build_mds_plot(metric=metric, mds_alg=mds_alg)
    script, div = get_bokeh_components(plot1)

    js_name = naming_rules.get_js_name(upload_obj)
    # IDS
    all_ids = mds.use_ids
    js_ids = mds.js_ids

    # INDICES IN JS_OBJECT OF PLOT POINTS. plot->js
    js_inds = []
    for i in all_ids:
        try:
            js_inds.append(js_ids.index(i))
        except ValueError:
            js_inds.append(-1)
    assert js_inds == range(len(all_ids))
    # INDICES IN PLOT OF JS_OBJECT ITEMS (A SUBSET)
    plot_inds = range(len(all_ids))
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
        'resources': resources,
        'metric': metric,
        'alg': mds_alg
        }


text_size = '14pt'
grid_line_width = 3
grid_line_alpha = 0.4


def prettify_axes(plot):
    """Remove ticks and tick labels."""
    axx, axy = plot.xaxis[0], plot.yaxis[0]
    for ax in [axx, axy]:
        ax.axis_label_text_font_size = text_size
        ax.major_label_text_font_size = text_size
#         ax.axis_line_color = 'None'
        ax.minor_tick_line_color = None
        ax.visible = False
        ax.major_label_text_font_size = '8pt'
        ax.major_tick_line_alpha = 0
    xgrid, ygrid = plot.xgrid[0], plot.ygrid[0]
    for grid in [xgrid, ygrid]:
        grid.grid_line_width = grid_line_width
        grid.grid_line_alpha = grid_line_alpha
    plot.outline_line_alpha = 0


def embolden_nonselected(plot):
    """For circles, set opacity of nonselected."""
    circle_renderers = [i for i in plot.select(GlyphRenderer)
                        if type(i.glyph) == Circle]
    for i in circle_renderers:
        nonselected = i.nonselection_glyph
        nonselected.fill_alpha = 0.2
        nonselected.line_alpha = 0.4
