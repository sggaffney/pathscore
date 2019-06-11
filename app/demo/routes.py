from flask import render_template, flash, redirect, url_for, abort,\
    request, current_app, send_file, send_from_directory, g, jsonify
from flask_login import current_user, login_user  # login_required
from werkzeug.utils import secure_filename
import os
import signal
from collections import OrderedDict
import pandas as pd
import numpy as np
import bokeh
from bokeh.plotting import figure, ColumnDataSource
from bokeh.layouts import widgetbox, gridplot, Spacer
from bokeh.models import Range1d
from bokeh.models.tools import HoverTool
from bokeh.models.widgets import TextInput, RadioGroup
from bokeh.models.callbacks import CustomJS

from . import demo  # FileTester, TempFile
from .helpers import get_all_demos, get_single_demo
from ..uploads import MutationFile, BmrFile
from ..errors import ValidationError
from ..models import UserFile, create_anonymous_user, initialize_project, \
    CustomBMR, BmrProcessor
from ..get_effective_pathways import run_analysis, load_pathway_list_from_file
from ..admin import zip_project
from .. import naming_rules
from .. import misc, db
from ..decorators import no_ssl, limit_user_uploads
from .. import plot_fns


DIM_COMP_SM = 60  # comparison peripheral plot small dimension (pixels)
DIM_COMP_H = 400  # comparison major plot width (pixels)
DIM_COMP_W = 550  # comparison major plot height (pixels)

SCATTER_KW = dict(size=10, color="red", alpha=0.1, line_color="firebrick",
                  line_alpha=0.5)


@demo.route('/')
# @login_required
def index():
    upload_list = get_all_demos()
    upload_list = upload_list[::-1]
    return render_template('pway/index2.html', projects=upload_list, is_demo=True)


@demo.route('/scatter')
# @login_required
@no_ssl
def scatter():
    # if proj among arguments, show this tree first.
    try:
        show_proj = int(request.args.get('proj', None))
    except (TypeError, ValueError):
        show_proj = None
    include = request.args.get('include', None)
    # list of projects (and proj_names) used to create dropdown project selector
    upload_list = get_all_demos()

    if upload_list:
        # Use specified project from args or highest file_id as CURRENT PROJECT
        current_proj = upload_list[-1]  # override if valid proj specified
        if show_proj:
            current_temp = [u for u in upload_list
                            if u.file_id == show_proj]
            # if not among user's finished projects, use highest file_id
            if len(current_temp) == 1:
                current_proj = current_temp[0]
            else:
                current_proj = upload_list[-1]
        scatter_dict = plot_fns.get_scatter_dict(current_proj)
        # includes js_name, js_inds, plot_inds, has_cnv, script, div
    else:  # no projects yet!
        flash("No demo results to show yet.", "warning")
        return redirect(url_for('.index'))

    return render_template('pway/scatter.html', current_proj=current_proj,
                           projects=upload_list, user_id=current_proj.user_id,
                           include_genes=include, **scatter_dict)


@demo.route('/mds')
# @login_required
@no_ssl
def mds():
    # if proj among arguments, show this tree first.
    try:
        show_proj = int(request.args.get('proj', None))
    except (TypeError, ValueError):
        show_proj = None
    metric = request.args.get('metric', 'jaccard')
    alg = request.args.get('alg', 'NMDS')
    include = request.args.get('include', None)
    # list of projects (and proj_names) used to create dropdown project selector

    upload_list = get_all_demos(mds_only=True)
    if not upload_list:
        flash("No project results to show yet.", "warning")
        return redirect(url_for('.index'))
    # Use specified project from args | highest file_id as CURRENT PROJECT
    current_proj = upload_list[-1]  # override if valid proj specified
    if show_proj:
        current_temp = [u for u in upload_list
                        if u.file_id == show_proj]
        # if not among user's finished projects, use highest file_id
        if len(current_temp) == 1:
            current_proj = current_temp[0]
        else:
            current_proj = upload_list[-1]

    scatter_dict = plot_fns.get_mds_dict(current_proj, metric=metric,
                                         mds_alg=alg)

    return render_template('pway/mds.html', current_proj=current_proj,
                           projects=upload_list, user_id=current_proj.user_id,
                           include_genes=include, **scatter_dict)



@demo.route('/compare')
# @login_required
@no_ssl
def compare():
    # if proj among arguments, show this tree first.
    try:
        proj_a = int(request.args.get('proj_a', None))
        proj_b = int(request.args.get('proj_b', None))
    except (TypeError, ValueError):
        proj_a = proj_b = None
    include = request.args.get('include', None)

    # list of projects (and proj_names) used to create dropdown project selector
    upload_list = get_all_demos(mds_only=False)

    if len(upload_list) > 1:
        # Use specified project from args or highest file_id as CURRENT PROJECT
        current_proj = upload_list[-1]  # override if valid proj specified
        if proj_a and proj_b:
            current_temp_a = [u for u in upload_list if u.file_id == proj_a]
            current_temp_b = [u for u in upload_list if u.file_id == proj_b]
            # if not among user's finished projects, use highest file_id
            if len(current_temp_a) == 1 and len(current_temp_b) == 1:
                current_proj_a = current_temp_a[0]
                current_proj_b = current_temp_b[0]
            else:
                current_proj_a = upload_list[-2]
                current_proj_b = upload_list[-1]
        else:
            current_proj_a = upload_list[-2]
            current_proj_b = upload_list[-1]
        user_ids = [current_proj_a.user_id, current_proj_b.user_id]
        detail_path1 = naming_rules.get_detailed_path(current_proj_a)
        detail_path2 = naming_rules.get_detailed_path(current_proj_b)
        js_name1 = naming_rules.get_js_name(current_proj_a)
        js_name2 = naming_rules.get_js_name(current_proj_b)
        xlabel = u"Effect size ({})".format(current_proj_a.get_fancy_filename())
        ylabel = u"Effect size ({})".format(current_proj_b.get_fancy_filename())

        # load pathways with 1+ mutation in 1+ patients,
        # ignoring ones with 'cancer' etc in name
        all_paths1 = load_pathway_list_from_file(detail_path1)
        all_paths2 = load_pathway_list_from_file(detail_path2)

        # IDs with p<0.05 and +ve effect
        sig_p = OrderedDict(
            [(i.path_id, i.nice_name) for i in all_paths1 if i.gene_set])
        sig_pids1 = [i for i in sig_p]
        sig_p2 = OrderedDict(
            [(i.path_id, i.nice_name) for i in all_paths2 if i.gene_set])
        sig_pids2 = [i for i in sig_p2]
        sig_p.update(sig_p2)  # ORDERED by proj1 effect size

        # BUILD DATAFRAME WITH ALL sig PATHWAYS, proj1 object order.
        pway_names = list(sig_p.values())  # order important
        columns = ['path_id', 'pname', 'ind1', 'ind2', 'e1', 'e2', 'e1_only',
                   'e2_only', 'q1', 'q2']
        df = pd.DataFrame(index=sig_p.keys(), data={'pname': pway_names},
                          columns=columns)
        for path_group, evar, qvar, ind, sigs in \
                [(all_paths1, 'e1', 'q1', 'ind1', sig_pids1),
                 (all_paths2, 'e2', 'q2', 'ind2', sig_pids2)]:
            for path in path_group:
                path_id = path.path_id
                if path_id not in sig_p:
                    continue
                df.loc[path_id, evar] = get_effect(path)
                df.loc[path_id, qvar] = get_q(path)
                temp_ind = sigs.index(path_id) if path_id in sigs else -1
                df.loc[path_id, ind] = temp_ind
        df.ind1.fillna(-1, inplace=True)
        df.ind2.fillna(-1, inplace=True)
        df.e1_only = df.where(df.e2.isnull())['e1']
        df.e2_only = df.where(df.e1.isnull())['e2']

        inds1 = list(df.ind1)
        inds2 = list(df.ind2)

        source = ColumnDataSource(data=df)
        source_full = ColumnDataSource(data=df)
        source.name, source_full.name = 'data_visible', 'data_full'
        # SET UP FIGURE
        minx = df.e1.min()
        minx *= 1 - minx / abs(minx) * 0.2
        miny = df.e2.min()
        miny *= 1 - miny/abs(miny) * 0.2
        maxx = df.e1.max() * 1.2
        maxy = df.e2.max() * 1.2
        TOOLS = "lasso_select,box_select,hover,crosshair,pan,wheel_zoom,"\
                "box_zoom,reset,tap,help" # poly_select,lasso_select, previewsave

        # SUPLOTS
        p = figure(plot_width=DIM_COMP_W, plot_height=DIM_COMP_H, tools=TOOLS,
                   title=None, toolbar_location="above",
                   x_range=Range1d(minx, maxx), y_range=Range1d(miny, maxy),
                   x_axis_type="log", y_axis_type="log"
                   )
        pb = figure(plot_width=DIM_COMP_SM, plot_height=DIM_COMP_H, tools=TOOLS,
                    y_range=p.y_range, x_axis_type="log", y_axis_type="log")
        pa = figure(plot_width=DIM_COMP_W, plot_height=DIM_COMP_SM, tools=TOOLS,
                    x_range=p.x_range, x_axis_type="log", y_axis_type="log")
        pp = figure(plot_width=DIM_COMP_SM, plot_height=DIM_COMP_SM,
                    tools=TOOLS, outline_line_color=None)

        # SPANS
        p.add_layout(plot_fns.get_span(1, 'height'))
        p.add_layout(plot_fns.get_span(1, 'width'))
        pa.add_layout(plot_fns.get_span(1, 'height'))
        pb.add_layout(plot_fns.get_span(1, 'width'))

        # STYLE
        for ax in [p, pa, pb]:
            ax.grid.visible = False
            ax.outline_line_width = 2
            ax.background_fill_color = 'whitesmoke'
        for ax in [pa, pb]:
            ax.xaxis.visible = False
            ax.yaxis.visible = False

        pa.title.text = xlabel
        pb.title.text = ylabel
        pa.title_location, pa.title.align = 'below', 'center'
        pb.title_location, pb.title.align = 'left', 'center'

        # WIDGETS
        q_input = TextInput(value='', title="P* cutoff",
                            placeholder='e.g. 0.05')
        gene_input = TextInput(value='', title="Gene list",
                               placeholder='e.g. TP53,BRAF')
        radio_include = RadioGroup(labels=["Include", "Exclude"], active=0)
        widgets = widgetbox(q_input, gene_input, radio_include, width=200,
                            css_classes=['widgets_sg'])

        grid = gridplot([[pb, p, widgets],
                         [Spacer(width=DIM_COMP_SM), pa, Spacer()]],
                        sizing_mode='fixed')

        cb_inclusion = CustomJS(args=dict(genes=gene_input), code="""
            var gene_str = genes.value
            if (!gene_str)
                return;
            var include = cb_obj.active == 0 ? true : false
            selectPathwaysByGenes(gene_str, include);
            """)
        cb_genes = CustomJS(args=dict(radio=radio_include), code="""
            var gene_str = cb_obj.value
            if (!gene_str)
                return;
            var include = radio.active == 0 ? true : false
            selectPathwaysByGenes(gene_str, include);
            """)
        radio_include.js_on_change('active', cb_inclusion)
        gene_input.js_on_change('value', cb_genes)

        # SCATTER
        p.circle("e1", "e2", source=source, **SCATTER_KW)
        pa.circle('e1_only', 1, source=source, **SCATTER_KW)
        pb.circle(1, 'e2_only', source=source, **SCATTER_KW)

        # HOVER
        for hover in grid.select(dict(type=HoverTool)):
            hover.tooltips = OrderedDict([
                ("name", "@pname"),
                ("effects", "(@e1, @e2)"),
                ("P*", ("(@q1, @q2)"))
            ])

        # ADD Q FILTERING CALLBACK
        callback = CustomJS(args=dict(source=source, full=source_full), code="""
            // get old selection indices, if any
            var prv_selected = source.selected['1d'].indices;
            var prv_select_full = []
            for(var i=0; i<prv_selected.length; i++){
                prv_select_full.push(scatter_array[prv_selected[i]])
            }
            var new_selected = []
            var q_val = cb_obj.value;
            if(q_val == '')
                q_val = 1
            var fullset = full.data;
            var n_total = fullset['e1'].length;
            // Convert float64arrays to array
            var col_names = %s ;
            col_names.forEach(function(col_name){
                source.data[col_name] = [].slice.call(source.data[col_name])
                source.data[col_name].length = 0
            })
            scatter_array.length = 0;
            var j = -1;  // new glyph indices
            for (i = 0; i < n_total; i++) {
                this_q1 = fullset['q1'][i];
                this_q2 = fullset['q2'][i];
                if(this_q1 <= q_val || this_q2 <= q_val){
                    j++; // preserve previous selection if still visible
                    col_names.forEach(function(col){
                        source.data[col].push(fullset[col][i]);
                    })
                    scatter_array.push(i)
                    if($.inArray(i, prv_select_full) > -1){
                        new_selected.push(j);
                    }
                }
            }
            source.selected['1d'].indices = new_selected;
            source.trigger('change');
            updateIfSelectionChange_afterWait();
            """ % columns)

        q_input.js_on_change('value', callback)
        script, div = plot_fns.get_bokeh_components(grid)

        proj_dir_a = naming_rules.get_project_folder(current_proj_a)
        proj_dir_b = naming_rules.get_project_folder(current_proj_b)
        if os.path.exists(os.path.join(proj_dir_a, 'matrix_svg_cnv')) and \
                os.path.exists(os.path.join(proj_dir_b, 'matrix_svg_cnv')):
            has_cnv = True
        else:
            has_cnv = False

    else:  # not enough projects yet!
        flash("Two completed projects are required for a comparison.",
              "warning")
        return redirect(url_for('.index'))

    return render_template('pway/compare.html',
                           current_projs=[current_proj_a, current_proj_b],
                           inds_use=[inds1, inds2],
                           has_cnv=has_cnv,
                           js_name_a=js_name1,
                           js_name_b=js_name2,
                           projects=upload_list,
                           bokeh_script=script,
                           bokeh_div=div, include_genes=include,
                           resources=plot_fns.resources)


def get_effect_at_index(path_list, index):
    if index < 0:
        return pd.np.nan
    else:
        temp_path = path_list[index]
        return get_effect(temp_path)


def get_effect(pway_obj):
    return float(pway_obj.n_effective) / pway_obj.n_actual


def get_q_at_index(path_list, index):
    if index < 0:
        return pd.np.nan
    else:
        temp_path = path_list[index]
        return min(1, float(temp_path.p_value) * g.n_pathways)


def get_q(pway_obj):
    return min(1, float(pway_obj.p_value) * g.n_pathways)



@demo.route('/tree')
# @login_required
def tree():
    # if proj among arguments, show this tree first.
    try:
        show_proj = int(request.args.get('proj', None))
    except (TypeError, ValueError):
        show_proj = None
    # list of projects (and proj_names) used to create dropdown project selector
    upload_list = get_all_demos()
    proj_names = {int(i.file_id): i.get_local_filename() for i in upload_list}
    if upload_list:
        # Use specified project from args or highest file_id as CURRENT PROJECT
        current_proj = upload_list[-1]  # override if valid proj specified
        if show_proj:
            current_temp = [u for u in upload_list if u.file_id == show_proj]
            # if not among user's finished projects, use highest file_id
            if len(current_temp) == 1:
                current_proj = current_temp[0]
            else:
                current_proj = upload_list[-1]
        names_odict = plot_fns.get_tree_data(current_proj)
    else:  # no projects yet!
        flash("No project results to show yet.", "warning")
        return redirect(url_for('.index'))
    return render_template('pway/tree.html', current_proj=current_proj,
                           projects=upload_list, user_id=current_proj.user_id,
                           proj_names=proj_names, names_odict=names_odict)


@demo.route('/archive/<int:proj>')
# @login_required
def archive(proj):
    upload_obj = get_single_demo(proj=proj)
    zip_path = zip_project(upload_obj)
    filename = os.path.basename(zip_path)
    return send_file(zip_path, mimetype='application/zip',
                     as_attachment=True, attachment_filename=filename)


@demo.route('/results')
# @login_required
def results():
    try:
        show_proj = int(request.args.get('proj', None))
    except (TypeError, ValueError):
        show_proj = None
    include = request.args.get('include', None)
    upload_list = get_all_demos()
    user_ids = [u.user_id for u in upload_list]
    if not upload_list:
        flash("No project results to show yet.", "warning")
        return redirect(url_for('.index'))
    if show_proj is None and len(upload_list):
        show_proj = upload_list[-1].file_id
    proj_names = {int(i.file_id): i.get_local_filename() for i in upload_list}
    user_dict = {int(i.file_id): int(i.user_id) for i in upload_list}
    return render_template('pway/show_pathways_template.html',
                           projects=upload_list, user_dict=user_dict,
                           show_proj=show_proj, proj_names=proj_names,
                           include_genes=include, is_demo=True)
