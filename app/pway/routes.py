from flask import render_template, flash, redirect, url_for, abort,\
    request, current_app, send_file, send_from_directory, g
from flask_login import login_required, current_user, login_user
from werkzeug.utils import secure_filename
import os
import signal
from collections import OrderedDict
import numpy as np
from bokeh.plotting import figure, ColumnDataSource, hplot
from bokeh.resources import Resources
from bokeh.embed import components
from bokeh.models.tools import HoverTool
from bokeh.models.renderers import GlyphRenderer
from bokeh.models.markers import Circle
from bokeh.models.widgets import TextInput
from bokeh.models import Callback
from bokeh.io import vform

from . import pway  # FileTester, TempFile
from ..maf import MutationFile
from ..errors import ValidationError
from .forms import UploadForm
from ..models import UserFile, create_anonymous_user, initialize_project
from ..get_effective_pathways import run_analysis, load_pathway_list_from_file
from ..admin import zip_project
from .. import naming_rules
from .. import misc
from ..decorators import no_ssl, limit_user_uploads


@pway.route('/')
@login_required
def index():
    upload_list = UserFile.query.filter_by(user_id=current_user.id).all()
    return render_template('pway/index2.html', projects=upload_list)


@pway.route('/demo')
@login_required
def demo():
    return render_template('pway/show_pathways_demo.html')


@pway.route('/restart')
@login_required
def reset():
    role_names = [r.name for r in current_user.roles]
    if 'vip' not in role_names:
        abort(404)
    if request.environ['mod_wsgi.process_group'] != '':
        os.kill(os.getpid(), signal.SIGINT)
        flash("Restarted.", "info")
        return redirect(url_for('.index'))


@pway.route('/scatter')
@login_required
@no_ssl
def scatter():
    # if proj among arguments, show this tree first.
    show_proj = request.args.get('proj', None)
    include = request.args.get('include', None)
    # show_proj = '57'
    # list of projects (and proj_names) used to create dropdown project selector
    upload_list = UserFile.query.filter_by(user_id=current_user.id).\
        filter_by(run_complete=True).order_by(UserFile.file_id).all()
    # proj_names = {int(i.file_id): i.get_local_filename() for i in upload_list}
    resources = Resources(mode="cdn")
    if upload_list:
        # Use specified project from args or highest file_id as CURRENT PROJECT
        current_proj = upload_list[-1]  # override if valid proj specified
        if show_proj:
            current_temp = [u for u in upload_list
                            if u.file_id == int(show_proj)]
            # if not among user's finished projects, use highest file_id
            if len(current_temp) == 1:
                current_proj = current_temp[0]
            else:
                current_proj = upload_list[-1]
        detail_path = naming_rules.get_detailed_path(current_proj)
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
        script, div = components(plot1, resources)
        js_name = naming_rules.get_js_name(current_proj)
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
        proj_dir = naming_rules.get_project_folder(current_proj)
        if os.path.exists(os.path.join(proj_dir, 'matrix_svg_cnv')):
            has_cnv = True
        else:
            has_cnv = False

    else:  # no projects yet!
        flash("No project results to show yet.", "warning")
        return redirect(url_for('.index'))

    return render_template('pway/scatter.html', current_proj=current_proj,
                           projects=upload_list, js_name=js_name,
                           js_inds=js_inds, plot_inds=plot_inds,
                           has_cnv=has_cnv, user_id=current_user.id,
                           bokeh_script=script, bokeh_div=div,
                           include_genes=include, resources=resources)


@pway.route('/compare')
@login_required
@no_ssl
def compare():
    # if proj among arguments, show this tree first.
    proj_a = request.args.get('proj_a', None)
    proj_b = request.args.get('proj_b', None)
    include = request.args.get('include', None)
    resources = Resources(mode="cdn")
    show_logged = False

    # list of projects (and proj_names) used to create dropdown project selector
    upload_list = UserFile.query.filter_by(user_id=current_user.id).\
        filter_by(run_complete=True).order_by(UserFile.file_id).all()

    if len(upload_list) > 1:
        # Use specified project from args or highest file_id as CURRENT PROJECT
        current_proj = upload_list[-1]  # override if valid proj specified
        if proj_a and proj_b:
            current_temp_a = [u for u in upload_list if u.file_id == int(proj_a)]
            current_temp_b = [u for u in upload_list if u.file_id == int(proj_b)]
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
        detail_path1 = naming_rules.get_detailed_path(current_proj_a)
        detail_path2 = naming_rules.get_detailed_path(current_proj_b)

        js_name1 = naming_rules.get_js_name(current_proj_a)
        js_name2 = naming_rules.get_js_name(current_proj_b)

        # load pathways with 1+ mutation in 1+ patients,
        # ignoring ones with 'cancer' etc in name
        all_paths1 = load_pathway_list_from_file(detail_path1)
        all_paths2 = load_pathway_list_from_file(detail_path2)
        # all_paths1 = [i for i in all_paths1 if i.gene_set]
        # all_paths2 = [i for i in all_paths2 if i.gene_set]

        sig_pids1 = [i.path_id for i in all_paths1 if i.gene_set]  # sig path ids in proj1
        sig_pids2 = [i.path_id for i in all_paths2 if i.gene_set]  # sig path ids in proj2
        all_pids1 = [i.path_id for i in all_paths1]  # all path ids in proj1
        all_pids2 = [i.path_id for i in all_paths2]  # all path ids in proj1
        use_pids1 = set.intersection(set(sig_pids1), set(all_pids2))
        use_pids2 = set.intersection(set(sig_pids2), set(all_pids1))
        good_paths = set.union(set(use_pids1), set(use_pids2))  # sig pids contained in both projects

        # get use_paths 1&2: pathway objects in proj1 order
        use_paths1 = [i for i in all_paths1 if i.path_id in good_paths]  # sig paths (either proj) in proj1
        use_path_ids = [i.path_id for i in use_paths1]  # sig path ids (either proj) ordered by proj1 index
        use_paths2 = list()
        for pid in use_path_ids:
            use_paths2.extend([i for i in all_paths2 if i.path_id == pid])  # paths from proj2

        # save indices of chosen pathways -- these are the indices in js file
        inds1 = list()
        inds2 = list()  # [0, 1, 4, 11, -1, 7, 10, 5, 6, 319, ...]
        for i in use_path_ids:
            try:
                inds1.append(sig_pids1.index(i))
            except ValueError:
                inds1.append(-1)
            try:
                inds2.append(sig_pids2.index(i))
            except ValueError:
                inds2.append(-1)

        if show_logged:
            effects1 = [np.log10(float(i.n_effective) / i.n_actual) for i in use_paths1]
            effects2 = [np.log10(float(i.n_effective) / i.n_actual) for i in use_paths2]
            xlabel = "Log10 effect size ({})".format(current_proj_a.proj_suffix)
            ylabel = "Log10 effect size ({})".format(current_proj_b.proj_suffix)
        else:
            effects1 = [float(i.n_effective) / i.n_actual for i in use_paths1]
            effects2 = [float(i.n_effective) / i.n_actual for i in use_paths2]
            xlabel = "Effect size ({})".format(current_proj_a.proj_suffix)
            ylabel = "Effect size ({})".format(current_proj_b.proj_suffix)
        q1 = [float(i.p_value)*g.n_pathways for i in use_paths1]
        q2 = [float(i.p_value)*g.n_pathways for i in use_paths2]
        pnames = [misc.strip_contributors(i.nice_name) for i in use_paths1]
        source = ColumnDataSource(data={'effects1': effects1,
                                        'effects2': effects2,
                                        'q1': q1, 'q2': q2,
                                        'pname': pnames})
        # SET UP FIGURE
        minx = min(effects1)
        minx *= 1 - minx/abs(minx)*0.2
        miny = min(effects2)
        miny *= 1 - miny/abs(miny)*0.2
        maxx = max(effects1)*1.2
        maxy = max(effects2)*1.2
        tools = "crosshair,pan,wheel_zoom,box_zoom,reset,tap," \
                "box_select,hover"  # poly_select,lasso_select, previewsave
        if show_logged:
            plot = figure(tools=tools, plot_height=400, plot_width=600,
                          title=None, logo=None, toolbar_location="above",
                          x_axis_label=xlabel,
                          y_axis_label=ylabel,
                          x_range=[0, maxx], y_range=[0, maxy])
        else:
            plot = figure(tools=tools, plot_height=400, plot_width=600,
                          title=None, logo=None, toolbar_location="above",
                          x_axis_label=xlabel, y_axis_label=ylabel,
                          x_range=[minx, maxx], y_range=[miny, maxy],
                          x_axis_type="log", y_axis_type="log")
        plot.min_border_top = 15
        plot.min_border_right = 15

        plot.xaxis.axis_label_text_font_size = "12pt"
        plot.yaxis.axis_label_text_font_size = "12pt"
        plot.line([1, 1], [miny, maxy], line_width=2, color="blue", alpha=1, line_dash=[6, 6])
        plot.line([minx, maxx], [1, 1], line_width=2, color="blue", alpha=1, line_dash=[6, 6])

        # radius=radii, fill_color=colors, fill_alpha=0.6, line_color=None
        plot.scatter("effects1", "effects2", source=source, size=10,
                     color="red", alpha=0.1, marker="circle",
                     line_color="firebrick", line_alpha=0.5)
        hover = plot.select(dict(type=HoverTool))
        hover.tooltips = OrderedDict([
            ("name", "@pname")
        ])
        renderers = [i for i in plot.renderers
                     if type(i) == GlyphRenderer and type(i._glyph) == Circle]
        hover[0].renderers.extend(renderers)

        # ADD TEXT INPUT OR SLIDER
        callback = Callback(args=dict(source=source), code="""
            // get old selection indices, if any
            var model = getDataModel();
            prv_selected = model.attributes.selected['1d'].indices;
            prv_select_full = []
            for(var i=0; i<prv_selected.length; i++){
                prv_select_full.push(scatter_array[prv_selected[i]])
            }
            new_selected = []

            var data = source.get('data');
            var q_val = q_widget.get('value');

            var n_total = fullset['effects1'].length;
            var e1_use = data['effects1'];
            var e2_use = data['effects2'];
            var n_use = data['pname'];
            e1_use.length = 0;
            e2_use.length = 0;
            n_use.length = 0;
            scatter_array = []
            var j = -1;  // new glyph indices
            for (i = 0; i < n_total; i++) {
                this_q1 = fullset['q1'][i];
                this_q2 = fullset['q2'][i];
                if(this_q1 <= q_val || this_q2 <= q_val){
                    j++;
                    e1_use.push(fullset['effects1'][i]);
                    e2_use.push(fullset['effects2'][i]);
                    n_use.push(fullset['pname'][i]);
                    scatter_array.push(i)
                    if($.inArray(i, prv_select_full) > -1){
                        new_selected.push(j);
                    }
                }
            }
            source.trigger('change');
            model.attributes.selected['1d'].indices = new_selected;
            model.trigger('change');
            updateIfSelectionChange_afterWait();
            """)
        text_input = TextInput(value='',
                               title="Q cutoff:",
                               callback=callback)
        callback.args["q_widget"] = text_input
        callback.args["fullset"] = source.clone().data
        layout = hplot(plot, vform(text_input))
        script, div = components(layout, resources)

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
                           user_id=current_user.id, bokeh_script=script,
                           bokeh_div=div, include_genes=include,
                           resources=resources)


@pway.route('/tree')
@login_required
def tree():
    # if proj among arguments, show this tree first.
    show_proj = request.args.get('proj', None)
    # list of projects (and proj_names) used to create dropdown project selector
    upload_list = UserFile.query.filter_by(user_id=current_user.id).\
        filter_by(run_complete=True).order_by(UserFile.file_id).all()
    proj_names = {int(i.file_id): i.get_local_filename() for i in upload_list}
    if upload_list:
        # Use specified project from args or highest file_id as CURRENT PROJECT
        current_proj = upload_list[-1]  # override if valid proj specified
        if show_proj:
            current_temp = [u for u in upload_list if u.file_id == int(show_proj)]
            # if not among user's finished projects, use highest file_id
            if len(current_temp) == 1:
                current_proj = current_temp[0]
            else:
                current_proj = upload_list[-1]
        # load ordered dictionary of path_ids : pathway_display_name
        names_path, tree_path = naming_rules.get_tree_score_paths(current_proj)[1:3]
        names_ordered_path = names_path + '.reorder'
        tree_path = naming_rules.get_apache_path(tree_path) + 'z'
        names_odict = OrderedDict()  # ordered dictionary of path_id: name
        with open(names_ordered_path, 'rU') as f:
            for line in f:
                vals = line.strip('\n').split('\t')
                if len(vals) != 2:
                    continue
                names_odict[vals[0]] = vals[1]
    else:  # no projects yet!
        flash("No project results to show yet.", "warning")
        return redirect(url_for('.index'))
    return render_template('pway/tree.html', current_proj=current_proj,
                           projects=upload_list, user_id=current_user.id,
                           proj_names=proj_names, names_odict=names_odict,
                           tree_path=tree_path)


@pway.route('/faq')
def faq():
    return render_template('pway/faq.html')


@pway.route('/archive/<int:proj>')
@login_required
def archive(proj):
    upload_obj = UserFile.query.\
        filter_by(user_id=current_user.id, file_id=proj).\
        first_or_404()
    zip_path = zip_project(upload_obj)
    filename = os.path.basename(zip_path)
    return send_file(zip_path, mimetype='application/zip',
                     as_attachment=True, attachment_filename=filename)


@pway.route('/demofile')
def demo_file():
    # return render_template('pway/show_pathways_template.html')
    return send_from_directory(current_app.config['DATA_ROOT'],
                               'skcm_ns_500.txt', as_attachment=True)


@pway.route('/filtered')
def get_filtered():

    proj_id = request.args.get('projId', None)
    filter_type = request.args.get('type', None)  # ignored or rejected
    if not proj_id.isdigit() or filter_type not in ('ignored', 'rejected'):
        abort(404)
    proj_id = int(proj_id)
    upload_obj = UserFile.query.filter_by(user_id=current_user.id).\
        filter_by(file_id=proj_id).all()[0]
    if not upload_obj:
        abort(404)
    if filter_type == 'ignored':
        file_path = naming_rules.get_unused_gene_path(upload_obj)
    else:
        file_path = naming_rules.get_rejected_gene_path(upload_obj)

    try:
        return send_file(file_path, as_attachment=True)
    except IOError:
        abort(404)


@pway.route('/results')
@login_required
def results():
    show_proj = request.args.get('proj', None)
    include = request.args.get('include', None)
    upload_list = UserFile.query.filter_by(user_id=current_user.id).\
        filter_by(run_complete=True).all()
    if not upload_list:
        flash("No project results to show yet.", "warning")
        return redirect(url_for('.index'))
    if show_proj is None and len(upload_list):
        show_proj = upload_list[-1].get_local_filename()
    proj_names = {int(i.file_id): i.get_local_filename() for i in upload_list}
    return render_template('pway/show_pathways_template.html',
                           projects=upload_list, user_id=current_user.id,
                           show_proj=int(show_proj), proj_names=proj_names,
                           include_genes=include)


@pway.route('/upload', methods=('GET', 'POST'))
@limit_user_uploads
def upload():
    """http://flask.pocoo.org/docs/0.10/patterns/fileuploads/"""

    form = UploadForm()
    if form.validate_on_submit():
        try:
            mut_file = MutationFile(form.mut_file.data)
        except ValidationError as e:
            flash(str(e), 'danger')
            return render_template('pway/upload.html', form=form)

        # CREATE NEW USER IF UNAUTHENTICATED (HERE, FILE IS VALID)
        if not current_user.is_authenticated():
            # create guest user
            temp_user, temp_pswd = create_anonymous_user()
            flash('Your temporary username is {} and password is {}. '.format(
                temp_user.email, temp_pswd) + 'This account will be deleted in '
                '{} days.'.format(current_app.config['ANONYMOUS_MAX_AGE_DAYS']),
                'info')
            login_user(temp_user, force=True, remember=True)

        # CREATE USERFILE FROM FORM
        mut_filename = secure_filename(form.mut_file.data.filename)
        user_upload = UserFile(filename=mut_filename, user_id=current_user.id)
        form.to_model(user_upload)
        # CREATE PROJECT FOLDER, TWEAK USERFILE, MOVE MUTATIONS
        out = initialize_project(user_upload=user_upload, mut_file=mut_file)
        user_upload, proj_folder, file_path = out

        flash('File accepted and validated. Analysis in progress.',
              'success')
        # run analysis asynchronously
        # http://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xi-email-support
        run_analysis(proj_folder, file_path, user_upload.file_id)
        return redirect(url_for('.index'))

    # MESSAGES FOR INITIAL UPLOAD PAGE ACCESS OR FAILED UPLOAD
    if current_user.is_authenticated():
        message = "You've run {} projects in the last week. ".format(g.n_week)
        message += "Your weekly limit is {}.".format(g.n_week_max)
        if g.incomplete:
            message += "\nYou have {} jobs still running."\
                .format(len(g.incomplete))
        flash(message, "info")
    else:
        flash("You can upload here as a new guest user, but with severe "
              "usage limits and no email alerts. Consider registering instead.",
              "danger")
    return render_template('pway/upload.html', form=form)
