from flask import render_template, flash, redirect, url_for, abort,\
    request, current_app, send_file, send_from_directory
from flask_login import login_required, current_user
from werkzeug.utils import secure_filename
import os
from datetime import datetime, timedelta
from collections import OrderedDict
import numpy as np
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.plotting import ColumnDataSource
from bokeh.models.tools import HoverTool

from . import pway, FileTester
from .forms import UploadForm
from ..models import UserFile
from .. import db
from ..get_effective_pathways import run_analysis, load_pathway_list_from_file
from ..admin import get_project_folder, get_user_folder, zip_project, \
    delete_project_folder
from .. import naming_rules
from .. import misc


@pway.route('/')
@login_required
def index():

    upload_list = UserFile.query.filter_by(user_id=current_user.id).all()
    return render_template('pway/index2.html', projects=upload_list)


@pway.route('/demo')
@login_required
def demo():
    return render_template('pway/show_pathways_demo.html')


@pway.route('/tree_test')
@login_required
def tree_test():
    return render_template('pway/tree_test.html')


@pway.route('/scatter')
@login_required
def scatter():
    # if proj among arguments, show this tree first.
    show_proj = request.args.get('proj', None)
    # show_proj = '57'
    # list of projects (and proj_names) used to create dropdown project selector
    upload_list = UserFile.query.filter_by(user_id=current_user.id).\
        filter_by(run_complete=True).order_by(UserFile.file_id).all()
    # proj_names = {int(i.file_id): i.get_local_filename() for i in upload_list}
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

        detail_path = naming_rules.get_detailed_path(current_proj)
        allPathways = load_pathway_list_from_file(detail_path)

        data_pways, data_pvals, data_effect = [], [], []
        for p in allPathways:
            if p.n_effective >= p.n_actual:
                pass
            pval = float(p.p_value)
            if pval >= 0.05:
                continue
            data_pvals.append(float(p.p_value))
            data_effect.append(float(p.n_effective) / p.n_actual)
            data_pways.append(p)

        x = np.log2(np.array(data_effect))  # effect size
        y = -np.log10(np.array(data_pvals))  # p-value
        # adjust zero pvalues. e-15.9 seems to be minimum.
        max_y = max([np.ceil(max([i for i in x if i != np.inf])),
                     np.float64(17)])
        y[y == np.inf] = max_y
        pnames = [misc.strip_contributors(p.name) for p in data_pways]
        source = ColumnDataSource(data={'x': x, 'y': y, 'pname': pnames})
        tools = "resize,crosshair,pan,wheel_zoom,box_zoom,reset,tap," \
                "box_select,hover"  # poly_select,lasso_select, previewsave
        plot = figure(tools=tools, plot_height=400, plot_width=600,
                      title='Effect size vs p-value', logo=None,
                      x_axis_label="log2 fold change",
                      y_axis_label="-log10 p-value",)
        plot.scatter("x", "y", source=source, size=10, color="red", alpha=0.1,
                     marker="circle", line_color="firebrick", line_alpha=0.5)
        hover = plot.select({'type': HoverTool})
        hover.tooltips = OrderedDict([
            ("name", "@pname")  # optional: ("index", "$index")
        ])
        script, div = components(plot, CDN)

    else:  # no projects yet!
        flash("No project results to show yet.", "info")
        current_proj = None
        script, div = None, None

    return render_template('pway/scatter.html', current_proj=current_proj,
                           projects=upload_list,
                           user_id=current_user.id, bokeh_script=script,
                           bokeh_div=div)


@pway.route('/compare')
@login_required
def compare():
    # if proj among arguments, show this tree first.
    proj_a = request.args.get('proj_a', None)
    proj_b = request.args.get('proj_b', None)

    # list of projects (and proj_names) used to create dropdown project selector
    upload_list = UserFile.query.filter_by(user_id=current_user.id).\
        filter_by(run_complete=True).order_by(UserFile.file_id).all()

    if len(upload_list) > 2:
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
        all_paths1 = [i for i in all_paths1 if i.gene_set]
        all_paths2 = [i for i in all_paths2 if i.gene_set]

        # temporary variables
        all_ids1 = [p.path_id for p in all_paths1]
        all_ids2 = [p.path_id for p in all_paths2]
        common_ids = set.intersection(set(all_ids1), set(all_ids2))
        ids_sorted = [i for i in all_ids1 if i in common_ids]  # uses project A
        inds1 = [all_ids1.index(i) for i in ids_sorted]  # indices into all_paths1
        inds2 = [all_ids2.index(i) for i in ids_sorted]  # indices into all_paths2

        # key data for export
        # inds_del1 = [i for i in xrange(len(all_paths1)) if i not in inds1]
        # inds_del2 = [i for i in xrange(len(all_paths2)) if i not in inds2]
        effects1 = [np.log10(float(all_paths1[i].ne_low)
                             / all_paths1[i].n_actual) for i in inds1]
        effects2 = [np.log10(float(all_paths2[i].ne_low)
                             / all_paths2[i].n_actual) for i in inds2]
        pnames = [misc.strip_contributors(all_paths1[i].name) for i in inds1]
        source = ColumnDataSource(data={'x': effects1, 'y': effects2,
                                        'pname': pnames})
        maxx = max(effects1)*1.1
        maxy = max(effects2)*1.1

        tools = "resize,crosshair,pan,wheel_zoom,box_zoom,reset,tap," \
                "box_select,hover"  # poly_select,lasso_select, previewsave
        xlabel = "Log10 effect size ({})".format(current_proj_a.proj_suffix)
        ylabel = "Log10 effect size ({})".format(current_proj_b.proj_suffix)
        plot = figure(tools=tools, plot_height=400, plot_width=600, title=None,
                      logo=None, toolbar_location="right",
                      x_axis_label=xlabel,
                      y_axis_label=ylabel,
                      x_range=[0, maxx], y_range=[0, maxy])

        # radius=radii, fill_color=colors, fill_alpha=0.6, line_color=None
        plot.scatter("x", "y", source=source, size=10, color="red", alpha=0.1,
                  marker="circle", line_color="firebrick", line_alpha=0.5)
        hover = plot.select(dict(type=HoverTool))
        hover.tooltips = OrderedDict([
            ("name", "@pname")
        ])

        script, div = components(plot, CDN)

    else:  # no projects yet!
        flash("No project results to show yet.", "info")
        current_proj_a = None
        current_proj_b = None
        script, div = None, None
        inds_del1 = None
        inds_del2 = None

    return render_template('pway/compare.html',
                           current_projs=[current_proj_a, current_proj_b],
                           inds_use=[inds1, inds2],
                           js_name_a=js_name1,
                           js_name_b=js_name2,
                           projects=upload_list,
                           user_id=current_user.id, bokeh_script=script,
                           bokeh_div=div)


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
        current_proj = upload_list[0]  # override if valid proj specified
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
        flash("No project results to show yet.", "info")
        current_proj = None
        names_odict = OrderedDict()
        tree_path = None
    return render_template('pway/tree.html', current_proj=current_proj,
                           projects=upload_list, user_id=current_user.id,
                           proj_names=proj_names, names_odict=names_odict,
                           tree_path=tree_path)


@pway.route('/faq')
@login_required
def faq():
    return render_template('pway/faq.html')


@pway.route('/archive')
@login_required
def archive():
    proj_id = int(request.args.get('proj', 0))
    upload_obj = None
    if proj_id:
        upload_obj = UserFile.query.get(proj_id)
    if not upload_obj or current_user.id != upload_obj.user_id:
        abort(404)
    zip_path = zip_project(upload_obj)
    filename = os.path.basename(zip_path)
    return send_file(zip_path, mimetype='application/zip',
                     as_attachment=True, attachment_filename=filename)

@pway.route('/demofile')
@login_required
def demo_file():
    # return render_template('pway/show_pathways_template.html')
    return send_from_directory(current_app.config['UPLOAD_FOLDER'],
                               'skcm_ns_500.txt', as_attachment=True)


@pway.route('/results')
@login_required
def results():
    show_proj = request.args.get('proj', None)
    upload_list = UserFile.query.filter_by(user_id=current_user.id).filter_by(run_complete=True).all()
    proj_names = {int(i.file_id): i.get_local_filename() for i in upload_list}
    if not(upload_list):
        flash("No project results to show yet.", "info")
    return render_template('pway/show_pathways_template.html',
                           projects=upload_list, user_id=current_user.id,
                           show_proj=show_proj, proj_names=proj_names)


@pway.route('/upload', methods=('GET', 'POST'))
@login_required
def upload():
    """http://flask.pocoo.org/docs/0.10/patterns/fileuploads/"""
    # only accept this file if user has no running jobs.
    user_id = current_user.id
    incomplete = UserFile.query.filter_by(user_id=user_id)\
        .filter_by(run_complete=0).all()
    role_names = [r.name for r in current_user.roles]
    if 'townsend' not in role_names and incomplete:
        flash("Sorry, you must wait until your currently running projects "
              "have finished.", "danger")
        return index()

    # ENFORCE WEEKLY LIMIT AND DISPLAY INFO MESSAGE
    week_ago = datetime.now() - timedelta(days=7)
    week_complete = UserFile.query.filter_by(user_id=user_id)\
        .filter_by(run_complete=True).filter(UserFile.upload_time > week_ago)\
        .all()
    n_week = len(week_complete)
    n_week_max = int(max([r.uploads_pw for r in current_user.roles]))
    # if len(week_complete>9):
    if n_week >= n_week_max:
        first_time = min([p.upload_time for p in week_complete])
        wait_time = first_time + timedelta(days=7) - datetime.now()
        diff_str = '{:.1f}'.format(wait_time.seconds/(60*60.))

        flash("Sorry, you've used your allotted runs for this week. You can "
              "try again in {} hours.".format(diff_str), "danger")
        return index()

    form = UploadForm()
    if form.validate_on_submit():
        # upload = Upload(uploader=current_user)
        mut_filename = secure_filename(form.mut_file.data.filename)

        # will create folder for job: <run_id>
        user_upload = UserFile(filename=mut_filename, user_id=user_id)
        form.to_model(user_upload)
        db.session.add(user_upload)
        db.session.commit()

        user_folder = get_user_folder(user_id)
        proj_folder = get_project_folder(user_upload)
        if not os.path.exists(user_folder):
            os.mkdir(user_folder)
        os.mkdir(proj_folder)
        file_path = os.path.join(proj_folder,
                                 user_upload.get_local_filename())
        form.mut_file.data.save(file_path)

        # FILE SHOULD BE IN UPLOADS FOLDER NOW.
        # RESPOND BASED ON FILE VALIDITY:
        #     VALID: START RUN, FLASH 'EXPECT EMAIL' on STATUS PAGE
        #     INVALID: FLASH 'BAD FILE', LIST ERRORS ON UPLOAD PAGE
        file_tester = FileTester(file_path)

        if file_tester.good_headers:
            if file_tester.good_headers and file_tester.data_present:
                user_upload.is_valid = True
                user_upload.run_complete = False
                upload_id = user_upload.file_id
                db.session.add(user_upload)
                db.session.commit()
                flash('File accepted and validated. Analysis in progress.',
                      'success')
                # want analysis to run asynchronously!
                # http://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xi-email-support
                run_analysis(proj_folder, file_path, upload_id)
                return redirect(url_for('.index'))
            else:  # good headers but no data
                flash('Your file seems to be missing data. Please try again.',
                      'danger')
        else:  # bad headers
            flash("Your headers don't look right. Expected headers are:\n{}"
                  .format('\t'.join(FileTester.want_headers)), 'danger')
            delete_project_folder(user_upload)
            db.session.delete(user_upload)
            db.session.commit()
    else:
        message = "You've run {} projects in the last week. ".format(n_week)
        message += "Your weekly limit is {}.".format(n_week_max)
        if incomplete:
            message += "\nYou have {} jobs still running."\
                .format(len(incomplete))
        flash(message, "info")
    return render_template('pway/upload.html', form=form)


"""

orig defaults:
 - output files saved to cwd
 - txt files: pway_pvalues_<PROJ_NAME>_<SUFFIX>(_pretty).txt
 - directories:
    - pathways_svg
    - matrix_txt
        - matrix_svg

UPLOADS FOLDER
==============

F 1_testfile.txt
F 2_melanoma_mk2014.txt
...
D output_1/
D output_2/
...


PROJECT RUN FOLDER
==================

 - create folder for project if text file is valid and user allowance not full

D pathways_svg/
D matrix_svg/
D matrix_txt/
F FILE_ID.txt
F FILE_ID_detail.txt
F FILE_ID.js

 - txt file name/path (default pathways_pvalues_PROJSTR) determined by
  GenericPathwayFileProcessor

root_svg_target = "pathways_svg/"
root_svg_matrix = "matrix_txt/matrix_svg/"
EXPECTED BY HTML: root_svg_matrix2 = "matrix_svg/"
proj_svg = 'proj_' + project_str + '/'

"""