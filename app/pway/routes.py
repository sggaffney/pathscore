from flask import render_template, flash, redirect, url_for, abort,\
    request, current_app, send_file, send_from_directory
from flask_login import login_required, current_user, login_user
from werkzeug.utils import secure_filename
import os
import signal
from datetime import datetime, timedelta
from collections import OrderedDict
import numpy as np
from bokeh.plotting import figure, gridplot, ColumnDataSource
from bokeh.resources import Resources
from bokeh.embed import components
from bokeh.models.tools import HoverTool
from bokeh.models.renderers import GlyphRenderer
from bokeh.models.markers import Circle

from . import pway, FileTester, TempFile
from .forms import UploadForm
from ..models import UserFile, User, Role
from .. import db
from ..get_effective_pathways import run_analysis, load_pathway_list_from_file
from ..admin import zip_project, \
    delete_project_folder
from .. import naming_rules
from .. import misc


@pway.route('/')
@login_required
def index():
    current_app.logger.info('Loaded home page, info.')
    current_app.logger.warn('Loaded home page, warn.')
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
            current_temp = [u for u in upload_list if u.file_id == int(show_proj)]
            # if not among user's finished projects, use highest file_id
            if len(current_temp) == 1:
                current_proj = current_temp[0]
            else:
                current_proj = upload_list[-1]
        detail_path = naming_rules.get_detailed_path(current_proj)
        all_pathways = load_pathway_list_from_file(detail_path)
        data_pways, data_pvals, data_effect, data_D = [], [], [], []
        for p in all_pathways:
            pval = float(p.p_value)
            if pval >= 0.05:
                continue
            data_pvals.append(pval)
            if p.n_actual < p.n_effective:
                data_effect.append(float(p.ne_low) / p.n_actual)
            else:  # undermutated
                if p.ne_high < p.n_actual:  # high-effect less than n_actual
                    data_effect.append(float(p.ne_high) / p.n_actual)
                else:
                    data_effect.append(float(p.n_effective) / p.n_actual)
            data_pways.append(p)
            data_D.append(p.D)
        x = np.log2(np.array(data_effect))  # effect size
        y = -np.log10(np.array(data_pvals))  # p-value
        D = np.array(data_D)  # alternatively: np.log2...
        # adjust zero pvalues. e-15.9 seems to be minimum.
        max_y = max([np.ceil(max([i for i in x if i != np.inf])),
                     np.float64(17)])
        y[y == np.inf] = max_y
        pnames = [misc.strip_contributors(p.nice_name) for p in data_pways]
        xyvalues = ColumnDataSource({'effect': x,
                                     'pvals': y,
                                     'D': D,
                                     'pname': pnames})
        tools = "resize,crosshair,pan,wheel_zoom,box_zoom,reset,tap," \
                "box_select,hover"  # poly_select,lasso_select, previewsave
        plot_config = dict(plot_height=400, plot_width=400, logo=None, tools=tools)

        plot1 = figure(title='Effect size vs p-value',
                       x_axis_label="log2 fold change",
                       y_axis_label="-log10 p-value",
                       **plot_config)
        plot1.scatter('effect', 'pvals', source=xyvalues, size=10, color="red", alpha=0.1,
                     marker="circle", line_color="firebrick", line_alpha=0.5)

        plot2 = figure(title='Effect size vs deviance',
                       x_axis_label="log2 fold change",
                       y_axis_label="D",
                       **plot_config)
        plot2.scatter('effect', 'D', source=xyvalues, size=10, color="red", alpha=0.1,
                      marker="circle", line_color="firebrick", line_alpha=0.5)

        for plot in [plot1, plot2]:
            hover = plot.select(dict(type=HoverTool))
            hover.tooltips = OrderedDict([
            #     ("index", "$index"),
                ("name", "@pname")
            ])
        plot = gridplot([[plot1, plot2]])
        script, div = components(plot, resources)
        js_name = naming_rules.get_js_name(current_proj)
        # IDS
        all_ids = [p.path_id for p in all_pathways]
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
                           js_inds=js_inds, plot_inds=plot_inds, has_cnv=has_cnv,
                           user_id=current_user.id, bokeh_script=script,
                           bokeh_div=div, include_genes=include,
                           resources=resources)


@pway.route('/compare')
@login_required
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

        use_paths1 = [i for i in all_paths1 if i.path_id in good_paths]  # sig paths (either proj) in proj1
        use_path_ids = [i.path_id for i in use_paths1]  # sig path ids (either proj) ordered by proj1 index

        use_paths2 = list()
        for pid in use_path_ids:
            use_paths2.extend([i for i in all_paths2 if i.path_id == pid])  # paths from proj2

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

        # key data for export
        # inds_del1 = [i for i in xrange(len(all_paths1)) if i not in inds1]
        # inds_del2 = [i for i in xrange(len(all_paths2)) if i not in inds2]
        if show_logged:
            effects1 = [np.log10(float(i.ne_low) / i.n_actual) for i in use_paths1]
            effects2 = [np.log10(float(i.ne_low) / i.n_actual) for i in use_paths2]
            xlabel = "Log10 effect size ({})".format(current_proj_a.proj_suffix)
            ylabel = "Log10 effect size ({})".format(current_proj_b.proj_suffix)
        else:
            effects1 = [float(i.ne_low) / i.n_actual for i in use_paths1]
            effects2 = [float(i.ne_low) / i.n_actual for i in use_paths2]
            xlabel = "Effect size ({})".format(current_proj_a.proj_suffix)
            ylabel = "Effect size ({})".format(current_proj_b.proj_suffix)
        pnames = [misc.strip_contributors(i.nice_name) for i in use_paths1]
        source = ColumnDataSource(data={'x': effects1, 'y': effects2,
                                        'pname': pnames})
        minx = min(effects1)
        minx *= 1 - minx/abs(minx)*0.2
        miny = min(effects2)
        miny *= 1 - miny/abs(miny)*0.2
        maxx = max(effects1)*1.2
        maxy = max(effects2)*1.2

        tools = "resize,crosshair,pan,wheel_zoom,box_zoom,reset,tap," \
                "box_select,hover"  # poly_select,lasso_select, previewsave

        if show_logged:
            plot = figure(tools=tools, plot_height=400, plot_width=600, title=None,
                          logo=None, toolbar_location="right",
                          x_axis_label=xlabel,
                          y_axis_label=ylabel,
                          x_range=[0, maxx], y_range=[0, maxy])
        else:
            plot = figure(tools=tools, plot_height=400, plot_width=600, title=None,
                          logo=None, toolbar_location="right",
                          x_axis_label=xlabel, y_axis_label=ylabel,
                          x_range=[minx, maxx], y_range=[miny, maxy],
                          x_axis_type="log", y_axis_type="log")

        plot.line([1,1], [miny, maxy], line_width=2, color="blue", alpha=1, line_dash=[6, 6])
        plot.line([minx, maxx], [1, 1], line_width=2, color="blue", alpha=1, line_dash=[6, 6])

        # radius=radii, fill_color=colors, fill_alpha=0.6, line_color=None
        plot.scatter("x", "y", source=source, size=10, color="red", alpha=0.1,
                  marker="circle", line_color="firebrick", line_alpha=0.5)
        hover = plot.select(dict(type=HoverTool))
        hover.tooltips = OrderedDict([
            ("name", "@pname")
        ])
        renderers = [i for i in plot.renderers
                     if type(i) == GlyphRenderer and type(i._glyph) == Circle]
        hover[0].renderers.extend(renderers)

        script, div = components(plot, resources)

        proj_dir_a = naming_rules.get_project_folder(current_proj_a)
        proj_dir_b = naming_rules.get_project_folder(current_proj_b)
        if os.path.exists(os.path.join(proj_dir_a, 'matrix_svg_cnv')) and \
                os.path.exists(os.path.join(proj_dir_b, 'matrix_svg_cnv')):
            has_cnv = True
        else:
            has_cnv = False

    else:  # not enough projects yet!
        flash("Two completed projects are required for a comparison.", "warning")
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
        flash("No project results to show yet.", "warning")
        return redirect(url_for('.index'))
    return render_template('pway/tree.html', current_proj=current_proj,
                           projects=upload_list, user_id=current_user.id,
                           proj_names=proj_names, names_odict=names_odict,
                           tree_path=tree_path)


@pway.route('/faq')
@login_required
def faq():
    current_app.logger.info('Loaded faq page, info.')
    current_app.logger.warn('Loaded faq page, warn.')
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
    return send_from_directory(current_app.config['DATA_ROOT'],
                               'skcm_ns_500.txt', as_attachment=True)


@pway.route('/results')
@login_required
def results():
    show_proj = request.args.get('proj', None)
    include = request.args.get('include', None)
    upload_list = UserFile.query.filter_by(user_id=current_user.id).filter_by(run_complete=True).all()
    proj_names = {int(i.file_id): i.get_local_filename() for i in upload_list}
    if not(upload_list):
        flash("No project results to show yet.", "warning")
        return redirect(url_for('.index'))
    return render_template('pway/show_pathways_template.html',
                           projects=upload_list, user_id=current_user.id,
                           show_proj=show_proj, proj_names=proj_names,
                           include_genes=include)


@pway.route('/upload', methods=('GET', 'POST'))
def upload():
    """http://flask.pocoo.org/docs/0.10/patterns/fileuploads/"""

    # CHECK IF RUNNING PROJECT COUNT IS WITHIN USER LIMITS
    if current_user.is_authenticated():
        incomplete = UserFile.query.filter_by(user_id=current_user.id)\
            .filter_by(run_complete=0).all()
        role_names = [r.name for r in current_user.roles]
        if 'vip' not in role_names and incomplete:
            flash("Sorry, you must wait until your currently running projects "
                  "have finished.", "danger")
            return index()

        # ENFORCE WEEKLY LIMIT AND DISPLAY INFO MESSAGE
        week_ago = datetime.utcnow() - timedelta(days=7)
        week_complete = UserFile.query.filter_by(user_id=current_user.id)\
            .filter_by(run_complete=True).filter(UserFile.upload_time > week_ago)\
            .all()
        n_week = len(week_complete)
        n_week_max = int(max([r.uploads_pw for r in current_user.roles]))
        # if len(week_complete>9):
        if n_week >= n_week_max:
            first_time = min([p.upload_time for p in week_complete])
            wait_time = first_time + timedelta(days=7) - datetime.utcnow()
            wait_str = misc.get_wait_time_string(wait_time)
            flash("Sorry, you've used your allotted runs for this week. "
                  "You can try again in {}.".format(wait_str), "danger")
            return index()

    form = UploadForm()
    if form.validate_on_submit():
        # upload = Upload(uploader=current_user)
        mut_filename = secure_filename(form.mut_file.data.filename)
        temp_file = TempFile(form.mut_file.data)
        file_tester = FileTester(temp_file.path)

        # ABORT IF FILE HAS ISSUES
        if file_tester.data_issues:
            for issue in file_tester.data_issues:
                flash(issue, 'danger')
            return render_template('pway/upload.html', form=form)

        # CREATE NEW USER IF UNAUTHENTICATED
        if not current_user.is_authenticated():
            # create guest user
            temp_pswd = User.generate_random_password()
            temp_user = User(password_raw=temp_pswd)
            db.session.add(temp_user)
            db.session.commit()
            temp_uname = User.get_guest_username(temp_user.id)
            temp_user.email = temp_uname
            temp_user.confirmed_at = datetime.utcnow()
            temp_user.roles.append(Role.query.filter_by(name='anonymous').one())
            db.session.commit()
            flash('Your temporary username is {} and password is {}. '.format(
                temp_uname, temp_pswd) + 'This account will be deleted in '
                '{} days.'.format(current_app.config['ANONYMOUS_MAX_AGE_DAYS']),
                'info')
            login_user(temp_user, force=True, remember=True)

        # CREATE USERFILE OBJECT
        user_id = current_user.id
        user_upload = UserFile(filename=mut_filename, user_id=user_id)
        form.to_model(user_upload)
        user_upload.is_valid = True
        user_upload.run_complete = False
        db.session.add(user_upload)
        db.session.commit()
        current_app.logger.info("Project {} uploaded by {}".format(
            user_upload.file_id, current_user.email))
        user_folder = naming_rules.get_user_folder(user_id)
        proj_folder = naming_rules.get_project_folder(user_upload)
        if not os.path.exists(user_folder):
            os.mkdir(user_folder)
        os.mkdir(proj_folder)
        file_path = os.path.join(proj_folder,
                                 user_upload.get_local_filename())
        # MOVE/COPY TEMP FILE
        if file_tester.line_endings == '\n':
            os.rename(temp_file.path, file_path)
        else:  # rewrite file with \n line endings
            with open(temp_file.path, 'rU') as file:
                with open(file_path, 'w') as out:
                    for line in file:
                        out.write(line)
            os.remove(temp_file.path)

        flash('File accepted and validated. Analysis in progress.',
              'success')
        # want analysis to run asynchronously!
        # http://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xi-email-support
        run_analysis(proj_folder, file_path, user_upload.file_id)
        return redirect(url_for('.index'))

    # MESSAGES FOR INITIAL UPLOAD PAGE ACCESS OR FAILED UPLOAD
    if current_user.is_authenticated():
        message = "You've run {} projects in the last week. ".format(n_week)
        message += "Your weekly limit is {}.".format(n_week_max)
        if incomplete:
            message += "\nYou have {} jobs still running."\
                .format(len(incomplete))
        flash(message, "info")
    else:
        flash("You can upload here as a new guest user, but with severe "
              "usage limits and no email alerts. Consider registering instead.",
              "danger")
    return render_template('pway/upload.html', form=form)
