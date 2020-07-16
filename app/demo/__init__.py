import os
from flask import Blueprint, g, flash, url_for, current_app, redirect
from ..get_effective_pathways import ref_info
from ..errors import LimitError

demo = Blueprint('demo', __name__)


@demo.context_processor
def inject_n_pathways():
    return dict(
        n_pathways=len(ref_info.path_info_dict),
        get_common_js=get_common_js,
    )


@demo.before_request
def get_n_pathways():
    g.n_pathways = len(ref_info.path_info_dict)


@demo.errorhandler(LimitError)
def validation_error(e):
    flash(str(e), "danger")
    return redirect(url_for('.index'), code=307)


def get_common_js(is_archive=False):
    """Get text of common_js or script tag with link."""
    if not is_archive:
        js_text = '<script src="/static/pway_display_common.js"></script>'
    else:
        js_path = os.path.join(current_app.static_folder,
                               'pway_display_common.js')
        with open(js_path, 'r') as infile:
            js_text = infile.read()
            js_text = "<script>\n" + js_text + "\n</script>"
            js_text = js_text.replace(
                """'<img src="static/data/' + user_id + '/' + projId + '/""",
                """'<img src=\"""")
    return js_text


from . import routes
