import os
import pathlib

from flask import current_app
from flask_login import current_user

from .decorators import make_async
from .errors import ValidationError
from .models import UserFile
from .naming_rules import get_all_pways_path, get_comparison_dir
from .admin import export_detail_path

from helpers import compare_projects as cp


def get_comparison_table(proj_ids: list, nice_names=None, caption=None):
    """Get pathway enrichment comparison for multiple pathways.

    Order of input Project IDs determines reference (1st ID), comparison (2nd),
    and other (3rd onwards, optional).

    Args:
        proj_ids (list): list of project IDs.
        nice_names (list): (optional) list of names to use for projects in
            output, otherwise uses proj_suffix specified at upload time.
        caption (str): (optional) Name for comparison shown in first cell of
            styled output.

    Returns:
        styler pd.DataFrame if get_styled else ordinary pd.DataFrame.

    """
    # requires app context
    is_vip = True in [r.name == 'vip' for r in current_user.roles]
    n_ids = len(proj_ids)
    if is_vip:
        query = UserFile.query.filter_by(run_complete=1)
    else:
        query = UserFile.query.filter_by(user_id=current_user.id, run_complete=1)
    projs = query.filter(UserFile.file_id.in_(proj_ids)).all()
    if len(projs) != n_ids:
        raise ValidationError("Invalid project(s) requested")
    if nice_names is None:
        categ_names = [p.proj_suffix for p in projs]
    else:
        if len(nice_names) != n_ids:
            raise ValidationError("IDs and names have different lengths.")
        categ_names = list(nice_names)
    tsv_dict = dict()
    n_dict = dict()
    for proj_name, proj in zip(categ_names, projs):
        proj_id = proj.file_id
        all_path = get_all_pways_path(proj)
        summary_exists = os.path.exists(all_path)
        if not summary_exists:
            print(f"Gathering pathway details for proj={proj_id}.")
            export_detail_path(proj)
        tsv_dict[proj_name] = all_path
        n_dict[proj_name] = proj.n_patients
    # app = current_app._get_current_object()
    out_dir = get_comparison_dir(current_user.id, proj_ids)
    os.makedirs(out_dir, exist_ok=True)
    build_comparison_tables(categ_names, tsv_dict, n_dict,
                            caption=caption, out_dir=out_dir)


@make_async
def build_comparison_tables(categ_names, tsv_dict, n_dict, caption=None,
                            out_dir=None):
    complete_file = os.path.join(out_dir, '.complete')
    if os.path.exists(complete_file):
        return  # Files already created
    df = cp.gather_pathways(categ_names, tsv_dict, n_dict)
    ref_name = categ_names[0]
    compared_name = categ_names[1]
    other_compared = categ_names[2:]
    comp = cp.identify_enrichment(ref_name, compared_name, df,
                                  other_compared=other_compared)

    caption = f"{ref_name} > {compared_name}" if caption is None else caption
    styled = cp.get_styled_comparison(comp, caption)
    # SAVE EXCEL
    excel_path = os.path.join(out_dir, 'comparison.xlsx')
    styled.to_excel(excel_path, sheet_name=f"{caption}")
    # SAVE HTML
    html = styled.render()
    html_path = os.path.join(out_dir, 'comparison.html')
    with open(html_path, 'w') as out:
        out.write(html)
    # SAVE 'COMPLETE' INDICATOR FILE
    pathlib.Path(complete_file).touch()


def proj_str_from_ids(proj_ids):
    proj_str = '/'.join([str(i) for i in proj_ids])
    return proj_str


def proj_ids_from_str(proj_id_str):
    proj_ids = proj_id_str.split('/')
    return proj_ids

