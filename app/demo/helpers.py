from ..models import UserFile


def get_all_demos(mds_only=False):
    """Get list of upload objects for demo projects."""
    if not mds_only:
        return UserFile.query.filter_by(is_demo=1, run_complete=True).\
            order_by(UserFile.file_id).all()
    else:
        return UserFile.query.filter_by(is_demo=1, run_complete=True).\
            filter_by(has_mds=True).order_by(UserFile.file_id).all()


def get_single_demo(proj=None):
    """Get single demo project, by id."""
    return UserFile.query.filter_by(is_demo=1, file_id=proj).first_or_404()
