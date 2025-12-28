from flask import abort
from sqlalchemy import select
from ..models import UserFile
from .. import db


def get_all_demos(mds_only=False):
    """Get list of upload objects for demo projects."""
    if not mds_only:
        stmt = select(UserFile).where(
            UserFile.is_demo == 1,
            UserFile.run_complete == True
        ).order_by(UserFile.file_id)
    else:
        stmt = select(UserFile).where(
            UserFile.is_demo == 1,
            UserFile.run_complete == True,
            UserFile.has_mds == True
        ).order_by(UserFile.file_id)
    return db.session.execute(stmt).scalars().all()


def get_single_demo(proj=None):
    """Get single demo project, by id."""
    stmt = select(UserFile).where(
        UserFile.is_demo == 1,
        UserFile.file_id == proj
    )
    result = db.session.execute(stmt).scalar_one_or_none()
    if result is None:
        abort(404)
    return result
