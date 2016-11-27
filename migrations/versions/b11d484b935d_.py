"""empty message

Revision ID: b11d484b935d
Revises: 63a13ecafbd
Create Date: 2016-11-28 10:55:35.698971

"""

# revision identifiers, used by Alembic.
revision = 'b11d484b935d'
down_revision = '63a13ecafbd'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import mysql

def upgrade():
    op.add_column('uploads', sa.Column('has_mds', sa.Boolean(), nullable=True))


def downgrade():
    op.drop_column('uploads', 'has_mds')

