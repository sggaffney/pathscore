"""empty message

Revision ID: 4a4f183423f3
Revises: eef96661db5
Create Date: 2016-02-01 16:52:09.648589

"""

# revision identifiers, used by Alembic.
revision = '4a4f183423f3'
down_revision = 'eef96661db5'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import mysql

def upgrade():
    op.alter_column('uploads', 'algorithm',
                    server_default=sa.text(u"'bmr_length'"))


def downgrade():
    op.alter_column('uploads', 'algorithm',
                    server_default=sa.text(u"'gene_length'"))
