"""empty message

Revision ID: eef96661db5
Revises: 42ff7701dc87
Create Date: 2016-01-05 16:08:00.663203

"""

# revision identifiers, used by Alembic.
revision = 'eef96661db5'
down_revision = '42ff7701dc87'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import mysql

def upgrade():
    op.add_column('uploads', sa.Column('is_queued', sa.Boolean(),
                                       nullable=False,
                                       server_default=sa.text('FALSE')))


def downgrade():
    op.drop_column('uploads', 'is_queued')

