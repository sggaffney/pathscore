"""empty message

Revision ID: 36d5d10ab490
Revises: 2e569b19ba7a
Create Date: 2015-12-17 11:50:17.708118

"""

# revision identifiers, used by Alembic.
revision = '36d5d10ab490'
down_revision = '2e569b19ba7a'

from alembic import op
import sqlalchemy as sa


def upgrade():
    op.alter_column('uploads', 'algorithm', server_default='gene_length')


def downgrade():
    op.alter_column('uploads', 'algorithm', server_default='gene_count')
