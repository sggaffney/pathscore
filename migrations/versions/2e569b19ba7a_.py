"""empty message

Revision ID: 2e569b19ba7a
Revises: 54225938b150
Create Date: 2015-09-29 10:49:45.881512

"""

# revision identifiers, used by Alembic.
revision = '2e569b19ba7a'
down_revision = '54225938b150'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import mysql

def upgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.add_column('uploads', sa.Column('n_loaded', sa.Integer(), nullable=True))
    ### end Alembic commands ###


def downgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.drop_column('uploads', 'n_loaded')
    ### end Alembic commands ###
