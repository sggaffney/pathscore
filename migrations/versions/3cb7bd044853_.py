"""empty message

Revision ID: 3cb7bd044853
Revises: 463a674df1c8
Create Date: 2015-04-22 14:15:46.696988

"""

# revision identifiers, used by Alembic.
revision = '3cb7bd044853'
down_revision = '463a674df1c8'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import mysql

def upgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.add_column('uploads', sa.Column('n_patients', sa.Integer(), nullable=True))
    ### end Alembic commands ###


def downgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.drop_column('uploads', 'n_patients')
    ### end Alembic commands ###
