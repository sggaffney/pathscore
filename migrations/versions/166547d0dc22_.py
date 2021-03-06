"""empty message

Revision ID: 166547d0dc22
Revises: 2e13e634399e
Create Date: 2014-12-15 15:31:52.845523

"""

# revision identifiers, used by Alembic.
revision = '166547d0dc22'
down_revision = '2e13e634399e'

from alembic import op
import sqlalchemy as sa


def upgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.add_column('uploads', sa.Column('filename', sa.String(length=255), nullable=True))
    op.add_column('uploads', sa.Column('genome_size', sa.String(length=255), nullable=True))
    op.add_column('uploads', sa.Column('ignore_genes', sa.Text(), nullable=True))
    op.add_column('uploads', sa.Column('n_cutoff', sa.Integer(), nullable=True))
    op.add_column('uploads', sa.Column('proj_suffix', sa.String(length=255), nullable=True))
    op.add_column('uploads', sa.Column('required_genes', sa.Text(), nullable=True))
    ### end Alembic commands ###


def downgrade():
    ### commands auto generated by Alembic - please adjust! ###
    op.drop_column('uploads', 'required_genes')
    op.drop_column('uploads', 'proj_suffix')
    op.drop_column('uploads', 'n_cutoff')
    op.drop_column('uploads', 'ignore_genes')
    op.drop_column('uploads', 'genome_size')
    op.drop_column('uploads', 'filename')
    ### end Alembic commands ###
