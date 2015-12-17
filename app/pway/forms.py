# -*- coding: utf-8 -*-

from flask.ext.wtf import Form
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import SubmitField, RadioField, TextAreaField, \
    StringField
from flask.ext.wtf.html5 import IntegerField
from wtforms.validators import Length, DataRequired, optional, length, Regexp
from ..misc import GeneListTester


def hugo_validator():
    return Regexp(GeneListTester.hugo_re_str, flags=0,
                  message=u'Comma-separated list of valid HUGO symbols please.')


class UploadForm(Form):
    mut_file = FileField('Patient-mutations file', validators=[
        FileRequired(),
        FileAllowed(['txt', 'tsv'], 'Use txt or tsv extension')
    ])
    algorithm = RadioField('Algorithm', validators=[
        DataRequired(), Length(1, 128)], choices=[
        ('gene_length', 'Gene length'), ('gene_count', 'Gene count')],
        default='gene_length')
    # genome_size = RadioField('Genome size (gene count alg)', validators=[
    #     DataRequired(), Length(1, 128)], choices=[
    #     ('protein-coding', 'Protein-coding (20462)'),
    #     ('no_pseudo', 'All minus pseudogenes (28795)'),
    #     ('all', 'All (45466)'),
    #     ('inc_misc_chr', "All plus 'misc' chr genes (46286)")],
    #     default='protein-coding')
    n_cutoff = IntegerField('Gene limit per patient', validators=[optional()],
                            default=500)
    required_genes = TextAreaField(
        'Required genes (optional)',
        validators=[optional(), length(max=600, message="600 character limit."),
                    hugo_validator()])
    ignore_genes = TextAreaField(
        'Ignore genes (optional)',
        validators=[optional(), length(max=600, message="600 character limit."),
                    hugo_validator()])
    proj_suffix = StringField('Project name',
                              validators=[DataRequired(), length(
                                  max=40, message="40 character limit.")])
    submit = SubmitField('Upload')

    def to_model(self, upload):
        # upload.mut_file = self.mut_file.data
        # upload.genome_size = self.genome_size.data
        upload.n_cutoff = self.n_cutoff.data
        upload.algorithm = self.algorithm.data
        if self.required_genes.data:
            upload.required_genes = str(self.required_genes.data)
        if self.ignore_genes.data:
            upload.ignore_genes = str(self.ignore_genes.data)
        if self.proj_suffix.data:
            upload.proj_suffix = self.proj_suffix.data
