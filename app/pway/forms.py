# -*- coding: utf-8 -*-

from flask_wtf import Form
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import SubmitField, RadioField, TextAreaField, \
    StringField, SelectField
from wtforms.fields.html5 import IntegerField
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
        ('bmr_length', 'Gene length, BMR-scaled'),
        ('gene_length', 'Gene length, unscaled'),
        ('gene_count', 'Gene count')],
        default='bmr_length')
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
    bmr = SelectField('Custom BMR', validators=[DataRequired()], coerce=int)
    submit = SubmitField('Upload')

    def to_model(self, upload):
        # upload.mut_file = self.mut_file.data
        # upload.genome_size = self.genome_size.data
        upload.n_cutoff = self.n_cutoff.data
        upload.algorithm = self.algorithm.data
        if self.bmr.data > 0:
            upload.bmr_id = self.bmr.data
        if self.required_genes.data:
            upload.required_genes = str(self.required_genes.data)
        if self.ignore_genes.data:
            upload.ignore_genes = str(self.ignore_genes.data)
        if self.proj_suffix.data:
            upload.proj_suffix = self.proj_suffix.data


class BmrForm(Form):
    bmr_file = FileField('Custom BMR file', validators=[
        FileRequired(),
        FileAllowed(['txt', 'tsv'], 'Use txt or tsv extension')
    ])
    # title, tissue, description
    title = StringField(u'Brief name', [DataRequired(), length(
        max=32, message="32 character limit.")])
    tissue = StringField(u'Tissue type (optional)', validators=
        [optional(), length(max=100, message="100 character limit.")])
    description = TextAreaField('Description', validators=
        [optional(), length(max=255, message="255 character limit.")])
    submit = SubmitField('Upload')

    def to_model(self, bmr):
        # upload.mut_file = self.mut_file.data
        # upload.genome_size = self.genome_size.data
        bmr.title = self.title.data
        bmr.tissue = self.tissue.data
        if self.description.data:
            bmr.description = self.description.data
