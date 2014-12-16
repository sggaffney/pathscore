# -*- coding: utf-8 -*-

from flask.ext.wtf import Form
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import SubmitField, RadioField, IntegerField, TextAreaField, \
    StringField
from wtforms.validators import Optional, Length, Required, DataRequired, \
    optional, length


class UploadForm(Form):
    mut_file = FileField('Patient-mutations file', validators=[
        FileRequired(),
        FileAllowed(['txt', 'tsv'], 'Use txt or tsv extension')
    ])
    genome_size = RadioField('Genome size', validators=[
        DataRequired(), Length(1, 128)], choices=[
        ('protein-coding', 'Protein-coding (20462)'),
        ('no_pseudo', 'All minus pseudogenes (28795)'),
        ('all', 'All (45466)'),
        ('inc_misc_chr', "All plus 'misc' chr genes (46286)")],
        default='client')  # TODO: check validation enforces choices
    n_cutoff = IntegerField('Gene limit per patient', [optional()])
    required_genes = TextAreaField('Required genes', validators=[optional(),
        length(max=600, message="600 character limit.")])
    ignore_genes = TextAreaField('Ignore genes', validators=[optional(),
        length(max=600, message="600 character limit.")])
    proj_suffix = StringField('Project suffix string', validators=[optional()])
    submit = SubmitField('Upload')

    def to_model(self, upload):
        upload.mut_file = self.mut_file.data
        upload.genome_size = self.genome_size.data
        upload.n_cutoff = self.n_cutoff.data
        upload.required_genes = self.required_genes.data
        upload.ignore_genes = self.ignore_genes.data
        upload.proj_suffix = self.proj_suffix.data
