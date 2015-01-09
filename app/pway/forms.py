# -*- coding: utf-8 -*-

from flask.ext.wtf import Form
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import SubmitField, RadioField, TextAreaField, \
    StringField
from flask.ext.wtf.html5 import IntegerField
from wtforms.validators import Optional, Length, DataRequired, DataRequired, \
    optional, length, Regexp


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
        default='protein-coding')  # TODO: check validation enforces choices
    n_cutoff = IntegerField('Gene limit per patient', validators=[optional()],
                            default=500)
    required_genes = TextAreaField('Required genes (optional)', validators=
        [optional(), length(max=600, message="600 character limit."),
         Regexp("([A-Z0-9-]{1,}[,]{0,})+$", flags=0,
                message=u'Comma-separated list of valid HUGO symbols please.')])
    ignore_genes = TextAreaField('Ignore genes (optional)', validators=
        [optional(), length(max=600, message="600 character limit."),
         Regexp("([A-Z0-9-]{1,}[,]{0,})+$", flags=0,
                message=u'Comma-separated list of valid HUGO symbols please.')])
    proj_suffix = StringField('Project name',
                              validators=[DataRequired(),
                                          length(max=40, message=
                                          "40 character limit.")])
    submit = SubmitField('Upload')

    def to_model(self, upload):
        upload.mut_file = self.mut_file.data
        upload.genome_size = self.genome_size.data
        upload.n_cutoff = self.n_cutoff.data
        if self.required_genes.data:
            upload.required_genes = str(self.required_genes.data)
        if self.ignore_genes.data:
            upload.ignore_genes = str(self.ignore_genes.data)
        if self.proj_suffix.data:
            upload.proj_suffix = self.proj_suffix.data
