from flask.ext.wtf import Form
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import SubmitField


class UploadForm(Form):
    mut_file = FileField('Patient-mutations file', validators=[
        FileRequired(),
        FileAllowed(['txt', 'tsv'], 'Use txt or tsv extension')
    ])
    submit = SubmitField('Upload')
