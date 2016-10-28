from flask import current_app, render_template, url_for
from flask_mail import Message
from . import mail, celery
from .models import UserFile

_email_thread = None


def get_notification_email(email, subject, body_text, body_html):
    msg = Message(subject, recipients=[email])
    msg.body = body_text
    msg.html = body_html
    return msg


def run_finished_notification(upload_id):
    upload = UserFile.query.get(upload_id)
    # Prevent email attempt if user is anonymous
    if 'anonymous' in [r.name for r in upload.uploader.roles]:
        return
    msg = get_notification_email(
        email=upload.uploader.email,
        subject='[pathscore] Run complete',
        body_text=render_template('email/notify.txt', project=upload),
        body_html=render_template('email/notify.html', project=upload)
        )
    run_finished_notification_async.delay(msg)


# @async
@celery.task
def run_finished_notification_async(msg):
    with mail.connect() as conn:
        conn.send(msg)
