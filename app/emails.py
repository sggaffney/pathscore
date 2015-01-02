from flask import current_app, render_template, url_for
from flask.ext.mail import Message
from . import mail
from .decorators import async

_email_thread = None


def get_notification_email(email, subject, body_text, body_html):
    msg = Message(subject, recipients=[email])
    msg.body = body_text
    msg.html = body_html
    return msg


def run_finished_notification(upload):
    app = current_app._get_current_object()
    run_finished_notification_async(app, upload)


@async
def run_finished_notification_async(app, upload):
    with app.app_context():
        msg = get_notification_email(
            email=upload.uploader.email,
            subject='[pathway search] Run complete',
            body_text=render_template('email/notify.txt', project=upload),
            body_html=render_template('email/notify.html', project=upload)
            )
        with mail.connect() as conn:
            conn.send(msg)
