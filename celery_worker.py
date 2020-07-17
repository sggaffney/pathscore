import os
from dotenv import load_dotenv
load_dotenv()


from app import celery, create_app

app = create_app(os.getenv('FLASK_CONFIG') or 'default')
app.app_context().push()
