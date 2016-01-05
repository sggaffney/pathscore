import os
if os.path.exists('.env'):
    print('Importing environment from .env...')
    for line in open('.env'):
        var = line.strip().split('=')
        if len(var) == 2:
            os.environ[var[0]] = var[1]


from app import celery, create_app
print 'flask config:'
print os.getenv('FLASK_CONFIG')


app = create_app(os.getenv('FLASK_CONFIG') or 'default')
app.app_context().push()
