# Phase 3: Celery 5.x Migration - Completed

**Date:** 2025-12-28

## Summary

Successfully migrated from Celery 4.x configuration to Celery 5.x patterns. The Celery worker now starts correctly and registers all tasks.

## Changes Made

### 1. app/config.py - Celery Configuration Namespace

Changed from flat `CELERY_*` config keys to nested `CELERY` dict with lowercase keys:

```python
# Before (Celery 4.x)
CELERY_BROKER_URL = 'redis://localhost:6379/0'
CELERY_RESULT_BACKEND = 'redis://localhost:6379/0'
CELERY_RESULT_DB_SHORT_LIVED_SESSIONS = True
CELERY_DEFAULT_QUEUE = 'default'
CELERY_CREATE_MISSING_QUEUES = True
CELERY_ROUTES = {...}

# After (Celery 5.x)
CELERY = {
    'broker_url': os.environ.get('CELERY_BROKER_URL', 'redis://localhost:6379/0'),
    'result_backend': os.environ.get('CELERY_RESULT_BACKEND', 'redis://localhost:6379/0'),
    'result_extended': True,
    'task_default_queue': 'default',
    'task_create_missing_queues': True,
    'task_routes': {...},
    'broker_connection_retry_on_startup': True,
}
```

### 2. app/__init__.py - Celery Factory Pattern

Added proper Celery 5.x initialization with Flask app context:

```python
from celery import Celery, Task

def celery_init_app(app: Flask) -> Celery:
    """Initialize Celery with Flask app context (Celery 5.x pattern)."""
    class FlaskTask(Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery_app = Celery(app.name, task_cls=FlaskTask)
    celery_app.config_from_object(app.config['CELERY'])
    celery_app.set_default()
    app.extensions['celery'] = celery_app
    return celery_app
```

Changed initialization from:
```python
celery = Celery(__name__, broker=Config.CELERY_BROKER_URL)
# ...
celery.conf.update(app.config)
```

To:
```python
celery = Celery()
# ...
global celery
celery = celery_init_app(app)
```

### 3. app/emails.py - SQLAlchemy 2.0 Update

Fixed missed SQLAlchemy 1.x pattern:
```python
# Before
upload = UserFile.query.get(upload_id)

# After
upload = db.session.get(UserFile, upload_id)
```

### 4. pathscore.py - Celery Export

Added celery import after app creation for proper command-line access:
```python
from app import create_app
app = create_app(os.getenv('FLASK_CONFIG') or 'default')

# Import celery AFTER create_app() since it's initialized there
from app import celery
```

### 5. docker-compose.yml - Fixes

- Removed `./app:/app/app:ro` volume mount from `flask` service (was overwriting compiled Cython)
- Added `DEV_DATABASE_URL` to `celery-worker` service (needed for development config)

## Verification

Celery worker starts successfully and registers all tasks:

```
-------------- celery@74413cd16562 v5.6.0 (recovery)
.> transport:   redis://redis:6379/0
.> results:     redis://redis:6379/0
.> concurrency: 2 (prefork)

[tasks]
  . app.decorators.wrapped
  . app.emails.run_finished_notification_async
  . app.get_effective_pathways.run_analysis_async

celery@74413cd16562 ready.
```

## Files Modified

| File | Changes |
|------|---------|
| [app/config.py](app/config.py) | Celery 5.x config namespace |
| [app/__init__.py](app/__init__.py) | `celery_init_app()` factory, FlaskTask class |
| [app/emails.py](app/emails.py) | SQLAlchemy 2.0 query pattern |
| [pathscore.py](pathscore.py) | Export celery after app creation |
| [docker-compose.yml](docker-compose.yml) | Volume mount fix, DEV_DATABASE_URL |

## Task Decorator Compatibility

The existing `@celery.task` decorator works with Celery 5.x without changes. The FlaskTask base class automatically provides Flask app context to all tasks.

## Next Steps

- Phase 4: Bokeh 3.x Migration
- Phase 5: Python Modernization
