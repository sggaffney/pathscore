# Phase 2: SQLAlchemy 2.0 Migration - Complete

**Date:** 2025-12-28

## Overview

Migrated all SQLAlchemy 1.x query patterns to SQLAlchemy 2.0 compatible syntax. The main changes involve:

1. Replacing `Model.query.get(id)` with `db.session.get(Model, id)`
2. Replacing `Model.query.filter_by()` chains with `select(Model).where()` statements
3. Wrapping raw SQL strings with `text()` for explicit textual SQL
4. Replacing `result.rowcount` iteration patterns with `result.all()` for reliable result handling
5. Updating deprecated Flask `attachment_filename` parameter to `download_name`

## Files Modified

### Core Application Files

| File | Changes |
|------|---------|
| [app/models.py](app/models.py) | Added `select, text` imports; replaced `query.get()` → `db.session.get()`; replaced `filter_by().one()` → `select().scalar_one()`; wrapped all raw SQL with `text()` in CustomBMR and BmrProcessor classes |
| [app/db_lookups.py](app/db_lookups.py) | Replaced `result.rowcount` + fetch loop patterns → `result.all()` iteration for reliable result handling in all query functions |
| [app/get_effective_pathways.py](app/get_effective_pathways.py) | Removed `six` import (Python 2 compat); added `text` import; replaced `query.get()` → `db.session.get()`; wrapped SQL with `text()`; modernized string handling with f-strings |
| [app/decorators.py](app/decorators.py) | Added `select` and `db` imports; replaced `filter_by()` chains → `select().where()` in `limit_user_uploads` decorator |

### Blueprint Route Files

| File | Changes |
|------|---------|
| [app/auth/routes.py](app/auth/routes.py) | Added `select` and `db` imports; replaced `query.filter_by().first()` → `select().scalar_one_or_none()` for login |
| [app/api/routes.py](app/api/routes.py) | Added `select` import; replaced all `query.filter_by()` and `first_or_404()` patterns; changed `attachment_filename` → `download_name` |
| [app/pway/routes.py](app/pway/routes.py) | Added `select` and `or_` imports; replaced all `query.filter_by()`, `query.get()`, and `first_or_404()` patterns; changed `attachment_filename` → `download_name` |
| [app/demo/routes.py](app/demo/routes.py) | Changed `attachment_filename` → `download_name` |
| [app/demo/helpers.py](app/demo/helpers.py) | Added `select`, `db`, `abort` imports; replaced `query.filter_by()` → `select().where()` patterns |

## Key Pattern Changes

### Query.get() → db.session.get()

```python
# Before (SQLAlchemy 1.x)
user = User.query.get(id)

# After (SQLAlchemy 2.0)
user = db.session.get(User, id)
```

### Query.filter_by() → select().where()

```python
# Before (SQLAlchemy 1.x)
uploads = UserFile.query.filter_by(user_id=current_user.id).all()

# After (SQLAlchemy 2.0)
stmt = select(UserFile).where(UserFile.user_id == current_user.id)
uploads = db.session.execute(stmt).scalars().all()
```

### filter_by().first_or_404() → Manual check

```python
# Before (SQLAlchemy 1.x)
project = UserFile.query.filter_by(file_id=id).first_or_404()

# After (SQLAlchemy 2.0)
stmt = select(UserFile).where(UserFile.file_id == id)
project = db.session.execute(stmt).scalar_one_or_none()
if project is None:
    abort(404)
```

### Raw SQL → text() wrapper

```python
# Before (SQLAlchemy 1.x)
db.session.execute("DROP TABLE {}".format(table_name))

# After (SQLAlchemy 2.0)
db.session.execute(text(f"DROP TABLE `{table_name}`"))
```

### Result iteration with rowcount → result.all()

```python
# Before (SQLAlchemy 1.x) - unreliable in 2.0
result = db.session.execute(cmd)
if result.rowcount:
    for row in result:
        process(row)

# After (SQLAlchemy 2.0)
result = db.session.execute(cmd)
rows = result.all()
if rows:
    for row in rows:
        process(row)
```

## Testing

- Docker containers rebuilt successfully
- Flask application starts without errors
- Demo page loads correctly at http://localhost:5001/demo/
- Database queries execute without SQLAlchemy deprecation warnings

## Notes

- The `rowcount` attribute is unreliable for SELECT statements in SQLAlchemy 2.0; use `len(result.all())` instead
- Flask's `send_file()` parameter `attachment_filename` was renamed to `download_name` in Flask 2.0+
- Python's `six` library (Python 2/3 compatibility) was removed as no longer needed for Python 3.12

## Next Phase

Phase 3: Celery 5.x Migration
- Update config namespace (`CELERY_BROKER_URL` → `broker_url`)
- Create proper `make_celery()` factory function
- Update task decorators as needed
