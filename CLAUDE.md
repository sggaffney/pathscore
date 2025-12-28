# CLAUDE.md - AI Agent Operational Reference

This file provides context for AI agents working on the PathScore codebase.

## Build, Test, and Run Commands

```bash
# Build containers
docker compose build flask
docker compose build mysql  # Only needed if Dockerfile.mysql changed

# Start services (development)
docker compose up -d mysql flask

# Start with Celery (production mode)
docker compose --profile celery up -d

# View logs
docker compose logs -f flask
docker compose logs flask --tail 50

# Rebuild and restart single service
docker compose build flask && docker compose up -d flask

# Run tests (not yet configured - Phase 6)
docker compose exec flask pytest

# Stop everything
docker compose down

# Nuclear reset (removes data volumes)
docker compose down -v
```

## Architecture Overview

### Request Flow

1. **Entry**: `pathscore.py` loads `.env`, creates Flask app via `app.create_app()`
2. **Routing**: Blueprints (`pway/`, `demo/`, `api/`, `auth/`) handle requests
3. **Authentication**: Flask-Security-Too manages users/roles via `models.py`
4. **Analysis Pipeline**: Upload triggers `get_effective_pathways.run_analysis()`
5. **Async Execution**: `decorators.make_async` routes to threads or Celery based on `PARALLEL_MODE`
6. **Database**: SQLAlchemy ORM for `pway` database; raw SQL for `refs` database queries

### Databases

- **`refs`** (read-only reference data): Gene lengths, pathway-gene links, NCBI data. Hardcoded database name throughout codebase.
- **`pway`** (application data): Users, roles, projects, BMR files. Name configurable via `MYSQLDB_DB`.

### Key Files by Function

| File | Purpose | Complexity |
|------|---------|------------|
| `app/__init__.py` | App factory, extension init, blueprint registration | Medium |
| `app/config.py` | Environment-based config classes | Low |
| `app/models.py` | SQLAlchemy models, user loader, BMR processor | High |
| `app/get_effective_pathways.py` | Core analysis pipeline (~1200 lines) | Very High |
| `app/db_lookups.py` | Raw SQL queries against `refs` database | Medium |
| `app/plot_fns.py` | Bokeh visualizations | High |
| `app/decorators.py` | `make_async`, upload limiting | Medium |
| `app/comb_functions.pyx` | Cython likelihood calculations | Low (don't modify) |

## Directory Structure

```
pway_app/
├── pathscore.py              # Entry point, loads .env
├── celery_worker.py          # Celery app creation (when PARALLEL_MODE=celery)
├── setup_cython.py           # Cython build script with NumPy includes
├── app/
│   ├── __init__.py           # create_app factory, celery init
│   ├── config.py             # Config, DevelopmentConfig, ProductionConfig
│   ├── models.py             # User, Role, UserFile, CustomBMR, BmrProcessor
│   ├── decorators.py         # make_async (off/threads/celery), limit_user_uploads
│   ├── get_effective_pathways.py  # run_analysis, MutationTable, pathway scoring
│   ├── db_lookups.py         # SQL queries for pathway/gene lookups
│   ├── plot_fns.py           # Bokeh scatter, MDS plots
│   ├── emails.py             # Celery mail tasks
│   ├── uploads.py            # File validation classes
│   ├── naming_rules.py       # Path generation for project files
│   ├── admin.py              # Cleanup threads, project deletion
│   ├── errors.py             # Custom exception classes
│   ├── misc.py               # Utility functions
│   ├── compare.py            # Project comparison features
│   ├── comb_functions.pyx    # Cython: get_pway_likelihood_cython
│   ├── plot/                 # matplotlib gene matrix plotting
│   ├── pway/                 # Main blueprint (authenticated users)
│   │   ├── routes.py         # /upload, /scatter, /results, /mds
│   │   └── forms.py          # WTForms for upload
│   ├── demo/                 # Public demo blueprint
│   │   ├── routes.py         # /demo/scatter, /demo/results
│   │   └── helpers.py        # get_all_demos, get_single_demo
│   ├── api/                  # REST API blueprint
│   │   ├── routes.py         # /api/projects/, /api/bmr/
│   │   ├── auth.py           # HTTP Basic Auth
│   │   └── decorators.py     # @json, @etag, @collection
│   └── auth/                 # Auth blueprint (legacy, mostly Flask-Security now)
│       └── routes.py         # /alogin, /alogout
├── data/                     # Reference SQL files (loaded into MySQL)
│   ├── refs_*.sql.gz         # Gene/pathway reference data
│   ├── create_dbs.sql        # Database creation script
│   └── refs_pathways.sql     # Minimal pathways for testing
├── helpers/                  # Compare projects utilities
│   └── compare_projects.py
├── docker-compose.yml        # Service orchestration
├── Dockerfile                # Flask app container
├── Dockerfile.mysql          # MySQL with ref data preloaded
└── requirements.txt          # Python dependencies
```

## Patterns and Conventions

### SQLAlchemy 2.0 Query Style

All database queries use SQLAlchemy 2.0 patterns:

```python
# Get by ID
user = db.session.get(User, id)

# Query with filters
stmt = select(UserFile).where(
    UserFile.user_id == current_user.id,
    UserFile.run_complete == True
)
results = db.session.execute(stmt).scalars().all()

# Single result
result = db.session.execute(stmt).scalar_one_or_none()
if result is None:
    abort(404)

# Raw SQL must use text()
db.session.execute(text("DROP TABLE `{}`".format(table_name)))
```

### Blueprint Organization

Each blueprint follows the pattern:
- `__init__.py` - Creates Blueprint, optional `before_request`
- `routes.py` - Route handlers
- `forms.py` - WTForms (if needed)
- `helpers.py` - Shared functions (if needed)

### PARALLEL_MODE Execution

```python
# In decorators.py
@make_async
def some_task():
    pass

# make_async checks PARALLEL_MODE:
# - 'off': runs synchronously
# - 'threads': spawns Thread
# - 'celery': queues as Celery task
```

### File Paths

All project files go through `naming_rules.py`:

```python
naming_rules.get_project_folder(upload_obj)  # /data/user_<id>/proj_<id>/
naming_rules.get_detailed_path(upload_obj)   # .../pathways_detailed.txt
naming_rules.get_js_name(upload_obj)         # JavaScript variable name
```

## Anti-patterns to Avoid

1. **Don't use `query.filter_by()`** - Use `select().where()` (SQLAlchemy 2.0)

2. **Don't use `result.rowcount` for SELECTs** - Use `len(result.all())` instead

3. **Don't hardcode database names except `refs`** - Use config for `pway` database

4. **Don't import `six`** - Python 2 compatibility removed

5. **Don't use `xrange`** - Use `range` (Python 3)

6. **Don't use `attachment_filename`** - Use `download_name` (Flask 2.0+)

7. **Don't iterate result objects multiple times** - Call `.all()` first, then iterate

8. **Don't skip `text()` wrapper for raw SQL** - SQLAlchemy 2.0 requires it

## File Dependencies

### If you change X, also update Y

| Changed File | Also Update |
|--------------|-------------|
| `app/models.py` (add column) | Database migration, possibly `naming_rules.py` |
| `app/config.py` (add config) | `docker-compose.yml` environment, `.env.template` |
| `requirements.txt` | Rebuild Docker image |
| `app/comb_functions.pyx` | Run `python setup_cython.py build_ext --inplace` |
| `data/*.sql` | Rebuild `Dockerfile.mysql` |
| Blueprint routes | Check both `pway/` and `demo/` for parallel routes |

### Shared Query Patterns

`pway/routes.py` and `demo/routes.py` have parallel implementations for many routes (scatter, mds, tree, compare). Changes to one often need mirroring to the other. Demo uses `get_all_demos()` helper instead of user-filtered queries.

## Known Gotchas

1. **MySQL initialization takes ~30 seconds** on first container start. Flask will fail to connect if started too early. The `depends_on: condition: service_healthy` in docker-compose.yml handles this.

2. **`refs` database is hardcoded** throughout `db_lookups.py` as `refs.tablename`. Cannot be configured.

3. **Cython extension must be pre-compiled** in Docker. The `setup_cython.py` script handles NumPy include paths. If you see `pyximport` warnings in logs, the pre-compiled extension isn't being found.

4. **`flask-security-too` not `flask-security`** - The package was forked; use the `-too` variant.

5. **Port 5000 conflict on macOS** - AirPlay Receiver uses port 5000. Either disable it or use `FLASK_PORT=5001` in `.env`.

6. **`LOCAL_MODE` not yet implemented** - References exist in config but the feature isn't complete.

7. **Mail features require SMTP config** - Without `MAIL_SERVER`, notification emails silently fail.

8. **`UserFile.user_id` vs `current_user.id`** - In `pway/routes.py` line 559, there's a bug using `current_user.user_id` instead of `current_user.id`. Fixed in Phase 2 but watch for similar issues.

## Technical Debt

1. **Bokeh 3.x migration incomplete** - `plot_fns.py` and route files still use old Bokeh API (`plot_width` instead of `width`)

2. **Celery 5.x migration incomplete** - Config namespace not updated (`CELERY_BROKER_URL` should be `broker_url`)

3. **No test suite** - Phase 6 planned but not implemented

4. **MDS visualization deprioritized** - Works but uses older patterns

5. **Compare feature complexity** - `compare.py` and comparison routes are complex and fragile

6. **Mixed raw SQL and ORM** - `db_lookups.py` uses raw SQL for `refs`, ORM for `pway`

## Common Tasks

### Add a new route

1. Add route in appropriate blueprint (`app/pway/routes.py` or `app/demo/routes.py`)
2. If authenticated, use `@login_required` decorator
3. Use SQLAlchemy 2.0 query patterns
4. Add template in `app/templates/<blueprint>/`

### Add a new model field

1. Add column in `app/models.py`
2. Run app once to trigger `db.create_all()` (or use migration)
3. Update any queries that should include the new field

### Debug database issues

```bash
# Connect to MySQL container
docker compose exec mysql mysql -u www -ppathscore_dev pway

# Check table structure
DESCRIBE user_files;

# Check refs database
USE refs;
SHOW TABLES;
```

### Check Celery task status

```bash
docker compose logs celery-worker --tail 100
```

### Force rebuild everything

```bash
docker compose down -v
docker compose build --no-cache
docker compose up -d mysql flask
```
