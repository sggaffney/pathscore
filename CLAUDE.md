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

# Run tests (62 unit tests, 5 integration stubs)
docker compose run --rm flask pytest -v

# Run only unit tests (skip integration tests that need refs database)
docker compose run --rm flask pytest -v -m "not integration"

# Run with coverage
docker compose run --rm flask pytest --cov=app

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
5. **Computation**: Main loop builds numpy arrays per pathway, calls `analyze_pathway()` from `app.computation`
6. **Async Execution**: `decorators.make_async` routes to threads or Celery based on `PARALLEL_MODE`
7. **Database**: SQLAlchemy ORM for `pway` database; raw SQL for `refs` database queries

### Layered Architecture (Phase 7/7b)

The codebase has a new layered architecture alongside the legacy code:

```
┌─────────────────────────────────────────────────────────────────┐
│                    DOMAIN LAYER (app/domain/)                    │
│  Gene, GeneSet, Pathway, PathwayMetadata, SizeAlgorithm         │
│  MutationDataset, PatientAnalysisData                           │
│  - Immutable frozen dataclasses                                 │
│  - No database or Flask dependencies                            │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│              PURE COMPUTATION LAYER (app/computation/)           │
│  analyze_pathway(), compute_effective_size()                     │
│  → PathwayAnalysisResult, PatientProbabilities                  │
│  - Pure functions: primitives in, dataclasses out               │
│  - No side effects, no database access                          │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│              DATA ACCESS LAYER (app/repositories/)              │
│  GeneRepository, PathwayRepository                              │
│  - Load domain objects from refs database                       │
│  - Currently use Flask-SQLAlchemy db.session                    │
└─────────────────────────────────────────────────────────────────┘
```

**Migration status**: `analyze_pathway()` is integrated into the main analysis loop (replaced `LCalculator`). Legacy classes `PathwaySummary`/`PathwaySummaryParsed` and `RefInfo` still exist but are slated for replacement.

### Databases

- **`refs`** (read-only reference data): Gene lengths, pathway-gene links, NCBI data. Hardcoded database name throughout codebase.
- **`pway`** (application data): Users, roles, projects, BMR files. Name configurable via `MYSQLDB_DB`.

### Key Files by Function

| File | Purpose | Complexity |
|------|---------|------------|
| `app/__init__.py` | App factory, extension init, blueprint registration | Medium |
| `app/config.py` | Environment-based config classes | Low |
| `app/models.py` | SQLAlchemy models, user loader, BMR processor | High |
| `app/get_effective_pathways.py` | Analysis pipeline: MutationTable, main loop, file I/O | Very High |
| `app/db_lookups.py` | Raw SQL queries against `refs` database | Medium |
| `app/plot_fns.py` | Bokeh visualizations | High |
| `app/decorators.py` | `make_async`, upload limiting | Medium |
| `app/comb_functions.pyx` | Cython likelihood calculations | Low (don't modify) |
| `app/domain/` | Immutable domain objects: Gene, GeneSet, Pathway, SizeAlgorithm | Medium |
| `app/computation/` | Pure functions: `analyze_pathway()`, PathwayAnalysisResult | Medium |
| `app/repositories/` | GeneRepository, PathwayRepository (data access layer) | Medium |

## Directory Structure

```
pway_app/
├── pathscore.py              # Entry point, loads .env
├── celery_worker.py          # Celery app creation (when PARALLEL_MODE=celery)
├── setup_cython.py           # Cython build script with NumPy includes
├── pytest.ini                # Pytest config with markers (integration)
├── app/
│   ├── __init__.py           # create_app factory, celery init
│   ├── config.py             # Config, DevelopmentConfig, ProductionConfig, TestingConfig
│   ├── models.py             # User, Role, UserFile, CustomBMR, BmrProcessor
│   ├── decorators.py         # make_async (off/threads/celery), limit_user_uploads
│   ├── get_effective_pathways.py  # run_analysis, MutationTable, main loop, file I/O
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
│   ├── domain/               # Immutable domain objects (Phase 7)
│   │   ├── gene.py           # Gene dataclass (frozen)
│   │   ├── geneset.py        # GeneSet with algorithm-aware size
│   │   ├── pathway.py        # Pathway, PathwayMetadata
│   │   ├── mutations.py      # MutationDataset, PatientAnalysisData
│   │   └── enums.py          # SizeAlgorithm enum
│   ├── computation/          # Pure computation functions (Phase 7)
│   │   ├── likelihood.py     # analyze_pathway(), compute_effective_size()
│   │   └── results.py        # PathwayAnalysisResult, PatientProbabilities
│   ├── repositories/         # Data access layer (Phase 7)
│   │   ├── gene_repository.py    # GeneRepository
│   │   └── pathway_repository.py # PathwayRepository
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
├── tests/                    # Pytest test suite (Phase 6+7)
│   ├── conftest.py           # Fixtures: app, db_session, test_user, sample files
│   ├── test_models.py        # SQLAlchemy model tests
│   ├── test_routes.py        # Route/blueprint tests
│   ├── test_uploads.py       # File validation tests
│   ├── test_domain.py        # Domain object tests (Gene, GeneSet, Pathway, etc.)
│   └── test_computation.py   # Pure computation tests (likelihood, analysis)
├── data/                     # Reference SQL files (loaded into MySQL)
│   ├── refs_*.sql.gz         # Gene/pathway reference data
│   ├── create_dbs.sql        # Database creation script
│   └── refs_pathways.sql     # Minimal pathways for testing
├── docs/                     # Architecture docs and proposals
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

### Domain Objects (frozen dataclasses)

New code should use immutable domain objects from `app.domain`:

```python
from app.domain import Gene, GeneSet, Pathway, SizeAlgorithm

# Gene with mutational properties
braf = Gene(entrez_id=673, symbol='BRAF', length_bp=2301, mutation_rate=5.0)

# GeneSet with algorithm-aware sizing
genes = GeneSet([braf, kras])
genes.get_size(SizeAlgorithm.GENE_COUNT)   # 2.0
genes.get_size(SizeAlgorithm.BMR_LENGTH)   # sum of effective_bp

# Pathway as named gene collection
pathway = Pathway(path_id=1, name='MAPK_SIGNALING', genes=genes)
```

### Pure Computation Functions

Use `analyze_pathway()` for likelihood calculations (replaced `LCalculator`):

```python
from app.computation import analyze_pathway

result = analyze_pathway(
    pathway_id=123,
    pathway_size=100,
    genome_size=18000,
    n_mutated_array=np.array([50, 30, 20, 40]),
    is_mutated_array=np.array([1, 0, 1, 0])
)
# Returns PathwayAnalysisResult with p_value, n_effective, effect_size, etc.
```

### Repository Pattern

Database access for domain objects goes through repositories:

```python
from app.repositories import GeneRepository, PathwayRepository

gene_repo = GeneRepository()
pathway_repo = PathwayRepository(gene_repo)
pathway = pathway_repo.get_pathway(path_id=123)
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

9. **Don't use `pd.np`** - Use `np` directly (pandas 2.x removed `pd.np`)

10. **Don't use `applymap`** - Use `map` for element-wise DataFrame operations (pandas 2.x)

11. **Don't use `iteritems()`** - Use `items()` for pandas Series iteration (pandas 2.x)

12. **Don't use `'rU'` file mode** - Use `'r'` with default newline handling (Python 3)

13. **Don't create new classes wrapping computation** - Use `analyze_pathway()` from `app.computation` directly. `LCalculator` was deleted in Phase 7b.

14. **Don't add mutable state to domain objects** - `Gene`, `GeneSet`, `Pathway` are `@dataclass(frozen=True)`. Keep them immutable.

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
| `app/domain/` classes | Update `tests/test_domain.py`, check `app/computation/` |
| `app/computation/` functions | Update `tests/test_computation.py`, check `get_effective_pathways.py` |
| `app/repositories/` | Check callers in routes and `get_effective_pathways.py` |

### Shared Query Patterns

`pway/routes.py` and `demo/routes.py` have parallel implementations for many routes (scatter, mds, tree, compare). Changes to one often need mirroring to the other. Demo uses `get_all_demos()` helper instead of user-filtered queries.

## Bokeh JavaScript Callbacks

The scatter plot in `pway/routes.py` and `demo/routes.py` uses Bokeh 3.x with custom JavaScript callbacks:

### Selection API (Bokeh 3.x)
```javascript
// Get/set selected indices
source.selected.indices    // Array of selected glyph indices
source.change.emit()       // Trigger reactive update
```

### Q-filtering Callback
The main callback filters pathways by q-value threshold:
1. Stores previously selected indices
2. Filters data where q1 or q2 ≤ cutoff
3. Updates `scatter_array` (global JS array mapping visible → original indices)
4. Preserves selection if item still visible after filter
5. Calls `updateIfSelectionChange_afterWait()` (app-level function that loads pathway images)

### App-defined JavaScript Functions
The Python-generated CustomJS callbacks interact with app-level JavaScript functions:
- `selectPathwaysByGenes()` - Called by gene inclusion/exclusion callbacks
- `updateIfSelectionChange_afterWait()` - Loads pathway images on selection change

These functions are defined in the HTML templates, not in Python code.

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

1. **Legacy classes coexist with new architecture** - `PathwaySummary`, `PathwaySummaryParsed`, and `RefInfo` still exist in `get_effective_pathways.py`. `LCalculator` has been deleted and replaced by `analyze_pathway()`. Next steps: replace `PathwaySummary` classes and `RefInfo` with domain/repository equivalents.

2. **Repositories not yet wired into main pipeline** - `GeneRepository` and `PathwayRepository` exist but `get_effective_pathways.py` still uses `RefInfo` and `db_lookups.py` for data access.

3. **MDS visualization deprioritized** - Works but uses older patterns

4. **Compare feature complexity** - `compare.py` and comparison routes are complex and fragile

5. **Mixed raw SQL and ORM** - `db_lookups.py` uses raw SQL for `refs`, ORM for `pway`

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
