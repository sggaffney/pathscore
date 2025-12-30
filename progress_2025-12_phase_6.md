# Phase 6: Testing Infrastructure - Completed

**Date:** 2025-12-29

## Summary

Successfully implemented pytest-based testing infrastructure with 19 passing unit tests and 5 integration test stubs.

## Files Created

| File | Purpose |
|------|---------|
| [pytest.ini](pytest.ini) | Pytest configuration with markers |
| [tests/__init__.py](tests/__init__.py) | Test package |
| [tests/conftest.py](tests/conftest.py) | Pytest fixtures for app, db, users |
| [tests/test_models.py](tests/test_models.py) | SQLAlchemy model unit tests |
| [tests/test_routes.py](tests/test_routes.py) | Route/blueprint tests |
| [tests/test_uploads.py](tests/test_uploads.py) | File upload validation tests |

## Configuration Changes

### app/config.py
Added `SKIP_REF_LOADING = True` to TestingConfig to allow unit tests to run without the MySQL refs database.

### app/__init__.py
Added conditional check to skip `set_refs()` when `SKIP_REF_LOADING` is True:
```python
if not app.config.get('SKIP_REF_LOADING'):
    get_effective_pathways.set_refs(app)
```

### Dockerfile
Added tests directory to container:
```dockerfile
COPY pytest.ini .
COPY tests tests/
```

## Test Results

```
======================== 19 passed, 5 skipped in 4.65s =========================
```

### Passing Tests (19)

**Model Tests:**
- TestUser::test_create_user
- TestUser::test_user_has_fs_uniquifier
- TestUser::test_user_roles
- TestRole::test_default_roles_exist
- TestUserFile::test_user_file_creation
- TestUserFile::test_get_local_filename
- TestUserFile::test_user_file_belongs_to_user

**Route Tests:**
- TestAuthRoutes::test_login_page_accessible
- TestAuthRoutes::test_register_page_accessible
- TestProtectedRoutes::test_api_projects_requires_auth
- TestAPIRoutes::test_api_returns_json

**Upload Validation Tests:**
- TestMutationFile::test_valid_mutation_file
- TestMutationFile::test_valid_mutation_file_with_annot
- TestMutationFile::test_invalid_headers_rejected
- TestMutationFile::test_empty_file_rejected
- TestMutationFile::test_header_only_rejected
- TestMutationFile::test_invalid_entrez_id_rejected
- TestBmrFile::test_valid_bmr_file
- TestBmrFile::test_bmr_headers_required

### Skipped Tests (5)

These require the MySQL refs database and are marked with `@pytest.mark.integration`:
- TestPublicRoutes::test_demo_index_redirects
- TestPublicRoutes::test_home_page
- TestProtectedRoutes::test_upload_requires_auth_or_guest
- TestFileUpload::test_upload_invalid_file_rejected
- TestFileUpload::test_upload_empty_file_rejected

## Running Tests

```bash
# Run all tests (unit tests will pass, integration tests will skip)
docker compose run --rm flask pytest -v

# Run only unit tests (exclude integration tests)
docker compose run --rm flask pytest -v -m "not integration"

# Run with coverage
docker compose run --rm flask pytest --cov=app --cov-report=term-missing
```

## Key Design Decisions

1. **Session-scoped app fixture**: Creates the app once per test session for efficiency.

2. **Function-scoped db_session fixture**: Creates fresh database tables per test for isolation.

3. **SKIP_REF_LOADING config option**: Allows unit tests to run without the refs database (which requires MySQL and populated reference data).

4. **Integration test marker**: Tests requiring refs database are marked and skipped in unit test mode.

5. **Flask-Security-Too compatibility**: User fixtures explicitly set `fs_uniquifier` field.

6. **Primary key awareness**: UserFile uses `file_id` not `id` as primary key.

## Fixture Usage Guide

```python
# Basic app and client
def test_something(client, app):
    with app.app_context():
        response = client.get('/some-route')

# Database tests
def test_model(db_session):
    obj = MyModel(...)
    db_session.add(obj)
    db_session.commit()

# User tests
def test_with_user(test_user, db_session):
    # test_user is a User with email='test@example.com'
    assert test_user.email == 'test@example.com'

# File upload tests
def test_file_upload(app, sample_mutation_file):
    # sample_mutation_file is bytes of a valid mutation file
    file_storage = FileStorage(stream=BytesIO(sample_mutation_file), ...)
```

## Next Steps

- Add more test coverage as features are developed
- Consider adding integration test suite that runs against MySQL
- Add CI/CD pipeline (GitHub Actions)
