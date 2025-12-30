"""Pytest fixtures for PathScore tests."""

import os
import tempfile
import pytest


@pytest.fixture(scope='session')
def app():
    """Create application for testing.

    Uses TestingConfig which sets SKIP_REF_LOADING=True to avoid
    needing the MySQL refs database for unit tests.
    """
    # Set minimal required environment variables for testing
    os.environ.setdefault('SECRET_KEY', 'test-secret-key-do-not-use-in-production')
    os.environ.setdefault('SECURITY_PASSWORD_SALT', 'test-salt-do-not-use-in-production')
    os.environ.setdefault('SECURITY_PASSWORD_HASH', 'bcrypt')
    os.environ.setdefault('PARALLEL_MODE', 'off')

    # Create temp directories for testing
    temp_dir = tempfile.mkdtemp()
    os.environ.setdefault('DATA_ROOT', temp_dir)
    os.environ.setdefault('TEMP_FOLDER', temp_dir)
    os.environ.setdefault('LOG_PATH', os.path.join(temp_dir, 'test.log'))

    from app import create_app
    app = create_app('testing')

    app.config['WTF_CSRF_ENABLED'] = False  # Disable CSRF for testing
    app.config['SECURITY_REGISTERABLE'] = True

    yield app

    # Cleanup temp directory
    import shutil
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture(scope='function')
def client(app):
    """Create test client."""
    return app.test_client()


@pytest.fixture(scope='function')
def db_session(app):
    """Create database tables and provide session."""
    from app import db
    from app.models import Role

    with app.app_context():
        db.create_all()

        # Create default roles if they don't exist
        if not db.session.get(Role, 1):
            general_role = Role(name='general', description='General user')
            admin_role = Role(name='admin', description='Administrator')
            db.session.add(general_role)
            db.session.add(admin_role)
            db.session.commit()

        yield db.session

        db.session.remove()
        db.drop_all()


@pytest.fixture(scope='function')
def test_user(db_session):
    """Create a test user."""
    import uuid
    from flask_security import hash_password
    from app.models import User, Role

    user = User(
        email='test@example.com',
        password=hash_password('testpassword'),
        active=True,
        confirmed_at=None,
        fs_uniquifier=str(uuid.uuid4())  # Required by Flask-Security-Too
    )
    db_session.add(user)
    db_session.commit()

    # Add general role
    general_role = db_session.query(Role).filter_by(name='general').first()
    if general_role:
        user.roles.append(general_role)
        db_session.commit()

    return user


@pytest.fixture(scope='function')
def authenticated_client(client, test_user, app):
    """Create an authenticated test client."""
    with app.test_request_context():
        with client.session_transaction() as sess:
            # Flask-Security uses fs_uniquifier for session identification
            sess['_user_id'] = test_user.fs_uniquifier
    return client


@pytest.fixture
def sample_mutation_file():
    """Create a sample mutation file for testing."""
    content = """hugo_symbol\tentrez_id\tpatient_id
TP53\t7157\tpatient1
BRCA1\t672\tpatient1
KRAS\t3845\tpatient2
EGFR\t1956\tpatient2
PIK3CA\t5290\tpatient3
"""
    return content.encode('utf-8')


@pytest.fixture
def sample_mutation_file_with_annot():
    """Create a sample mutation file with annotation column."""
    content = """hugo_symbol\tentrez_id\tpatient_id\tannot
TP53\t7157\tpatient1\tmissense
BRCA1\t672\tpatient1\tnonsense
KRAS\t3845\tpatient2\tmissense
"""
    return content.encode('utf-8')
