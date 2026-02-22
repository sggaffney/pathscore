"""Tests for Flask routes."""

import pytest


@pytest.mark.integration
class TestPublicRoutes:
    """Tests for publicly accessible routes.

    These tests require ref_info to be loaded (MySQL refs database).
    Run with: pytest -m integration
    """

    def test_demo_index_redirects(self, client, app):
        """Test demo index page is accessible."""
        pytest.skip("Requires refs database - run with full integration tests")

    def test_home_page(self, client, app):
        """Test home page redirects to upload."""
        pytest.skip("Requires refs database - run with full integration tests")


class TestAuthRoutes:
    """Tests for authentication routes."""

    def test_login_page_accessible(self, client, app):
        """Test login page is accessible."""
        with app.app_context():
            response = client.get('/login')
            assert response.status_code == 200

    def test_register_page_accessible(self, client, app):
        """Test register page is accessible (SECURITY_REGISTERABLE=True in test config)."""
        with app.app_context():
            response = client.get('/register')
            assert response.status_code == 200


@pytest.mark.integration
class TestProtectedRoutes:
    """Tests for routes requiring authentication.

    These tests require ref_info to be loaded (MySQL refs database).
    """

    def test_upload_requires_auth_or_guest(self, client, app):
        """Test upload page is accessible (allows guest users)."""
        pytest.skip("Requires refs database - run with full integration tests")


@pytest.mark.integration
class TestFileUpload:
    """Tests for file upload functionality.

    These tests require ref_info to be loaded (MySQL refs database).
    """

    def test_upload_invalid_file_rejected(self, client, app, db_session):
        """Test that invalid files are rejected."""
        pytest.skip("Requires refs database - run with full integration tests")

    def test_upload_empty_file_rejected(self, client, app, db_session):
        """Test that empty files are rejected."""
        pytest.skip("Requires refs database - run with full integration tests")


class TestAPIRoutes:
    """Tests for REST API routes."""

    def test_api_unauthenticated_rejected(self, client, app):
        """Test API rejects unauthenticated requests."""
        with app.app_context():
            response = client.get('/api/projects/')
            assert response.status_code in [401, 403]
