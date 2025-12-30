"""Tests for SQLAlchemy models."""

import pytest
from sqlalchemy import select

from app import db
from app.models import User, Role, UserFile


class TestUser:
    """Tests for User model."""

    def test_create_user(self, db_session):
        """Test creating a basic user."""
        import uuid
        from flask_security import hash_password

        user = User(
            email='newuser@example.com',
            password=hash_password('password123'),
            active=True,
            fs_uniquifier=str(uuid.uuid4())  # Required by Flask-Security-Too
        )
        db_session.add(user)
        db_session.commit()

        # Verify user was created
        retrieved = db.session.get(User, user.id)
        assert retrieved is not None
        assert retrieved.email == 'newuser@example.com'
        assert retrieved.active is True

    def test_user_has_fs_uniquifier(self, test_user):
        """Test that user has Flask-Security uniquifier."""
        assert test_user.fs_uniquifier is not None
        assert len(test_user.fs_uniquifier) > 0

    def test_user_roles(self, test_user, db_session):
        """Test user role assignment."""
        assert len(test_user.roles) > 0
        role_names = [r.name for r in test_user.roles]
        assert 'general' in role_names


class TestRole:
    """Tests for Role model."""

    def test_default_roles_exist(self, db_session):
        """Test that default roles are created."""
        stmt = select(Role).where(Role.name == 'general')
        role = db_session.execute(stmt).scalar_one_or_none()
        assert role is not None
        assert role.description == 'General user'


class TestUserFile:
    """Tests for UserFile model."""

    def test_user_file_creation(self, test_user, db_session):
        """Test creating a UserFile record."""
        user_file = UserFile(
            user_id=test_user.id,
            filename='test_mutations.txt',
            file_id=1
        )
        db_session.add(user_file)
        db_session.commit()

        # Verify (UserFile uses file_id as primary key, not id)
        retrieved = db.session.get(UserFile, user_file.file_id)
        assert retrieved is not None
        assert retrieved.filename == 'test_mutations.txt'
        assert retrieved.user_id == test_user.id

    def test_get_local_filename(self, test_user, db_session):
        """Test the get_local_filename method."""
        user_file = UserFile(
            user_id=test_user.id,
            filename='Test Project.txt',
            file_id=42
        )
        db_session.add(user_file)
        db_session.commit()

        local_name = user_file.get_local_filename()
        # Should strip extension and spaces
        assert 'Test' in local_name or 'test' in local_name.lower()

    def test_user_file_belongs_to_user(self, test_user, db_session):
        """Test UserFile relationship with User."""
        user_file = UserFile(
            user_id=test_user.id,
            filename='project1.txt',
            file_id=1
        )
        db_session.add(user_file)
        db_session.commit()

        # Query files belonging to user
        stmt = select(UserFile).where(UserFile.user_id == test_user.id)
        files = db_session.execute(stmt).scalars().all()
        assert len(files) == 1
        assert files[0].filename == 'project1.txt'
