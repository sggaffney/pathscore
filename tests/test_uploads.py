"""Tests for file upload validation."""

import pytest
from io import BytesIO
from werkzeug.datastructures import FileStorage

from app.uploads import MutationFile, BmrFile, BasicFile
from app.errors import ValidationError


class TestMutationFile:
    """Tests for MutationFile validation."""

    def test_valid_mutation_file(self, app, sample_mutation_file):
        """Test valid mutation file is accepted."""
        with app.app_context():
            file_storage = FileStorage(
                stream=BytesIO(sample_mutation_file),
                filename='mutations.txt'
            )
            # Should not raise
            mut_file = MutationFile(file_storage)
            assert mut_file.has_annot is False

    def test_valid_mutation_file_with_annot(self, app, sample_mutation_file_with_annot):
        """Test valid mutation file with annotation is accepted."""
        with app.app_context():
            file_storage = FileStorage(
                stream=BytesIO(sample_mutation_file_with_annot),
                filename='mutations.txt'
            )
            mut_file = MutationFile(file_storage)
            assert mut_file.has_annot is True

    def test_invalid_headers_rejected(self, app):
        """Test file with wrong headers is rejected."""
        content = b"wrong\theaders\there\ndata\t1\tpatient1\n"
        with app.app_context():
            file_storage = FileStorage(
                stream=BytesIO(content),
                filename='mutations.txt'
            )
            with pytest.raises(ValidationError) as exc_info:
                MutationFile(file_storage)
            assert "headers" in str(exc_info.value).lower()

    def test_empty_file_rejected(self, app):
        """Test empty file is rejected."""
        with app.app_context():
            file_storage = FileStorage(
                stream=BytesIO(b''),
                filename='empty.txt'
            )
            with pytest.raises(ValidationError):
                MutationFile(file_storage)

    def test_header_only_rejected(self, app):
        """Test file with only headers (no data) is rejected."""
        content = b"hugo_symbol\tentrez_id\tpatient_id\n"
        with app.app_context():
            file_storage = FileStorage(
                stream=BytesIO(content),
                filename='header_only.txt'
            )
            with pytest.raises(ValidationError) as exc_info:
                MutationFile(file_storage)
            assert "data" in str(exc_info.value).lower()

    def test_invalid_entrez_id_rejected(self, app):
        """Test non-integer entrez_id is rejected."""
        content = b"hugo_symbol\tentrez_id\tpatient_id\nTP53\tnotanumber\tpatient1\n"
        with app.app_context():
            file_storage = FileStorage(
                stream=BytesIO(content),
                filename='invalid.txt'
            )
            with pytest.raises(ValidationError) as exc_info:
                MutationFile(file_storage)
            assert "invalid" in str(exc_info.value).lower()


class TestBmrFile:
    """Tests for BmrFile (background mutation rate) validation."""

    def test_valid_bmr_file(self, app):
        """Test valid BMR file is accepted."""
        content = b"hugo_symbol\tentrez_id\tper_mb\nTP53\t7157\t1.5\nBRCA1\t672\t0.8\n"
        with app.app_context():
            file_storage = FileStorage(
                stream=BytesIO(content),
                filename='bmr.txt'
            )
            # Should not raise
            bmr_file = BmrFile(file_storage)
            assert bmr_file is not None

    def test_bmr_headers_required(self, app):
        """Test BMR file requires specific headers."""
        content = b"gene\tid\trate\nTP53\t7157\t1.5\n"
        with app.app_context():
            file_storage = FileStorage(
                stream=BytesIO(content),
                filename='bmr.txt'
            )
            with pytest.raises(ValidationError):
                BmrFile(file_storage)
