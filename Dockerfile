# PathScore Flask Application Dockerfile
# Python 3.12 with all dependencies for Flask web app and Celery workers

FROM python:3.12-slim-bookworm

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Build tools for Cython and native extensions
    gcc \
    g++ \
    # MySQL client libraries
    libmariadb-dev \
    libmariadb-dev-compat \
    pkg-config \
    # For matplotlib and other visualization
    libfreetype6-dev \
    libpng-dev \
    # Cleanup
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy and install Python dependencies first (better layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Create data directories
RUN mkdir -p /data/temp /data/logs /data/.matplotlib \
    && chmod -R 777 /data

# Copy application code
COPY pathscore.py .
COPY celery_worker.py .
COPY setup_cython.py .
COPY app app/
COPY helpers helpers/

# Compile Cython extension if present (using setup script for proper numpy includes)
RUN if [ -f app/comb_functions.pyx ]; then \
        python setup_cython.py build_ext --inplace; \
    fi

# Set environment variables
ENV FLASK_APP=pathscore.py
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV MPLCONFIGDIR=/data/.matplotlib

# Expose Flask port
EXPOSE 5000

# Default command (can be overridden in docker-compose)
CMD ["gunicorn", "-b", "0.0.0.0:5000", "-w", "2", "--timeout", "120", "pathscore:app"]
