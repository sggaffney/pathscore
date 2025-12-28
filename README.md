# PathScore

A pathway burden analysis web application for cancer genomics research.

PathScore analyzes patient-gene mutation data against pathway databases to identify pathways with statistically significant mutation burden, accounting for gene length and background mutation rates.

**Public instance:** [http://pathscore.publichealth.yale.edu](http://pathscore.publichealth.yale.edu)

## Quick Start (Docker)

```bash
# Clone and enter directory
git clone <repository-url>
cd pway_app

# Start MySQL and Flask
docker compose up -d mysql flask

# Wait for MySQL to initialize (~30 seconds on first run), then access:
# http://localhost:5000/demo/
```

## Requirements

- Docker and Docker Compose (recommended)
- OR: Python 3.12+, MySQL 8.x, Redis (optional)

## Running with Docker Compose

### Basic Usage (No Celery)

```bash
# Start services (threads mode - no Celery required)
docker compose up -d mysql flask

# View logs
docker compose logs -f flask

# Stop services
docker compose down

# Stop and remove data volumes
docker compose down -v
```

### With Celery (Production)

```bash
# Start all services including Redis and Celery worker
docker compose --profile celery up -d

# Check worker status
docker compose logs celery-worker
```

### Development Mode

```bash
# Start with Flask's dev server (hot reload)
docker compose --profile dev up -d mysql flask-dev
```

### Configuration

Create a `.env` file to override defaults:

```bash
# Required for production
SECRET_KEY=your-secure-secret-key
SECURITY_PASSWORD_SALT=your-password-salt

# Database credentials (defaults provided for dev)
MYSQL_USER=www
MYSQL_PASSWORD=your-db-password
MYSQL_ROOT_PASSWORD=your-root-password

# Ports
FLASK_PORT=5000
MYSQL_PORT=3307

# Execution mode: off, threads, or celery
PARALLEL_MODE=threads

# Mail (optional, for notifications)
MAIL_SERVER=smtp.example.com
MAIL_PORT=587
MAIL_USERNAME=user
MAIL_PASSWORD=pass
```

## Architecture

```
pway_app/
├── pathscore.py           # Application entry point
├── app/
│   ├── __init__.py        # Flask app factory
│   ├── config.py          # Environment-based configuration
│   ├── models.py          # SQLAlchemy models
│   ├── get_effective_pathways.py  # Core analysis pipeline
│   ├── db_lookups.py      # Database queries
│   ├── plot_fns.py        # Bokeh visualizations
│   ├── comb_functions.pyx # Cython likelihood calculations
│   ├── pway/              # Main authenticated routes
│   ├── demo/              # Public demo routes
│   ├── api/               # REST API
│   └── auth/              # Authentication
├── data/                  # Reference data SQL files
├── docker-compose.yml     # Container orchestration
├── Dockerfile             # Flask application image
└── Dockerfile.mysql       # MySQL with reference data
```

### Execution Modes

| Mode | `PARALLEL_MODE` | Services | Use Case |
|------|-----------------|----------|----------|
| Synchronous | `off` | MySQL, Flask | Debugging |
| Threaded | `threads` | MySQL, Flask | Development, light use |
| Celery | `celery` | MySQL, Redis, Flask, Celery | Production |

## Manual Installation

For development without Docker, see [plan_2025-12.md](plan_2025-12.md) for detailed requirements. Key steps:

1. Install Python 3.12+ with pip
2. Install MySQL 8.x server
3. Create databases: `refs` (reference data) and `pway` (application data)
4. Load reference data from `data/*.sql.gz`
5. Install Python dependencies: `pip install -r requirements.txt`
6. Create `.env` file with database credentials
7. Run: `flask run` or `gunicorn pathscore:app`

## Testing

```bash
# With Docker
docker compose exec flask pytest

# Local
pytest
```

Sample mutation file for testing: `skcm_ns_500_sm.txt`

## Documentation

- [plan_2025-12.md](plan_2025-12.md) - Modernization plan and architecture details
- [progress_2025-12_phase_1.md](progress_2025-12_phase_1.md) - Phase 1 completion notes
- [progress_2025-12_phase_2.md](progress_2025-12_phase_2.md) - Phase 2 completion notes
- [CLAUDE.md](CLAUDE.md) - AI/developer operational reference

## License

Copyright (C) 2015-2025 Stephen Gaffney

PathScore is free software under the GNU General Public License v3.0. See LICENSE for details.
