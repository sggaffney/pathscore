PathScore
=========

A pathway burden analysis web app.

Hosted at [http://pathscore.publichealth.yale.edu](http://pathscore.publichealth.yale.edu).


## Installation

Installation is unnecessary for most users, as the app is available for anyone to
use at [http://pathscore.publichealth.yale.edu](http://pathscore.publichealth.yale.edu).
But if you wish to install your own local copy of PathScore, the dependencies and
main installation steps are outlined below.

**Software dependencies**:

- [MySQL server]
- [Apache server]
- [nodejs], for [svgo]
- [Redis] and [mod_wsgi] (both installed automatically with Python virtual
  environment setup.) 


**Main installation steps**:

1. create required MySQL databases and set up users with suitable permissions;
1. (optional) install Apache. This is required for serving the app via 
mod_wsgi, which is more stable than flask's built-in server. Apache must be in 
place before installing mod_wsgi.
1. set up a Python environment with all required dependencies;
1. create an .env file with environment variables tailored to your server;
1. decide on one or more user 'roles' (which determine project upload frequency);
1. set up the directories used for data storage;
1. set up and start Redis daemon;
1. try serving the app using flask, with celery running in shell;
1. (optional) serve the app using mod_wsgi-express and Apache.
   1. test mod_wsgi-express can successfully serve its demo app;
   2. set up and start daemons to run mod_wsgi-express and celery.


Details on individual steps are below. This document is a work in progress, 
and is not yet comprehensive.


### Homebrew

If you're a Mac user, you can use [Homebrew] to install all of the non-python 
dependencies: MySQL server, Apache, and svgo. Homebrew can also set up the daemon 
that starts the MySQL server at boot time. (It can do the same 
for Redis if you're happy to install it twice -- it will also be 
installed in your Python environment when you set that up.)

```bash
brew install apache2
brew install mysql
brew install svgo  # also installs Node.js
brew install redis  # optional
```

You can use `brew services` to create the necessary daemon processes.
```bash
brew services mysql start  # start mysql server
brew services redis start  # start redis-server
```


### Background processes (daemons)

To have MySQL and Redis available as a background process, you can set them up
to run at boot or login time. On a Mac this is done using `.plist` files. The 
paid software [Lingon X] can help you create these files.

For Redis, the plist file might include the following command: 
```bash
/path/to/redis-server /optional/path/to/redis.conf "--daemonize no"
```

For MySQL the appropriate command would have the form:
```bash
/path/to/mysqld_safe --datadir=/path/to/mysql/dir
```


### MySQL databases setup

PathScore expects the following databases to be in place:
1. a *reference* database called `refs` (this is currently hardcoded) which 
    contains gene and pathway data.
2. a database to track projects, users, roles and background mutation rate files.
    You may wish to set up three versions of this database: development, testing,
    and production.    

The tables required by the `refs` database are provided as gzipped SQL files in 
the `data/` directory:
- refs_pathway_gene_link.sql.gz
- refs_entrez_length.sql.gz
- refs_ncbi_entrez.sql.gz

The tables in the projects database will be created automatically by the app
when it is loaded for the first time. You will need to manually add a role to 
the roles table after the table is created -- a future update will automate this
process.

```sql
INSERT INTO `roles` (`name`, `description`, `uploads_pw`)
VALUES ('general', 'The general public', 10);
```


### Setting up the Python environment

I recommend *pyenv* for easy management of multiple Python distributions. To 
create an appropriate virtual environment within a *conda* Python distribution, 
run the following in your shell:  
 
```bash
# Installation with shared libraries is required by mod_wsgi
PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install miniconda2-latest

# Create virtual environment using metadata yaml file
conda env create -f environment.lock.yaml

# Activate the environment
conda activate pway

# To deactivate, use `conda deactivate`
```

### .env environment file

Below is an example environment file that requires customization for your server.
Use this to create a file called `.env` in the repository root directory.

```bash
FLASK_CONFIG=development
SECRET_KEY=some-secure-password
MAIL_SERVER=smtp.example.com
MAIL_SENDER=your-email@example.com
MAIL_DEFAULT_SENDER=your-email@example.com
SECURITY_EMAIL_SENDER=your-email@example.com
MAIL_PORT=465
MAIL_PASSWORD=your-email-password
MAIL_USE_SSL=True
MAIL_USE_TLS=False
MAIL_ERROR_RECIPIENT=your-email@example.com
DATABASE_URL=mysql+pymysql://db-username:db-username-password@127.0.0.1:3306/pathscore
DEV_DATABASE_URL=mysql+pymysql://db-username:db-username-password@127.0.0.1:3306/pathscore_dev
TEST_DATABASE_URL=mysql+pymysql://db-username:db-username-password@127.0.0.1:3306/pathscore_test
SECURITY_CONFIRMABLE=True
SECURITY_REGISTERABLE=True
SECURITY_CHANGEABLE=True
SECURITY_PASSWORD_HASH=bcrypt
SECURITY_PASSWORD_SALT=some-secret-word
DATA_ROOT=/www/pway/
LOG_PATH=/www/pway/pathscore.log
TEMP_FOLDER=/www/pway/temp/
MYSQLDB_DB_DEV=pathscore_dev
MYSQLDB_DB_TEST=pathscore_test
MYSQLDB_DB=pathscore
MYSQLDB_CNF=/path/to/.my.cnf.pathscore
MYSQLDB_HOST=localhost
SQLALCHEMY_ECHO=False
SERVER_NAME=localhost:8000
MPLCONFIGDIR=/www/.matplotlib
SVGO_PATH=/path/to/svgo
ANALYTICS_ID=UA-123456789-1
```


## License

Copyright (C) 2015 Stephen Gaffney

PathScore is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PathScore is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PathScore.  If not, see <http://www.gnu.org/licenses/>.


[MySQL server]: https://dev.mysql.com/downloads/mysql/
[Apache server]: https://httpd.apache.org/
[Redis]: https://redis.io/
[svgo]: https://www.npmjs.com/package/svgo
[nodejs]: https://nodejs.org/
[Homebrew]: https://brew.sh/
[mod_wsgi]: https://pypi.org/project/mod-wsgi/