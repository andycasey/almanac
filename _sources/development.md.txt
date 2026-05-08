# Development Guide

This guide covers development setup, contributing guidelines, and testing for the `almanac` project.

## Development Setup

### Prerequisites

- **Python**: 3.8 or higher
- **Git**: For version control
- **uv**: Recommended for dependency management
- **Access**: SDSS database access for full functionality

### Installation for Development

1. **Clone the repository**:
```bash
git clone https://github.com/sdss/almanac.git
cd almanac
```

2. **Install development dependencies**:
```bash
# Using uv (recommended)
uv pip install -e ".[dev]"

# Or using pip
pip install -e ".[dev]"
```

3. **Verify installation**:
```bash
almanac --help
almanac config show
```

## Project Structure

```
almanac/
├── src/almanac/         # Main package source
│   ├── __init__.py     # Version information
│   ├── cli.py          # Command-line interface
│   ├── config.py       # Configuration management
│   ├── database.py     # Database connectivity
│   ├── apogee.py       # APOGEE data processing
│   ├── io.py           # HDF5 I/O operations
│   ├── display.py      # Output formatting
│   ├── logger.py       # Logging configuration
│   └── utils.py        # Utility functions
├── docs/               # Documentation
├── tests/              # Test suite (if present)
├── bin/                # Utility scripts
├── pyproject.toml      # Project configuration
├── CHANGELOG.rst       # Change log
├── LICENSE.md          # License
└── README.md           # Main documentation
```

## Development Workflow

### Code Style

The project follows these style conventions:

#### Python Code Style
- **Line length**: 99 characters (configured in `pyproject.toml`)
- **Import sorting**: Uses `isort` with custom sections for SDSS packages
- **Linting**: Uses `flake8` with custom rules

#### Import Organization
```python
# Standard library imports
import os
import sys

# Third-party imports
import numpy as np
import click

# SDSS-specific imports
import sdssdb
from sdsstools import get_logger

# Local imports
from almanac.config import config
from almanac.logger import logger
```

### Code Quality Tools

#### Linting
```bash
# Run flake8
flake8 src/almanac/

# Run with specific config
flake8 --config pyproject.toml src/
```

#### Import Sorting
```bash
# Check import order
isort --check-only src/almanac/

# Fix import order
isort src/almanac/
```

#### Code Coverage
```bash
# Run tests with coverage
pytest --cov almanac --cov-report html

# View coverage report
open htmlcov/index.html
```

## Testing

### Test Framework

The project uses `pytest` for testing with additional plugins:

- **pytest-cov**: Coverage reporting
- **pytest-mock**: Mocking support  
- **pytest-asyncio**: Async testing
- **pytest-sugar**: Enhanced output

### Running Tests

```bash
# Run all tests
pytest

# Run with verbose output
pytest -v

# Run with coverage
pytest --cov almanac --cov-report html

# Run specific test file
pytest tests/test_config.py

# Run specific test
pytest tests/test_config.py::test_load_config
```

### Test Configuration

Tests are configured in `pyproject.toml`:

```toml
[tool.pytest.ini_options]
addopts = "--cov almanac --cov-report html -W ignore"
```

### Writing Tests

#### Test Structure
```python
import pytest
from almanac.config import Config, load_config_file

def test_default_config():
    """Test default configuration values."""
    config = Config()
    assert config.logging_level == 20
    assert config.sdssdb.port == 5432

@pytest.fixture
def temp_config_file(tmp_path):
    """Create temporary config file for testing."""
    config_file = tmp_path / "test_config.yaml"
    config_file.write_text("""
logging_level: 10
sdssdb:
  host: test-host.org
""")
    return config_file

def test_load_config_file(temp_config_file):
    """Test loading configuration from file."""
    config_data = load_config_file(str(temp_config_file))
    assert config_data['logging_level'] == 10
    assert config_data['sdssdb']['host'] == 'test-host.org'
```

#### Mocking Database Connections
```python
@pytest.fixture
def mock_database(mocker):
    """Mock database connection for testing."""
    mock_db = mocker.patch('almanac.database.database')
    mock_db.set_profile.return_value = True
    return mock_db

def test_database_functionality(mock_database):
    """Test database-dependent functionality."""
    from almanac.database import is_database_available
    assert is_database_available is True
```

## Contributing

### Contribution Process

1. **Fork** the repository on GitHub
2. **Clone** your fork locally
3. **Create** a feature branch: `git checkout -b feature/new-feature`
4. **Make** your changes following the code style guidelines
5. **Add** tests for new functionality
6. **Run** the test suite: `pytest`
7. **Commit** your changes: `git commit -m "Add new feature"`
8. **Push** to your fork: `git push origin feature/new-feature`
9. **Submit** a pull request

### Pull Request Guidelines

#### Before Submitting
- [ ] All tests pass: `pytest`
- [ ] Code follows style guidelines: `flake8 src/` and `isort --check-only src/`
- [ ] Documentation updated for new features
- [ ] Changelog updated (`CHANGELOG.rst`)
- [ ] Version number updated if appropriate

#### Pull Request Description
Include in your PR description:
- **Summary**: Brief description of changes
- **Motivation**: Why this change is needed
- **Testing**: How the change was tested
- **Breaking Changes**: Any backwards compatibility issues

### Code Review Process

1. **Automated checks**: GitHub Actions run tests and style checks
2. **Maintainer review**: Code review by project maintainers
3. **Feedback incorporation**: Address review comments
4. **Approval and merge**: After approval, changes are merged

## Building and Distribution

### Local Development Build

```bash
# Build the package
python -m build

# Install locally
pip install -e .
```

### Version Management

Version is managed in `src/almanac/__init__.py`:

```python
__version__ = "0.1.11"
```

Update version for releases following [semantic versioning](https://semver.org/).

### Release Process

1. **Update version** in `__init__.py`
2. **Update CHANGELOG.rst** with release notes
3. **Commit changes**: `git commit -m "Release version X.Y.Z"`
4. **Tag release**: `git tag vX.Y.Z`
5. **Push changes**: `git push origin main --tags`

## Documentation

### Building Documentation

Documentation is in Markdown format. For Sphinx documentation (if configured):

```bash
# Install docs dependencies
pip install -e ".[docs]"

# Build HTML documentation
sphinx-build -b html docs/ docs/_build/html/

# View documentation
open docs/_build/html/index.html
```

### Documentation Standards

- **API documentation**: Docstrings for all public functions and classes
- **User guides**: Step-by-step instructions with examples
- **Code examples**: Working examples that can be copy-pasted
- **Changelog**: Document all user-facing changes

## Debugging

### Development Debugging

```bash
# Enable debug logging
almanac config set logging_level 10

# Run with verbose output
almanac -vv --mjd -1

# Use debugger
python -m pdb -m almanac --mjd -1
```

### Common Development Issues

#### Database Connection
```bash
# Test database connectivity
almanac config show
python -c "from almanac.database import is_database_available; print(is_database_available)"
```

#### Import Issues
```bash
# Check package installation
python -c "import almanac; print(almanac.__version__)"

# Check module imports
python -c "from almanac import config, database, apogee"
```

#### Performance Profiling
```bash
# Profile execution
python -m cProfile -o profile.stats -m almanac --mjd -1
python -c "import pstats; pstats.Stats('profile.stats').sort_stats('cumulative').print_stats(10)"
```

## Continuous Integration

### GitHub Actions

The project uses GitHub Actions for:
- **Testing**: Run test suite on multiple Python versions
- **Code quality**: Check code style and import order  
- **Coverage**: Generate and report test coverage

### Configuration

CI configuration in `.github/workflows/`:
- `test.yml`: Run tests on push and pull requests
- `lint.yml`: Check code style and formatting
- `coverage.yml`: Generate coverage reports

## Security

### Sensitive Information

- **Never commit**: Database passwords, API keys, or credentials
- **Use environment variables**: For sensitive configuration
- **Review dependencies**: Regularly update and audit dependencies

### Security Testing

```bash
# Check for known vulnerabilities
safety check

# Audit dependencies
pip-audit
```

## Performance Considerations

### Optimization Guidelines

- **Database queries**: Minimize database roundtrips
- **Parallel processing**: Use multiprocessing for I/O-bound operations
- **Memory usage**: Process data in chunks for large datasets
- **Caching**: Cache expensive computations when appropriate

### Profiling

```bash
# Memory profiling
mprof run almanac --mjd-start -30 --fibers
mprof plot

# Time profiling
python -m cProfile -s cumulative -m almanac --mjd -7
```

## Getting Help

### Development Support

- **Issues**: Open GitHub issues for bugs and feature requests
- **Discussions**: Use GitHub Discussions for questions
- **Email**: Contact maintainers for sensitive issues

### Resources

- [Python Packaging Guide](https://packaging.python.org/)
- [pytest Documentation](https://docs.pytest.org/)
- [Click Documentation](https://click.palletsprojects.com/)
- [SDSS Developer Documentation](https://sdss-python-dev.readthedocs.io/)