# Installation Guide

This guide covers installation of `almanac` in various environments.

## Requirements

- **Python**: 3.8 or higher
- **System Requirements**: Access to raw APOGEE data frames (local disk or mounted filesystem)
- **Database Access**: SDSS database access (automatic at Utah, requires whitelisting elsewhere)

## Installation Methods

### 1. At Utah (SDSS Computing Environment)

The simplest installation method if you're working at Utah:

```bash
module load almanac
```

**Setting up Python Environment at Utah** (if needed):
```bash
module load miniconda/3.8.5_astra
```

### 2. Using uv (Recommended)

Install using `uv` for dependency management:

```bash
uv pip install git+https://github.com/sdss/almanac
```

### 3. Using pip

Standard pip installation:

```bash
pip install git+https://github.com/sdss/almanac
```

### 4. Development Installation

For development work:

```bash
git clone https://github.com/sdss/almanac.git
cd almanac
uv pip install -e ".[dev]"
```

## Post-Installation Setup

### Database Configuration

If installing outside Utah, you'll need to:

1. **Get database access**: Ensure your IP address is whitelisted for SDSS database access
2. **Configure connection**: Use `almanac config` to set database parameters if needed

```bash
almanac config show  # View current configuration
almanac config set sdssdb.host your-host-address
```

### Data Access Setup

**Outside Utah**: Set up Globus transfer for raw APOGEE data frames to ensure local disk access.

### Verification

Verify installation:

```bash
almanac --help
almanac config show
```

Test basic functionality:

```bash
almanac --mjd -1  # Query yesterday's observations
```

## Dependencies

Core dependencies are automatically installed:

- **Data Processing**: `numpy`, `astropy`, `h5py`
- **Database**: `sdssdb`, `pydl`
- **CLI/Display**: `click`, `tqdm`, `rich`, `colorlog`

Development dependencies (with `[dev]` install):

- **Testing**: `pytest`, `pytest-cov`, `pytest-mock`
- **Documentation**: `sphinx`, `sphinx_bootstrap_theme`
- **Code Quality**: `flake8`, `isort`, `codecov`

## Troubleshooting

### Common Issues

**Database Connection Timeout**:
- Check network connectivity to SDSS database
- Verify IP whitelisting
- Adjust `database_connect_time_warning` config if needed

**Missing Raw Data**:
- Ensure local access to APOGEE raw data frames
- Set up proper filesystem mounts or Globus transfers

**Permission Issues**:
- Check file permissions for data directories
- Ensure write access for output files

### Getting Help

- Check [FAQ](faq.md) for common issues
- Open an issue on [GitHub](https://github.com/sdss/almanac/issues)
- Contact the development team