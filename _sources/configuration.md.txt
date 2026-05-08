# Configuration Guide

This guide covers all configuration options available in `almanac`.

## Configuration Overview

`almanac` uses a hierarchical configuration system with:
- **Default values** built into the code
- **Configuration file** for persistent settings
- **Command-line options** for one-time overrides

## Configuration File

### Location

The configuration file is automatically created at:
- **Linux/macOS**: `~/.config/almanac/config.yaml`
- **Windows**: `%APPDATA%/almanac/config.yaml`

### Viewing Configuration

```bash
# Show all current settings and config file path
almanac config show

# Get a specific setting
almanac config get logging_level
```

### Setting Values

```bash
# Set a configuration value
almanac config set logging_level 10

# Set database connection parameters
almanac config set sdssdb.host your-database-host.org
almanac config set sdssdb.port 5432
```

## Configuration Reference

### Database Settings

#### `sdssdb.user`
- **Type**: String
- **Default**: `"sdss_user"`
- **Description**: Database username for SDSS database connection

#### `sdssdb.host`
- **Type**: String  
- **Default**: `"operations.sdss.org"`
- **Description**: Database host address

#### `sdssdb.port`
- **Type**: Integer
- **Default**: `5432`
- **Description**: Database port number

#### `sdssdb.domain`
- **Type**: String
- **Default**: `"operations.sdss.*"`
- **Description**: Database domain pattern for connection

### Performance Settings

#### `database_connect_time_warning`
- **Type**: Integer
- **Default**: `3`
- **Units**: Seconds
- **Description**: Time threshold for database connection warnings

### Data Processing Settings

#### `sdssdb_exposure_min_mjd.apo`
- **Type**: Integer
- **Default**: `59558`
- **Description**: Minimum MJD for Apache Point Observatory exposure queries

#### `sdssdb_exposure_min_mjd.lco`
- **Type**: Integer
- **Default**: `59558` 
- **Description**: Minimum MJD for Las Campanas Observatory exposure queries

### Logging Settings

#### `logging_level`
- **Type**: Integer
- **Default**: `20` (INFO level)
- **Description**: Python logging level
- **Values**:
  - `10`: DEBUG (most verbose)
  - `20`: INFO (default)
  - `30`: WARNING
  - `40`: ERROR
  - `50`: CRITICAL (least verbose)

## Configuration File Format

The configuration file uses YAML format:

```yaml
# Example almanac configuration
database_connect_time_warning: 3
logging_level: 20

sdssdb:
  user: sdss_user
  host: operations.sdss.org
  port: 5432
  domain: operations.sdss.*

sdssdb_exposure_min_mjd:
  apo: 59558
  lco: 59558
```

## Environment-Specific Configuration

### Utah Computing Environment

At Utah, the default configuration typically works without modification:
```bash
# Verify default settings work
almanac config show
almanac --mjd -1  # Test query
```

### External Installations

For installations outside Utah, configure database access:

```bash
# Set custom database host (if different)
almanac config set sdssdb.host your-database-host.edu

# Adjust connection timeout warning if needed
almanac config set database_connect_time_warning 10
```

### Development Environment

For development work, you might want more verbose logging:

```bash
# Enable DEBUG-level logging
almanac config set logging_level 10

# Reduce connection time warning for local development
almanac config set database_connect_time_warning 1
```

## Configuration Management

### Backup Configuration

```bash
# Find config file location
almanac config show | grep "Config file"

# Backup the file
cp ~/.config/almanac/config.yaml ~/almanac-config-backup.yaml
```

### Reset to Defaults

To reset all settings to defaults, simply delete the configuration file:

```bash
rm ~/.config/almanac/config.yaml
```

The file will be recreated with default values on the next run.

### Multiple Configurations

For different environments, you can manage multiple config files:

```bash
# Save current config
cp ~/.config/almanac/config.yaml ~/configs/almanac-utah.yaml

# Create test config
cp ~/configs/almanac-test.yaml ~/.config/almanac/config.yaml
```

## Command-Line Overrides

Many settings can be overridden on the command line:

```bash
# Override verbosity (affects logging output)
almanac -vv --mjd -1

# Override parallel processing
almanac --processes 8 --mjd-start -7
```

Note that command-line overrides don't modify the configuration file.

## Troubleshooting Configuration

### Common Issues

**Configuration Not Taking Effect**:
- Check configuration file path with `almanac config show`
- Verify YAML syntax is correct
- Ensure proper data types (numbers vs. strings)

**Database Connection Issues**:
- Verify `sdssdb.host` and `sdssdb.port` settings
- Check network connectivity
- Confirm IP address is whitelisted for database access

**Permission Errors**:
- Check write permissions for config directory
- Verify file ownership

### Validation

```bash
# Test configuration with a simple query
almanac config show
almanac --mjd -1 -v

# Test database connectivity
almanac config get sdssdb.host
```

## Advanced Configuration

### Custom Logging

For more advanced logging configuration, you can modify the logging level:

```bash
# Very verbose (DEBUG level)
almanac config set logging_level 10

# Quiet (WARNING level only)
almanac config set logging_level 30
```

### Performance Tuning

```bash
# Adjust database timeout warning for slower connections
almanac config set database_connect_time_warning 30

# Set historical data limits if needed
almanac config set sdssdb_exposure_min_mjd.apo 60000
```

## Configuration Schema

The configuration follows this schema structure:

```python
Config:
  sdssdb: DatabaseConfig
    user: str
    host: str  
    port: int
    domain: str
  database_connect_time_warning: int
  sdssdb_exposure_min_mjd: ObservatoryMJD
    apo: int
    lco: int  
  logging_level: int
```

All settings are optional and will fall back to defaults if not specified.