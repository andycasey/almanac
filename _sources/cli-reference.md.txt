# CLI Reference

Complete command-line interface reference for `almanac`.

## Main Command

```bash
almanac [OPTIONS] [COMMAND]
```

The main `almanac` command queries APOGEE observational data. When run without subcommands, it performs data queries based on the specified options.

## Global Options

### Date and Time Selection

#### `--mjd <integer>`
Query specific Modified Julian Date.
- Positive values: Absolute MJD (e.g., `59300`)
- Negative values: Relative to current MJD (e.g., `-1` for yesterday)
- **Example**: `almanac --mjd 59300`, `almanac --mjd -7`

#### `--mjd-start <integer>`
Start of MJD range for queries.
- **Example**: `almanac --mjd-start 59300 --mjd-end 59310`

#### `--mjd-end <integer>`
End of MJD range for queries (inclusive).
- **Example**: `almanac --mjd-start -30 --mjd-end -1`

#### `--date <YYYY-MM-DD>`
Query specific calendar date (UTC).
- **Format**: ISO date format (YYYY-MM-DD)
- **Example**: `almanac --date 2021-01-01`

#### `--date-start <YYYY-MM-DD>`
Start of calendar date range.
- **Example**: `almanac --date-start 2021-01-01 --date-end 2021-01-31`

#### `--date-end <YYYY-MM-DD>`
End of calendar date range (inclusive).
- **Example**: `almanac --date-start 2024-01-01 --date-end 2024-12-31`

### Observatory Selection

#### `--apo`
Query Apache Point Observatory data only.
- **Example**: `almanac --apo --mjd -1`

#### `--lco`
Query Las Campanas Observatory data only.
- **Example**: `almanac --lco --date 2024-01-01`

### Data Options

#### `--fibers`, `--fibres`
Include fiber-to-target mappings in output.
- **Example**: `almanac --mjd 60000 --fibers`

#### `--no-x-match`
Skip cross-matching targets with SDSS database.
- Only effective when combined with `--fibers`
- Faster processing but less complete target information
- **Example**: `almanac --fibers --no-x-match --mjd-start -7`

### Output Control

#### `--output <path>`, `-O <path>`
Write output to HDF5 file at specified path.
- **Incremental**: Appends to existing files, preserves existing data
- **Example**: `almanac --output results.h5 --mjd-start -30`

#### `-v`, `--verbosity`
Control output verbosity (stackable).
- No flag: Minimal output
- `-v`: Show progress display
- `-vv`: Show progress display and exposure metadata
- **Example**: `almanac -vv --mjd -1`

### Performance Options

#### `--processes <integer>`, `-p <integer>`
Number of parallel processes for data processing.
- **Default**: Automatic based on available CPU cores
- **Example**: `almanac --processes 4 --mjd-start -30`


## Configuration Commands

### `almanac config show`
Display all current configuration settings and config file location.

```bash
almanac config show
```

**Output includes**:
- All configuration values  
- Configuration file path
- Database connection status

### `almanac config get <key>`
Retrieve specific configuration value.

```bash
almanac config get logging_level
almanac config get sdssdb.host
```

**Nested keys**: Use dot notation for nested configuration values.

### `almanac config set <key> <value>`
Set configuration value persistently.

```bash
almanac config set logging_level 10
almanac config set sdssdb.host custom-host.org
almanac config set database_connect_time_warning 5
```

**Data types**: Values are automatically converted to appropriate types.

## Usage Examples

### Basic Queries

```bash
# Today's observations from both observatories
almanac

# Yesterday's observations with details
almanac --mjd -1 -vv

# Specific observatory
almanac --apo --mjd -1
```

### Date Range Queries  

```bash
# Last week's data
almanac --mjd-start -7 --mjd-end -1

# Specific month
almanac --date-start 2024-01-01 --date-end 2024-01-31

# Single historical date
almanac --date 2021-06-15
```

### Fiber Mapping Analysis

```bash
# Include fiber mappings
almanac --mjd 60000 --fibers

# Fast fiber query (no cross-matching)
almanac --mjd 60000 --fibers --no-x-match

# Get fibers for specific MJD
almanac --mjd 60000 --fibers
```

### Output and Performance

```bash
# Save to file with progress
almanac --output survey.h5 --mjd-start -30 -v

# Parallel processing
almanac --processes 8 --mjd-start -7 --output week.h5

# Verbose output
almanac -vv --mjd-start -7
```

### Configuration Management

```bash
# View all settings
almanac config show

# Set debug logging
almanac config set logging_level 10

# Configure database
almanac config set sdssdb.host your-host.edu
almanac config set sdssdb.port 5432
```

## Advanced Usage Patterns

### Survey Monitoring

```bash
# Daily monitoring script
almanac --mjd -1 --output daily_$(date +%Y%m%d).h5 -v

# Weekly summary with fibers
almanac --mjd-start -7 --fibers --output weekly.h5 -vv
```

### Data Extraction

```bash
# Export specific time period for analysis
almanac --date-start 2024-01-01 --date-end 2024-03-31 \\
        --fibers --output q1_2024.h5 -v

# Export recent observations
almanac --mjd-start -30 --output exposures.h5
```

### Troubleshooting

```bash
# Debug mode with full verbosity
almanac config set logging_level 10
almanac -vv --mjd -1

# Test database connectivity
almanac config show
almanac --mjd -1 --no-x-match
```

## Exit Codes

- **0**: Success
- **1**: General error (invalid arguments, file errors)  
- **2**: Database connection error
- **3**: Data processing error

## Environment Variables

No special environment variables are required. Configuration is managed through the config file and command-line options.

## Shell Completion

For bash completion support, consider adding alias:

```bash
# Add to ~/.bashrc or ~/.bash_profile  
alias almanac='almanac'
complete -W '--mjd --mjd-start --mjd-end --date --date-start --date-end --apo --lco --fibers --fibres --no-x-match --output --processes --verbosity --help' almanac
```

## Tips

1. **Performance**: Use `--no-x-match` for faster queries when SDSS identifiers aren't needed
2. **Storage**: HDF5 files support incremental updates - reuse the same output file
3. **Debugging**: Use `-vv` and set `logging_level` to `10` for maximum detail
4. **Large queries**: Break large date ranges into smaller chunks for better performance
5. **Column selection**: Use custom column lists to reduce output size and processing time