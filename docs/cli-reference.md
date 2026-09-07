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

#### `--mjds <list>`
Comma-separated list of explicit MJDs to query. The list does not need to be contiguous.
- Useful for testbeds, reprocessing selected nights, or filling gaps
- Also accepted by `almanac add metadata`
- **Example**: `almanac --mjds 58011,59337,60000`

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

#### `--skip-existing`
Skip observatory/MJD pairs that already exist in the output file (resume mode).
- Requires `--output`; pairs with an existing `raw/{observatory}/{mjd}/exposures` group are skipped
- An interrupted or partially failed run can simply be re-run with the same command: only the missing pairs are processed
- Existing `/missing_exposures` entries for skipped pairs are preserved
- **Example**: `almanac --mjd-start 59300 --mjd-end 60300 --output big.h5 --skip-existing`

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

### Fault Tolerance

Long queries are fault-isolated per observatory/MJD pair:

- **Retry with backoff**: Transient database connection errors are retried automatically (see the `database_retry_attempts` and `database_retry_backoff` [configuration settings](configuration.md))
- **Per-night isolation**: A night that fails outright is recorded and skipped rather than aborting the whole run; failed pairs are listed at the end
- **Resume**: Re-run the same command with `--skip-existing` to retry only the failed pairs
- **Missing-exposure report**: Every detected missing exposure is recorded in the `/missing_exposures` table of the output file with a reason code (see [Data Formats](data-formats.md))

## Adding to an Existing File

### `almanac add fibers <file>`
Add fiber-to-target mappings to an existing almanac HDF5 file that was created without `--fibers`. The exposures already in the file are used to locate the confSummary (FPS era) or plugmap (plate era) files, so the raw exposure headers are not read again. Mappings are written to `raw/<observatory>/<mjd>/fibers/<config_id or plate_id>`; a night that already has mappings has them replaced. Also accepted as `almanac add fibres`.

```bash
almanac add fibers results.h5
almanac add fibers results.h5 --mjd 60000 --apo
almanac add fibers results.h5 --no-x-match
almanac add fibers results.h5 --processes 8
```

**Options**:
- Date/MJD selection: `--mjd`, `--mjds`, `--mjd-start`/`--mjd-end`, `--date`, `--date-start`/`--date-end` (defaults to every MJD present in the file)
- Observatory selection: `--apo`, `--lco`
- `--no-x-match`: Do not cross-match targets with the SDSS-V database (no `sdss_id` values are assigned)
- `--processes`, `-p <integer>`: Number of processes to use (default: serial; negative values use every CPU)

A night that fails (for example, a missing confSummary file) is recorded and skipped rather than aborting the run; failed nights are listed at the end.

### `almanac add metadata <file>`
Decorate an existing almanac HDF5 file with astrometry, photometry, and targeting flags for every cross-matched target. Results are written to a `meta/` group in the same file, queried once per unique `sdss_id`. The file must already contain fiber mappings (created with `--fibers`, or added with `almanac add fibers`); if none are found for the selected nights, the command warns and exits with a non-zero status without touching the file.

```bash
almanac add metadata results.h5
almanac add metadata results.h5 --mjds 58011,59337 --apo
almanac add metadata results.h5 --query-workers 4
```

**Options**:
- Date/MJD selection: `--mjd`, `--mjds`, `--mjd-start`/`--mjd-end`, `--date`, `--date-start`/`--date-end` (defaults to every MJD present in the file)
- Observatory selection: `--apo`, `--lco`
- `--p <integer>`: Number of workers used to read `sdss_id` values from the file
- `--query-workers <integer>`: Number of parallel database query workers (default: the `catalog_query_max_workers` configuration setting). Each worker opens its own database connection; keep this low when connecting through a single SSH tunnel.

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

# Non-contiguous set of nights
almanac --mjds 58011,59337,60000 --fibers --output testbed.h5
```

### Resumable Bulk Generation

```bash
# Large multi-year run; safe to interrupt and re-run
almanac --mjd-start 59300 --mjd-end 60300 --fibers \
        --output survey.h5 --skip-existing -v

# Then decorate with catalog metadata
almanac add metadata survey.h5 --query-workers 4

# Review what was missing before downstream processing
# (see the /missing_exposures table in survey.h5)
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
complete -W '--mjd --mjds --mjd-start --mjd-end --date --date-start --date-end --apo --lco --fibers --fibres --no-x-match --output --skip-existing --processes --verbosity --help' almanac
```

## Tips

1. **Performance**: Use `--no-x-match` for faster queries when SDSS identifiers aren't needed
2. **Storage**: HDF5 files support incremental updates - reuse the same output file
3. **Debugging**: Use `-vv` and set `logging_level` to `10` for maximum detail
4. **Large queries**: Break large date ranges into smaller chunks for better performance
5. **Column selection**: Use custom column lists to reduce output size and processing time