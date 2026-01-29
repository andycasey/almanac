# User Guide

This comprehensive guide covers all aspects of using `almanac` to query and analyze APOGEE observational data.

## Overview

`almanac` scrapes headers from raw APOGEE image files and cross-matches them against the SDSS database to provide comprehensive summaries of observational data. It supports both SDSS-IV (plate-based) and SDSS-V (fiber positioner-based) observations.

## Basic Usage

### Default Behavior

Run `almanac` without arguments to see today's observations from both observatories:

```bash
almanac
```

### Observatory Selection

Query specific observatories:

```bash
almanac --apo    # Apache Point Observatory only
almanac --lco    # Las Campanas Observatory only
```

## Date and Time Queries

### Single Date Queries

**Using MJD (Modified Julian Date)**:
```bash
almanac --mjd 59300        # Specific MJD
almanac --mjd -1           # Yesterday
almanac --mjd -7           # One week ago
```

**Using Calendar Dates** (UTC):
```bash
almanac --date 2021-01-01  # Specific date
almanac --date 2024-12-25  # Christmas Day 2024
```

### Date Range Queries

**MJD Ranges**:
```bash
almanac --mjd-start 59300 --mjd-end 59310    # 10-day range
almanac --mjd-start -30                       # Last 30 days
```

**Calendar Date Ranges**:
```bash
almanac --date-start 2021-01-01 --date-end 2021-01-31  # January 2021
almanac --date-start 2021-01-01 --date-end 2021-12-31  # All of 2021
```

## Fiber Mapping Analysis

### Basic Fiber Information

Include fiber-to-target mappings:

```bash
almanac --mjd 60000 --fibers    # or --fibres
```

This provides:
- Target assignments to fibers
- Target properties and coordinates
- SDSS identifiers (cross-matched)

### Cross-matching Control

Skip database cross-matching for faster results:

```bash
almanac --mjd 60000 --fibers --no-x-match
```

## Output and Verbosity

### Verbosity Levels

Control output detail:

```bash
almanac                    # Minimal output
almanac -v                 # Show progress display
almanac -vv                # Show progress + exposure metadata
```

### Output Files

Save results to structured HDF5 files:

```bash
almanac --output results.h5 --mjd-start -7  # Last week to HDF5
almanac -O monthly.h5 --date-start 2024-01-01 --date-end 2024-01-31
```

**Incremental Updates**: Running `almanac` multiple times with the same output file appends new data while preserving existing entries.

## Advanced Features

### Data Access

The CLI provides access to all available data columns automatically. Column selection and filtering are handled through the package's data model, providing a complete view of survey operations.

### Performance Tuning

**Parallel Processing**:
```bash
almanac --processes 4      # Use 4 CPU cores
almanac -p 8               # Use 8 CPU cores
```

## Data Structure Examples

### Exposure Data

Typical exposure information includes:
- **Basic**: Observatory, MJD, exposure number, type
- **Configuration**: Design ID, field ID, cart ID, config ID
- **Instrumentation**: Number of reads, lamp states, dither pixels
- **Conditions**: Seeing, focus, comments

### Fiber Mapping Data

**FPS (SDSS-V) Data**:
- SDSS identifier and catalog ID
- Program and carton information
- Target coordinates (RA, Dec)
- Fiber assignments

**Plate (SDSS-IV) Data**:
- Target identifiers and coordinates
- Target and source types
- Fiber assignments

## Common Workflows

### Daily Operations Check

```bash
# Check today's observations with details
almanac -vv

# Check specific observatory
almanac --apo -vv
```

### Survey Planning and Review

```bash
# Review recent survey progress
almanac --mjd-start -30 --fibers -O recent_survey.h5

# Analyze specific survey period
almanac --date-start 2024-01-01 --date-end 2024-03-31 --fibers
```

### Data Quality Assessment

```bash
# Check recent exposures with full output
almanac --mjd-start -7 -vv

# Focus on yesterday's data
almanac --mjd-start -1 -v
```

### Export for Analysis

```bash
# Create comprehensive dataset for analysis
almanac --date-start 2024-01-01 --date-end 2024-12-31 \\
        --fibers --output 2024_survey.h5 -v

# Quick daily summary
almanac --mjd -1 --output daily_$(date +%Y%m%d).h5
```

## Tips and Best Practices

### Performance

- Use `--no-x-match` when cross-matching isn't needed
- Adjust `--processes` based on available CPU cores
- For large date ranges, consider breaking into smaller chunks

### Data Management

- Use descriptive output filenames with dates
- Leverage incremental updates for ongoing monitoring
- Regularly backup important datasets

### Troubleshooting

- Use `-vv` for detailed progress information
- Check database connectivity with `almanac config show`
- Verify data paths and permissions

## Next Steps

- See [Configuration Guide](configuration.md) for customization options
- Check [CLI Reference](cli-reference.md) for complete command documentation
- Visit [API Reference](api-reference.md) for programmatic usage
