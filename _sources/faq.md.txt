# Frequently Asked Questions

Common questions and troubleshooting for `almanac`.

## Installation and Setup

### Q: I'm getting "command not found" after installation

**A:** This usually means the installation directory isn't in your PATH, or there was an installation issue.

**Solutions:**
```bash
# Check if almanac is installed
python -c "import almanac; print(almanac.__version__)"

# If installed but not in PATH, find it:
python -c "import almanac; print(almanac.__file__)"

# Try reinstalling:
uv pip install --force-reinstall git+https://github.com/sdss/almanac

# Or add to PATH (example path):
export PATH="$HOME/.local/bin:$PATH"
```

### Q: How do I install at Utah?

**A:** At Utah, use the module system:
```bash
module load almanac
```

If you need a specific Python environment:
```bash
module load miniconda/3.8.5_astra
module load almanac
```

### Q: The installation fails with dependency conflicts

**A:** Try using `uv` for better dependency resolution:
```bash
# Install uv first
pip install uv

# Use uv to install almanac
uv pip install git+https://github.com/sdss/almanac
```

## Database Connection Issues

### Q: "Unable to connect to SDSS database" error

**A:** This is a common issue outside Utah. Check these solutions:

1. **Verify your IP is whitelisted** for SDSS database access
2. **Check network connectivity**:
   ```bash
   ping operations.sdss.org
   ```
3. **Test with reduced functionality**:
   ```bash
   almanac --mjd -1 --no-x-match
   ```
4. **Check configuration**:
   ```bash
   almanac config show
   almanac config get sdssdb.host
   ```

### Q: Database connection is very slow

**A:** Adjust the timeout warning:
```bash
almanac config set database_connect_time_warning 10
```

Or check if you're using the correct database host:
```bash
almanac config set sdssdb.host closer-database-mirror.org
```

### Q: Can I use almanac without database access?

**A:** Limited functionality is available without database access:
```bash
# This will work but with limited target information
almanac --mjd -1 --no-x-match

# Raw exposure data doesn't require database
almanac --mjd -1
```

## Data and File Access

### Q: "No raw data found" or similar errors

**A:** `almanac` needs access to raw APOGEE data files:

1. **At Utah**: Data should be automatically available
2. **Elsewhere**: Set up Globus transfer or mount the data filesystem
3. **Check data paths**: Verify the data directories are accessible

### Q: How do I access data from many years ago?

**A:** Historical data availability depends on:
- **Data archive access**: Older data may need special access
- **Database coverage**: Very old data may not be in the current database
- **MJD limits**: Check minimum MJD settings:
  ```bash
  almanac config get sdssdb_exposure_min_mjd.apo
  almanac config get sdssdb_exposure_min_mjd.lco
  ```

### Q: Output files are very large

**A:** Control file size with:
```bash
# Save to smaller file
almanac --mjd -1 --output smaller.h5

# Skip fiber data if not needed
almanac --mjd-start -30 --output no-fibers.h5

# Use date ranges instead of large ranges
almanac --date-start 2024-01-01 --date-end 2024-01-31 --output jan.h5
```

## Usage and Commands

### Q: What's the difference between --mjd and --date?

**A:** 
- **--mjd**: Uses Modified Julian Date (astronomical standard)
  - Example: `--mjd 59300`
  - Supports relative values: `--mjd -1` (yesterday)
- **--date**: Uses calendar dates in UTC
  - Example: `--date 2021-01-01`
  - More intuitive for specific dates

### Q: How do I query a specific month or year?

**A:**
```bash
# Specific month
almanac --date-start 2024-01-01 --date-end 2024-01-31

# Entire year
almanac --date-start 2024-01-01 --date-end 2024-12-31

# Using MJD (if you know them)
almanac --mjd-start 59300 --mjd-end 59330
```

### Q: What do the verbosity levels do?

**A:**
- **No flag**: Minimal output (just summary)
- **-v**: Shows progress bars and basic information
- **-vv**: Shows progress bars plus detailed exposure tables

```bash
almanac --mjd -1      # Just summary
almanac --mjd -1 -v   # With progress
almanac --mjd -1 -vv  # With full details
```

### Q: How do I save data to a file?

**A:**
```bash
# Save to HDF5 file
almanac --mjd -1 --output results.h5

# Files are incremental - you can append:
almanac --mjd -2 --output results.h5  # Adds to same file
```

## Data Interpretation

### Q: What does "missing exposures" mean?

**A:** `almanac` detects gaps in exposure sequences. This could indicate:
- **Technical issues**: Instrument problems during observations
- **Weather**: Observations stopped due to poor conditions  
- **Scheduled breaks**: Intentional gaps in observing
- **Data processing**: Exposures not yet processed or archived

### Q: What are "sequences"?

**A:** Sequences are groups of consecutive exposures that form logical observational units, typically:
- **Science sequences**: Multiple exposures of the same target/field
- **Calibration sequences**: Sets of calibration exposures
- **Standard star sequences**: Observations of standard stars

### Q: What's the difference between plates and FPS?

**A:**
- **Plates (SDSS-IV)**: Fixed fiber positions drilled into metal plates
- **FPS (SDSS-V)**: Robotic fiber positioners that can be reconfigured

Use different column options for each:
```bash
# Get all available data for plate-based and FPS observations
almanac --fibers --mjd -1
```

## Performance and Optimization

### Q: almanac is running slowly

**A:** Try these optimizations:

1. **Use parallel processing**:
   ```bash
   almanac --processes 4 --mjd-start -7
   ```

2. **Skip cross-matching when not needed**:
   ```bash
   almanac --mjd -1 --fibers --no-x-match
   ```

3. **Limit date ranges**:
   ```bash
   # Instead of a full year, try monthly chunks
   almanac --date-start 2024-01-01 --date-end 2024-01-31
   ```

4. **Select only needed columns**:
   ```bash
   almanac --mjd-start -7
   ```

### Q: How many processes should I use?

**A:** Generally:
- **Local machine**: Number of CPU cores (check with `nproc`)
- **Shared systems**: Be conservative (2-4 processes)
- **High I/O load**: Fewer processes may be faster

```bash
# Check CPU cores
nproc

# Use half your cores
almanac --processes $(expr $(nproc) / 2) --mjd-start -7
```

## Configuration Issues

### Q: How do I reset my configuration?

**A:**
```bash
# Find config file location
almanac config show

# Remove config file (it will be recreated with defaults)
rm ~/.config/almanac/config.yaml
```

### Q: My configuration changes aren't taking effect

**A:**
1. **Check the config file location**:
   ```bash
   almanac config show
   ```

2. **Verify the setting was saved**:
   ```bash
   almanac config get your_setting_name
   ```

3. **Check for YAML syntax errors** in the config file

4. **Restart your shell** after major configuration changes

### Q: Can I have different configurations for different projects?

**A:** Yes, by managing config files manually:
```bash
# Save current config
cp ~/.config/almanac/config.yaml ~/project1-config.yaml

# Switch configurations
cp ~/project2-config.yaml ~/.config/almanac/config.yaml
```

## Error Messages

### Q: "Permission denied" errors

**A:**
- **File permissions**: Check write permissions for output directories
- **Data access**: Verify access to raw data directories
- **Config file**: Check write permissions for config directory

### Q: "Module not found" errors

**A:** Usually indicates incomplete installation:
```bash
# Reinstall in development mode
pip install -e .

# Or regular install
pip install git+https://github.com/sdss/almanac
```

### Q: Memory errors with large queries

**A:**
- **Reduce date range**: Query smaller time periods
- **Skip fiber data**: Remove `--fibers` flag if not needed
- **Use streaming**: Process data in smaller chunks
- **Increase swap space**: If on your own system

## Data Quality and Validation

### Q: How do I verify my results are correct?

**A:**
1. **Cross-check with known observations**:
   ```bash
   almanac --mjd 60000 -vv  # Use a known good MJD
   ```

2. **Compare observatories**:
   ```bash
   almanac --apo --mjd -1
   almanac --lco --mjd -1
   ```

3. **Check sequence patterns**:
   ```bash
   almanac --mjd -1
   ```

### Q: Some targets show "N/A" or missing information

**A:** This is normal when:
- **Cross-matching disabled**: Use without `--no-x-match`
- **Database incomplete**: Some targets may not be in the database
- **Historical data**: Older observations may have less complete information

## Getting Help

### Q: None of these solutions work for my issue

**A:** Get help through:

1. **GitHub Issues**: [Report bugs or request features](https://github.com/sdss/almanac/issues)
2. **Verbose output**: Run with `-vv` and include output in your report
3. **System information**: Include your OS, Python version, and installation method
4. **Configuration**: Include output from `almanac config show`

### Q: How do I report a bug effectively?

**A:** Include in your bug report:
```bash
# System information
almanac --version  # (if available)
python --version
uname -a

# Configuration
almanac config show

# Error reproduction
almanac --mjd -1 -vv  # Whatever command fails
```

### Q: Is there a user community or forum?

**A:** 
- **GitHub Discussions**: For questions and community discussions
- **SDSS Collaboration**: Internal SDSS communication channels
- **Documentation**: This documentation and the main README