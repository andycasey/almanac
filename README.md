## almanac
`almanac` scrapes headers from raw image files and cross-matches those against the SDSS database to create a comprehensive summary of everything ever observed with an APOGEE instrument.

## Getting Started

Here are a few example cases of how `almanac` might be helpful:

List all exposures taken yesterday from either telescope:
```bash
almanac --mjd -1 -vv
```

Or just from Apache Point Observatory:
```bash
almanac --mjd -1 -vv --apo
```

Write out all exposures taken in the last month to `january.h5`:
```bash
almanac -o january.h5 --mjd-start -30
```

Write out all fiber observations during 2021, where we switched from plates to robotic fiber positioners:
```bash
almanac -vv --date-start 2021-01-01 --date-end 2021-31-12 --fibers
```

And it looks pretty, even when it warns you about missing exposures:

![](https://github.com/sdss/almanac/blob/83159e03632e3edbb45bb0c8de9810dec2dc49f1/docs/almanac-example-1.gif)


## Installation


### At Utah

If you want to use this at Utah, you can install it with:

```bash
module load almanac
```

> [!TIP]
> We recommend you manage your own Python environment, but if you don't have one set up at Utah then you can use `module load miniconda/3.8.5_astra`. 

### Anywhere else

`almanac` needs local disk access to raw APOGEE data frames. If you are going to run it somewhere else, you should set up a Globus transfer of raw APOGEE frames, and ensure your internet address is whitelisted to remotely access the SDSS database.

We recommend using `uv` to manage Python environments. Using `uv`, you can install `almanac` with:
```bash
uv add sdss-almanac
```

## Usage

Use `almanac` to see details on data taken today from both observatories, or specify the observatory:

```bash
almanac
almanac --apo # Apache Point Observatory
almanac --lco # Las Campanas Observatory
```

### Specifying a date

If you want a particular day, either use the ``--mjd`` or ``--date`` (UTC) flags:

```bash
almanac --mjd 59300
almanac --date 2021-01-01
```

You can use negative MJD values to indicate days relative to today:

```bash
almanac --mjd -1 # Yesterday
almanac --mjd -7 # Last week
```

You can also specify a range of days:

```bash
almanac --mjd-start 59300 --mjd-end 59310 # Give me these 10 days
almanac --date-start 2021-01-01 --date-end 2021-01-31 # Give me all of January 2021
```

### Fiber mappings

You can also use `almanac` to see the fiber mappings for a given plate (SDSS-IV) or FPS pointing (SDSS-V) by specifing the ``--fibers`` (or ``--fibres``) flag. This will give you the mapping of fibers to targets, and the target properties. 

```bash
almanac --mjd 60000 --fibres
```

The fiber mapping tables are cross-matched to the SDSS database to include the SDSS identifiers for each target. If you don't want to do this cross-match, you can use the ``--no-x-match`` flag. The ``--no-x-match`` flag is ignored if ``--fibers`` is not used.

### Verbosity

By default there is minimal output to the terminal. You can adjust the verbosity level using `-v`:
- `-v`: show progress display only
- `-vv`: show progress display and exposure metadata

In verbose mode you can see exposure information in the terminal, and additional per-fiber metadata is stored in the HDF5 files that `almanac` creates.

![](https://github.com/sdss/almanac/blob/e3f46c8ce66b88843de943ca31eec88d12be8f06/docs/almanac-example-2.gif)

### Outputs

You can write the outputs to a structured HDF5 file by specifying an output path with the ``--output`` (or ``-O``) flag. If the output path already exists, the default behaviour is to overwrite existing entries *only*. So if you run `almanac` once for MJD 60000 and output to a file, and then run it again for MJD 60001 and output to the same file, your file will have data for both MJDs. 

```bash
almanac --output /path/to/file.h5 # Append today's data to existing file
```

An example structure of the HDF5 file is below:

```
raw/apo/59300/exposures        # a data table of exposures
raw/apo/59300/sequences        # a Nx2 array of exposure numbers (inclusive) that form a sequence
raw/apo/59300/fibers/fps/1     # a data table of fiber mappings for FPS configuration id 1
raw/apo/59300/fibers/plates/2  # a data table of fiber mappings for plate id 2
```

## Configuration

You can view and change the `almanac` configuration settings through the `almanac config` interface. To view all current settings and to see the configuration file path:

```bash
almanac config show
```

### To get a single configuration value
```bash
almanac config get logging_level
```

### To set a configuration value
```bash
almanac config set logging_level 10
```
