# Data Formats

This document describes the data models and formats used by `almanac`, which are built using Pydantic for type validation and structured data handling.

## Data Model Architecture

`almanac` uses three main Pydantic data models to represent SDSS survey data:

- **`Exposure`**: Metadata for individual exposures/observations
- **`FPSTarget`**: Target information for SDSS-V Fiber Positioning System observations
- **`PlateTarget`**: Target information for SDSS-IV plate-based observations

These models provide type validation, automatic conversion, and structured access to survey data.

## Exposure Model

The `Exposure` model (`almanac.data_models.Exposure`) represents the metadata for a single astronomical exposure.

### Basic Information

| Field | Type | Description |
|-------|------|-------------|
| `observatory` | `Observatory` | Observatory name ("apo" or "lco") |
| `mjd` | `int` | Modified Julian Date |
| `exposure` | `int` | Exposure number (≥ 1) |
| `prefix` | `Optional[Prefix]` | Raw exposure basename prefix ("apR", "asR") |

### Exposure Metadata

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `name` | `str` | | Plugged plate name or configuration identifier |
| `n_read` | `int` | `nread` | Number of detector reads (≥ 0) |
| `image_type` | `ImageType` | `imagetyp` | Image type (see ImageType literals) |
| `observer_comment` | `str` | `obscmnt` | Observer comments |

### Survey Identifiers

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `map_id` | `int` | `mapid` | Map identifier |
| `cart_id` | `int` | `cartid` | Cartridge ID |
| `plate_id` | `int` | `plateid` | Plate ID (SDSS-IV) |
| `field_id` | `int` | `fieldid` | Field identifier |
| `design_id` | `int` | `designid` | Design identifier |
| `config_id` | `int` | `configid` | Configuration ID (SDSS-V) |

### Observing Conditions

| Field | Type | Description |
|-------|------|-------------|
| `seeing` | `float` | Seeing in arcseconds |

### Instrument State

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `focus` | `float` | | Focus position |
| `collpist` | `float` | | Collimator piston |
| `colpitch` | `float` | | Collimator pitch |
| `dithered_pixels` | `float` | `dithpix` | Dither offset in pixels |
| `lamp_quartz` | `int` | `lampqrtz` | Quartz lamp state (-1, 0, 1) |
| `lamp_thar` | `int` | `lampthar` | ThAr lamp state (-1, 0, 1) |
| `lamp_une` | `int` | `lampune` | UNe lamp state (-1, 0, 1) |

### Computed Properties

| Property | Type | Description |
|----------|------|-------------|
| `exposure_string` | `str` | Formatted exposure string for file paths |
| `fps` | `bool` | Whether this is from the FPS era |
| `flagged_bad` | `bool` | Whether exposure is flagged as bad |
| `chip_flags` | `int` | Bitmask indicating available chips |

### Type Definitions

**`Observatory`**: `Literal["apo", "lco"]`

**`Prefix`**: `Literal["apR", "asR"]`

**`ImageType`**: `Literal["blackbody", "dark", "object", "domeflat", "arclamp", "twilightflat", "internalflat", "quartzflat", "missing"]`

## FPSTarget Model

The `FPSTarget` model (`almanac.data_models.FPSTarget`) represents targets observed with the SDSS-V Fiber Positioning System.

### Target Identification

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `sdss_id` | `Int64` | | SDSS unique identifier |
| `catalogid` | `Int64` | | Catalog identifier |
| `twomass_designation` | `str` | `tmass_id` | 2MASS designation |
| `category` | `Category` | | Target category |
| `cadence` | `str` | | Observing cadence |
| `firstcarton` | `str` | | Primary carton name |
| `program` | `str` | | Survey program |

### Positioner Information

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `positioner_id` | `int` | `positionerId` | Positioner ID |
| `hole_id` | `str` | `holeId` | Hole identifier |
| `hole_type` | `HoleType` | `holeType` | Hole type |
| `planned_hole_type` | `HoleType` | `holetype` | Planned hole type |
| `fiber_type` | `str` | `fiberType` | Fiber type |
| `assigned` | `bool` | | Target assigned to fiber |

### Status Flags

| Field | Type | Description |
|-------|------|-------------|
| `on_target` | `bool` | Fiber positioned on target |
| `disabled` | `bool` | Fiber is disabled |
| `valid` | `bool` | Valid coordinate conversion |
| `decollided` | `bool` | Positioner was decollided |

### Coordinates

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `ra` | `float` | `racat` | Right Ascension (degrees) |
| `dec` | `float` | `deccat` | Declination (degrees) |
| `alt` | `float` | `alt_observed` | Observed altitude (degrees) |
| `az` | `float` | `az_observed` | Observed azimuth (degrees) |

### Position Coordinates

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `x_wok` | `float` | `xwok` | X coordinate in wok frame |
| `y_wok` | `float` | `ywok` | Y coordinate in wok frame |
| `z_wok` | `float` | `zwok` | Z coordinate in wok frame |
| `x_focal` | `float` | `xFocal` | X coordinate in focal plane |
| `y_focal` | `float` | `yFocal` | Y coordinate in focal plane |

### Positioner Angles

| Field | Type | Description |
|-------|------|-------------|
| `alpha` | `float` | Alpha angle of positioner arm |
| `beta` | `float` | Beta angle of positioner arm |

### Wavelength Information

| Field | Type | Description |
|-------|------|-------------|
| `lambda_design` | `float` | Design wavelength |
| `lambda_eff` | `float` | Effective wavelength |
| `coord_epoch` | `float` | Coordinate epoch |

### Instrument Identifiers

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `spectrograph_id` | `int` | `spectrographId` | Spectrograph ID |
| `fiber_id` | `int` | `fiberId` | Fiber ID |

### Position Offsets

| Field | Type | Description |
|-------|------|-------------|
| `delta_ra` | `float` | RA offset applied to fiber |
| `delta_dec` | `float` | Dec offset applied to fiber |

### Target of Opportunity

| Field | Type | Description |
|-------|------|-------------|
| `too` | `bool` | Target of opportunity flag |
| `too_id` | `int` | TOO identifier |
| `too_program` | `str` | TOO program name |

## PlateTarget Model

The `PlateTarget` model (`almanac.data_models.PlateTarget`) represents targets observed with SDSS-IV plates.

### Target Identification

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `sdss_id` | `Int64` | | SDSS unique identifier |
| `catalogid` | `Int64` | | Catalog identifier |
| `twomass_designation` | `str` | | Computed 2MASS designation |
| `twomass_id` | `str` | `tmass_id` | 2MASS ID |
| `target_ids` | `str` | `targetids` | Target ID string |
| `category` | `Category` | `targettype` | Target category |

### Survey Information

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `observatory` | `Observatory` | | Observatory for mapping fixes |
| `hole_type` | `HoleType` | `holeType` | Hole type |
| `planned_hole_type` | `HoleType` | `holetype` | Planned hole type |
| `obj_type` | `ObjType` | `objType` | Object type |
| `assigned` | `bool` | | Assignment flag |

### Status Flags

| Field | Type | Description |
|-------|------|-------------|
| `conflicted` | `bool` | Fiber conflict flag |
| `ranout` | `bool` | Ran out of targets flag |
| `outside` | `bool` | Outside survey area flag |

### Coordinates

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `ra` | `float` | | Right Ascension (degrees) |
| `dec` | `float` | | Declination (degrees) |
| `x_focal` | `float` | `xfocal` | X focal plane coordinate |
| `y_focal` | `float` | `yfocal` | Y focal plane coordinate |
| `xf_default` | `float` | | Default X focal coordinate |
| `yf_default` | `float` | | Default Y focal coordinate |

### Wavelength Information

| Field | Type | Description |
|-------|------|-------------|
| `lambda_eff` | `float` | Effective wavelength |
| `zoffset` | `float` | Z offset |

### Instrument Identifiers

| Field | Type | Alias | Description |
|-------|------|-------|-------------|
| `spectrograph_id` | `int` | `spectrographId` | Spectrograph ID |
| `fiber_id` | `int` | `fiberId` | Fiber ID |
| `planned_fiber_id` | `int` | `fiberid` | Planned fiber ID |
| `throughput` | `int` | | Throughput value |

### Plate-Specific Fields

| Field | Type | Description |
|-------|------|-------------|
| `iplateinput` | `int` | Plate input ID |
| `pointing` | `int` | Pointing number |
| `offset` | `int` | Offset value |
| `block` | `int` | Block number |
| `iguide` | `int` | Guide flag |
| `bluefiber` | `int` | Blue fiber flag |
| `chunk` | `int` | Chunk number |
| `ifinal` | `int` | Final flag |
| `plugged_mjd` | `int` | MJD when plate was plugged |
| `fix_fiber_flag` | `int` | Fiber mapping fix indicator |

### Physical Properties

| Field | Type | Description |
|-------|------|-------------|
| `diameter` | `float` | Target diameter |
| `buffer` | `float` | Buffer size |
| `priority` | `int` | Target priority |

## Type Definitions

### Category Types

**`Category`**: `Literal["", "science", "sky_apogee", "sky_boss", "standard_apogee", "standard_boss", "open_fiber"]`

### Hole Types

**`HoleType`**: `Literal["object", "coherent_sky", "guide", "light_trap", "alignment", "quality", "manga", "manga_single", "manga_alignment", "acquisition_center", "acquisition_offaxis", "apogee", "center", "trap", "boss", "apogee_shared", "apogee_south", "bosshalf", "boss_shared", "fps"]`

### Object Types

**`ObjType`**: `Literal["galaxy", "qso", "star_bhb", "star_carbon", "star_brown_dwarf", "star_sub_dwarf", "star_caty_var", "star_red_dwarf", "star_white_dwarf", "redden_std", "spectrophoto_std", "hot_std", "rosat_a", "rosat_b", "rosat_c", "rosat_d", "serendipity_blue", "serendipity_first", "serendipity_red", "serendipity_distant", "serendipity_manual", "qa", "sky", "na"]`

## Data Validation and Processing

### Automatic Data Cleaning

The Pydantic models include automatic validation and cleaning:

- **Empty strings to integers**: Convert empty strings to default integer values
- **Lamp states**: Convert 'F'/'T' strings to 0/1 integers
- **Float validation**: Handle invalid values by setting to NaN
- **Case normalization**: Automatic lowercase conversion for categorical fields
- **Comment sanitization**: Clean observer comments and handle None values

### Model Validators

#### Exposure Model
- Validates prefixes based on observatory
- Sanitizes observer comments
- Converts lamp states from boolean strings
- Automatically detects twilight flats from comments

#### FPSTarget Model
- Fixes early fiber mapping duplicates for spectrograph 2
- Standardizes 2MASS designations

#### PlateTarget Model
- Applies historical fiber mapping fixes for specific MJD ranges
- Translates plate-era categories to FPS-era equivalents
- Handles observatory-specific fiber corrections

### Computed Properties

Models include computed properties that provide additional derived information:

- **Exposure strings**: Formatted identifiers for file paths
- **Era detection**: Automatic FPS/plate era classification
- **Quality flags**: Bad exposure detection
- **2MASS designations**: Standardized target identifiers

## HDF5 Output Structure

When saving data, `almanac` organizes the Pydantic model data into HDF5 format:

```
filename.h5
├── {observatory}/           # Observatory (apo/lco)
│   └── {mjd}/              # Modified Julian Date
│       ├── exposures       # Exposure model data
│       └── targets/        # Target model data
│           ├── fps/        # FPSTarget data by config_id
│           └── plates/     # PlateTarget data by plate_id
```

### Data Access

**Python Example**:
```python
from almanac.data_models import Exposure, FPSTarget, PlateTarget

# Create exposure from file headers
exposure = Exposure.from_path("/path/to/exposure/file")

# Access computed properties
print(f"FPS era: {exposure.fps}")
print(f"Exposure string: {exposure.exposure_string}")

# Access targets (automatically loaded)
for target in exposure.targets:
    if isinstance(target, FPSTarget):
        print(f"SDSS ID: {target.sdss_id}, Fiber: {target.fiber_id}")
    elif isinstance(target, PlateTarget):
        print(f"SDSS ID: {target.sdss_id}, Fiber: {target.fiber_id}")
```

## Data Quality and Validation

### Missing Value Handling

- **Numeric fields**: Default to appropriate values (NaN for floats, -1 for IDs)
- **String fields**: Default to empty strings
- **Boolean fields**: Default to False
- **Optional fields**: Can be None

### Type Safety

Pydantic ensures type safety throughout the data pipeline:
- Automatic type conversion where possible
- Validation errors for incompatible data
- Consistent data structures across the application

### Quality Assurance

- **Bad exposure flagging**: Automatic detection of known problematic exposures
- **Fiber mapping corrections**: Historical fixes for known fiber mapping errors
- **Coordinate validation**: Ensures valid sky coordinates
- **Identifier consistency**: Validates survey identifiers and relationships