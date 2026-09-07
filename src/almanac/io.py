import h5py as h5
import numpy as np
from tqdm import tqdm
from typing import List, Literal, Dict, Any, Tuple, Union, get_origin, get_args
from pydantic import BaseModel
from pydantic.fields import (FieldInfo, ComputedFieldInfo)
from pydantic_core import PydanticUndefined
import datetime
from enum import Enum

from almanac.data_models import Exposure

def write_almanac(
    output: str,
    results: List[Tuple[str, int, List[Exposure], Dict[str, List[Any]]]],
    fibers: bool = False,
    verbose: bool = False,
    compression: bool = True
):
    """
    Write the results of an Almanac query to an HDF5 file.

    :param output:
        Path to the output HDF5 file.

    :param results:
        List of tuples containing (observatory, mjd, exposures, sequences).
        - observatory: str, e.g., "apo" or "lco"
        - mjd: int, Modified Julian Date
        - exposures: List[Exposure], list of Exposure models
        - sequences: Dict[str, List[Any]], dictionary of sequences by image type

    :param fibers:
        Whether to include fiber data in the output.

    :param verbose:
        Whether to print progress information.

    :param compression:
        Compression algorithm to use for datasets. If True, uses 'gzip'.
    """

    kwds = dict(fibers=fibers, verbose=verbose, compression=compression)
    with h5.File(output, "a") as fp:
        for args in sorted(results, key=lambda x: (x[0], x[1])):
            update(fp, *args, **kwds)

def update(
    fp,
    observatory,
    mjd,
    exposures,
    sequences,
    fibers: bool = False,
    verbose: bool = False,
    compression: Union[bool, str] = True
):
    _print = print if verbose else lambda *args, **kwargs: None

    group = get_or_create_group(fp, f"raw/{observatory}/{mjd}")
    _print(f"\traw/{observatory}/{mjd}")

    delete_hdf5_entry(group, "exposures")
    write_models_to_hdf5_group(
        exposures,
        group.create_group("exposures", track_order=True)
    )

    _print(f"\traw/{observatory}/{mjd}/exposures")

    if len(sequences) > 0:
        delete_hdf5_entry(group, "sequences")
        sequences_group = group.create_group("sequences")
        for image_type, entries in sequences.items():
            # Guard empty sequences: an empty np.array() defaults to float64
            # with shape (0,); all sequences are (start, end) integer pairs.
            data = (
                np.array(entries, dtype=np.int64)
                if len(entries) else
                np.zeros((0, 2), dtype=np.int64)
            )
            sequences_group.create_dataset(image_type, data=data)
            _print(f"\traw/{observatory}/{mjd}/sequences/{image_type}")

    if fibers:
        write_fibers(fp, observatory, mjd, exposures, verbose=verbose)


def write_fibers(fp, observatory, mjd, exposures, verbose: bool = False):
    """
    Write the fiber-to-target mappings for one observatory and MJD to
    `raw/{observatory}/{mjd}/fibers/{config_id or plate_id}`, one group per
    configuration (FPS era) or plate (plate era). Existing groups for the same
    configuration are replaced.

    Only the first exposure of each configuration is consulted, and only
    exposures whose `targets` have been populated (see
    `almanac.apogee.cross_match_targets`) are written.
    """
    _print = print if verbose else lambda *args, **kwargs: None

    fibers_group = get_or_create_group(fp, f"raw/{observatory}/{mjd}/fibers")
    done = set()
    for exposure in exposures:
        reference_id_string = str(
            exposure.config_id if exposure.fps else exposure.plate_id
        )
        if reference_id_string in done:
            # Check `done` BEFORE touching `exposure.targets`: `targets`
            # is a lazy property, and for exposures whose configuration
            # has already been written, accessing it would re-parse the
            # confSummary/plugmap yanny files only to discard the result.
            # Workers pre-compute (and pickle) targets for the exposures
            # that are actually written, so this ordering keeps the
            # parent process from redoing ~seconds of parsing per night.
            continue

        if not exposure.targets:
            continue

        delete_hdf5_entry(fibers_group, reference_id_string)
        write_models_to_hdf5_group(
            exposure.targets,
            fibers_group.create_group(reference_id_string, track_order=True)
        )
        done.add(reference_id_string)
        _print(f"\traw/{observatory}/{mjd}/fibers/{reference_id_string}")


def read_exposures(fp, observatory, mjd) -> List[Exposure]:
    """
    Reconstruct the `Exposure` models for one observatory and MJD from an
    existing almanac file, in exposure-number order.

    This reads only the `raw/{observatory}/{mjd}/exposures` group; no raw
    exposure headers are touched. Returns an empty list if the group does not
    exist or is empty.
    """
    path = f"raw/{observatory}/{mjd}/exposures"
    if path not in fp:
        return []
    group = fp[path]

    columns = {
        name: group[name][:]
        for name in Exposure.model_fields
        if name in group
    }
    if not columns:
        return []
    num_records = len(next(iter(columns.values())))

    exposures = []
    for i in range(num_records):
        kwds = {}
        for name, values in columns.items():
            value = values[i]
            if isinstance(value, bytes):
                value = value.decode("utf-8")
                if value in ("", "None"):
                    # Unset optional string (stored as an empty or "None"
                    # byte string): let the model default apply.
                    continue
            elif isinstance(value, np.generic):
                value = value.item()
            kwds[name] = value
        kwds.setdefault("prefix", dict(apo="apR", lco="asR").get(kwds["observatory"]))
        exposures.append(Exposure(**kwds))
    return exposures


def read_sequences(fp, observatory, mjd) -> Dict[str, List[Tuple[int, int]]]:
    """
    Read the exposure sequences (e.g. "objects", "arclamps", "missing") for
    one observatory and MJD from an existing almanac file, as a dictionary of
    (start, end) pairs of 1-indexed exposure numbers. Returns an empty
    dictionary if the sequences group does not exist.
    """
    path = f"raw/{observatory}/{mjd}/sequences"
    sequences = {}
    if path in fp:
        for image_type in fp[path]:
            sequences[image_type] = [
                (int(start), int(end)) for start, end in fp[path][image_type][:]
            ]
    return sequences



MISSING_EXPOSURES_GROUP = "missing_exposures"

_MISSING_EXPOSURES_DESCRIPTIONS = {
    "observatory": "Observatory name",
    "mjd": "MJD of the night",
    "exposure": (
        "Within-night exposure number (1-indexed); "
        "-1 for observatory/MJD-level records (e.g. reason=db_unavailable)"
    ),
    "expected_max_db": (
        "Maximum within-night exposure number expected from the operations "
        "database (-1 if unknown)"
    ),
    "reason": (
        "Why the exposure is flagged: one of "
        "hole | trailing | db_unavailable | db_no_file | file_no_db"
    ),
}


def read_missing_exposures(fp) -> List[Dict[str, Any]]:
    """
    Read the run-level /missing_exposures table from an almanac HDF5 file.

    :param fp:
        An open h5py File (or Group).

    :returns:
        A list of records (dicts) with keys: observatory, mjd, exposure,
        expected_max_db, reason. Empty if the table is absent.
    """
    if MISSING_EXPOSURES_GROUP not in fp:
        return []
    group = fp[MISSING_EXPOSURES_GROUP]
    rows = []
    for observatory, mjd, exposure, expected_max_db, reason in zip(
        group["observatory"][:],
        group["mjd"][:],
        group["exposure"][:],
        group["expected_max_db"][:],
        group["reason"][:],
    ):
        rows.append(
            dict(
                observatory=observatory.decode(),
                mjd=int(mjd),
                exposure=int(exposure),
                expected_max_db=int(expected_max_db),
                reason=reason.decode(),
            )
        )
    return rows


def write_missing_exposures(fp, rows, replace_keys=None):
    """
    Write the run-level /missing_exposures table to an almanac HDF5 file.

    :param fp:
        An open h5py File (or Group).

    :param rows:
        List of records (dicts) with keys: observatory, mjd, exposure,
        expected_max_db, reason (see `apogee.classify_missing_exposures`).

    :param replace_keys: [optional]
        A set of (observatory, mjd) keys that were processed in this run.
        Existing table entries for other keys are preserved (this supports
        partial re-runs, e.g. with --skip-existing). If `None`, the whole
        table is replaced by `rows`.
    """
    if replace_keys is not None:
        preserved = [
            r for r in read_missing_exposures(fp)
            if (r["observatory"], r["mjd"]) not in replace_keys
        ]
        rows = preserved + list(rows)

    rows = sorted(rows, key=lambda r: (r["observatory"], r["mjd"], r["exposure"]))

    delete_hdf5_entry(fp, MISSING_EXPOSURES_GROUP)
    group = fp.create_group(MISSING_EXPOSURES_GROUP, track_order=True)
    datasets = dict(
        observatory=np.array([r["observatory"] for r in rows], dtype="S3"),
        mjd=np.array([r["mjd"] for r in rows], dtype=np.int64),
        exposure=np.array([r["exposure"] for r in rows], dtype=np.int64),
        expected_max_db=np.array(
            [r["expected_max_db"] for r in rows], dtype=np.int64
        ),
        reason=np.array([r["reason"] for r in rows], dtype="S16"),
    )
    for name, data in datasets.items():
        dataset = group.create_dataset(name, data=data)
        dataset.attrs["description"] = _MISSING_EXPOSURES_DESCRIPTIONS[name]
    return group


def get_or_create_group(fp, group_name):
    try:
        group = fp[group_name]
    except KeyError:
        group = fp.create_group(group_name)
    finally:
        return group


def delete_hdf5_entry(fp, group_name):
    try:
        del fp[group_name]
    except KeyError:
        pass


def get_hdf5_dtype(pydantic_type, sample_value=None):
    """
    Map Pydantic field types to appropriate HDF5/NumPy dtypes.

    Args:
        pydantic_type: The Pydantic field type annotation
        sample_value: A sample value to help determine string lengths, etc.

    Returns:
        Appropriate NumPy dtype for HDF5
    """
    # Handle Union types (including Optional)
    if get_origin(pydantic_type) is Union:
        # For Optional[T] (Union[T, None]), use the non-None type
        args = get_args(pydantic_type)
        non_none_types = [arg for arg in args if arg is not type(None)]
        if non_none_types:
            pydantic_type = non_none_types[0]

    # Handle List types
    if get_origin(pydantic_type) is list:
        inner_type = get_args(pydantic_type)[0]
        return get_hdf5_dtype(inner_type, sample_value)

    # Basic type mappings
    type_mapping = {
        np.int64: np.int64,
        int: np.int64,
        float: np.float64,
        bool: np.bool_,
        str: 'S',  # Will be handled specially for variable length
        bytes: np.bytes_,
        datetime.datetime: 'S19',  # ISO format YYYY-MM-DDTHH:MM:SS
        datetime.date: 'S10',      # ISO format YYYY-MM-DD
        datetime.time: 'S8',       # Format HH:MM:SS
        np.ndarray: np.ndarray,
    }

    # Direct type mapping
    if pydantic_type in type_mapping:
        dtype = type_mapping[pydantic_type]

        # Handle string length determination
        if dtype == 'S' and sample_value is not None:
            if isinstance(sample_value, (list, tuple, np.ndarray)):
                max_len = max((len(str(v)) for v in sample_value), default=1)
            else:
                max_len = len(str(sample_value)) if sample_value else 1
            # Guard against 'S0' (all-empty strings), which h5py rejects
            return f'S{max(max_len, 1)}'
        elif dtype == 'S':
            return 'S100'  # Default string length

        return dtype

    # Handle Enum types
    if isinstance(pydantic_type, type) and issubclass(pydantic_type, Enum):
        # Store enum values as strings
        if sample_value is not None:
            if isinstance(sample_value, (list, tuple)):
                max_len = max(len(str(v.value)) for v in sample_value) if sample_value else 1
            else:
                max_len = len(str(sample_value.value)) if sample_value else 1
            return f'S{max_len}'
        return 'S50'

    # Handle Literal types
    if get_origin(pydantic_type) is Literal:
        args = get_args(pydantic_type)
        if all(isinstance(arg, str) for arg in args):
            max_len = max(len(arg) for arg in args) if args else 1
            return f'S{max_len}'
        elif all(isinstance(arg, int) for arg in args):
            return np.int64
        elif all(isinstance(arg, float) for arg in args):
            return np.float64
        elif all(isinstance(arg, bool) for arg in args):
            return np.bool_

    # Default fallback - try to convert to string
    return 'S100'

def extract_field_data(models: List[BaseModel], field_name: str) -> List[Any]:
    """Extract data for a specific field from all models."""
    return [getattr(model, field_name) for model in models]

def convert_value_for_hdf5(value, target_dtype):
    """Convert a Python value to be compatible with HDF5 storage."""
    if value is None:
        if target_dtype.char == 'S':
            return b''
        elif target_dtype == np.bool_:
            return False
        else:
            return 0  # or np.nan for float types

    if isinstance(value, Enum):
        return str(value.value).encode('utf-8') if target_dtype.char == 'S' else str(value.value)

    if isinstance(value, datetime.datetime):
        return value.isoformat().encode('utf-8')

    if isinstance(value, datetime.date):
        return value.isoformat().encode('utf-8')

    if isinstance(value, datetime.time):
        return value.isoformat().encode('utf-8')

    if isinstance(value, str) and target_dtype.char == 'S':
        return value.encode('utf-8')

    if isinstance(value, list):
        # Handle lists by converting each element
        return [convert_value_for_hdf5(v, target_dtype) for v in value]

    return value


def write_models_to_hdf5_group(
    models: List[BaseModel],
    hdf5_group: h5.Group,
    chunk_size: int = 1000,
    compression: str = None
):
    """
    Write a list of Pydantic models to an HDF5 group as separate datasets per field.

    Args:
        models: List of Pydantic model instances (all same type)
        hdf5_group: HDF5 group to write datasets to
        chunk_size: Chunk size for HDF5 datasets (for performance)
        compression: Compression algorithm ('gzip', 'lzf', 'szip', None)
    """
    if not models:
        # Nothing to write (e.g. a night with zero cross-matched sources);
        # leave the group empty rather than crash on models[0].
        return
    model_type = type(models[0])

    fields = { **model_type.model_fields, **model_type.model_computed_fields }

    data = {
        field_name: extract_field_data(models, field_name) for field_name in fields.keys()
    }
    return _write_models_to_hdf5_group(
        fields,
        data,
        hdf5_group,
        chunk_size=chunk_size,
        compression=compression
    )


def get_default_array(field_spec, num_records: int):
    if isinstance(field_spec, FieldInfo):
        field_type = field_spec.annotation
    else:
        field_type = field_spec.return_type
    # Get the HDF5 dtype
    hdf5_dtype = get_hdf5_dtype(field_type)

    # Get default from field_spec, with fallback for required fields
    default = getattr(field_spec, 'default', PydanticUndefined)
    if default is PydanticUndefined or default is None:
        if np.issubdtype(np.dtype(hdf5_dtype), np.floating):
            default = np.nan
        elif np.issubdtype(np.dtype(hdf5_dtype), np.integer):
            default = -1
        elif np.dtype(hdf5_dtype).kind == 'S':
            default = b""
        elif np.dtype(hdf5_dtype) == np.bool_:
            default = False
        else:
            default = 0

    # Create array filled with default value
    return np.full(num_records, default, dtype=hdf5_dtype)


def _write_models_to_hdf5_group(
    fields,
    data,
    hdf5_group,
    chunk_size: int = None,
    compression: str = None,
    callback: Any = None,
    num_records: int = None,
):
    if callback is None:
        callback = lambda *_, **__: None

    # Determine num_records from existing data if not provided
    if num_records is None:
        for field_name in fields:
            if field_name in data:
                num_records = len(data[field_name])
                break

    for field_name, field_spec in fields.items():

        # Extract data for this field from all models
        try:
            field_data = data[field_name]
        except KeyError:
            # Field is missing - create array with appropriate default value
            if num_records is None:
                print(f"Warning: missing {field_name} and cannot determine array size")
                callback(field_name)
                continue

            # Get the field type
            if isinstance(field_spec, FieldInfo):
                field_type = field_spec.annotation
            else:
                field_type = field_spec.return_type

            # Get the HDF5 dtype
            hdf5_dtype = get_hdf5_dtype(field_type)

            # Get default from field_spec, with fallback for required fields
            default = getattr(field_spec, 'default', PydanticUndefined)
            if default is PydanticUndefined or default is None:
                if np.issubdtype(np.dtype(hdf5_dtype), np.floating):
                    default = np.nan
                elif np.issubdtype(np.dtype(hdf5_dtype), np.integer):
                    default = -1
                elif np.dtype(hdf5_dtype).kind == 'S':
                    default = b""
                elif np.dtype(hdf5_dtype) == np.bool_:
                    default = False
                else:
                    default = 0

            # Create array filled with default value
            field_data = np.full(num_records, default, dtype=hdf5_dtype)
        else:
            field_data = np.atleast_1d(field_data)

        if num_records is None:
            num_records = len(field_data)


        # Determine the appropriate HDF5 dtype
        if isinstance(field_spec, FieldInfo):
            field_type = field_spec.annotation
        else:
            field_type = field_spec.return_type

        hdf5_dtype = get_hdf5_dtype(field_type, field_data)

        chunks = (min(chunk_size, num_records),) if chunk_size is not None and num_records > chunk_size else None
        if field_data.ndim > 1 and chunks is not None:
            chunks = (chunks[0], *field_data.shape[1:])
        compression_setting = compression if chunk_size is not None and num_records > chunk_size else None

        if getattr(field_spec, "annotation", None) is np.ndarray:
            dataset = hdf5_group.create_dataset(
                field_name,
                data=field_data,
                chunks=chunks,
                compression=compression_setting
            )

        elif isinstance(field_data, (np.ndarray, )):
            dataset = hdf5_group.create_dataset(
                field_name,
                data=field_data.astype(np.dtype(hdf5_dtype)),
                chunks=chunks,
                compression=compression_setting
            )

        else:
            # Convert values for HDF5 storage
            converted_data = [convert_value_for_hdf5(value, np.dtype(hdf5_dtype))
                            for value in field_data]

            # Handle variable-length data (like lists)
            if any(isinstance(value, list) for value in converted_data):
                # Create variable-length dataset
                dt = h5.special_dtype(vlen=np.dtype(hdf5_dtype))
                dataset = hdf5_group.create_dataset(
                    field_name,
                    (num_records,),
                    dtype=dt,
                    chunks=True if num_records > chunk_size else None,
                    compression=compression if num_records > chunk_size else None
                )
                dataset[:] = converted_data
            else:
                # Create regular dataset
                np_array = np.array(converted_data, dtype=hdf5_dtype)

                dataset = hdf5_group.create_dataset(
                    field_name,
                    data=np_array,
                    chunks=chunks,
                    compression=compression_setting
                )

        # Add description, even if it is empty string.
        dataset.attrs["description"] = field_spec.description or ""
        callback(field_name)
