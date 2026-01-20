import numpy as np
from time import time
from datetime import datetime
from itertools import cycle, groupby
from typing import Tuple, Union, Optional, List, Dict, Any
from astropy.table import Table
from astropy.time import Time
from astropy.io.registry import register_identifier, register_reader, register_writer
from pydl.pydlutils.yanny import is_yanny, read_table_yanny, write_table_yanny


register_identifier("yanny", Table, is_yanny)
register_reader("yanny", Table, read_table_yanny)
register_writer("yanny", Table, write_table_yanny)

def adjusted_fiber_index_to_fiber_id(adjusted_fiber_index):
    """
    Adjusted fiber index runs from 1 to 600 (inclusive) to account for APO and LCO.

    Separate from whether it is APO or LCO, the fiber index runs in the opposite
    direction as what we call `fiber_id`.
    """
    adjusted_fiber_index = np.atleast_1d(adjusted_fiber_index)
    adjusted_fiber_index = np.where(
        adjusted_fiber_index > 300,
        adjusted_fiber_index - 300,
        adjusted_fiber_index
    )
    return 301 - adjusted_fiber_index

def group_contiguous(v):
    groups = []
    for k, g in groupby(enumerate(sorted(v)), lambda x: x[1] - x[0]):
        group = list(map(lambda x: x[1], g))
        groups.append((group[0], group[-1]))
    return groups


def get_observatories(apo: bool, lco: bool) -> Tuple[str, ...]:
    """Get observatory names based on boolean flags.

    Args:
        apo: Whether to include APO observatory
        lco: Whether to include LCO observatory

    Returns:
        Tuple of observatory names ("apo", "lco", or both)
    """
    if apo and not lco:
        return ("apo",)
    elif lco and not apo:
        return ("lco",)
    else:
        return ("apo", "lco")


def timestamp_to_mjd(v: float) -> float:
    """Convert Unix timestamp to Modified Julian Date (MJD).

    Args:
        v: Unix timestamp in seconds

    Returns:
        Modified Julian Date as float
    """
    return (v / 86400.0) + 40587.5


def get_current_mjd() -> int:
    """Get current Modified Julian Date as integer.

    Returns:
        Current MJD as integer
    """
    return int(timestamp_to_mjd(time()))


def datetime_to_mjd(date: str) -> int:
    """Convert date string to Modified Julian Date.

    Args:
        date: Date string in format "YYYY-MM-DD"

    Returns:
        Modified Julian Date as integer
    """
    return int(timestamp_to_mjd(datetime.strptime(date, "%Y-%m-%d").timestamp()))

def mjd_to_datetime(mjd: float) -> datetime:
    """Convert Modified Julian Date to datetime object.

    Args:
        mjd: Modified Julian Date

    Returns:
        Datetime object
    """
    return Time(mjd, format='mjd').datetime

def parse_mjds(
    mjd: Optional[int],
    mjd_start: Optional[int],
    mjd_end: Optional[int],
    date: Optional[str],
    date_start: Optional[str],
    date_end: Optional[str],
    earliest_mjd: int = 0,
    return_nones: bool = False
) -> Tuple[Union[int, range, Tuple[int, ...]], int, int]:
    """Parse MJD and date parameters to determine observation date range.

    Args:
        mjd: Single MJD value (can be negative for relative to current)
        mjd_start: Start MJD for range (can be negative for relative to current)
        mjd_end: End MJD for range (can be negative for relative to current)
        date: Single date string in "YYYY-MM-DD" format
        date_start: Start date string in "YYYY-MM-DD" format
        date_end: End date string in "YYYY-MM-DD" format
        earliest_mjd: Earliest allowed MJD value (default: 0)
        return_nones: If True, return None for mjd range if no dates are specified (default: False)

    Returns:
        Tuple containing:
            - MJD values (single int, range, or tuple)
            - Start MJD (int)
            - End MJD (int)

    Raises:
        ValueError: If more than one time specification method is provided
        RuntimeError: If no valid time specification is found
    """
    has_mjd_range = mjd_start is not None or mjd_end is not None
    has_date_range = date_start is not None or date_end is not None

    current_mjd = get_current_mjd()
    n_given = sum([has_mjd_range, has_date_range, mjd is not None, date is not None])
    if n_given > 1:
        raise ValueError(
            "Cannot specify more than one of --mjd, --mjd-start/--mjd-end, --date, --date-start/--date-end"
        )
    if n_given == 0:
        if return_nones:
            return (None, None, None)
        else:
            return ((current_mjd, ), current_mjd, current_mjd)
    if mjd is not None:
        if mjd < 0:
            mjd += current_mjd
        return ((mjd, ), mjd, mjd)
    if has_mjd_range:
        mjd_start = mjd_start or earliest_mjd
        if mjd_start < 0:
            mjd_start += current_mjd
        mjd_end = mjd_end or current_mjd
        if mjd_end < 0:
            mjd_end += current_mjd
        return (range(mjd_start, 1 + mjd_end), mjd_start, mjd_end)
    if date is not None:
        mjd = datetime_to_mjd(date)
        return ((mjd, ), mjd, mjd)
    if has_date_range:
        mjd_start = earliest_mjd if date_start is None else datetime_to_mjd(date_start)
        mjd_end = current_mjd if date_end is None else datetime_to_mjd(date_end)
        return (range(mjd_start, 1 + mjd_end), mjd_start, mjd_end)

    raise RuntimeError("Should not be able to get here")
