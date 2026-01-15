"""
Targeting flags computation for SDSS-V observations.

This module provides functions for aggregating targeting carton
assignments into SDSS-V target flags arrays using sdss_semaphore.
"""

import numpy as np
from typing import List, Optional, Tuple


def compute_targeting_flags(
    carton_pks_list: List[Optional[List[int]]],
    n_total: int,
) -> Tuple[np.ndarray, int]:
    """
    Aggregate targeting carton primary keys into SDSS-V target flags.

    Parameters
    ----------
    carton_pks_list : list
        List of carton primary key lists for each observation.
        Each element can be None, an empty list, or a list of integers.
    n_total : int
        Total number of observations (used to initialize the flags array).

    Returns
    -------
    tuple
        (flags_array, n_unknown) where:
        - flags_array: numpy array of targeting flags with shape (n_total, 1)
        - n_unknown: count of carton assignments that couldn't be mapped

    Examples
    --------
    >>> carton_pks = [[1, 2, 3], None, [4, 5]]
    >>> flags, n_unknown = compute_targeting_flags(carton_pks, 3)
    """
    from sdss_semaphore.targeting import TargetingFlags

    flags = TargetingFlags(np.zeros((n_total, 1)))
    n_unknown = 0

    for i, carton_pks in enumerate(carton_pks_list):
        for carton_pk in (carton_pks or []):
            try:
                flags.set_bit_by_carton_pk(i, carton_pk)
            except KeyError:
                n_unknown += 1

    return flags.array, n_unknown
