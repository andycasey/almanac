"""
Utility functions for efficient array operations in postprocessing.

This module provides numpy-based alternatives to Python loops for
grouping and indexing operations.
"""

import numpy as np
from typing import Dict, List, Tuple, Any


def group_indices_by_keys(*arrays) -> Dict[Tuple, np.ndarray]:
    """
    Group row indices by unique combinations of key arrays.

    This is a numpy-based alternative to building groupby dictionaries
    with Python loops. It's significantly faster for large arrays.

    Parameters
    ----------
    *arrays : np.ndarray
        One or more arrays of the same length to use as grouping keys.

    Returns
    -------
    dict
        Dictionary mapping unique key tuples to arrays of row indices.

    Examples
    --------
    >>> tele = np.array(['apo', 'lco', 'apo', 'apo'])
    >>> mjd = np.array([60000, 60000, 60001, 60000])
    >>> groups = group_indices_by_keys(tele, mjd)
    >>> groups[('apo', 60000)]
    array([0, 3])
    """
    n = len(arrays[0])

    # Stack arrays and create a structured view for lexsort
    # Convert all to string representation for consistent sorting
    if len(arrays) == 1:
        arr = np.asarray(arrays[0])
        sort_idx = np.argsort(arr, kind='stable')
        sorted_arr = arr[sort_idx]

        # Find where values change
        change_mask = np.concatenate([[True], sorted_arr[1:] != sorted_arr[:-1]])
        change_indices = np.where(change_mask)[0]

        # Build result dictionary
        result = {}
        for i, start in enumerate(change_indices):
            end = change_indices[i + 1] if i + 1 < len(change_indices) else n
            key = sorted_arr[start]
            # Convert numpy types to Python types for hashable keys
            if hasattr(key, 'item'):
                key = key.item()
            result[key] = sort_idx[start:end]

        return result

    # For multiple arrays, create a record array for sorting
    # Convert arrays to consistent types
    processed = []
    for arr in arrays:
        arr = np.asarray(arr)
        if arr.dtype.kind in ('U', 'S', 'O'):
            # String array - convert to fixed-width string
            processed.append(np.char.encode(arr.astype(str), 'utf-8'))
        else:
            processed.append(arr)

    # Use lexsort (sorts by last key first, so reverse the order)
    sort_idx = np.lexsort(processed[::-1])

    # Sort all arrays by the combined key
    sorted_arrays = [np.asarray(arr)[sort_idx] for arr in arrays]

    # Find where any key changes
    change_mask = np.ones(n, dtype=bool)
    change_mask[0] = True
    for sarr in sorted_arrays:
        if sarr.dtype.kind in ('U', 'S'):
            change_mask[1:] |= (sarr[1:] != sarr[:-1])
        else:
            change_mask[1:] |= (sarr[1:] != sarr[:-1])

    change_indices = np.where(change_mask)[0]

    # Build result dictionary
    result = {}
    for i, start in enumerate(change_indices):
        end = change_indices[i + 1] if i + 1 < len(change_indices) else n
        key = tuple(
            sarr[start].item() if hasattr(sarr[start], 'item') else sarr[start]
            for sarr in sorted_arrays
        )
        result[key] = sort_idx[start:end]

    return result


def group_indices_by_array(keys: np.ndarray) -> Dict[Any, np.ndarray]:
    """
    Group row indices by unique values in a single array.

    Optimized version of group_indices_by_keys for single-array case.

    Parameters
    ----------
    keys : np.ndarray
        Array of keys to group by.

    Returns
    -------
    dict
        Dictionary mapping unique keys to arrays of row indices.
    """
    keys = np.asarray(keys)
    sort_idx = np.argsort(keys, kind='stable')
    sorted_keys = keys[sort_idx]

    # Find boundaries between groups
    if sorted_keys.dtype.kind in ('U', 'S', 'O'):
        change_mask = np.concatenate([[True], sorted_keys[1:] != sorted_keys[:-1]])
    else:
        change_mask = np.concatenate([[True], sorted_keys[1:] != sorted_keys[:-1]])

    group_starts = np.where(change_mask)[0]
    group_ends = np.concatenate([group_starts[1:], [len(keys)]])

    # Build result
    result = {}
    for start, end in zip(group_starts, group_ends):
        key = sorted_keys[start]
        if hasattr(key, 'item'):
            key = key.item()
        result[key] = sort_idx[start:end]

    return result


def unique_indices(keys: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get unique values, their first occurrence indices, and inverse mapping.

    Similar to np.unique with return_index and return_inverse, but returns
    the indices in a more useful format for vectorized assignment.

    Parameters
    ----------
    keys : np.ndarray
        Array of keys.

    Returns
    -------
    tuple
        (unique_keys, first_indices, inverse_indices) where:
        - unique_keys: array of unique values
        - first_indices: index of first occurrence of each unique value
        - inverse_indices: mapping from each row to its unique key index
    """
    return np.unique(keys, return_index=True, return_inverse=True)


def vectorized_assign_by_group(
    target: np.ndarray,
    source_dict: Dict[Any, Any],
    key_to_indices: Dict[Any, np.ndarray],
) -> None:
    """
    Vectorized assignment from a dictionary to an array using pre-computed indices.

    Parameters
    ----------
    target : np.ndarray
        Target array to assign values to (modified in place).
    source_dict : dict
        Dictionary mapping keys to values.
    key_to_indices : dict
        Dictionary mapping keys to arrays of target indices.
    """
    for key, indices in key_to_indices.items():
        if key in source_dict:
            target[indices] = source_dict[key]
