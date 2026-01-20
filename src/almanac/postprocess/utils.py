"""
Utility functions for efficient array operations in postprocessing.

This module provides numpy-based alternatives to Python loops for
grouping and indexing operations.
"""

import numpy as np
from typing import Dict, List, Tuple, Any


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

    for key, indices in result.items():
        assert np.all(keys[indices] == key)

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
