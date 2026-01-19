"""
Data loading utilities for postprocessing.

This module provides functions for loading data from various sources
including almanac HDF5 files, arMADGICS outputs, ar1Dunical files,
and metadata pickle files.
"""

import os
import pickle
import h5py as h5
from typing import Tuple, Dict, List, Any


# Keys to extract from ar1Dunical metadata
AR1D_KEYS = (
    "bitmsk_relFluxFile",
    "cartid",
    "extraction_method",
    "mjd_mid_exposure",
    "ndiff_used",
    "git_commit",
    "git_branch",
    "git_clean",
    "nread_total",
    "trace_found_match",
    "trace_orig_param_fname",
    "trace_type",
    "wavecal_type",
)


def load_ar1d_unical_meta(outdir: str, obs: str, mjd: int, expnum: int) -> dict:
    """
    Read metadata from a single ar1Dunical file.

    Parameters
    ----------
    outdir : str
        Base output directory containing the apred folder.
    obs : str
        Observatory identifier ('apo' or 'lco').
    mjd : int
        Modified Julian Date.
    expnum : int
        Exposure number.

    Returns
    -------
    dict
        Dictionary containing ar1Dunical metadata fields, or
        {'flag_missing_ar1dunical': True} if file doesn't exist.
    """
    ar1dunical_path = f"{outdir}/apred/{mjd}/ar1Dunical_{obs}_{mjd}_{expnum:04d}_object.h5"

    # Check existence first (cheaper than catching exception)
    if not os.path.exists(ar1dunical_path):
        print(f"Could not find path {ar1dunical_path}")
        return {"flag_missing_ar1dunical": True}

    d = {"flag_missing_ar1dunical": False}
    with h5.File(ar1dunical_path, "r") as ar1d:
        meta = ar1d["metadata"]
        for key in AR1D_KEYS:
            d[key] = meta[key][()]
    return d


def load_ar1d_unical_meta_batch(args: Tuple[str, int, List[Tuple[str, int]]]) -> dict:
    """
    Process all exposures for a single MJD in one batch.

    This reduces ProcessPoolExecutor overhead by having fewer, larger tasks.
    Files for the same MJD are co-located on disk, improving cache performance.

    Parameters
    ----------
    args : tuple
        Tuple of (outdir, mjd, tele_exp_list) where tele_exp_list is a list
        of (telescope, exposure_number) tuples.

    Returns
    -------
    dict
        Dictionary mapping (tele, mjd, exp) tuples to metadata dictionaries.
    """
    outdir, mjd, tele_exp_list = args
    results = {}
    for tele, exp in tele_exp_list:
        results[(tele, mjd, exp)] = load_ar1d_unical_meta(outdir, tele, mjd, exp)
    return results


def load_almanac_file(
    input_path: str,
    mjd_min: int = None,
    mjd_max: int = None,
) -> Tuple[Dict, Dict, List[Tuple[str, int]]]:
    """
    Load all almanac exposure and fiber data into memory.

    Parameters
    ----------
    input_path : str
        Path to the almanac HDF5 file.
    mjd_min : int, optional
        Minimum MJD to load (inclusive). If None, no lower bound.
    mjd_max : int, optional
        Maximum MJD to load (inclusive). If None, no upper bound.

    Returns
    -------
    tuple
        (almanac_exposures, almanac_fibers, obs_mjd_keys) where:
        - almanac_exposures: dict mapping (obs, mjd) to exposure data
        - almanac_fibers: dict mapping (obs, mjd) to fiber data by config
        - obs_mjd_keys: list of (obs, mjd) tuples that were loaded
    """
    almanac_exposures = {}
    almanac_fibers = {}
    obs_mjd_keys = []

    with h5.File(input_path, "r") as fp:
        all_keys = [(obs, int(mjd)) for obs in fp.keys() for mjd in fp[obs].keys()]

        for obs, mjd in all_keys:
            # Apply MJD filtering
            if mjd_min is not None and mjd < mjd_min:
                continue
            if mjd_max is not None and mjd > mjd_max:
                continue

            key = (obs, mjd)
            obs_mjd_keys.append(key)

            g = fp[f"{obs}/{mjd}"]
            almanac_exposures[key] = {k: v[:] for k, v in g["exposures"].items()}
            almanac_fibers[key] = {
                int(cfg): {k: v[:] for k, v in g[f"fibers/{cfg}"].items()}
                for cfg in g.get("fibers", {}).keys()
            }

    return almanac_exposures, almanac_fibers, obs_mjd_keys


def load_armadgics_files(armadgics_paths: List[str]) -> dict:
    """
    Load arMADGICS scalar arrays from HDF5 files.

    Parameters
    ----------
    armadgics_paths : list of str
        List of paths to arMADGICS output HDF5 files.

    Returns
    -------
    dict
        Dictionary mapping field names to numpy arrays.
    """
    result = {}
    for path in armadgics_paths:
        basename = os.path.basename(path)
        # Extract key from filename: arMADGICS_out_<key>.h5
        key = basename[14:-3]
        with h5.File(path, "r") as fp:
            if fp[key].ndim == 1:
                result[key] = fp[key][:]
    return result


def load_metadata_pickle(pickle_path: str) -> dict:
    """
    Load metadata from pickle file.

    Parameters
    ----------
    pickle_path : str
        Path to the pickle file containing source metadata.

    Returns
    -------
    dict
        Dictionary mapping sdss_id to metadata dictionaries.
    """
    with open(pickle_path, "rb") as fp:
        return pickle.load(fp)
