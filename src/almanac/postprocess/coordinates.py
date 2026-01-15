"""
Compute astronomical coordinate quantities for observations.

This module provides efficient batched calculations for:
- Moon phase (illumination fraction)
- Moon separation (angular distance from observation to moon)
- Airmass (optical path length through atmosphere)
"""

import numpy as np
import warnings
from astropy import units as u
from astropy.coordinates import (
    EarthLocation,
    SkyCoord,
    AltAz,
    get_body,
)
from astropy.time import Time
from erfa import ErfaWarning
from typing import Literal

# Observatory definitions
OBSERVATORIES = {
    "apo": EarthLocation.of_site("Apache Point Observatory"),
    "lco": EarthLocation.of_site("Las Campanas Observatory"),
}

# Minimum valid MJD: roughly MJD 40000 (around 1968) to avoid ERFA warnings
# about "dubious year" for dates before UTC was well-defined
MIN_VALID_MJD = 40000


def compute_moon_phase(time: Time) -> np.ndarray:
    """
    Compute the moon illumination fraction at given times.

    The moon phase is computed as the fraction of the moon's disk that
    is illuminated, ranging from 0 (new moon) to 1 (full moon).

    Parameters
    ----------
    time : astropy.time.Time
        Observation times.

    Returns
    -------
    np.ndarray
        Moon illumination fraction (0 to 1).
    """
    # Get Sun and Moon positions in GCRS (geocentric) frame
    # Both are in the same frame so no transformation warning
    sun = get_body("sun", time)
    moon = get_body("moon", time)

    # Elongation: angular separation between Sun and Moon as seen from Earth
    elongation = sun.separation(moon)

    # Moon phase (illumination fraction) using the approximation:
    # phase = (1 - cos(elongation)) / 2
    # This gives 0 at new moon (elongation=0) and 1 at full moon (elongation=180)
    phase = (1 - np.cos(elongation.rad)) / 2

    return phase


def compute_moon_separation(
    ra: np.ndarray,
    dec: np.ndarray,
    time: Time,
) -> np.ndarray:
    """
    Compute the angular separation between observations and the Moon.

    Parameters
    ----------
    ra : np.ndarray
        Right ascension in degrees.
    dec : np.ndarray
        Declination in degrees.
    time : astropy.time.Time
        Observation times.

    Returns
    -------
    np.ndarray
        Angular separation in degrees.
    """
    # Get Moon position at each observation time and transform to ICRS
    # to avoid NonRotationTransformationWarning when computing separation
    moon = get_body("moon", time).transform_to("icrs")

    # Create SkyCoord for observation positions
    obs_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

    # Compute separation
    separation = obs_coord.separation(moon)

    return separation.deg


def compute_airmass(
    ra: np.ndarray,
    dec: np.ndarray,
    time: Time,
    location: EarthLocation,
) -> np.ndarray:
    """
    Compute the airmass for observations at a given location.

    Airmass is computed using the secant of the zenith angle, with
    a correction for atmospheric refraction at low altitudes.

    Parameters
    ----------
    ra : np.ndarray
        Right ascension in degrees.
    dec : np.ndarray
        Declination in degrees.
    time : astropy.time.Time
        Observation times.
    location : astropy.coordinates.EarthLocation
        Observatory location.

    Returns
    -------
    np.ndarray
        Airmass values. NaN for objects below the horizon.
    """
    # Create SkyCoord for observation positions
    obs_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

    # Transform to AltAz frame
    altaz_frame = AltAz(obstime=time, location=location)
    altaz = obs_coord.transform_to(altaz_frame)

    # Compute airmass using secant approximation with Pickering (2002) correction
    # This is more accurate at high zenith angles than simple sec(z)
    altitude = altaz.alt.deg

    # Initialize with NaN
    airmass = np.full(len(altitude), np.nan)

    # Only compute for objects above horizon
    above_horizon = altitude > 0
    if np.any(above_horizon):
        alt = altitude[above_horizon]
        # Pickering (2002) formula for airmass
        # X = 1 / sin(h + 244/(165 + 47*h^1.1))
        # where h is altitude in degrees
        arg = alt + 244.0 / (165.0 + 47.0 * alt**1.1)
        airmass[above_horizon] = 1.0 / np.sin(np.deg2rad(arg))

    return airmass


def compute_altaz(
    ra: np.ndarray,
    dec: np.ndarray,
    time: Time,
    location: EarthLocation,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute altitude and azimuth for observations at a given location.

    Parameters
    ----------
    ra : np.ndarray
        Right ascension in degrees.
    dec : np.ndarray
        Declination in degrees.
    time : astropy.time.Time
        Observation times.
    location : astropy.coordinates.EarthLocation
        Observatory location.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Altitude and azimuth in degrees.
    """
    # Create SkyCoord for observation positions
    obs_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

    # Transform to AltAz frame
    altaz_frame = AltAz(obstime=time, location=location)
    altaz = obs_coord.transform_to(altaz_frame)

    return altaz.alt.deg, altaz.az.deg


# =============================================================================
# Chunk worker functions for parallel processing
# These are designed to be called by ProcessPoolExecutor
# =============================================================================


def compute_observation_metadata_chunk(ra, dec, mjd, obs_name):
    """
    Compute moon phase, moon separation, airmass, alt, and az for a chunk.

    This function is designed to be called by a ProcessPoolExecutor.

    Parameters
    ----------
    ra : np.ndarray
        Right ascension in degrees.
    dec : np.ndarray
        Declination in degrees.
    mjd : np.ndarray
        Modified Julian Date of observations.
    obs_name : str
        Observatory identifier ('apo' or 'lco').

    Returns
    -------
    dict
        Dictionary with 'moon_phase', 'moon_separation', 'airmass', 'alt', 'az'.
    """
    location = OBSERVATORIES[obs_name]

    # Suppress ERFA warnings about "dubious year"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ErfaWarning)

        t = Time(mjd, format="mjd")

        # Moon phase
        moon_phase = compute_moon_phase(t)

        # Moon separation
        moon_separation = compute_moon_separation(ra, dec, t)

        # Airmass
        airmass = compute_airmass(ra, dec, t, location)

        # Alt/Az
        alt, az = compute_altaz(ra, dec, t, location)

    return {
        "moon_phase": moon_phase,
        "moon_separation": moon_separation,
        "airmass": airmass,
        "alt": alt,
        "az": az,
    }


def fill_missing_altaz(
    ra: np.ndarray,
    dec: np.ndarray,
    mjd: np.ndarray,
    observatory: np.ndarray,
    alt: np.ndarray,
    az: np.ndarray,
    mjd_precision: float = 0.01,
    progress_callback: callable = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Fill in missing (NaN) altitude and azimuth values.

    This function computes alt/az only for entries where the existing
    values are NaN, preserving values that were already computed by
    the pipeline.

    Parameters
    ----------
    ra : np.ndarray
        Right ascension in degrees.
    dec : np.ndarray
        Declination in degrees.
    mjd : np.ndarray
        Modified Julian Dates.
    observatory : np.ndarray
        Observatory identifiers ('apo' or 'lco') for each observation.
    alt : np.ndarray
        Existing altitude values (may contain NaN).
    az : np.ndarray
        Existing azimuth values (may contain NaN).
    mjd_precision : float, optional
        Precision for grouping MJDs (default: 0.01 days ~ 15 min).
    progress_callback : callable, optional
        Function to call after each batch with number of items processed.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Updated altitude and azimuth arrays with NaN values filled in.
    """
    ra = np.atleast_1d(ra).astype(float)
    dec = np.atleast_1d(dec).astype(float)
    mjd = np.atleast_1d(mjd).astype(float)
    observatory = np.atleast_1d(observatory)
    alt = np.atleast_1d(alt).astype(float).copy()
    az = np.atleast_1d(az).astype(float).copy()

    # Handle string arrays that may be bytes
    if observatory.dtype.kind == "S":
        observatory = observatory.astype(str)

    n = len(ra)

    # Find entries that need computation (either alt or az is NaN)
    needs_compute = ~np.isfinite(alt) | ~np.isfinite(az)

    # Also need valid coordinates and dates
    valid_coords = (
        np.isfinite(ra)
        & np.isfinite(dec)
        & np.isfinite(mjd)
        & (ra >= 0) & (ra <= 360)
        & (dec >= -90) & (dec <= 90)
        & (mjd >= MIN_VALID_MJD)
    )

    to_compute = needs_compute & valid_coords

    if not np.any(to_compute):
        return alt, az

    # Round MJDs to group similar times
    mjd_rounded = np.round(mjd / mjd_precision) * mjd_precision

    obs_lower = np.char.lower(np.char.strip(observatory.astype(str)))
    unique_obs = np.unique(obs_lower[to_compute])

    for obs in unique_obs:
        obs_str = str(obs).lower().strip()
        if obs_str not in OBSERVATORIES:
            continue

        location = OBSERVATORIES[obs_str]
        obs_mask = to_compute & (obs_lower == obs_str)

        if not np.any(obs_mask):
            continue

        # Group by unique rounded MJD for this observatory
        obs_mjds_rounded = mjd_rounded[obs_mask]
        unique_obs_mjds = np.unique(obs_mjds_rounded)

        for mjd_val in unique_obs_mjds:
            mjd_mask = obs_mask & (mjd_rounded == mjd_val)
            if not np.any(mjd_mask):
                continue

            # Suppress ERFA warnings about "dubious year" for dates at the edge
            # of ERFA's leap second predictions - these are benign for our use case
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=ErfaWarning)
                t = Time(mjd_val, format="mjd")
                computed_alt, computed_az = compute_altaz(
                    ra[mjd_mask], dec[mjd_mask], t, location
                )
            alt[mjd_mask] = computed_alt
            az[mjd_mask] = computed_az

            if progress_callback:
                progress_callback(np.sum(mjd_mask))

    return alt, az


def compute_observation_metadata(
    ra: np.ndarray,
    dec: np.ndarray,
    mjd: np.ndarray,
    observatory: np.ndarray,
    mjd_precision: float = 0.01,
    progress_callback: callable = None,
) -> dict:
    """
    Compute moon phase, moon separation, and airmass for batched observations.

    This is the main entry point for computing observation metadata for a batch
    of observations. It handles observations from different observatories and
    at different times efficiently by grouping computations.

    Parameters
    ----------
    ra : np.ndarray
        Right ascension in degrees.
    dec : np.ndarray
        Declination in degrees.
    mjd : np.ndarray
        Modified Julian Dates.
    observatory : np.ndarray
        Observatory identifiers ('apo' or 'lco') for each observation.
    mjd_precision : float, optional
        Precision for grouping MJDs (default: 0.01 days ~ 15 min).
        Observations within this range use the same ephemeris calculation.
    progress_callback : callable, optional
        Function to call after each batch with number of items processed.

    Returns
    -------
    dict
        Dictionary with keys 'moon_phase', 'moon_separation', 'airmass',
        each containing an ndarray of values.

    Examples
    --------
    >>> import numpy as np
    >>> ra = np.array([120.0, 180.0, 240.0])
    >>> dec = np.array([-30.0, -45.0, -60.0])
    >>> mjd = np.array([60000.5, 60000.5, 60001.5])
    >>> obs = np.array(['lco', 'lco', 'apo'])
    >>> result = compute_observation_metadata(ra, dec, mjd, obs)
    """
    ra = np.atleast_1d(ra).astype(float)
    dec = np.atleast_1d(dec).astype(float)
    mjd = np.atleast_1d(mjd).astype(float)
    observatory = np.atleast_1d(observatory)

    # Handle string arrays that may be bytes
    if observatory.dtype.kind == "S":
        observatory = observatory.astype(str)

    n = len(ra)
    if not (len(dec) == len(mjd) == len(observatory) == n):
        raise ValueError(
            f"All input arrays must have the same length. "
            f"Got ra={len(ra)}, dec={len(dec)}, mjd={len(mjd)}, "
            f"observatory={len(observatory)}"
        )

    # Initialize output arrays
    moon_phase = np.full(n, np.nan)
    moon_separation = np.full(n, np.nan)
    airmass = np.full(n, np.nan)

    # Create mask for valid coordinates and dates
    valid_coords = (
        np.isfinite(ra)
        & np.isfinite(dec)
        & np.isfinite(mjd)
        & (ra >= 0) & (ra <= 360)
        & (dec >= -90) & (dec <= 90)
        & (mjd >= MIN_VALID_MJD)
    )

    if not np.any(valid_coords):
        return {
            "moon_phase": moon_phase,
            "moon_separation": moon_separation,
            "airmass": airmass,
        }

    # Round MJDs to group similar times (reduces ephemeris calculations)
    mjd_rounded = np.round(mjd / mjd_precision) * mjd_precision

    # Process by observatory for airmass calculation
    # Moon phase and separation don't depend on observatory location
    obs_lower = np.char.lower(np.char.strip(observatory.astype(str)))
    unique_obs = np.unique(obs_lower[valid_coords])

    # Get unique rounded MJDs for efficiency
    valid_mjds_rounded = mjd_rounded[valid_coords]
    unique_mjds = np.unique(valid_mjds_rounded)

    # Suppress ERFA warnings about "dubious year" for dates at the edge
    # of ERFA's leap second predictions - these are benign for our use case
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ErfaWarning)

        # For moon phase, we only need to compute once per unique rounded MJD
        mjd_to_phase = {}
        for mjd_val in unique_mjds:
            t = Time(mjd_val, format="mjd")
            mjd_to_phase[mjd_val] = compute_moon_phase(t)

        # Assign moon phase values using rounded MJD lookup
        valid_indices = np.where(valid_coords)[0]
        for i in valid_indices:
            moon_phase[i] = mjd_to_phase[mjd_rounded[i]]

        # For moon separation, we need (ra, dec, mjd) - group by rounded MJD
        for mjd_val in unique_mjds:
            mjd_mask = valid_coords & (mjd_rounded == mjd_val)
            if not np.any(mjd_mask):
                continue

            t = Time(mjd_val, format="mjd")
            moon_separation[mjd_mask] = compute_moon_separation(
                ra[mjd_mask], dec[mjd_mask], t
            )

            if progress_callback:
                progress_callback(np.sum(mjd_mask))

        # For airmass, we need (ra, dec, mjd, observatory)
        for obs in unique_obs:
            obs_str = str(obs).lower().strip()
            if obs_str not in OBSERVATORIES:
                continue

            location = OBSERVATORIES[obs_str]
            obs_mask = valid_coords & (obs_lower == obs_str)

            if not np.any(obs_mask):
                continue

            # Group by unique rounded MJD for this observatory
            obs_mjds_rounded = mjd_rounded[obs_mask]
            unique_obs_mjds = np.unique(obs_mjds_rounded)

            for mjd_val in unique_obs_mjds:
                mjd_mask = obs_mask & (mjd_rounded == mjd_val)
                if not np.any(mjd_mask):
                    continue

                t = Time(mjd_val, format="mjd")
                airmass[mjd_mask] = compute_airmass(
                    ra[mjd_mask], dec[mjd_mask], t, location
                )

    return {
        "moon_phase": moon_phase,
        "moon_separation": moon_separation,
        "airmass": airmass,
    }
