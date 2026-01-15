"""
Compute radial velocities for APOGEE observations.

This module provides functions for computing heliocentric and barycentric
radial velocity corrections from pixel offsets measured by the arMADGICS
pipeline.
"""

import numpy as np
import warnings
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
from astropy.time import Time
from erfa import ErfaWarning


# Speed of light in km/s
C_KM_S = 299792.458

# Minimum valid MJD: roughly MJD 40000 (around 1968) to avoid ERFA warnings
# about "dubious year" for dates before UTC was well-defined
MIN_VALID_MJD = 40000


def propagate_pixels_to_z(p, delta_lambda=6e-6):
    """
    Propagate uncertainty from pixel offset to redshift.

    Parameters
    ----------
    p : np.ndarray
        Pixel offset values.
    delta_lambda : float, optional
        Wavelength step per pixel in log-lambda (default: 6e-6).

    Returns
    -------
    np.ndarray
        Derivative dz/dp for error propagation.
    """
    return delta_lambda * np.log(10) * 10 ** (p * delta_lambda)


def propagate_z_to_v(z):
    """
    Propagate uncertainty from redshift to velocity.

    Parameters
    ----------
    z : np.ndarray
        Redshift values.

    Returns
    -------
    np.ndarray
        Derivative dv/dz for error propagation.
    """
    return np.abs(4 * (z + 1) / (((z + 1) ** 2 + 1) ** 2)) * C_KM_S


def pixels_to_z(x, delta=6e-6):
    """
    Convert pixel offset to redshift.

    Parameters
    ----------
    x : np.ndarray
        Pixel offset values.
    delta : float, optional
        Wavelength step per pixel in log-lambda (default: 6e-6).

    Returns
    -------
    np.ndarray
        Redshift values.
    """
    return 10 ** (x * delta) - 1


def z_to_v(z):
    """
    Convert redshift to velocity using the relativistic formula.

    Parameters
    ----------
    z : np.ndarray
        Redshift values.

    Returns
    -------
    np.ndarray
        Velocity in km/s.
    """
    return ((z + 1) ** 2 - 1) / ((z + 1) ** 2 + 1) * C_KM_S


def v_to_z(v):
    """
    Convert velocity to redshift.

    Parameters
    ----------
    v : np.ndarray
        Velocity in km/s.

    Returns
    -------
    np.ndarray
        Redshift values.
    """
    return v / C_KM_S


def compute_barycentric_correction_chunk(ra, dec, mjd, obs_name):
    """
    Compute barycentric velocity correction for a chunk of observations.

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
    np.ndarray
        Barycentric velocity correction in km/s.
    """
    location = EarthLocation.of_site(obs_name)

    # Suppress ERFA warnings about "dubious year"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ErfaWarning)
        result = (
            SkyCoord(
                ra=ra * u.deg,
                dec=dec * u.deg
            )
            .radial_velocity_correction(
                kind="barycentric",
                obstime=Time(mjd, format="mjd"),
                location=location,
            )
            .to(u.km / u.s)
            .value
        )

    return result


def finalize_radial_velocities(
    rv_pixoff_final: np.ndarray,
    rv_pix_var: np.ndarray,
    v_bary_corr: np.ndarray,
) -> dict:
    """
    Finalize radial velocity computation from pixel offsets and barycentric corrections.

    This function takes pre-computed barycentric corrections and combines them
    with pixel offset measurements to produce final radial velocities.

    Parameters
    ----------
    rv_pixoff_final : np.ndarray
        Final pixel offset from arMADGICS.
    rv_pix_var : np.ndarray
        Pixel variance from arMADGICS.
    v_bary_corr : np.ndarray
        Pre-computed barycentric velocity corrections in km/s.

    Returns
    -------
    dict
        Dictionary containing:
        - 'v_barycentric_correction': Barycentric correction in km/s
        - 'v_rel': Relative velocity (no barycentric correction) in km/s
        - 'v_rad': Barycentric-corrected radial velocity in km/s
        - 'e_v_rad': Uncertainty in radial velocity in km/s
    """
    rv_pixoff_final = np.atleast_1d(rv_pixoff_final).astype(float)
    rv_pix_var = np.atleast_1d(rv_pix_var).astype(float)

    # Compute redshift and its uncertainty from pixel offsets
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        e_z = propagate_pixels_to_z(rv_pixoff_final) * np.sqrt(rv_pix_var)

    z = pixels_to_z(rv_pixoff_final)

    # Convert barycentric correction to redshift
    z_corr = v_to_z(v_bary_corr)

    # Compute velocities
    v_rel = z_to_v(z)
    v_rad = z_to_v(z + z_corr + z_corr * z)
    e_v_rad = propagate_z_to_v(z) * e_z

    return {
        "v_barycentric_correction": v_bary_corr,
        "v_rel": v_rel,
        "v_rad": v_rad,
        "e_v_rad": e_v_rad,
    }
