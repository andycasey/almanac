"""
Calculate the shadow height for astronomical observations.

Shadow height is the altitude above Earth's surface where a line-of-sight
from an observatory (at given RA, Dec coordinates) intersects the Earth's
umbral shadow cone. This is relevant for observations of Earth satellites
and other objects that may be in shadow.

The algorithm finds the intersection of a ray (from observatory pointing
at celestial coordinates) with the Earth's shadow cone, cast by the Sun.
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import (
    EarthLocation,
    GCRS,
    get_body_barycentric,
    SkyCoord,
    AltAz,
)
from astropy.time import Time
from typing import Literal, Optional

# Observatory definitions using astropy's EarthLocation
OBSERVATORIES = {
    "apo": EarthLocation.of_site("Apache Point Observatory"),
    "lco": EarthLocation.of_site("Las Campanas Observatory"),
}

# Physical constants
EARTH_RADIUS_M = 6.357e6
SUN_RADIUS_M = 695.700e6
AU_TO_M = 1.495978707e11


class ShadowHeightCalculator:
    """
    Calculator for shadow heights at a given observatory.

    This class computes where lines-of-sight from an observatory intersect
    the Earth's umbral shadow cone. The calculation uses the geometry of
    the Sun-Earth system and the ray-cone intersection algorithm.

    Parameters
    ----------
    observatory : str
        Observatory identifier ('apo' or 'lco').

    Attributes
    ----------
    observatory : str
        Observatory identifier.
    shadow_cone_cos_theta_sqr : float
        Squared cosine of the shadow cone opening angle.
    d_ec_au : float
        Distance from Earth to the shadow cone tip in AU.
    """

    def __init__(
        self,
        observatory: Literal["apo", "lco"],
    ):
        if observatory not in OBSERVATORIES:
            raise ValueError(
                f"Unknown observatory '{observatory}'. "
                f"Available: {list(OBSERVATORIES.keys())}"
            )

        self.observatory = observatory
        self._location = OBSERVATORIES[observatory]

        # Physical constants
        earth_radius = EARTH_RADIUS_M * u.m
        sun_radius = SUN_RADIUS_M * u.m

        # Distance from Sun to Earth (1 AU typical)
        d_se = 1 * u.au

        # Distance from Earth to shadow cone tip
        d_ec = d_se * (earth_radius / sun_radius)
        self.d_ec_au = d_ec.to(u.au).value

        # Shadow cone opening angle
        shadow_cone_theta = np.arctan((earth_radius / d_ec).decompose().value)
        self.shadow_cone_cos_theta_sqr = np.cos(shadow_cone_theta) ** 2

        # Cache for time-dependent quantities
        self._cached_jd = None
        self._xyz_earth = None  # Earth barycentric position in AU
        self._xyz_sun = None  # Sun barycentric position in AU
        self._xyz_observatory = None  # Observatory barycentric position in AU
        self._xyz_cone_tip = None  # Cone tip position in AU
        self._cone_axis = None  # Unit vector along cone axis
        self._co = None  # Vector from cone tip to observatory

    def _update_positions(self, jd: float) -> None:
        """
        Update Earth, Sun, and observatory positions for a given Julian Date.

        Parameters
        ----------
        jd : float
            Julian Date (TT scale).
        """
        if self._cached_jd == jd:
            return

        t = Time(jd, format="jd", scale="tt")

        # Get barycentric positions in AU
        earth_bary = get_body_barycentric("earth", t)
        sun_bary = get_body_barycentric("sun", t)

        self._xyz_earth = np.array([
            earth_bary.x.to(u.au).value,
            earth_bary.y.to(u.au).value,
            earth_bary.z.to(u.au).value,
        ])
        self._xyz_sun = np.array([
            sun_bary.x.to(u.au).value,
            sun_bary.y.to(u.au).value,
            sun_bary.z.to(u.au).value,
        ])

        # Observatory position relative to Earth center in AU
        obs_gcrs = self._location.get_gcrs(t)
        obs_offset = np.array([
            obs_gcrs.cartesian.x.to(u.au).value,
            obs_gcrs.cartesian.y.to(u.au).value,
            obs_gcrs.cartesian.z.to(u.au).value,
        ])
        self._xyz_observatory = self._xyz_earth + obs_offset

        # Shadow cone axis unit vector (pointing from Earth toward Sun)
        earth_to_sun = self._xyz_sun - self._xyz_earth
        self._cone_axis = earth_to_sun / np.linalg.norm(earth_to_sun)

        # Shadow cone tip position (on opposite side of Earth from Sun)
        self._xyz_cone_tip = self._xyz_earth - self._cone_axis * self.d_ec_au

        # Vector from cone tip to observatory
        self._co = self._xyz_cone_tip - self._xyz_observatory

        self._cached_jd = jd

    def _ra_dec_to_unit_vectors(
        self, ra: np.ndarray, dec: np.ndarray
    ) -> np.ndarray:
        """
        Convert RA/Dec to unit vectors in the ICRF frame.

        Parameters
        ----------
        ra : np.ndarray
            Right ascension in degrees.
        dec : np.ndarray
            Declination in degrees.

        Returns
        -------
        np.ndarray
            Unit vectors with shape (N, 3).
        """
        ra_rad = np.deg2rad(ra)
        dec_rad = np.deg2rad(dec)

        cos_dec = np.cos(dec_rad)

        unit_vectors = np.column_stack([
            np.cos(ra_rad) * cos_dec,
            np.sin(ra_rad) * cos_dec,
            np.sin(dec_rad),
        ])

        return unit_vectors

    def compute(
        self,
        ra: np.ndarray,
        dec: np.ndarray,
        jd: float,
        unit: str = "km",
    ) -> np.ndarray:
        """
        Compute shadow heights for given celestial coordinates and time.

        Parameters
        ----------
        ra : np.ndarray
            Right ascension in degrees.
        dec : np.ndarray
            Declination in degrees.
        jd : float
            Julian Date (TT scale).
        unit : str, optional
            Output unit for heights (default: 'km').

        Returns
        -------
        np.ndarray
            Shadow heights in the specified unit. NaN for coordinates
            that do not intersect the shadow cone or have negative distances.
        """
        ra = np.atleast_1d(ra).astype(float)
        dec = np.atleast_1d(dec).astype(float)

        # Update positions for this time
        self._update_positions(jd)

        # Convert RA/Dec to unit vectors (pointing directions)
        pointing = self._ra_dec_to_unit_vectors(ra, dec)
        n = len(ra)

        # Solve ray-cone intersection using quadratic formula
        # Ray: P(t) = observatory + t * pointing
        # Cone: ((P - cone_tip) · axis)^2 = cos^2(theta) * |P - cone_tip|^2
        #
        # This leads to: a*t^2 + b*t + c = 0
        # where:
        #   a = (pointing · axis)^2 - cos^2(theta)
        #   b = 2 * [(pointing · axis)(co · axis) - (pointing · co) * cos^2(theta)]
        #   c = (co · axis)^2 - |co|^2 * cos^2(theta)

        v = self._cone_axis
        cos2 = self.shadow_cone_cos_theta_sqr

        # Dot products (vectorized)
        pointing_dot_v = np.sum(pointing * v, axis=1)  # (N,)
        co_dot_v = np.dot(self._co, v)  # scalar
        pointing_dot_co = np.sum(pointing * self._co, axis=1)  # (N,)
        co_dot_co = np.dot(self._co, self._co)  # scalar

        # Quadratic coefficients
        a = pointing_dot_v ** 2 - cos2
        b = 2 * (pointing_dot_v * co_dot_v - pointing_dot_co * cos2)
        c = co_dot_v ** 2 - co_dot_co * cos2

        # Discriminant
        delta = b ** 2 - 4 * a * c

        # Initialize distances to NaN
        dist = np.full(n, np.nan)

        # Solve for positive delta
        valid = delta >= 0
        if np.any(valid):
            sqrt_delta = np.sqrt(delta[valid])
            a_valid = a[valid]
            b_valid = b[valid]

            # Two solutions
            t1 = (-b_valid + sqrt_delta) / (2 * a_valid)
            t2 = (-b_valid - sqrt_delta) / (2 * a_valid)

            # We want the nearest positive intersection
            # (negative means behind the observer)
            t1[t1 < 0] = np.nan
            t2[t2 < 0] = np.nan

            dist[valid] = np.nanmin(np.column_stack([t1, t2]), axis=1)

        # Calculate intersection points in AU
        intersection_xyz = (
            self._xyz_observatory + dist[:, np.newaxis] * pointing
        )

        # Calculate height above Earth's surface
        dist_from_earth_center = np.linalg.norm(
            intersection_xyz - self._xyz_earth, axis=1
        )

        # Convert to specified unit
        heights_au = dist_from_earth_center - (EARTH_RADIUS_M / AU_TO_M)
        heights = (heights_au * u.au).to(unit).value

        return heights


def compute_shadow_heights(
    ra: np.ndarray,
    dec: np.ndarray,
    jd: np.ndarray,
    observatory: np.ndarray,
    unit: str = "km",
) -> np.ndarray:
    """
    Compute shadow heights for batched observations at multiple observatories.

    This is the main entry point for computing shadow heights for a batch of
    observations. It handles observations from different observatories and
    at different times efficiently.

    Parameters
    ----------
    ra : np.ndarray
        Right ascension in degrees.
    dec : np.ndarray
        Declination in degrees.
    jd : np.ndarray
        Julian Dates (TT scale).
    observatory : np.ndarray
        Observatory identifiers ('apo' or 'lco') for each observation.
    unit : str, optional
        Output unit for heights (default: 'km').

    Returns
    -------
    np.ndarray
        Shadow heights in the specified unit. NaN for coordinates
        that do not intersect the shadow cone.

    Examples
    --------
    >>> import numpy as np
    >>> ra = np.array([120.0, 180.0, 240.0])
    >>> dec = np.array([-30.0, -45.0, -60.0])
    >>> jd = np.array([2460000.5, 2460000.5, 2460001.5])
    >>> obs = np.array(['lco', 'lco', 'apo'])
    >>> heights = compute_shadow_heights(ra, dec, jd, obs)
    """
    ra = np.atleast_1d(ra).astype(float)
    dec = np.atleast_1d(dec).astype(float)
    jd = np.atleast_1d(jd).astype(float)
    observatory = np.atleast_1d(observatory)

    # Handle string arrays that may be bytes
    if observatory.dtype.kind == "S":
        observatory = observatory.astype(str)

    n = len(ra)
    if not (len(dec) == len(jd) == len(observatory) == n):
        raise ValueError(
            f"All input arrays must have the same length. "
            f"Got ra={len(ra)}, dec={len(dec)}, jd={len(jd)}, "
            f"observatory={len(observatory)}"
        )

    # Initialize calculators for each observatory
    calculators = {}
    for obs in OBSERVATORIES:
        calculators[obs] = ShadowHeightCalculator(obs)

    # Initialize output array
    heights = np.full(n, np.nan)

    # Group by observatory and JD for efficient batch processing
    unique_obs = np.unique(observatory)

    for obs in unique_obs:
        obs_str = str(obs).lower().strip()
        if obs_str not in calculators:
            continue

        calc = calculators[obs_str]
        obs_mask = np.char.lower(np.char.strip(observatory.astype(str))) == obs_str

        # Get unique JDs for this observatory
        obs_jds = jd[obs_mask]
        unique_jds = np.unique(obs_jds)

        for jd_val in unique_jds:
            # Get indices for this observatory + JD combination
            jd_mask = obs_mask & (jd == jd_val)

            if not np.any(jd_mask):
                continue

            # Compute heights for this batch
            heights[jd_mask] = calc.compute(
                ra[jd_mask],
                dec[jd_mask],
                jd_val,
                unit=unit,
            )

    return heights
