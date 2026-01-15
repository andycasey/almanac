from typing import Optional
from pydantic import BaseModel, Field

from almanac.data_models.types import *
from almanac.data_models.source import Source

class ReducedSpectrum(Source):
    """
    A data model combining exposure-level, fiber-level, and target-level information.

    This model merges fields from Exposure, FPSTarget, and PlateTarget classes, and
    adds metadata from ARMADGICS.
    """

    #> Exposure-level Information
    observatory: Observatory = Field(description="Observatory name")
    mjd: int = Field(description="MJD of the exposure")
    exposure: int = Field(description="Exposure number", ge=1)
    prefix: Optional[Prefix] = Field(description="Raw exposure basename prefix", default=None)

    #> Exposure Metadata
    name: Optional[str] = Field(
        default="",
        description=(
            "The `name` field in the exposure header often refers to the plugged "
            "plate name, which describes which targets were observed."
        )
    )
    n_read: int = Field(default=0, alias="nread", ge=0)
    image_type: ImageType = Field(alias="imagetyp")
    observer_comment: Optional[str] = Field(default="", alias="obscmnt")

    #> Exposure Identifiers
    map_id: int = Field(default=-1, alias="mapid")
    cart_id: int = Field(default=-1, alias="cartid")
    plate_id: int = Field(default=-1, alias="plateid")
    field_id: int = Field(default=-1, alias="fieldid")
    design_id: int = Field(default=-1, alias="designid")
    config_id: int = Field(default=-1, alias="configid")

    #> Observing Conditions
    seeing: float = Field(default=float('NaN'))

    #> Instrument State
    focus: float = Field(default=float('NaN'))
    collpist: float = Field(default=float('NaN'))
    colpitch: float = Field(default=float('NaN'))
    dithered_pixels: float = Field(default=float('NaN'), alias="dithpix")
    lamp_quartz: int = Field(default=-1, alias="lampqrtz", ge=-1, le=1)
    lamp_thar: int = Field(default=-1, alias="lampthar", ge=-1, le=1)
    lamp_une: int = Field(default=-1, alias="lampune", ge=-1, le=1)

    #> Instrument Identifiers
    spectrograph_id: int = Field(description="Spectrograph identifier", alias='spectrographId', default=-1)
    fiber_id: int = Field(description="Fiber identifier", alias='fiberId', default=-1)
    planned_fiber_id: int = Field(alias="fiberid", description="Planned fiber ID (plate era)", default=-1)
    throughput: int = Field(description="Throughput value", default=-1)

    #> Wavelength Information
    lambda_design: float = Field(description="Design wavelength (FPS era)", default=0.0)
    lambda_eff: float = Field(description="Effective wavelength", default=0.0)
    coord_epoch: float = Field(description="Coordinate epoch (FPS era)", default=0.0)
    zoffset: float = Field(description="Z offset (plate era)", default=0.0)

    #> Target Identification
    sdss_id: Int64 = Field(default=-1)
    catalogid: Int64 = Field(default=-1)
    twomass_designation: Annotated[str, Field(default="", alias="tmass_id"), BeforeValidator(validate_str)]
    target_ids: str = Field(alias="targetids", default="", description="Legacy target IDs (plate era)")
    category: Category = Field(description="Category of the target")
    cadence: str = Field(description="Cadence identifier", default="")
    firstcarton: str = Field(description="Main carton from which this target was drawn", default="")
    program: str = Field(description="Program for 'firstcarton'", default="")

    #> Positioner and Hole Identifiers
    positioner_id: int = Field(alias='positionerId', description="Positioner identifier", default=-1)
    hole_id: str = Field(alias='holeId', description="Hole ID in which the positioner is sitting", default="")
    hole_type: HoleType = Field(alias="holeType", description="Type of hole", default="fps")
    planned_hole_type: HoleType = Field(alias="holetype", description="Hole type string", default="fps")
    obj_type: ObjType = Field(alias="objType", description="Object type (plate era)", default="na")
    fiber_type: str = Field(alias='fiberType', description="Type of fiber", default="")
    assigned: bool = Field(
        default=False,
        description=(
            "Target is assigned to this fiber. If False, no target assigned for this fiber, "
            "and no targeting information available"
        )
    )

    #> Position Coordinates (Common)
    x_focal: float = Field(description="x-coordinate in the focal plane", default=float('NaN'), alias='xFocal')
    y_focal: float = Field(description="y-coordinate in the focal plane", default=float('NaN'), alias='yFocal')

    #> Plate-specific Information
    iplateinput: int = Field(description="Plate input ID", default=-1)
    pointing: int = Field(description="Pointing number", default=-1)
    offset: int = Field(description="Offset value", default=-1)
    block: int = Field(description="Block number", default=-1)
    iguide: int = Field(description="Guide flag", default=-1)
    bluefiber: int = Field(description="Blue fiber flag", default=-1)
    chunk: int = Field(description="Chunk number", default=-1)
    ifinal: int = Field(description="Final flag", default=-1)
    plugged_mjd: int = Field(description="MJD when this plate was plugged", default=-1)
    fix_fiber_flag: int = Field(default=0, description="Whether this fiber mapping was fixed in software")
    diameter: float = Field(default=-1, description="Diameter")
    buffer: float = Field(default=-1, description="Buffer size")
    priority: int = Field(default=-1, description="Target priority")

    #> Position Coordinates (Plate)
    xf_default: float = Field(description="Default X focal coordinate (plate era)", default=float('NaN'))
    yf_default: float = Field(description="Default Y focal coordinate (plate era)", default=float('NaN'))

    #> Status Flags (Plate)
    conflicted: bool = Field(description="Conflicted flag", default=False)
    ranout: bool = Field(description="Ran out flag", default=False)
    outside: bool = Field(description="Outside flag", default=False)

    #> Position Coordinates (FPS)
    x_wok: float = Field(description="x-coordinate in the wok frame (FPS era)", default=float('NaN'), alias="xwok")
    y_wok: float = Field(description="y-coordinate in the wok frame (FPS era)", default=float('NaN'), alias="ywok")
    z_wok: float = Field(description="z-coordinate in the wok frame (FPS era)", default=float('NaN'), alias="zwok")

    #> Positioner Angles (FPS)
    alpha: float = Field(description="Alpha angle of the positioner arm (FPS era)", default=float('NaN'))
    beta: float = Field(description="Beta angle of the positioner arm (FPS era)", default=float('NaN'))

    #> Status Flags (FPS)
    on_target: bool = Field(description="Fiber placed on target", default=False)
    disabled: bool = Field(description="Fiber is disabled", default=False)
    valid: bool = Field(description="Converted on-sky coordinates to robot (α,β)", default=False)
    decollided: bool = Field(description="Positioner had to be moved to decollide it", default=False)

    #> Position Deltas (FPS)
    delta_ra: float = Field(description="The amount in RA this fiber has been offset (FPS era)", default=float('NaN'))
    delta_dec: float = Field(description="The amount in Dec this fiber has been offset (FPS era)", default=float('NaN'))

    #> Target Coordinates
    ra: float = Field(alias="racat", description="Right Ascension [deg]")
    dec: float = Field(alias="deccat", description="Declination [deg]")
    alt: float = Field(description="Altitude of the fiber on the sky [deg]", default=float('NaN'), alias="alt_observed")
    az: float = Field(description="Azimuth of the fiber on the sky [deg]", default=float('NaN'), alias="az_observed")

    #> Target of Opportunity (FPS)
    too: bool = Field(default=False, description="Target of opportunity (FPS era)")
    too_id: int = Field(default=-1, description="Target of opportunity ID (FPS era)")
    too_program: str = Field(default="", description="Target of opportunity program (FPS era)")

    #> ARMADGICS
    RV_flag: int = Field(description="Grid search flag value from apMADGICS from the search for stellar radial velocity. See https://github.com/andrew-saydjari/apMADGICS.jl for bit interpretations.", default=-1)
    RV_minchi2_final: float = Field(description="Value of the minimum on the delta chi2 surface for the stellar radial velocity determinination step in apMADGICS", default=float('NaN'))
    RV_pix_var: float = Field(description="Stellar radial velocity uncertainty expressed as a variance in the pixel offset [pixels]", default=float('NaN'))
    RV_pixoff_disc_final: float = Field(description="Pixel offset at nearest stellar radial velocity optimum [pixels]", default=-9999) # TODO
    RV_pixoff_final: float = Field(description="Pixel offset at nearest stellar radial velocity optimum", default=float('NaN'))
    RVchi2_residuals: float = Field(description="Chi2 value for the residual component after the RV fitting step for the apMADGICS component separation", default=float('NaN'))
    adjfiberindx: float = Field(description="A unique fiber identifier running 1 to 600 to handle fibers at both APO and LCO", default=float('NaN'))
    avg_flux_conservation: float = Field(description="Median fractional flux conservation of MADGICS component separation across the visit spectrum", default=float('NaN'))
    data_pix_cnt: float = Field(description="Number of unmasked pixels in the input spectrum to apMADGICS component separation", default=float('NaN'))
    final_pix_cnt: float = Field(description="Final number of unmasked pixels used for modeling the visit spectrum", default=float('NaN'))
    flux_nans: float = Field(description="Number of nan pixels in the reinterpolated flux spectrum", default=float('NaN'))
    fluxerr2_nans: float = Field(description="Number of nan pixels in the reinterpolated flux variance spectrum", default=float('NaN'))
    skyscale0: float = Field(description="First pass median sky flux in apMADGICS sky component model of the visit [ADU]", default=float('NaN'))
    skyscale1: float = Field(description="Second pass median sky flux in apMADGICS sky component model of the visit [ADU]", default=float('NaN'))
    starscale0: float = Field(description="First pass median star flux in apMADGICS star continuum component model of the visit [ADU]", default=float('NaN'), alias="starscale")
    starscale1: float = Field(description="Second pass median star flux in apMADGICS star continuum component model of the visit [ADU]", default=float('NaN'))
    tot_p5chi2_v0: float = Field(description="Total chi2 for the MADGICS model including only the sky continuum, (faint) sky lines, stellar continuum, stellar lines, and residual components", default=float('NaN'))

    #> Radial Velocities
    v_rad: float = Field(description="Barycentric rest frame radial velocity [km/s]", default=float('NaN'))
    e_v_rad: float = Field(description="Uncertainty on barycentric rest frame radial velocity [km/s]", default=float('NaN'))
    v_barycentric_correction: float = Field(description="Barycentric velocity correction applied [km/s]", default=float('NaN'))


    class Config:
        validate_by_name = True
        validate_assignment = True
        arbitrary_types_allowed = True
