import numpy as np
from typing import Literal
from pydantic import BaseModel, Field, validator, model_validator, computed_field

from almanac.data_models.types import *
from almanac.data_models.utils import sanitise_twomass_designation

class PlateTarget(BaseModel):

    """ A target that was observed with plates. """

    # Target information
    sdss_id: Int64 = Field(default=-1)
    catalogid: Int64 = Field(default=-1)

    @computed_field
    def twomass_designation(self) -> str:
        """ Convert a target ID to a standardized designation format. """
        return sanitise_twomass_designation(self.twomass_id or self.target_ids)

    twomass_id: str = Field(alias="tmass_id", default="")
    target_ids: str = Field(alias="targetids", default="")
    category: Category = Field(description="Category of the target", alias="targettype")

    # Positioner and hole identifiers
    observatory: Observatory = Field(description="Observatory") # necessary for fiber mapping fixes
    hole_type: HoleType = Field(alias="holeType", description="Type of hole")
    planned_hole_type: HoleType = Field(alias="holetype", description="Hole type string")
    obj_type: ObjType = Field(alias="objType", description="Object type", default="na")
    assigned: bool = Field(description="Assigned flag", default=False)

    # Status flags
    conflicted: bool = Field(description="Conflicted flag", default=False)
    ranout: bool = Field(description="Ran out flag", default=False)
    outside: bool = Field(description="Outside flag", default=False)

    # Position coordinates
    x_focal: float = Field(description="X focal plane coordinate", alias="xfocal", default=float('NaN'))
    y_focal: float = Field(description="Y focal plane coordinate", alias="yfocal", default=float('NaN'))
    xf_default: float = Field(description="Default X focal coordinate", default=float('NaN'))
    yf_default: float = Field(description="Default Y focal coordinate", default=float('NaN'))

    # Target coordinates
    ra: float = Field(description="Right ascension [deg]")
    dec: float = Field(description="Declination [deg]")

    # Wavelength information
    lambda_eff: float = Field(description="Effective wavelength", default=0.0)
    zoffset: float = Field(description="Z offset", default=0.0)

    # Instrument identifiers
    spectrograph_id: int = Field(alias="spectrographId", description="Spectrograph ID", default=-1)
    fiber_id: int = Field(alias="fiberId", description="Fiber ID", ge=1, le=300)
    planned_fiber_id: int = Field(alias="fiberid", description="Fiber ID", default=-1)
    throughput: int = Field(description="Throughput value", default=-1)

    # Plate-specific information
    iplateinput: int = Field(description="Plate input ID", default=-1)
    pointing: int = Field(description="Pointing number", default=-1)
    offset: int = Field(description="Offset value", default=-1)
    block: int = Field(description="Block number", default=-1)
    iguide: int = Field(description="Guide flag", default=-1)
    bluefiber: int = Field(description="Blue fiber flag", default=-1)
    chunk: int = Field(description="Chunk number", default=-1)
    ifinal: int = Field(description="Final flag", default=-1)
    plugged_mjd: int = Field(description="MJD when this plate was plugged", default=-1) # necessary for fiber mapping fixes
    fix_fiber_flag: int = Field(default=0, description="Whether this fiber mapping was fixed in software")

    # Physical properties
    diameter: float = Field(default=-1, description="Diameter")
    buffer: float = Field(default=-1, description="Buffer size")
    priority: int = Field(default=-1, description="Target priority")


    @validator('hole_type', 'planned_hole_type', 'obj_type', pre=True)
    def enforce_lower_case(cls, v):
        return v.lower()

    @validator('category', pre=True)
    def validate_category(cls, v):
        # Make consistent across FPS and plate era.
        translate_from_plate_to_fps = {
            "sky": "sky_apogee",
            "standard": "standard_apogee",
            "na": ""
        }
        return translate_from_plate_to_fps.get(v.lower(), v)

    class Config:
        validate_by_name = True
        validate_assignment = True
        arbitrary_types_allowed = True

    @property
    def expected_to_be_assigned_sdss_id(self) -> bool:
        """ A helper function so we don't try to cross-match sky targets for SDSS IDs. """
        return (
            (self.catalogid > 0) | (self.twomass_designation != "")

            and not self.category.startswith("sky_")
            and self.category != ""
        )


    @model_validator(mode="after")
    def fix_fiber_mappings(self):
        # Dates from /uufs/chpc.utah.edu/common/home/sdss09/software/apogee/Linux/apogee/trunk/data/cal/apogee-n.par
        # The Python mapping logic originates from:
        #   https://github.com/sdss/apogee_drp/blob/4ab6a04e448b279f2514550802b6732693e9847a/python/apogee_drp/utils/plugmap.py#L210-L238
        # but that code has a bug in it (see below)

        # Check against fix_fiber_flag to avoid recursively fixing things.
        if self.observatory == "apo" and self.fix_fiber_flag == 0:
            if 56764 <= self.plugged_mjd <= 56773 and self.fiber_id >= 0:
                self.fix_fiber_flag = 1
                sub_id = (self.fiber_id - 1) % 30
                bundle_id = (self.fiber_id - sub_id) // 30
                self.fiber_id = (9 - bundle_id) * 30 + sub_id + 1
            if (
                (58034 <= self.plugged_mjd <= 58046)
            and (self.spectrograph_id == 2)
            #and (self.hole_type == "object")
            ):
                self.fix_fiber_flag = 2
                # Note that the DRP code has this in a way where it ONLY changes
                # fibers 31, 37, 45, and 54, but their expressions are written in
                # a way that you would think they are meant to be ranges.

                # This is because the DRP code has a transcription error from the
                # original IDL code:
                #   https://github.com/sdss/apogee/blob/master/pro/apogeereduce/aploadplugmap.pro#L210-L221

                # There is a thinko here you should be aware of.
                # > Exposure(observatory=apo, mjd=58037, exposure=11, image_type=object)

                # When this is loaded there are 298 fibers. The `Exposure` code ensures
                # there are 300, so it adds two more and marks them as unplugged.

                # In an earlier version of this code I would check for "hole_type" == "object" to
                # match the DRP logic. However, those two extra fibers are marked as "unplugged",
                # which means some of the fibers get the offset applied and some do not.
                # This results in two fibers being mapped to fiber ID 36, and one fiber is now missing.
                # That obviously causes problems downstream.

                # For now I am going to remove the hole_type check, but we should look at this
                # problem more carefully with fresh eyes.

                # TODO: as above.



                offset_ranges = [
                    (31, 36, +23),
                    (37, 44, +8),
                    (45, 52, -8),
                    (54, 59, -23),
                ]
                for lower, upper, offset in offset_ranges:
                    if (lower <= self.fiber_id) & (self.fiber_id <= upper):
                        self.fiber_id += offset
                        break

                if self.fiber_id in (53, 60):
                    # The DRP indicates that these are unpopulated fibers.
                    self.hole_type = "unplugged"
                    self.category = "unplugged"
                    self.planned_hole_type = "unplugged"
                    self.ra = np.nan
                    self.dec = np.nan
                    self.obj_type = "na"

        return self
