import os
import numpy as np
from astropy.table import Table
from functools import partial, cached_property
from pydantic import BaseModel, Field, computed_field, field_validator, model_validator
from typing import Optional, Tuple, Union

from almanac import logger, config, utils

from almanac.data_models.fps import FPSTarget
from almanac.data_models.plate import PlateTarget
from almanac.data_models.types import *
from almanac.data_models.utils import (
    get_headers, match_planned_to_plugged, get_exposure_path,
    get_exposure_string, mjd_to_exposure_prefix
)

from almanac.qa import lookup_bad_exposures


class Exposure(BaseModel):

    #> Basic Information
    observatory: Observatory = Field(description="Observatory name")
    mjd: int = Field(description="MJD of the exposure")
    exposure: int = Field(description="Exposure number", ge=1)
    prefix: Optional[Prefix] = Field(description="Raw exposure basename prefix", default=None)

    #> Exposure Information
    name: Optional[str] = Field(
        default="",
        description=(
            "The `name` field in the exposure header often refers to the plugged"
            "plate name, which describes which targets were observed."
        )
    )
    n_read: int = Field(default=0, alias="nread", ge=0)
    image_type: ImageType = Field(alias="imagetyp")
    observer_comment: Optional[str] = Field(default="", alias="obscmnt")

    #> Identifiers
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

    #> File Metadata
    file_size: int = Field(default=-1, description="File size in bytes")

    _targets: Optional[Tuple[Union[FPSTarget, PlateTarget]]] = None

    @computed_field(description="Exposure string used in path")
    def exposure_string(self) -> str:
        return get_exposure_string(self.mjd, self.exposure)

    @computed_field(description="Whether this exposure is from the FPS era")
    def fps(self) -> bool:
        start = dict(apo=59423, lco=59809)[self.observatory]
        return self.mjd >= start

    #@computed_field(description="FPI")
    #def fpi(self) -> bool:
    #    return "fpi" in self.observer_comment.lower()

    #@computed_field(description="Sparse pak mode")
    #def sparse_pak(self) -> bool:
    #    return "sparse" in self.observer_comment.lower()

    def __str__(self):
        return (
            f"<Exposure(observatory={self.observatory},"
            f" mjd={self.mjd}, exposure={self.exposure},"
            f" image_type={self.image_type})>"
        )

    @computed_field
    def flagged_bad(self) -> bool:
        marked_bad = (self.observatory, self.mjd, self.exposure) in lookup_bad_exposures
        missing_plate_info = (
            (self.image_type == "object")
        &   ((self.fps and self.config_id <= 0) or (not self.fps and self.plate_id <= 0))
        )
        return marked_bad or missing_plate_info


    @computed_field
    def chip_flags(self) -> int:
        return int(np.sum(2**np.where(list(map(os.path.exists, self.paths)))[0]))

    # Validations
    @field_validator('prefix', mode="before")
    def validate_prefix(cls, v, values):
        if v is None:
            return dict(apo="apR", lco="asR").get(values.get("observatory"))
        return v

    @field_validator("observer_comment", mode="before")
    def sanitise_observer_comment(cls, v) -> str:
        if v is None or str(v).lower() == "none":
            return ""
        return v.strip()

    @field_validator('image_type', mode="before")
    def validate_descriptive_type(cls, v):
        return v.lower()

    @field_validator('cart_id', mode="before")
    def validate_cart_id(cls, v):
        if isinstance(v, str) and v.strip().upper() == 'FPS':
            return 0
        return empty_string_to_int(v, -1)

    @field_validator('map_id', 'plate_id', 'field_id', 'design_id', 'config_id', mode="before")
    def validate_identifiers(cls, v):
        return empty_string_to_int(v, -1)

    @field_validator('seeing', 'focus', 'collpist', 'colpitch', 'dithered_pixels', mode="before")
    def validate_floats(cls, v):
        try:
            return float(v)
        except:
            return float('NaN')

    @field_validator('lamp_quartz', 'lamp_thar', 'lamp_une', mode="before")
    def validate_lamps(cls, v):
        return {'F': 0, 'T': 1, '0': 0, '1': 1}.get(str(v).strip().upper(), -1)

    @model_validator(mode="after")
    def check_fields(self):
        if self.observer_comment is not None and self.image_type != "twilightflat":
            sanitised = self.observer_comment.lower().replace(' ', '')
            if 'skyflat' in sanitised or 'twilight' in sanitised:
                self.image_type = 'twilightflat'
        return self

    # TODO: we may want to change this to be way more descriptive, particularly
    #       when we start doing QA to make sure exposures look like they should
    @property
    def qa_metadata(self) -> Optional[dict]:
        print("Warning: The `qa_metadata` property will change in the future")
        return lookup_bad_exposures.get((self.observatory, self.mjd, self.exposure), None)

    @property
    def plugged_mjd(self) -> int:
        try:
            return int(self.name.split("-")[1])
        except:
            return -1

    @property
    def plugged_iteration(self) -> str:
        try:
            return self.name.split("-")[2]
        except:
            return ""

    @property
    def paths(self) -> Tuple[str]:
        return tuple(
            map(
                partial(
                    get_exposure_path,
                    self.observatory,
                    self.mjd,
                    self.prefix,
                    self.exposure
                ),
                "abc"
            )
        )

    @property
    def plate_hole_path(self):
        return (
            f"{config.platelist_dir}/"
            f"{str(self.plate_id)[:-2].zfill(4)}XX/"
            f"{self.plate_id:0>6.0f}/"
            f"plateHoles-{self.plate_id:0>6.0f}.par"
        )

    @property
    def plug_map_path(self):
        return (
            f"{config.mapper_dir}/"
            f"{self.observatory}/"
            f"{self.plugged_mjd}/"
            f"plPlugMapM-{self.plate_id}-{self.plugged_mjd}-{self.plugged_iteration}.par"
        )

    @property
    def config_summary_path(self):
        directory = (
            f"{config.sdsscore_dir}/"
            f"{self.observatory}/"
            f"summary_files/"
            f"{str(self.config_id)[:-3].zfill(3)}XXX/"
            f"{str(self.config_id)[:-2].zfill(4)}XX/"
        )
        # fall back to confSummary if confSummaryFS does not exist
        for flavor in ("FS", ""):
            path = f"{directory}/confSummary{flavor}-{self.config_id}.par"
            if os.path.exists(path):
                return path

        raise FileNotFoundError(f"Could not find confSummary file for config {self.config_id} in {directory}")

    @classmethod
    def from_keys(
        cls,
        mjd: int,
        observatory: str,
        exposure: int,
        prefix: Optional[Prefix] = None,
        chip: Optional[Chip] = None
    ) -> "Exposure":
        """
        Create an Exposure instance from basic identifying keys.

        :param mjd:
            MJD of the exposure.

        :param observatory:
            Observatory name (e.g., 'apo', 'lco').

        :param exposure:
            Exposure number.

        :param prefix: [optional]
            Prefix for the exposure file (e.g., 'apR', 'asR'). If not provided,
            defaults to 'apR' for APO and 'asR' for LCO.

        :param chip: [optional]
            Chip identifier ('a', 'b', or 'c'). If not provided, will check
            all chips in order.

        :returns:
            An instance of the Exposure class populated with data extracted from
            the file headers.
        """

        if prefix is None:
            prefix = dict(apo="apR", lco="asR").get(observatory, "apR")

        for chip in (chip or "abc"):
            path = get_exposure_path(observatory, mjd, prefix, exposure, chip)
            if os.path.exists(path):
                headers = get_headers(path)
                file_size = os.path.getsize(path)
                return cls(
                    observatory=observatory,
                    mjd=mjd,
                    exposure=exposure,
                    prefix=prefix,
                    file_size=file_size,
                    **headers
                )
        raise FileNotFoundError(f"No exposure files found for {observatory} {mjd} {exposure} {prefix}")

    @classmethod
    def from_path(cls, path: str) -> "Exposure":
        """
        Create an Exposure instance from a given file path.

        :param path:
            Full path to the exposure file.

        :returns:
            An instance of the Exposure class populated with data extracted from
            the file path and headers.
        """
        *_, observatory, mjd, basename = path.split("/")
        prefix, chip, cumulative_exposure = basename.split("-")
        exposure = int(cumulative_exposure.split(".")[0]) - mjd_to_exposure_prefix(mjd)
        headers = get_headers(path)
        file_size = os.path.getsize(path)
        return cls(
            observatory=observatory,
            mjd=mjd,
            exposure=exposure,
            prefix=prefix,
            file_size=file_size,
            **headers
        )

    @cached_property
    def headers(self) -> dict:
        for path in self.paths:
            if os.path.exists(path):
                return get_headers(path)
        raise FileNotFoundError(f"No exposure files found for {self.observatory} {self.mjd} {self.exposure} {self.prefix}")

    @cached_property
    def fiber_map(self):
        return Table.read(self.config_summary_path, format="yanny", tablename="FIBERMAP")

    @cached_property
    def plug_map(self):
        return Table.read(self.plug_map_path, format="yanny", tablename="PLUGMAPOBJ")

    @cached_property
    def plate_hole_map(self):
        return Table.read(self.plate_hole_path, format="yanny", tablename="STRUCT1")


    @cached_property
    def targets(self) -> Tuple[Union[FPSTarget, PlateTarget]]:
        if self._targets is None:
            if (
                (self.image_type == "object")
            &   (
                    (self.fps and self.config_id > 0)
                |   (not self.fps and self.plate_id > 0)
            )
            ):
                if self.fps:
                    factory = FPSTarget
                    rows = self.fiber_map
                    keep = (rows["fiberType"] == "APOGEE") & (rows["fiberId"] > 0)
                    rows = rows[keep]
                    missing = set(range(1, 301)) - set(rows["fiberId"])
                    for fiber_id in missing:
                        rows.add_row({
                            "fiberId": fiber_id,
                            "racat": np.nan,
                            "deccat": np.nan,
                            "category": "bonus",
                        })
                else:
                    factory = PlateTarget
                    bad_exposure_notes = (
                        lookup_bad_exposures
                        .get((self.observatory, self.mjd, self.exposure), {})
                        .get("notes", None)
                    )
                    if bad_exposure_notes == "missing_plug_map_file":
                        rows = []
                    else:
                        rows = match_planned_to_plugged(self.plate_hole_map, self.plug_map)
                        if rows:
                            # Plugged MJD is necessary to understand where the fiber mapping
                            # went wrong in early plate era.
                            rows["plugged_mjd"] = self.plugged_mjd
                            rows["observatory"] = self.observatory
                try:
                    rows.sort("fiberId")
                except:
                    None

                try:
                    self._targets = tuple([factory(**r) for r in rows])
                except Exception as e:
                    e.add_note(f"Originated from {self}")
                    raise
            else:
                self._targets = tuple()
        return self._targets

    def __repr__(self):
        return f"{self.__repr_name__()}(observatory={self.observatory}, mjd={self.mjd}, exposure={self.exposure}, image_type={self.image_type})"

    class Config:
        validate_by_name = True
        validate_assignment = True


def empty_string_to_int(v, default) -> int:
    if isinstance(v, str) and v.strip() == '':
        return default
    elif v is None:
        return default
    return int(v)
