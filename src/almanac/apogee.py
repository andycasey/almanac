import os
import numpy as np
from glob import glob
from subprocess import check_output
from astropy.table import Table, hstack, unique
from itertools import groupby
from typing import Optional, Tuple, Dict, List, Set, Generator, Any, Union

from scipy.spatial.distance import cdist

from almanac import config, logger, utils
from almanac.data_models import Exposure
from almanac.data_models.types import ImageType
from almanac.data_models.utils import mjd_to_exposure_prefix, get_exposure_path


APOGEE_ID_LOOKUP = None

def create_apogee_id_lookup() -> Dict[str, int]:
    """
    Create a lookup dictionary mapping APOGEE IDs to SDSS IDs.
    """
    global APOGEE_ID_LOOKUP
    """
    import pickle
    with open("apogee_id_lookup.pkl", "rb") as f:
        foo = pickle.load(f)
    APOGEE_ID_LOOKUP = {}
    for key, value in foo.items():
        APOGEE_ID_LOOKUP[key] = value
        for prefix in ("2M", "J", "AP"):
            if key.startswith(prefix):
                APOGEE_ID_LOOKUP[key[len(prefix):]] = value

    return APOGEE_ID_LOOKUP
    """

    from almanac.database import is_database_available, catalogdb

    from peewee import (ForeignKeyField, BooleanField, FloatField, TextField, IntegerField)

    class CatalogToAllStar_DR17_synspec_rev1(catalogdb.CatalogdbModel):
        catalog = ForeignKeyField(catalogdb.Catalog,
                                column_name='catalogid',
                                backref='+')
        target = ForeignKeyField(catalogdb.AllStar_DR17_synspec_rev1,
                                column_name='target_id',
                                backref='+')
        version = ForeignKeyField(catalogdb.Version,
                                column_name='version_id',
                                backref='+')
        best = BooleanField()
        distance = FloatField()
        plan_id = TextField(null=True)

        class Meta:
            table_name = "catalog_to_allstar_dr17_synspec_rev1"
            primary_key = False
            use_reflection = False


    q = (
        catalogdb.SDSS_ID_flat
        .select(
            catalogdb.SDSS_ID_flat.sdss_id,
            catalogdb.AllStar_DR17_synspec_rev1.apogee_id
        )
        .join(
            catalogdb.SDSS_ID_To_Catalog,
            on=(
                catalogdb.SDSS_ID_flat.catalogid
                == catalogdb.SDSS_ID_To_Catalog.catalogid
            )
        )
        .join(
            CatalogToAllStar_DR17_synspec_rev1,
            on=(
                catalogdb.SDSS_ID_To_Catalog.catalogid
                == CatalogToAllStar_DR17_synspec_rev1.catalogid
            ),
        )
        .join(
            catalogdb.AllStar_DR17_synspec_rev1,
            on=(
                CatalogToAllStar_DR17_synspec_rev1.target_id
                == catalogdb.AllStar_DR17_synspec_rev1.apstar_id
            ),
        )
        .where(catalogdb.SDSS_ID_flat.rank == 1)
        .tuples()
    )
    APOGEE_ID_LOOKUP = {}
    for sdss_id, apogee_id in q.iterator():
        if len(apogee_id.strip()) <= 2:
            # There is one problematic entry of "AP" which caused 180,000
            # fiber-exposures to be assigned to sdss_id 57991461
            continue
        APOGEE_ID_LOOKUP[apogee_id] = sdss_id
        for prefix in ("2M", "J", "AP"):
            if apogee_id.startswith(prefix):
                APOGEE_ID_LOOKUP[apogee_id[len(prefix):]] = sdss_id
    return APOGEE_ID_LOOKUP


def get_unique_exposure_paths(paths: List[str]) -> List[str]:
    """
    Process a list of file paths to find unique exposures and determine which chips are available.

    :param paths:
        List of file paths to APOGEE exposure files.

    :returns:
        List of exposure paths.
    """

    chip_mapping = {}
    for path in paths:
        _, observatory, mjd, basename = path.rsplit("/", 3)
        prefix, chip, exposure_apz = basename.split("-")

        key = (observatory, mjd, exposure_apz)
        chip_mapping.setdefault(key, [prefix, [False, False, False]])
        index = "abc".index(chip)
        chip_mapping[key][1][index] = True

    unique_exposure_paths = []
    for (observatory, mjd, exposure_apz), (prefix, chips) in chip_mapping.items():
        chip = "abc"[chips.index(True)]
        path = f"{config.apogee_dir}/{observatory}/{mjd}/{prefix}-{chip}-{exposure_apz}"
        unique_exposure_paths.append(path)

    return unique_exposure_paths


def get_exposures(observatory: str, mjd: int) -> Generator[Exposure, None, None]:
    """
    Generate exposures taken from a given observatory on a given MJD.

    :param observatory:
        The observatory name (e.g. "apo").

    :param mjd:
        The Modified Julian Date.

    :yields:
        Exposure instances for each unique exposure found on disk.
    """
    paths = glob(get_exposure_path(observatory, mjd, "a?R", "*", "*"))
    return organize_exposures(map(Exposure.from_path, get_unique_exposure_paths(paths)))


def get_expected_number_of_exposures(observatory: str, mjd: int) -> int:
    """
    Query the SDSS database to get the expected exposures for a given observatory and MJD.
    This is useful for identifying missing exposures.
    """

    if mjd < int(getattr(config.sdssdb_exposure_min_mjd, observatory)):
        return -1

    from almanac.database import opsdb
    from peewee import fn

    for model in (opsdb.Exposure, opsdb.ExposureFlavor):
        model._meta.schema = f"opsdb_{observatory}"

    start, end = map(mjd_to_exposure_prefix, (mjd, mjd + 1))

    q = (
        opsdb.Exposure.select(
            fn.max(opsdb.Exposure.exposure_no)
        )
        .where(
            (opsdb.Exposure.exposure_no > start)
        &   (opsdb.Exposure.exposure_no < end)
        )
        .join(
            opsdb.ExposureFlavor,
            on=(opsdb.ExposureFlavor.pk == opsdb.Exposure.exposure_flavor),
        )
    )
    try:
        return q.scalar() - start
    except:
        return -1


def organize_exposures(exposures: List[Exposure]) -> List[Exposure]:
    """
    Identify any missing exposures (based on non-contiguous exposure numbers)
    and fill them with missing image types.

    :param exposures:
        A list of `Exposure` instances.

    :returns:
        A list of organized `Exposure` instances.
    """

    exposures = sorted(exposures, key=lambda x: x.exposure)

    if len(exposures) == 0:
        return []

    observatory, mjd = (exposures[0].observatory, exposures[0].mjd)

    n_expected = get_expected_number_of_exposures(observatory, mjd)
    max_exposure = max(exposures[-1].exposure, n_expected)

    organized = []
    for i in range(1, max_exposure + 1):
        if exposures and exposures[0].exposure == i:
            organized.append(exposures.pop(0))
        else:
            organized.append(
                Exposure(
                    observatory=observatory,
                    exposure=i,
                    mjd=mjd,
                    image_type="missing"
                )
            )
    return organized


def get_sequences(
    exposures: List[Exposure],
    image_type: ImageType,
    fields: Tuple[str, ...],
    ignore_interrupting_darks: bool = False
) -> List[Tuple[int, int]]:
    """
    Get exposure number ranges for sequences of a specific image type.

    :param exposures:
        Astropy Table containing exposure metadata.
    :param image_type:
        The image type to search for (e.g., "Object", "ArcLamp").
    :param fields:
        Tuple of column names to group exposures by.
    :param ignore_interrupting_darks: [optional]
        Whether to allow (and ignore) dark exposures interrupting sequences.

    :returns:
        List of tuples containing (start_exposure, end_exposure) for each sequence.
    """
    if ignore_interrupting_darks:
        f = lambda x: x.image_type in (image_type, "dark")
    else:
        f = lambda x: x.image_type == image_type

    s = list(filter(f, exposures))

    sequence_exposure_numbers = []
    for v, group in groupby(s, key=lambda x: tuple(getattr(x, f) for f in fields)):
        for si, ei in utils.group_contiguous([e.exposure for e in group]):
            sequence_exposure_numbers.append((si, ei))
    return sequence_exposure_numbers


def get_arclamp_sequences(exposures: List[Exposure]) -> List[Tuple[int, int]]:
    """
    Return a list of tuples indicating the start and end exposure numbers for
    a sequence of arc lamp exposures.

    :param exposures:
        A list of `Exposure` instances.

    :returns:
        List of tuples containing (start_exposure, end_exposure) for each arc lamp sequence.
    """
    return get_sequences(
        exposures,
        "arclamp",
        ("dithered_pixels", ),
        ignore_interrupting_darks=True
    )


def get_science_sequences(exposures: List[Exposure]) -> List[Tuple[int, int]]:
    """
    Return a list of tuples indicating the start and end exposure numbers for
    a sequence of science exposures.

    :param exposures:
        A list of `Exposure` instances.

    :returns:
        List of tuples containing (start_exposure, end_exposure) for each science sequence.
    """
    return get_sequences(exposures, "object", ("field_id", "plate_id", "config_id", "image_type"))


def parse_target_identifier(target):

    # catalogids prior to 4204682057 don't have a direct sdss_id mapping.
    if target.catalogid >= 4204682057:
        return ("catalogid", target.catalogid)
    else:
        if target.twomass_designation.startswith("GAIA_DR2"):
            return ("gaia_dr2", np.int64(target.twomass_designation.split(" ")[1]))
        else:
            designation = target.twomass_designation
            for prefix in ("2MASS-J", "2MASS", "2M", "J", "APG-J", "AP"):
                if designation.startswith(prefix):
                    designation = designation[len(prefix):]
                    break
            return ("coordinate", designation)


def get_almanac_data(
    observatory: str,
    mjd: int,
    fibers: bool = False,
    meta: bool = False,
    only_assign_sdss_id_to_first_exposure_per_configuration: bool = True
):
    """
    Return comprehensive almanac data for all exposures taken from a given observatory on a given MJD.

    :param observatory:
        The observatory name (e.g. "apo").
    :param mjd:
        The Modified Julian Date.
    :param fibers:
        Whether to include fiber mapping information.
    :param xmatch:
        Whether to perform cross-matching with catalog database.
    :param kwargs:
        Additional keyword arguments passed to other functions.

    :returns:
        Tuple containing:
        - observatory name
        - MJD
        - A list of exposures
        - Table of exposure data
        - dictionary of sequence indices
        - dictionary of fiber mappings
    """

    exposures = get_exposures(observatory, mjd)
    sequences = {
        "objects": get_science_sequences(exposures),
        "arclamps": get_arclamp_sequences(exposures),
    }
    if fibers:
        identifiers = {}
        # We only need to get targets for one exposure in each science sequence.
        for si, ei in sequences["objects"]:
            exposure = exposures[si - 1]
            for target in exposure.targets:
                key, identifier = parse_target_identifier(target)
                identifiers.setdefault(key, set()).add(identifier)

        if meta:
            # We will often run `get_almanac_data` in parallel (through multiple processes),
            # so here we are avoiding opening a database connection until the child process starts.
            from almanac.database import is_database_available, catalogdb

            lookup = {}
            if "catalogid" in identifiers and is_database_available:
                q = (
                    catalogdb.SDSS_ID_flat
                    .select(
                        catalogdb.SDSS_ID_flat.sdss_id,
                        catalogdb.SDSS_ID_flat.catalogid
                    )
                    .where(
                        catalogdb.SDSS_ID_flat.catalogid.in_(tuple(identifiers["catalogid"]))
                    &   (catalogdb.SDSS_ID_flat.rank == 1)
                    )
                    .tuples()
                )
                lookup["catalogid"] = {}
                for sdss_id, catalogid in q:
                    lookup["catalogid"][catalogid] = sdss_id

            if "gaia_dr2" in identifiers and is_database_available:
                q = (
                    catalogdb.SDSS_ID_flat
                    .select(
                        catalogdb.SDSS_ID_flat.sdss_id,
                        catalogdb.CatalogToGaia_DR2.target
                    )
                    .join(
                        catalogdb.CatalogToGaia_DR2,
                        on=(
                            catalogdb.SDSS_ID_flat.catalogid == catalogdb.CatalogToGaia_DR2.catalog
                        )
                    )
                    .where(
                        catalogdb.CatalogToGaia_DR2.target.in_(tuple(identifiers["gaia_dr2"]))
                    &   (catalogdb.SDSS_ID_flat.rank == 1)
                    )
                    .tuples()
                )
                lookup["gaia_dr2"] = {}
                for sdss_id, gaia_dr2_source_id in q:
                    lookup["gaia_dr2"][gaia_dr2_source_id] = sdss_id


            if "coordinate" in identifiers and is_database_available:
                global APOGEE_ID_LOOKUP
                if APOGEE_ID_LOOKUP is None:
                    APOGEE_ID_LOOKUP = create_apogee_id_lookup()

                lookup["apogee_id"] = APOGEE_ID_LOOKUP
                remaining = set(identifiers["coordinate"]).difference(lookup["apogee_id"].keys())

                if remaining:
                    # For remaining things, check 2MASS.
                    q = (
                        catalogdb.SDSS_ID_flat
                        .select(
                            catalogdb.TwoMassPSC.designation,
                            catalogdb.SDSS_ID_flat.sdss_id
                        )
                        .join(
                            catalogdb.CatalogToTwoMassPSC,
                            on=(
                                catalogdb.SDSS_ID_flat.catalogid
                                == catalogdb.CatalogToTwoMassPSC.catalogid
                            ),
                        )
                        .join(
                            catalogdb.TwoMassPSC,
                            on=(
                                catalogdb.CatalogToTwoMassPSC.target_id
                                == catalogdb.TwoMassPSC.pts_key
                            ),
                        )
                        .where(
                            catalogdb.TwoMassPSC.designation.in_(tuple(remaining))
                        &   (catalogdb.SDSS_ID_flat.rank == 1)
                        )
                        .tuples()
                    )
                    lookup["twomass_psc"] = { d: s for d, s in q }

            # Add sdss_id to targets
            for si, ei in sequences["objects"]:
                iterable = [si - 1] if only_assign_sdss_id_to_first_exposure_per_configuration else range(si - 1, ei)
                for i in iterable:
                    exposure = exposures[i]
                    for target in exposure.targets:
                        if target.category == "unplugged":
                            # Don't assign sdss_id to unplugged targets.
                            continue

                        key, identifier = parse_target_identifier(target)
                        if key == "coordinate":
                            for lookup_key in ("apogee_id", "twomass_psc"):
                                target.sdss_id = lookup[lookup_key].get(identifier, -1)
                                if target.sdss_id != -1:
                                    break
                        else:
                            target.sdss_id = lookup[key].get(identifier, -1)

                        #if target.sdss_id == -1 and (target.category not in ("sky_apogee", "unplugged")):
                        #    print(f"{observatory} {mjd} Exposure {exposure.exposure}: Could not find SDSS ID for target {identifier} ({key}) {target}")

    return (observatory, mjd, exposures, sequences)
