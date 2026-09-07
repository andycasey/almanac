"""Shared helpers for the `almanac add` command tests."""

import h5py as h5
from click.testing import CliRunner

from almanac import io
from almanac.cli import main
from almanac.data_models import Exposure
from almanac.data_models.fps import FPSTarget

# FPS-era MJD (so `config_id` keys the fibers group). No raw data exists for
# it in the test environment, and nothing below touches the database.
MJD = 60000


def run_cli(args):
    runner = CliRunner()
    return runner.invoke(main, args, catch_exceptions=False)


def make_exposures(observatory="apo", mjd=MJD):
    """Two science exposures on one configuration, followed by an arc lamp."""
    return [
        Exposure(
            observatory=observatory,
            mjd=mjd,
            exposure=1,
            image_type="object",
            config_id=123,
            field_id=456,
            name="first",
            observer_comment="a comment",
            seeing=1.5,
        ),
        Exposure(
            observatory=observatory,
            mjd=mjd,
            exposure=2,
            image_type="object",
            config_id=123,
            field_id=456,
        ),
        Exposure(
            observatory=observatory,
            mjd=mjd,
            exposure=3,
            image_type="arclamp",
            lamp_thar=1,
        ),
    ]


SEQUENCES = {"objects": [(1, 2)], "arclamps": [(3, 3)], "missing": []}


def write_night(path, observatory="apo", mjd=MJD):
    exposures = make_exposures(observatory, mjd)
    with h5.File(path, "a") as fp:
        io.update(fp, observatory, mjd, exposures, SEQUENCES, fibers=False)
    return exposures


def fake_targets(exposure):
    return (
        FPSTarget(category="science", racat=10.0, deccat=-20.0, fiberId=1, catalogid=5),
        FPSTarget(category="sky_apogee", racat=11.0, deccat=-21.0, fiberId=2),
    )
