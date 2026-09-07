"""Tests for reading exposures and sequences back from an almanac file."""

import h5py as h5
import numpy as np

from almanac import io
from tests.helpers import MJD, write_night


class TestReadBack:
    def test_read_exposures_round_trip(self, tmp_path):
        path = str(tmp_path / "night.h5")
        written = write_night(path)
        with h5.File(path, "r") as fp:
            read = io.read_exposures(fp, "apo", MJD)

        assert len(read) == len(written)
        for a, b in zip(read, written):
            assert a.observatory == b.observatory
            assert a.mjd == b.mjd
            assert a.exposure == b.exposure
            assert a.image_type == b.image_type
            assert a.config_id == b.config_id
            assert a.plate_id == b.plate_id
            assert a.field_id == b.field_id
            assert a.name == b.name
            assert a.observer_comment == b.observer_comment
            assert a.lamp_thar == b.lamp_thar
            assert a.fps == b.fps
            # Unset on write; filled in from the observatory on read.
            assert a.prefix == "apR"
            # Floats: NaN-safe comparison
            assert np.isclose(a.seeing, b.seeing, equal_nan=True)

    def test_read_exposures_missing_group(self, tmp_path):
        path = str(tmp_path / "empty.h5")
        with h5.File(path, "a") as fp:
            assert io.read_exposures(fp, "apo", MJD) == []
            assert io.read_sequences(fp, "apo", MJD) == {}

    def test_read_sequences_round_trip(self, tmp_path):
        path = str(tmp_path / "night.h5")
        write_night(path)
        with h5.File(path, "r") as fp:
            sequences = io.read_sequences(fp, "apo", MJD)
        assert sequences["objects"] == [(1, 2)]
        assert sequences["arclamps"] == [(3, 3)]
        assert sequences["missing"] == []
