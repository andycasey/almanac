"""Tests for `almanac add fibers` and the `cross_match_targets` refactor."""

import h5py as h5
import pytest

import almanac.apogee as apogee
from tests.helpers import (
    MJD, SEQUENCES, fake_targets, make_exposures, run_cli, write_night
)


class TestAddFibers:
    def test_adds_fibers_group_from_existing_exposures(self, tmp_path, monkeypatch):
        path = str(tmp_path / "night.h5")
        write_night(path)

        calls = []

        def fake_cross_match(exposures, science_sequences, xmatch=True, **kwargs):
            calls.append((list(science_sequences), xmatch))
            for si, _ in science_sequences:
                exposures[si - 1]._targets = fake_targets(exposures[si - 1])

        monkeypatch.setattr(apogee, "cross_match_targets", fake_cross_match)
        # The raw exposure scan must NOT be used: exposures come from the file.
        monkeypatch.setattr(
            apogee,
            "get_almanac_data",
            lambda *a, **k: pytest.fail("raw exposures should not be re-read"),
        )

        result = run_cli(["add", "fibers", path])
        assert result.exit_code == 0
        assert calls == [([(1, 2)], True)]

        with h5.File(path, "r") as fp:
            group = fp[f"raw/apo/{MJD}/fibers"]
            assert list(group) == ["123"]
            assert list(group["123/fiber_id"][:]) == [1, 2]
            assert list(group["123/sdss_id"][:]) == [-1, -1]
            # Exposures and sequences are left untouched.
            assert len(fp[f"raw/apo/{MJD}/exposures/exposure"]) == 3
            assert f"raw/apo/{MJD}/sequences/objects" in fp

    def test_no_x_match_is_passed_through(self, tmp_path, monkeypatch):
        path = str(tmp_path / "night.h5")
        write_night(path)

        calls = []

        def fake_cross_match(exposures, science_sequences, xmatch=True, **kwargs):
            calls.append(xmatch)
            for exposure in exposures:
                exposure._targets = ()

        monkeypatch.setattr(apogee, "cross_match_targets", fake_cross_match)
        result = run_cli(["add", "fibers", path, "--no-x-match"])
        assert result.exit_code == 0
        assert calls == [False]

    def test_selects_nights_by_mjd_and_observatory(self, tmp_path, monkeypatch):
        path = str(tmp_path / "nights.h5")
        write_night(path, "apo", MJD)
        write_night(path, "apo", MJD + 1)
        write_night(path, "lco", MJD)

        seen = []

        def fake_cross_match(exposures, science_sequences, xmatch=True, **kwargs):
            seen.append((exposures[0].observatory, exposures[0].mjd))
            for exposure in exposures:
                exposure._targets = ()

        monkeypatch.setattr(apogee, "cross_match_targets", fake_cross_match)

        result = run_cli(["add", "fibers", path, "--mjd", str(MJD), "--apo"])
        assert result.exit_code == 0
        assert seen == [("apo", MJD)]

        seen.clear()
        result = run_cli(["add", "fibers", path])
        assert result.exit_code == 0
        assert sorted(seen) == [("apo", MJD), ("apo", MJD + 1), ("lco", MJD)]

    def test_failed_night_does_not_kill_the_run(self, tmp_path, monkeypatch):
        path = str(tmp_path / "nights.h5")
        write_night(path, "apo", MJD)
        write_night(path, "apo", MJD + 1)

        def fake_cross_match(exposures, science_sequences, xmatch=True, **kwargs):
            if exposures[0].mjd == MJD:
                raise FileNotFoundError("no confSummary file")
            for si, _ in science_sequences:
                exposures[si - 1]._targets = fake_targets(exposures[si - 1])

        monkeypatch.setattr(apogee, "cross_match_targets", fake_cross_match)
        result = run_cli(["add", "fibers", path])
        assert result.exit_code == 0
        with h5.File(path, "r") as fp:
            assert f"raw/apo/{MJD + 1}/fibers/123" in fp
            assert f"raw/apo/{MJD}/fibers/123" not in fp

    def test_empty_file_warns_and_exits_cleanly(self, tmp_path):
        path = str(tmp_path / "empty.h5")
        with h5.File(path, "a"):
            pass
        result = run_cli(["add", "fibers", path])
        assert result.exit_code == 0


class TestCrossMatchRefactor:
    def test_get_almanac_data_with_fibers_and_no_exposures(self):
        # MJD far before sdssdb coverage: no raw data and no database access.
        observatory, mjd, exposures, sequences, missing = apogee.get_almanac_data(
            "apo", 51000, fibers=True, meta=False
        )
        assert exposures == []
        assert sequences["objects"] == []

    def test_cross_match_targets_without_xmatch_loads_targets_only(self, monkeypatch):
        exposures = make_exposures()
        for exposure in exposures:
            exposure._targets = fake_targets(exposure)

        def boom(*args, **kwargs):
            raise AssertionError("database must not be touched when xmatch=False")

        monkeypatch.setattr(apogee, "_execute_query", boom)
        apogee.cross_match_targets(exposures, SEQUENCES["objects"], xmatch=False)
        assert all(t.sdss_id == -1 for t in exposures[0].targets)
