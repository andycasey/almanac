"""Tests for the no-fibers warning in `almanac add metadata`."""

import h5py as h5
import pytest

from tests.helpers import write_night, run_cli


class TestAddMetadataWarning:
    def test_no_fibers_group_warns_and_does_not_write_meta(self, tmp_path, monkeypatch):
        path = str(tmp_path / "night.h5")
        write_night(path)

        import almanac.catalog as catalog

        monkeypatch.setattr(
            catalog,
            "query",
            lambda *a, **k: pytest.fail("database should not be queried"),
        )

        result = run_cli(["add", "metadata", path])
        assert result.exit_code == 1
        with h5.File(path, "r") as fp:
            assert "meta" not in fp

    def test_existing_meta_group_is_preserved_when_no_fibers(self, tmp_path):
        path = str(tmp_path / "night.h5")
        write_night(path)
        with h5.File(path, "a") as fp:
            fp.create_group("meta").create_dataset("sdss_id", data=[1, 2, 3])

        result = run_cli(["add", "metadata", path])
        assert result.exit_code == 1
        with h5.File(path, "r") as fp:
            assert list(fp["meta/sdss_id"][:]) == [1, 2, 3]
