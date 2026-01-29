import os
import yaml
from typing import List, Dict
from dataclasses import dataclass, field, is_dataclass, asdict
from pathlib import Path


@dataclass
class DatabaseConfig:
    user: str = "sdss_user"
    host: str = "operations.sdss.org"
    port: int = 5432
    domain: str = "operations.sdss.*"


@dataclass
class ObservatoryMJD:
    apo: int = 59_558
    lco: int = 59_558


@dataclass
class Config:
    sdssdb: DatabaseConfig = field(default_factory=DatabaseConfig)
    database_connect_time_warning: int = 3  # seconds

    sdssdb_exposure_min_mjd: ObservatoryMJD = field(default_factory=ObservatoryMJD)
    logging_level: int = 20  # logging.INFO

    # Paths
    platelist_dir: str = "/uufs/chpc.utah.edu/common/home/sdss09/software/svn.sdss.org/data/sdss/platelist/trunk/plates/"
    sdsscore_dir: str = "/uufs/chpc.utah.edu/common/home/sdss50/software/git/sdss/sdsscore/main/"
    apogee_dir: str = "/uufs/chpc.utah.edu/common/home/sdss/sdsswork/data/apogee/"
    mapper_dir: str = "/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/data/mapper/"

    display_field_names: List[str] = field(
        default_factory=lambda: [
            "exposure",
            "image_type",
            "n_read",
            "field_id",
            "plate_id",
            "config_id",
            "design_id",
            "dithered_pixels",
            "lamp_quartz",
            "lamp_thar",
            "lamp_une",
            "name",
            "observer_comment"
        ]
    )

def get_config_path():
    config_dir = Path.home() / ".almanac"
    config_dir.mkdir(exist_ok=True)
    return config_dir / "config.yaml"


class ConfigManager:
    """A utility class to save and load dataclass configurations using YAML."""

    @staticmethod
    def save(config: object, file_path: str):
        """Saves a dataclass object to a YAML file."""
        if not is_dataclass(config):
            raise TypeError("Provided object is not a dataclass.")

        data = asdict(config)
        with open(file_path, "w") as f:
            yaml.dump(data, f, sort_keys=False)

    @staticmethod
    def load(cls, file_path: str):
        """Loads a dataclass object from a YAML file."""
        if not is_dataclass(cls):
            raise TypeError("Provided class is not a dataclass.")

        with open(file_path, "r") as f:
            data = yaml.safe_load(f)

        # Recursively create nested dataclasses from the dictionary
        def _load_recursive(cls, data):
            if not is_dataclass(cls):
                return data

            fields = {f.name: f.type for f in cls.__dataclass_fields__.values()}
            kwargs = {}
            for key, value in (data or {}).items():
                # Skip keys that don't exist in the dataclass (removed fields)
                if key not in fields:
                    continue
                field_type = fields.get(key)
                if is_dataclass(field_type):
                    kwargs[key] = _load_recursive(field_type, value)
                else:
                    kwargs[key] = value
            # Missing keys will use their default values from the dataclass
            return cls(**kwargs)

        if data:
            config = _load_recursive(cls, data)
        else:
            config = cls()
        return config


config_path = get_config_path()
if not os.path.exists(config_path):
    config = Config()
    ConfigManager.save(config, config_path)
else:
    config = ConfigManager.load(Config, config_path)
