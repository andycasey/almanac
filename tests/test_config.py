"""Tests for config backwards compatibility."""

import tempfile
import os
from dataclasses import dataclass, field

import pytest
import yaml

from almanac.config import ConfigManager


@dataclass
class SimpleConfig:
    """A simple test config."""
    name: str = "default"
    value: int = 42


@dataclass
class NestedInner:
    """Nested config inner class."""
    host: str = "localhost"
    port: int = 8080


@dataclass 
class NestedConfig:
    """Config with nested dataclass."""
    inner: NestedInner = field(default_factory=NestedInner)
    enabled: bool = True


class TestConfigBackwardsCompatibility:
    """Test that config loading handles new and removed keys gracefully."""

    def test_load_with_missing_key_uses_default(self):
        """New config keys should use their default values when not in YAML."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            # YAML file is missing the 'value' key
            yaml.dump({"name": "custom"}, f)
            f.flush()
            
            try:
                config = ConfigManager.load(SimpleConfig, f.name)
                assert config.name == "custom"
                assert config.value == 42  # Should use default
            finally:
                os.unlink(f.name)

    def test_load_with_extra_key_ignores_it(self):
        """Removed config keys in YAML should be ignored."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            # YAML file has an extra key that doesn't exist in dataclass
            yaml.dump({
                "name": "custom",
                "value": 100,
                "removed_field": "should be ignored"
            }, f)
            f.flush()
            
            try:
                config = ConfigManager.load(SimpleConfig, f.name)
                assert config.name == "custom"
                assert config.value == 100
                assert not hasattr(config, "removed_field")
            finally:
                os.unlink(f.name)

    def test_load_with_empty_yaml(self):
        """Empty YAML file should return config with all defaults."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("")  # Empty file
            f.flush()
            
            try:
                config = ConfigManager.load(SimpleConfig, f.name)
                assert config.name == "default"
                assert config.value == 42
            finally:
                os.unlink(f.name)

    def test_load_with_none_yaml(self):
        """YAML file with just 'null' should return config with all defaults."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("null")
            f.flush()
            
            try:
                config = ConfigManager.load(SimpleConfig, f.name)
                assert config.name == "default"
                assert config.value == 42
            finally:
                os.unlink(f.name)

    def test_load_nested_with_missing_inner_key(self):
        """Nested configs should also handle missing keys."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            # Inner config is missing 'port' key
            yaml.dump({
                "inner": {"host": "example.com"},
                "enabled": False
            }, f)
            f.flush()
            
            try:
                config = ConfigManager.load(NestedConfig, f.name)
                assert config.inner.host == "example.com"
                assert config.inner.port == 8080  # Should use default
                assert config.enabled == False
            finally:
                os.unlink(f.name)

    def test_load_nested_with_extra_inner_key(self):
        """Nested configs should ignore extra keys."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump({
                "inner": {
                    "host": "example.com",
                    "port": 9000,
                    "old_setting": "ignored"
                },
                "enabled": True
            }, f)
            f.flush()
            
            try:
                config = ConfigManager.load(NestedConfig, f.name)
                assert config.inner.host == "example.com"
                assert config.inner.port == 9000
                assert not hasattr(config.inner, "old_setting")
            finally:
                os.unlink(f.name)
