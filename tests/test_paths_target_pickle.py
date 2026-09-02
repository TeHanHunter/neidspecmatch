import configparser
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from neidspecmatch.paths import safe_filename_component
from neidspecmatch.target import Target
from neidspecmatch import utils


TARGET_DATA = {
    "ra": 1.0, "dec": 2.0, "pmra": 3.0, "pmdec": 4.0,
    "px": 5.0, "epoch": 2451545.0, "rv": 0.0,
}


class PathTargetAndPickleTests(unittest.TestCase):
    def test_safe_component_is_traversal_safe_and_collision_resistant(self):
        first = safe_filename_component("../A/B")
        second = safe_filename_component("../A B")
        self.assertNotIn("/", first)
        self.assertNotIn("..", first)
        self.assertNotEqual(first, second)
        self.assertEqual(first, safe_filename_component("../A/B"))

    def test_target_cache_name_is_safe_private_and_atomic(self):
        with tempfile.TemporaryDirectory() as temp, mock.patch.object(
                Target, "query_tic", return_value=TARGET_DATA):
            target = Target("TIC ../../escape", config_folder=temp,
                            allow_network=True, verbose=False)
            cache = Path(target.config_filename)
            self.assertEqual(cache.parent, Path(temp))
            self.assertTrue(cache.is_file())
            self.assertEqual(cache.stat().st_mode & 0o777, 0o600)
            self.assertFalse((Path(temp).parent / "escape.config").exists())

    def test_malformed_cache_does_not_trigger_network_fallback(self):
        with tempfile.TemporaryDirectory() as temp:
            cache = Path(temp) / f"{safe_filename_component('TIC 1')}.config"
            cache.write_text("[wrong]\nra=1\n", encoding="utf-8")
            with mock.patch.object(Target, "query_tic") as query:
                with self.assertRaises(configparser.NoSectionError):
                    Target("TIC 1", config_folder=temp, allow_network=True)
            query.assert_not_called()

    def test_network_is_disabled_by_default(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(FileNotFoundError):
                Target("TIC 1", config_folder=temp)

    def test_pickle_loading_requires_explicit_trust(self):
        with tempfile.TemporaryDirectory() as temp:
            filename = Path(temp) / "legacy.pkl"
            with self.assertWarns(DeprecationWarning):
                utils.pickle_dump(filename, {"value": 3})
            with self.assertRaises(ValueError):
                utils.pickle_load(filename)
            with self.assertWarns(RuntimeWarning):
                result = utils.pickle_load(filename, trusted=True)
            self.assertEqual(result, {"value": 3})


if __name__ == "__main__":
    unittest.main()
