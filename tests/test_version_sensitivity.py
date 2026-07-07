import tempfile
import unittest
from pathlib import Path

from PhenGO.predict.version_sensitivity import (
    discover_arff_files,
    expand_models,
    natural_key,
)


class VersionSensitivityHelperTests(unittest.TestCase):
    def test_expand_models_deduplicates_all(self):
        self.assertEqual(expand_models(["lr", "all", "lr"]), ["lr", "dt", "rf", "gb", "svm"])

    def test_natural_key_sorts_year_names(self):
        names = ["worm_2025", "worm_2017", "worm_2021"]
        self.assertEqual(sorted(names, key=natural_key), ["worm_2017", "worm_2021", "worm_2025"])

    def test_discover_arff_files_from_output_subdirs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "worm_2021").mkdir()
            (root / "worm_2017").mkdir()
            (root / "worm_2021" / "worm_PhenGO.arff").write_text("", encoding="utf-8")
            (root / "worm_2017" / "worm_PhenGO.arff").write_text("", encoding="utf-8")

            paths, names = discover_arff_files(str(root))

        self.assertEqual(names, ["worm_2017", "worm_2021"])
        self.assertTrue(paths[0].endswith("worm_2017/worm_PhenGO.arff"))


if __name__ == "__main__":
    unittest.main()
