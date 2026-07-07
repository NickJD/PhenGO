import tempfile
import unittest
from pathlib import Path

from PhenGO.core.PhenGO import write_arff_output
from PhenGO.scripts.compare_arff_genes import parse_arff_with_terms


class ARFFIOTests(unittest.TestCase):
    def test_writer_and_parser_preserve_comma_in_gene_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "quoted.arff"
            write_arff_output(
                {
                    "gene,with,comma": {
                        "status": "lethal",
                        "bin_vec": [1],
                    }
                },
                ["GO:0000001"],
                str(output),
            )

            genes, terms = parse_arff_with_terms(str(output))

        self.assertEqual(terms, ["GO:0000001"])
        self.assertIn("gene,with,comma", genes)
        self.assertEqual(genes["gene,with,comma"]["label"], "lethal")
        self.assertEqual(genes["gene,with,comma"]["features"]["GO:0000001"], "1")


if __name__ == "__main__":
    unittest.main()
