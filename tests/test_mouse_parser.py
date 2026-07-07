import csv
import gzip
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from PhenGO.core.go_handling import get_viability_go_data_mouse
from PhenGO.core.phenotype_handling import get_viable_inviable_mouse


def write_gzip_tsv(path, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerows(rows)


class MouseParserTests(unittest.TestCase):
    def test_mgi_genepheno_rows_use_marker_accession(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            lethal_terms = tmp / "mouse_terms.tsv.gz"
            phenotype = tmp / "MGI_GenePheno.rpt.gz"
            write_gzip_tsv(lethal_terms, [["ID", "Def"], ["MP:0001", "lethal"]])
            write_gzip_tsv(
                phenotype,
                [["allelic", "allele", "allele_id", "bg", "MP:0001", "pmid", "MGI:123", "GeneA", "geno"]],
            )

            result = get_viable_inviable_mouse(
                SimpleNamespace(mouse_phenotypes=str(lethal_terms), filter_mixed_terms=False),
                str(phenotype),
            )

        self.assertEqual(result, {"MGI:123": "lethal"})

    def test_mgi_phenogenomp_rows_use_marker_accession(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            lethal_terms = tmp / "mouse_terms.tsv.gz"
            phenotype = tmp / "MGI_PhenoGenoMP.rpt.gz"
            write_gzip_tsv(lethal_terms, [["ID", "Def"], ["MP:0001", "lethal"]])
            write_gzip_tsv(
                phenotype,
                [["allelic", "allele", "bg", "MP:0002", "pmid", "MGI:456", "geno"]],
            )

            result = get_viable_inviable_mouse(
                SimpleNamespace(mouse_phenotypes=str(lethal_terms), filter_mixed_terms=False),
                str(phenotype),
            )

        self.assertEqual(result, {"MGI:456": "viable"})

    def test_mouse_go_join_prefers_marker_accession(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gaf = Path(tmpdir) / "mgi.gaf"
            with open(gaf, "w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["MGI", "MGI:123", "GeneA", "", "GO:0000001"])

            result = get_viability_go_data_mouse(str(gaf), {"MGI:123": "lethal"})

        self.assertEqual(
            result,
            {"MGI:123": {"status": "lethal", "go_list": ["GO:0000001"]}},
        )


if __name__ == "__main__":
    unittest.main()
