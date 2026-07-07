import unittest

from PhenGO.core.PhenGO import removed_unused_gos


class GOFilteringTests(unittest.TestCase):
    def test_unused_filter_preserves_expanded_ancestor_terms(self):
        genes = {
            "gene_a": {
                "status": "lethal",
                "go_list": ["GO:child"],
                "expanded_go_list": ["GO:parent", "GO:child"],
                "bin_vec": [1, 1, 0],
            },
            "gene_b": {
                "status": "viable",
                "go_list": ["GO:parent"],
                "expanded_go_list": ["GO:parent"],
                "bin_vec": [1, 0, 0],
            },
        }

        filtered, terms = removed_unused_gos(genes, ["GO:parent", "GO:child", "GO:unused"])

        self.assertEqual(terms, ["GO:parent", "GO:child"])
        self.assertEqual(filtered["gene_a"]["expanded_go_list"], ["GO:parent", "GO:child"])
        self.assertEqual(filtered["gene_a"]["bin_vec"], [1, 1])
        self.assertEqual(filtered["gene_a"]["go_list"], ["GO:child"])
        self.assertEqual(filtered["gene_b"]["bin_vec"], [1, 0])


if __name__ == "__main__":
    unittest.main()
