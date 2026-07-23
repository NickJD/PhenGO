import math

from PhenGO.scripts.GO_Compare import calculate_go_term_stats
from PhenGO.scripts.GO_temporal_analysis import calculate_enrichment_stats


def _genes():
    return {
        "l1": {"label": "lethal", "features": {"GO:1": 1, "GO:2": 0}},
        "l2": {"label": "lethal", "features": {"GO:1": 1, "GO:2": 0}},
        "v1": {"label": "viable", "features": {"GO:1": 0, "GO:2": 1}},
        "v2": {"label": "viable", "features": {"GO:1": 0, "GO:2": 1}},
    }


def test_go_statistics_use_finite_effects_and_fdr():
    pairwise, _ = calculate_go_term_stats(_genes(), ["GO:1", "GO:2"])
    temporal, _ = calculate_enrichment_stats(_genes(), ["GO:1", "GO:2"])

    for collection in (pairwise, temporal):
        for values in collection.values():
            assert math.isfinite(values["enrichment_ratio"])
            assert 0 <= values["p_value"] <= 1
            assert 0 <= values["fdr"] <= 1
