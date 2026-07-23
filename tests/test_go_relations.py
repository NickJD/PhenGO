from pathlib import Path

import networkx as nx
import pytest

from PhenGO.core.PhenGO import assign_go_to_vector, validate_go_terms_against_graph
from PhenGO.core.obo_to_graph import obo_to_graph


OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: root
namespace: biological_process

[Term]
id: GO:0000002
name: is-a child
namespace: biological_process
alt_id: GO:0099992
is_a: GO:0000001 ! root

[Term]
id: GO:0000003
name: part-of child
namespace: biological_process
relationship: part_of GO:0000002 ! is-a child

[Term]
id: GO:0000004
name: obsolete term
namespace: biological_process
alt_id: GO:0099994
is_obsolete: true
replaced_by: GO:0000003
"""


def test_relation_selection_alt_ids_and_obsolete_replacement(tmp_path: Path):
    obo = tmp_path / "go.obo"
    obo.write_text(OBO, encoding="utf-8")

    graph, terms, obsolete = obo_to_graph(
        str(tmp_path), str(obo), relations=["is_a", "part_of"],
    )

    assert terms == ["GO:0000001", "GO:0000002", "GO:0000003"]
    assert set(nx.ancestors(graph, "GO:0000003")) == {"GO:0000001", "GO:0000002"}
    assert graph.graph["alt_id_to_id"]["GO:0099994"] == "GO:0000004"
    assert graph.graph["obsolete_replacements"]["GO:0000004"] == "GO:0000003"

    genes = {
        "gene-a": {"status": "lethal", "go_list": ["GO:0099992"]},
        "gene-b": {"status": "viable", "go_list": ["GO:0099994"]},
    }
    cleaned, report = validate_go_terms_against_graph(
        genes, graph, obsolete, str(tmp_path),
    )
    assert cleaned["gene-a"]["go_list"] == ["GO:0000002"]
    assert cleaned["gene-b"]["go_list"] == ["GO:0000003"]
    assert report["alternate_ids_mapped"] == 2
    assert report["obsolete_terms_replaced"] == 1

    expanded, _, _ = assign_go_to_vector(cleaned, graph, terms)
    assert expanded["gene-b"]["expanded_go_list"] == terms


def test_default_relations_traverse_is_a_and_part_of(tmp_path: Path):
    obo = tmp_path / "go.obo"
    obo.write_text(OBO, encoding="utf-8")
    graph, terms, _ = obo_to_graph(str(tmp_path), str(obo))

    assert nx.ancestors(graph, "GO:0000003") == {"GO:0000001", "GO:0000002"}
    genes = {"gene": {"status": "lethal", "go_list": ["GO:0000003"]}}
    expanded, _, _ = assign_go_to_vector(genes, graph, terms)
    assert expanded["gene"]["expanded_go_list"] == terms


def test_cross_namespace_propagation_requires_explicit_override(tmp_path: Path):
    obo = tmp_path / "go.obo"
    obo.write_text(
        """format-version: 1.2

[Term]
id: GO:0000001
name: process
namespace: biological_process

[Term]
id: GO:0000002
name: activity
namespace: molecular_function
relationship: part_of GO:0000001 ! process
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="cross-namespace"):
        obo_to_graph(str(tmp_path), str(obo))

    graph, _, _ = obo_to_graph(
        str(tmp_path), str(obo), allow_cross_namespace=True,
    )
    assert graph.graph["cross_namespace_edge_count"] == 1
