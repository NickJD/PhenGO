"""Parse GO OBO snapshots into a relation-aware canonical graph."""
from __future__ import annotations

import gzip
import json
import logging
import os
from collections import defaultdict

import networkx as nx

logger = logging.getLogger(__name__)

SUPPORTED_GO_RELATIONS = (
    "is_a",
    "part_of",
    "regulates",
    "positively_regulates",
    "negatively_regulates",
    "occurs_in",
    "capable_of",
    "capable_of_part_of",
)

GO_ROOT_TERMS = {"GO:0008150", "GO:0003674", "GO:0005575"}


def define_graph_from_file(edges_input):
    with open(edges_input, encoding="utf-8") as handle:
        data = json.load(handle)
    graph = nx.DiGraph()
    graph.add_nodes_from(data)
    for child, values in data.items():
        for parent in values.get("p", []):
            if parent in graph:
                graph.add_edge(parent, child, relations=["selected"])
    return graph


def getDescendents(goid, terms):
    found = {goid}
    for child in terms.get(goid, {}).get("c", []):
        found.update(getDescendents(child, terms))
    return found


def getAncestors(goid, terms):
    found = {goid}
    for parent in terms.get(goid, {}).get("p", []):
        found.update(getAncestors(parent, terms))
    return found


def _open_obo(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    return opener(path, "rt", encoding="utf-8", errors="replace")


def _parse_obo_header(path):
    header = {}
    with _open_obo(path) as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("["):
                break
            if not line or line.startswith("!") or ": " not in line:
                continue
            tag, value = line.split(": ", 1)
            if tag in {"format-version", "data-version", "ontology"}:
                header[tag.replace("-", "_")] = value
    return header


def _iter_term_stanzas(path):
    current = None
    with _open_obo(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                if current:
                    yield current
                current = defaultdict(list)
                continue
            if line.startswith("["):
                if current:
                    yield current
                current = None
                continue
            if current is None or not line or line.startswith("!") or ": " not in line:
                continue
            tag, value = line.split(": ", 1)
            current[tag].append(value)
    if current:
        yield current


def _parse_terms(path):
    records = {}
    for stanza in _iter_term_stanzas(path):
        if "id" not in stanza:
            continue
        term_id = stanza["id"][0]
        relations = defaultdict(list)
        for value in stanza.get("is_a", []):
            relations["is_a"].append(value.split()[0])
        for value in stanza.get("relationship", []):
            parts = value.split()
            if len(parts) >= 2:
                relations[parts[0]].append(parts[1])
        records[term_id] = {
            "id": term_id,
            "name": stanza.get("name", [""])[0],
            "namespace": stanza.get("namespace", [""])[0],
            "alt_ids": list(stanza.get("alt_id", [])),
            "obsolete": stanza.get("is_obsolete", ["false"])[0] == "true",
            "replaced_by": [value.split()[0] for value in stanza.get("replaced_by", [])],
            "consider": [value.split()[0] for value in stanza.get("consider", [])],
            "relations": {key: list(dict.fromkeys(values))
                          for key, values in relations.items()},
        }
    return records


def obo_to_graph(output_dir, go_obo_file, relations=("is_a", "part_of"), namespaces=None,
                 allow_cross_namespace=False):
    """Create a canonical GO graph using explicitly selected relation types.

    Alternate IDs are mapped to their canonical term and never become duplicate
    model features. Obsolete terms are retained only as replacement metadata.
    The returned graph stores ``alt_id_to_id``, ``obsolete_replacements`` and
    ``selected_relations`` in ``graph.graph``.
    """
    selected_relations = tuple(dict.fromkeys(relations or ("is_a",)))
    unknown = set(selected_relations) - set(SUPPORTED_GO_RELATIONS)
    if unknown:
        raise ValueError(f"Unsupported GO relation(s): {sorted(unknown)}")
    namespace_set = set(namespaces or ())

    logger.info("Parsing GO OBO file with relation(s): %s", ", ".join(selected_relations))
    obo_header = _parse_obo_header(go_obo_file)
    records = _parse_terms(go_obo_file)
    obsolete = {term_id for term_id, rec in records.items() if rec["obsolete"]}
    active = {
        term_id: rec for term_id, rec in records.items()
        if not rec["obsolete"] and (not namespace_set or rec["namespace"] in namespace_set)
    }

    alt_id_to_id = {}
    for term_id, rec in records.items():
        for alt_id in rec["alt_ids"]:
            if alt_id in alt_id_to_id and alt_id_to_id[alt_id] != term_id:
                raise ValueError(
                    f"Alternate GO ID {alt_id} maps to multiple terms: "
                    f"{alt_id_to_id[alt_id]} and {term_id}"
                )
            alt_id_to_id[alt_id] = term_id

    def resolve_active_replacement(term_id, seen=None):
        seen = set(seen or ())
        canonical = alt_id_to_id.get(term_id, term_id)
        if canonical in active:
            return canonical
        if canonical in seen or canonical not in records:
            return None
        seen.add(canonical)
        candidates = {
            replacement
            for raw in records[canonical]["replaced_by"]
            if (replacement := resolve_active_replacement(raw, seen)) is not None
        }
        return next(iter(candidates)) if len(candidates) == 1 else None

    obsolete_replacements = {}
    for term_id in obsolete:
        replacement = resolve_active_replacement(term_id)
        if replacement:
            obsolete_replacements[term_id] = replacement

    graph = nx.DiGraph()
    for term_id, rec in active.items():
        graph.add_node(term_id, name=rec["name"], namespace=rec["namespace"])

    parents_by_term = defaultdict(lambda: defaultdict(list))
    children_by_term = defaultdict(lambda: defaultdict(list))
    relation_edge_counts = defaultdict(int)
    cross_namespace_edges = []
    for child, rec in active.items():
        for relation in selected_relations:
            for raw_parent in rec["relations"].get(relation, []):
                parent = alt_id_to_id.get(raw_parent, raw_parent)
                parent = obsolete_replacements.get(parent, parent)
                if parent not in active or parent == child:
                    continue
                if graph.has_edge(parent, child):
                    edge_relations = graph[parent][child].setdefault("relations", [])
                    if relation not in edge_relations:
                        edge_relations.append(relation)
                else:
                    graph.add_edge(parent, child, relations=[relation])
                if active[parent]["namespace"] != rec["namespace"]:
                    cross_namespace_edges.append({
                        "parent": parent,
                        "parent_namespace": active[parent]["namespace"],
                        "child": child,
                        "child_namespace": rec["namespace"],
                        "relation": relation,
                    })
                parents_by_term[child][relation].append(parent)
                children_by_term[parent][relation].append(child)
                relation_edge_counts[relation] += 1

    if cross_namespace_edges and not allow_cross_namespace:
        example = cross_namespace_edges[0]
        raise ValueError(
            "Selected GO propagation graph contains cross-namespace edges; use a "
            "release-matched go-basic.obo or explicitly pass "
            "-allow_cross_namespace_go_edges. Example: "
            f"{example['child']} ({example['child_namespace']}) "
            f"-{example['relation']}-> {example['parent']} "
            f"({example['parent_namespace']})"
        )
    if not nx.is_directed_acyclic_graph(graph):
        cycle = nx.find_cycle(graph)
        raise ValueError(
            "Selected GO propagation relations form a cycle and cannot define "
            f"an ancestor feature hierarchy. Example cycle: {cycle[:5]}"
        )
    if cross_namespace_edges:
        logger.warning(
            "Retaining %d explicitly allowed cross-namespace GO edges",
            len(cross_namespace_edges),
        )

    graph.graph.update({
        "alt_id_to_id": alt_id_to_id,
        "obsolete_replacements": obsolete_replacements,
        "selected_relations": list(selected_relations),
        "namespaces": sorted(namespace_set),
        "obo_header": obo_header,
        "cross_namespace_edge_count": len(cross_namespace_edges),
    })

    serialised = {}
    for term_id, rec in active.items():
        relation_parents = {
            rel: sorted(set(parents_by_term[term_id].get(rel, [])))
            for rel in selected_relations
            if parents_by_term[term_id].get(rel)
        }
        relation_children = {
            rel: sorted(set(children_by_term[term_id].get(rel, [])))
            for rel in selected_relations
            if children_by_term[term_id].get(rel)
        }
        serialised[term_id] = {
            "name": rec["name"],
            "namespace": rec["namespace"],
            "p": sorted({p for values in relation_parents.values() for p in values}),
            "c": sorted({c for values in relation_children.values() for c in values}),
            "parents_by_relation": relation_parents,
            "children_by_relation": relation_children,
        }

    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, "GO_Children&Parents.json"), "w", encoding="utf-8") as handle:
        json.dump(serialised, handle, indent=2, sort_keys=True)
    unique_nodes = sorted(active)
    with open(os.path.join(output_dir, "Unique_GO_Nodes.txt"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(unique_nodes) + ("\n" if unique_nodes else ""))

    logger.info(
        "Parsed %d canonical active GO terms, %d alternate IDs and %d obsolete terms",
        len(active), len(alt_id_to_id), len(obsolete),
    )
    logger.info("Selected relation edges: %s", dict(relation_edge_counts))
    return graph, unique_nodes, obsolete
