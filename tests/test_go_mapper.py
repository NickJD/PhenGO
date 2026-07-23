import gzip
from pathlib import Path

from PhenGO.scripts.go_mapper import map_go_ids, parse_obo


def test_go_mapper_canonicalises_only_unambiguous_ids(tmp_path: Path):
    obo = tmp_path / "go.obo.gz"
    with gzip.open(obo, "wt", encoding="utf-8") as handle:
        handle.write(
            "format-version: 1.2\n\n"
            "[Term]\n"
            "id: GO:0001\n"
            "name: active\n"
            "namespace: biological_process\n"
            "alt_id: GO:ALT1\n\n"
            "[Term]\n"
            "id: GO:0002\n"
            "name: old\n"
            "is_obsolete: true\n"
            "replaced_by: GO:0001\n\n"
            "[Term]\n"
            "id: GO:0003\n"
            "name: ambiguous old\n"
            "is_obsolete: true\n"
            "replaced_by: GO:0001\n"
            "replaced_by: GO:0004\n\n"
            "[Term]\n"
            "id: GO:0004\n"
            "name: other active\n"
            "namespace: biological_process\n"
        )

    terms = parse_obo(obo)
    rows = map_go_ids(["GO:0001", "GO:ALT1", "GO:0002", "GO:0003", "GO:NOPE"], terms)

    assert [(row[0], row[1], row[2]) for row in rows] == [
        ("GO:0001", "GO:0001", "exact"),
        ("GO:ALT1", "GO:0001", "alt_id"),
        ("GO:0002", "GO:0001", "obsolete_replaced"),
        ("GO:0003", "", "obsolete_unresolved"),
        ("GO:NOPE", "", "missing"),
    ]
