"""Join phenotype-labelled genes to GAF annotations using stable identifiers."""
from __future__ import annotations

import csv
import gzip
import logging
from collections import Counter, defaultdict

logger = logging.getLogger(__name__)


def _open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def _evidence_allowed(row, options=None):
    if options is None:
        return True
    evidence = row[6].strip() if len(row) > 6 else ""
    include = set(getattr(options, "go_evidence_codes", None) or [])
    exclude = set(getattr(options, "exclude_go_evidence_codes", None) or [])
    return (not include or evidence in include) and evidence not in exclude


def _collect_gaf(path, labelled_genes, database, candidate_columns, options=None,
                 status_override=None):
    """Return records keyed by stable GAF DB Object ID (column 2).

    The complete filtered GAF is indexed before labels are assigned. Candidate
    columns are ordered by identifier authority: stable accession, canonical
    symbol, then historical synonym. Non-stable identifiers are accepted only
    when they resolve to exactly one stable accession. Strict snapshots disable
    synonym fallback entirely.
    """
    strict = bool(options is not None and getattr(options, "strict_snapshot", False))
    allowed_columns = tuple(
        index for index in candidate_columns if not (strict and index == 10)
    )
    identifier_names = {1: "stable_id", 2: "canonical_symbol", 10: "synonym"}
    counters = {
        "processed": 0,
        "not": 0,
        "evidence": 0,
        "assigned": 0,
        "malformed": 0,
        "stable_id": 0,
        "canonical_symbol": 0,
        "synonym": 0,
        "unmatched_identifiers": 0,
        "ambiguous_identifiers": 0,
        "contradictory_stable_ids": 0,
        "lower_priority_conflict_ids": 0,
    }
    genes = {}
    indices = {
        index: defaultdict(set) for index in allowed_columns
    }

    with _open_text(path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 5 or row[0] != database:
                continue
            counters["processed"] += 1
            qualifiers = {value for value in row[3].split("|") if value}
            if "NOT" in qualifiers:
                counters["not"] += 1
                continue
            if not _evidence_allowed(row, options):
                counters["evidence"] += 1
                continue

            stable_id = row[1].strip()
            if not stable_id:
                counters["malformed"] += 1
                continue
            gene = genes.setdefault(stable_id, {
                "symbols": Counter(),
                "go_list": [],
            })
            symbol = row[2].strip() if len(row) > 2 else ""
            if symbol:
                gene["symbols"][symbol] += 1
            go_id = row[4].strip()
            if go_id and go_id not in gene["go_list"]:
                gene["go_list"].append(go_id)

            for index in allowed_columns:
                if index >= len(row):
                    continue
                for candidate in row[index].split("|"):
                    candidate = candidate.strip()
                    if candidate:
                        indices[index][candidate].add(stable_id)

    audit_rows = []
    matches_by_stable = defaultdict(list)
    for identifier, raw_status in sorted(labelled_genes.items()):
        status = str(raw_status)
        selected = None
        for priority, index in enumerate(allowed_columns):
            stable_ids = indices[index].get(identifier, set())
            if len(stable_ids) == 1:
                selected = (priority, index, next(iter(stable_ids)))
                break
            if len(stable_ids) > 1:
                counters["ambiguous_identifiers"] += 1
                audit_rows.append({
                    "database": database,
                    "input_identifier": identifier,
                    "input_status": status,
                    "outcome": "excluded_ambiguous_identifier",
                    "identifier_type": identifier_names.get(
                        index, f"column_{index + 1}"
                    ),
                    "stable_gene_ids": "|".join(sorted(stable_ids)),
                    "details": "identifier maps to multiple stable GAF accessions",
                })
                selected = False
                break
        if selected is False:
            continue
        if selected is None:
            counters["unmatched_identifiers"] += 1
            audit_rows.append({
                "database": database,
                "input_identifier": identifier,
                "input_status": status,
                "outcome": "excluded_unmatched_identifier",
                "identifier_type": "",
                "stable_gene_ids": "",
                "details": "no permitted exact identifier match in filtered GAF",
            })
            continue
        priority, index, stable_id = selected
        identifier_kind = identifier_names.get(index, f"column_{index + 1}")
        counters[identifier_kind] = counters.get(identifier_kind, 0) + 1
        matches_by_stable[stable_id].append({
            "priority": priority,
            "identifier_type": identifier_kind,
            "identifier": identifier,
            "status": status,
        })

    records = {}
    for stable_id, matches in sorted(matches_by_stable.items()):
        best_priority = min(match["priority"] for match in matches)
        selected = [
            match for match in matches if match["priority"] == best_priority
        ]
        selected_statuses = {match["status"] for match in selected}
        if len(selected_statuses) != 1:
            counters["contradictory_stable_ids"] += 1
            audit_rows.append({
                "database": database,
                "input_identifier": "|".join(sorted(
                    match["identifier"] for match in selected
                )),
                "input_status": "|".join(sorted(selected_statuses)),
                "outcome": "excluded_contradictory_labels",
                "identifier_type": selected[0]["identifier_type"],
                "stable_gene_ids": stable_id,
                "details": "contradictory labels at the strongest identifier priority",
            })
            continue

        status = next(iter(selected_statuses))
        lower_conflicts = [
            match for match in matches
            if match["priority"] > best_priority and match["status"] != status
        ]
        if lower_conflicts:
            counters["lower_priority_conflict_ids"] += 1
            audit_rows.append({
                "database": database,
                "input_identifier": "|".join(sorted(
                    match["identifier"] for match in lower_conflicts
                )),
                "input_status": "|".join(sorted({
                    match["status"] for match in lower_conflicts
                })),
                "outcome": "ignored_lower_priority_conflict",
                "identifier_type": "|".join(sorted({
                    match["identifier_type"] for match in lower_conflicts
                })),
                "stable_gene_ids": stable_id,
                "details": "stronger exact identifier label retained",
            })

        gene = genes[stable_id]
        symbol = gene["symbols"].most_common(1)[0][0] if gene["symbols"] else ""
        if status_override is not None:
            status = str(status_override(stable_id, symbol, status))
        record = {
            "status": status,
            "go_list": gene["go_list"],
        }
        if options is not None:
            record["gene_symbol"] = symbol
        records[stable_id] = record

    counters["assigned"] = sum(len(record["go_list"]) for record in records.values())
    logger.info(
        "Processed %d %s rows; filtered %d NOT and %d evidence-code rows; assigned %d annotations",
        counters["processed"], database, counters["not"], counters["evidence"], counters["assigned"],
    )
    logger.info(
        "Label identifiers matched: stable ID=%d, canonical symbol=%d, synonym fallback=%d",
        counters["stable_id"], counters["canonical_symbol"], counters["synonym"],
    )
    if strict and 10 in candidate_columns:
        logger.info("Strict identifier matching: GAF synonym fallback disabled")
    if counters["lower_priority_conflict_ids"]:
        logger.warning(
            "Resolved %d stable-ID joins by identifier precedence despite conflicting "
            "lower-priority identifiers",
            counters["lower_priority_conflict_ids"],
        )
    excluded = (
        counters["ambiguous_identifiers"] +
        counters["contradictory_stable_ids"]
    )
    if excluded:
        logger.warning(
            "Excluded %d unsafe identifier assignments (%d ambiguous identifiers; "
            "%d stable IDs with contradictory strongest labels)",
            excluded, counters["ambiguous_identifiers"],
            counters["contradictory_stable_ids"],
        )
    if counters["unmatched_identifiers"]:
        logger.info(
            "Excluded %d phenotype identifiers with no permitted exact GAF match",
            counters["unmatched_identifiers"],
        )
    if options is not None:
        join_runs = list(getattr(options, "_identifier_join_runs", []))
        join_runs.append({
            "database": database,
            **counters,
        })
        options._identifier_join_runs = join_runs
        options._identifier_join_audit_rows = list(getattr(
            options, "_identifier_join_audit_rows", []
        )) + audit_rows
    return records


def _log_summary(species, labelled_genes, records):
    logger.info("Species: %s", species)
    logger.info("Lethal genes with GO terms: %d", sum(v["status"] == "lethal" for v in records.values()))
    logger.info("Viable genes with GO terms: %d", sum(v["status"] == "viable" for v in records.values()))
    logger.info("Phenotype-labelled source identifiers: %d", len(labelled_genes))
    logger.info("Canonical stable gene IDs with GO joins: %d", len(records))


def get_viability_go_data_yeast(gene_association_file, vi_inviable_genes, options=None):
    # SGD phenotype names are usually systematic names present first in GAF synonyms.
    records = _collect_gaf(
        gene_association_file, vi_inviable_genes, "SGD", (1, 2, 10), options,
    )
    _log_summary("yeast", vi_inviable_genes, records)
    return records


def _read_first_column(path):
    values = set()
    if not path:
        return values
    with _open_text(path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if row and row[0] and not row[0].startswith("#"):
                value = next((cell.strip() for cell in row[:2] if cell.strip()), "")
                if value:
                    values.add(value)
    return values


def _load_fly_assignment_ids(path):
    """Return unique release-matched FlyBase symbol-to-accession assignments."""
    symbol_ids = defaultdict(set)
    valid_symbols = set()
    with _open_text(path) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) <= 2 or row[0].startswith("#"):
                continue
            if row[1] != "melanogaster" or row[2] == "Withdrawn":
                continue
            symbol = row[0].strip()
            if not symbol:
                continue
            valid_symbols.add(symbol)
            if len(row) > 3 and row[3].strip().startswith("FBgn"):
                symbol_ids[symbol].add(row[3].strip())
    unique = {
        symbol: next(iter(stable_ids))
        for symbol, stable_ids in symbol_ids.items()
        if len(stable_ids) == 1
    }
    ambiguous = {
        symbol: stable_ids
        for symbol, stable_ids in symbol_ids.items()
        if len(stable_ids) > 1
    }
    return unique, ambiguous, valid_symbols


def get_viability_go_data_fly(options, gene_association_file, vi_inviable_genes,
                              apply_legacy_lethal_override=True):
    lethal_ids = _read_first_column(getattr(options, "fly_lethal_genes", None)) \
        if apply_legacy_lethal_override else set()

    def override(stable_id, symbol, status):
        return "lethal" if stable_id in lethal_ids or symbol in lethal_ids else status

    join_labels = vi_inviable_genes
    candidate_columns = (1, 2, 10)
    valid_symbols = set()
    if getattr(options, "fly_assignments", None):
        assignments, ambiguous_symbols, valid_symbols = _load_fly_assignment_ids(
            options.fly_assignments
        )
        if ambiguous_symbols:
            logger.warning(
                "Fly assignment file contains %d symbols mapped to multiple FBgn IDs; "
                "those symbols are excluded",
                len(ambiguous_symbols),
            )
        if assignments:
            translated = defaultdict(set)
            translation_audit = []
            for identifier, status in sorted(vi_inviable_genes.items()):
                if identifier in ambiguous_symbols:
                    translation_audit.append({
                        "database": "FB",
                        "input_identifier": identifier,
                        "input_status": str(status),
                        "outcome": "excluded_ambiguous_assignment",
                        "identifier_type": "fly_assignment_symbol",
                        "stable_gene_ids": "|".join(sorted(ambiguous_symbols[identifier])),
                        "details": "release assignment maps symbol to multiple FBgn IDs",
                    })
                    continue
                stable_id = identifier if identifier.startswith("FBgn") else assignments.get(identifier)
                if not stable_id:
                    translation_audit.append({
                        "database": "FB",
                        "input_identifier": identifier,
                        "input_status": str(status),
                        "outcome": "excluded_unmatched_assignment",
                        "identifier_type": "fly_assignment_symbol",
                        "stable_gene_ids": "",
                        "details": "symbol absent from release-matched D. melanogaster assignments",
                    })
                    continue
                translated[stable_id].add(str(status))
            join_labels = {}
            for stable_id, statuses in sorted(translated.items()):
                if len(statuses) == 1:
                    join_labels[stable_id] = next(iter(statuses))
                else:
                    translation_audit.append({
                        "database": "FB",
                        "input_identifier": stable_id,
                        "input_status": "|".join(sorted(statuses)),
                        "outcome": "excluded_contradictory_assignment_labels",
                        "identifier_type": "stable_id",
                        "stable_gene_ids": stable_id,
                        "details": "multiple release symbols assign contradictory labels to one FBgn ID",
                    })
            options._identifier_join_audit_rows = list(getattr(
                options, "_identifier_join_audit_rows", []
            )) + translation_audit
            candidate_columns = (1,)
            logger.info(
                "Translated %d fly phenotype symbols to %d release-matched FBgn IDs",
                len(vi_inviable_genes), len(join_labels),
            )
        elif getattr(options, "strict_snapshot", False):
            raise ValueError(
                "Strict fly snapshots require an assignment file with release-matched "
                "FBgn IDs in column 4"
            )

    records = _collect_gaf(
        gene_association_file, join_labels, "FB", candidate_columns, options, override,
    )

    if valid_symbols and candidate_columns != (1,):
        before = len(records)
        records = {
            gene: values for gene, values in records.items()
            if values.get("gene_symbol") in valid_symbols
        }
        logger.info("Fly species/status filter removed %d genes", before - len(records))
    _log_summary("fly", vi_inviable_genes, records)
    return records


def get_viability_go_data_fish(gene_association_file, vi_inviable_genes, options=None):
    records = _collect_gaf(
        gene_association_file, vi_inviable_genes, "ZFIN", (1, 2, 10), options,
    )
    _log_summary("fish", vi_inviable_genes, records)
    return records


def get_viability_go_data_worm(gene_association_file, vi_inviable_genes, options=None):
    records = _collect_gaf(
        gene_association_file, vi_inviable_genes, "WB", (1, 2, 10), options,
    )
    _log_summary("worm", vi_inviable_genes, records)
    return records


def get_viability_go_data_mouse(gene_association_file, vi_inviable_genes, options=None):
    records = _collect_gaf(
        gene_association_file, vi_inviable_genes, "MGI", (1, 2, 10), options,
    )
    _log_summary("mouse", vi_inviable_genes, records)
    return records
