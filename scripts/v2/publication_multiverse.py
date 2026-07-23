#!/usr/bin/env python3
"""Build a self-contained publication output from completed PhenGO V1 and V2 runs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ATTRIBUTE_RE = re.compile(
    r"^@attribute\s+(?:'([^']*)'|\"([^\"]*)\"|(\S+))\s+(.+)$", re.IGNORECASE
)
POSITIVE_LABELS = {"lethal", "inviable", "essential"}
NEGATIVE_LABELS = {"viable", "non-essential", "nonessential", "non_essential"}
PERFORMANCE_METRICS = (
    "average_precision", "balanced_accuracy", "mcc", "roc_auc", "brier_score"
)
COPY_RESERVE_BYTES = 512 * 1024 * 1024
CLONE_RESERVE_BYTES = 1024 * 1024 * 1024


def sha256_file(path: str | os.PathLike[str]) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: str | os.PathLike[str]) -> dict:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_tsv(path: Path, rows: list[dict], preferred: tuple[str, ...] = ()) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    fields = [key for key in preferred if key in fields] + [key for key in fields if key not in preferred]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        if not fields:
            handle.write("status\nno_rows\n")
            return
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_csv(path: Path, rows: list[dict], preferred: tuple[str, ...] = ()) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    fields = [key for key in preferred if key in fields] + [key for key in fields if key not in preferred]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        if not fields:
            handle.write("status\nno_rows\n")
            return
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def tree_file_inventory(root: Path) -> tuple[list[dict], dict]:
    rows = []
    aggregate = hashlib.sha256()
    total_bytes = 0
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = str(path.relative_to(root))
        digest = sha256_file(path)
        size = path.stat().st_size
        aggregate.update(relative.encode("utf-8"))
        aggregate.update(b"\0")
        aggregate.update(digest.encode("ascii"))
        aggregate.update(b"\n")
        total_bytes += size
        rows.append({"relative_path": relative, "size_bytes": size, "sha256": digest})
    return rows, {
        "file_count": len(rows),
        "total_bytes": total_bytes,
        "tree_sha256": aggregate.hexdigest(),
    }


def copy_verified_run(
    source: Path, destination: Path, source_run: str, copy_strategy: str = "copy2"
) -> tuple[list[dict], dict]:
    """Copy a completed run and verify every destination file byte-for-byte."""
    symlinks = [path for path in source.rglob("*") if path.is_symlink()]
    if symlinks:
        raise ValueError(
            f"Source run contains symbolic links and is not self-contained: {symlinks[:5]}"
        )
    if copy_strategy == "apfs_clone":
        completed = subprocess.run(
            ["cp", "-cR", str(source), str(destination)],
            check=False, capture_output=True, text=True,
        )
        if completed.returncode:
            raise OSError(
                f"APFS clone failed for {source}: {completed.stderr.strip()}"
            )
    elif copy_strategy == "copy2":
        shutil.copytree(source, destination, symlinks=False, copy_function=shutil.copy2)
    else:
        raise ValueError(f"Unknown copy strategy: {copy_strategy}")
    source_rows, source_summary = tree_file_inventory(source)
    copied_rows, copied_summary = tree_file_inventory(destination)
    source_map = {row["relative_path"]: row for row in source_rows}
    copied_map = {row["relative_path"]: row for row in copied_rows}
    if source_map != copied_map or source_summary != copied_summary:
        missing = sorted(set(source_map) - set(copied_map))[:5]
        extra = sorted(set(copied_map) - set(source_map))[:5]
        changed = sorted(
            path for path in set(source_map) & set(copied_map)
            if source_map[path] != copied_map[path]
        )[:5]
        raise ValueError(
            f"Copied run verification failed for {source_run}: "
            f"missing={missing}, extra={extra}, changed={changed}"
        )
    inventory = [
        {"source_run": source_run, "copied_root": str(destination), **row}
        for row in copied_rows
    ]
    return inventory, {
        "source_root": str(source),
        "copied_root": str(destination),
        "copy_strategy": copy_strategy,
        **copied_summary,
    }


def collect_phase_source_hashes(root: Path, source_run: str) -> list[dict]:
    """Inventory the exact source revision declared by each scientific phase."""
    rows = []
    phase_patterns = (
        ("arff_generation", "02_arff/**/PhenGO_manifest.json"),
        ("version_sensitivity", "04_ml/**/version_sensitivity_manifest.json"),
        ("temporal_analysis", "05_temporal/**/temporal_analysis_manifest.json"),
    )
    for phase, pattern in phase_patterns:
        for path in sorted(root.glob(pattern)):
            manifest = load_json(path)
            rows.append({
                "source_run": source_run,
                "phase": phase,
                "manifest_path": str(path.relative_to(root)),
                "schema_version": manifest.get("schema_version", ""),
                "source_tree_sha256": manifest.get("source_tree_sha256", ""),
                "git_commit": manifest.get("git_commit", ""),
                "tool_version": manifest.get("tool_version", ""),
            })
    return rows


def audit_scientific_run(root: Path, expected_source_hash: str | None = None) -> dict:
    """Fail closed when a completed run mixes code or unverifiable ARFFs."""
    root = root.resolve()
    errors = []
    warnings = []
    phase_rows = []
    arff_manifests = sorted(root.glob("02_arff/*/*/*/PhenGO_manifest.json"))
    arff_files = sorted(root.glob("02_arff/*/*/*/*_PhenGO.arff"))
    manifest_arffs = set()
    snapshot_ids = Counter()

    if not (root / "00_run_metadata/run.complete").is_file():
        errors.append("missing completion marker: 00_run_metadata/run.complete")
    if not arff_manifests:
        errors.append("no strict ARFF manifests were found")

    for manifest_path in arff_manifests:
        relative = manifest_path.relative_to(root)
        try:
            manifest = load_json(manifest_path)
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"unreadable manifest {relative}: {exc}")
            continue

        schema = manifest.get("schema_version")
        if not isinstance(schema, int) or schema < 3:
            errors.append(f"{relative}: manifest schema {schema!r} is older than 3")

        source_hash = manifest.get("source_tree_sha256")
        if not source_hash:
            errors.append(f"{relative}: missing source_tree_sha256")
        phase_rows.append({
            "phase": "arff_generation",
            "manifest_path": str(relative),
            "schema_version": schema,
            "source_tree_sha256": source_hash or "",
        })

        snapshot_id = manifest.get("snapshot_id")
        if snapshot_id:
            snapshot_ids[str(snapshot_id)] += 1
        else:
            errors.append(f"{relative}: missing snapshot_id")

        context = manifest.get("resource_context") or {}
        if context.get("snapshot_semantics") != "declared_composite_cross_section":
            errors.append(
                f"{relative}: snapshot semantics must be "
                "declared_composite_cross_section"
            )
        for field in (
            "phenotype_availability", "go_annotation_availability",
            "go_ontology_availability", "retrieval_route",
        ):
            value = str(context.get(field) or "").strip().lower()
            if value in {"", "unknown", "not_recorded", "unspecified"}:
                errors.append(f"{relative}: missing explicit resource_context.{field}")

        organism = relative.parts[2]
        expected_arff = manifest_path.parent / f"{organism}_PhenGO.arff"
        candidates = sorted(manifest_path.parent.glob("*_PhenGO.arff"))
        if expected_arff.is_file():
            arff_path = expected_arff
        elif len(candidates) == 1:
            arff_path = candidates[0]
        else:
            errors.append(f"{relative}: expected exactly one ARFF beside the manifest")
            continue
        manifest_arffs.add(arff_path.resolve())
        recorded = ((manifest.get("outputs") or {}).get("arff") or {}).get("sha256")
        if not recorded:
            errors.append(f"{relative}: missing recorded ARFF SHA-256")
            continue
        observed = sha256_file(arff_path)
        if observed != recorded:
            errors.append(
                f"{arff_path.relative_to(root)}: ARFF SHA-256 does not match its manifest"
            )

    for phase, pattern in (
        ("version_sensitivity", "04_ml/**/version_sensitivity_manifest.json"),
        ("temporal_analysis", "05_temporal/**/temporal_analysis_manifest.json"),
    ):
        for manifest_path in sorted(root.glob(pattern)):
            relative = manifest_path.relative_to(root)
            try:
                manifest = load_json(manifest_path)
            except (OSError, json.JSONDecodeError) as exc:
                errors.append(f"unreadable manifest {relative}: {exc}")
                continue
            source_hash = manifest.get("source_tree_sha256")
            if not source_hash:
                errors.append(f"{relative}: missing source_tree_sha256")
            phase_rows.append({
                "phase": phase,
                "manifest_path": str(relative),
                "schema_version": manifest.get("schema_version"),
                "source_tree_sha256": source_hash or "",
            })

    orphan_arffs = [
        str(path.relative_to(root)) for path in arff_files
        if path.resolve() not in manifest_arffs
    ]
    if orphan_arffs:
        errors.append(f"ARFF files without matching manifests: {orphan_arffs[:5]}")
    duplicates = sorted(key for key, count in snapshot_ids.items() if count > 1)
    if duplicates:
        errors.append(f"duplicate snapshot IDs: {duplicates[:5]}")

    source_hashes = sorted({
        row["source_tree_sha256"] for row in phase_rows
        if row["source_tree_sha256"]
    })
    if len(source_hashes) != 1:
        errors.append(
            "scientific manifests must declare exactly one source-tree hash; "
            f"observed {len(source_hashes)}"
        )
    elif expected_source_hash and source_hashes[0] != expected_source_hash:
        errors.append(
            "scientific run source-tree hash does not match the source revision "
            "launching this workflow"
        )
    if not any(row["phase"] == "version_sensitivity" for row in phase_rows):
        warnings.append("no version-sensitivity manifests were found")
    if not any(row["phase"] == "temporal_analysis" for row in phase_rows):
        warnings.append("no temporal-analysis manifests were found")

    return {
        "schema_version": 1,
        "audited_utc": datetime.now(timezone.utc).isoformat(),
        "run_root": str(root),
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "counts": {
            "arff_manifests": len(arff_manifests),
            "arff_files": len(arff_files),
            "scientific_manifests": len(phase_rows),
            "source_tree_hashes": len(source_hashes),
            "snapshot_ids": len(snapshot_ids),
        },
        "source_tree_sha256": source_hashes[0] if len(source_hashes) == 1 else None,
        "phase_manifests": phase_rows,
    }


def require_valid_scientific_run(
    root: Path, expected_source_hash: str | None = None
) -> dict:
    report = audit_scientific_run(root, expected_source_hash)
    if not report["valid"]:
        detail = "; ".join(report["errors"][:8])
        raise ValueError(f"Scientific run integrity failed for {root}: {detail}")
    return report


def available_bytes(path: Path) -> int:
    candidate = path
    while not candidate.exists() and candidate != candidate.parent:
        candidate = candidate.parent
    return shutil.disk_usage(candidate).free


def supports_apfs_clone(sources: tuple[Path, ...], destination_parent: Path) -> bool:
    """Probe macOS clonefile support on the exact source/destination devices."""
    if sys.platform != "darwin":
        return False
    destination_parent.mkdir(parents=True, exist_ok=True)
    try:
        destination_device = destination_parent.stat().st_dev
        if any(source.stat().st_dev != destination_device for source in sources):
            return False
        probe_source = next(
            path for source in sources for path in source.rglob("*") if path.is_file()
        )
    except (OSError, StopIteration):
        return False
    try:
        with tempfile.TemporaryDirectory(
            prefix=".phengo_clone_probe_", dir=destination_parent
        ) as directory:
            probe_destination = Path(directory) / "probe"
            completed = subprocess.run(
                ["cp", "-c", str(probe_source), str(probe_destination)],
                check=False, capture_output=True, text=True,
            )
            return (
                completed.returncode == 0
                and probe_destination.is_file()
                and sha256_file(probe_source) == sha256_file(probe_destination)
            )
    except OSError:
        return False


def availability_context(inventory: list[dict], events: list[dict]) -> None:
    for snapshot in inventory:
        year_text = snapshot.get("calendar_year", "")
        try:
            year = int(year_text)
        except (TypeError, ValueError):
            year = None
        matches = []
        for event in events:
            if event.get("organism") not in {snapshot.get("organism"), "all"}:
                continue
            try:
                start = int(event.get("start_year") or -10**9)
                end = int(event.get("end_year") or 10**9)
            except ValueError:
                continue
            if year is not None and start <= year <= end:
                matches.append(event.get("context_id", ""))
        snapshot["snapshot_semantics"] = (
            snapshot.get("snapshot_semantics") or "declared_composite_cross_section"
        )
        snapshot["temporal_alignment"] = "independently_released_components"
        snapshot["resource_context_ids"] = "|".join(sorted(item for item in matches if item))


def discover_snapshots(run_root: Path, source_run: str) -> list[dict]:
    snapshots = []
    arff_root = run_root / "02_arff"
    if not arff_root.is_dir():
        return snapshots
    for manifest_path in sorted(arff_root.glob("*/*/*/PhenGO_manifest.json")):
        relative = manifest_path.relative_to(arff_root)
        track, organism, timepoint = relative.parts[:3]
        manifest = load_json(manifest_path)
        arff_path = manifest_path.parent / f"{organism}_PhenGO.arff"
        if not arff_path.is_file():
            candidates = sorted(manifest_path.parent.glob("*_PhenGO.arff"))
            if len(candidates) != 1:
                raise ValueError(f"Expected one ARFF beside {manifest_path}")
            arff_path = candidates[0]
        recorded = ((manifest.get("outputs") or {}).get("arff") or {}).get("sha256")
        observed = sha256_file(arff_path)
        if recorded and recorded != observed:
            raise ValueError(f"ARFF hash mismatch: {arff_path}")
        snapshots.append({
            "source_run": source_run,
            "run_root": str(run_root),
            "track": track,
            "organism": organism,
            "timepoint": timepoint,
            "calendar_year": calendar_year(timepoint),
            "snapshot_id": manifest.get("snapshot_id", ""),
            "arff_path": str(arff_path),
            "arff_sha256": observed,
            "manifest_path": str(manifest_path),
            "manifest": manifest,
        })
    return snapshots


def calendar_year(value: str) -> str:
    match = re.search(r"(?:19|20)\d{2}", value or "")
    return match.group(0) if match else ""


def arff_contract(path: str) -> dict:
    attributes = []
    genes = {}
    in_data = False
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            if line.lower().startswith("@attribute"):
                match = ATTRIBUTE_RE.match(line)
                if not match:
                    raise ValueError(f"Malformed ARFF attribute in {path}: {line[:120]}")
                attributes.append(next(value for value in match.groups()[:3] if value is not None))
                continue
            if line.lower().startswith("@data"):
                in_data = True
                continue
            if not in_data:
                continue
            values = next(csv.reader([line]))
            if len(values) != len(attributes):
                raise ValueError(f"ARFF row width mismatch in {path}")
            gene = values[0].strip().strip("'\"")
            label = values[-1].strip().strip("'\"").lower()
            if label in POSITIVE_LABELS:
                encoded = 1
            elif label in NEGATIVE_LABELS:
                encoded = 0
            else:
                raise ValueError(f"Unsupported ARFF label {label!r} in {path}")
            if gene in genes and genes[gene] != encoded:
                raise ValueError(f"Conflicting duplicate label for {gene} in {path}")
            genes[gene] = encoded
    if len(attributes) < 3 or not genes:
        raise ValueError(f"Incomplete ARFF contract: {path}")
    return {"genes": genes, "features": set(attributes[1:-1])}


def inventory_row(snapshot: dict) -> dict:
    manifest = snapshot["manifest"]
    counts = manifest.get("counts") or {}
    policies = manifest.get("policies") or {}
    releases = manifest.get("releases") or {}
    inputs = manifest.get("inputs") or {}
    resource_context = manifest.get("resource_context") or {}
    return {
        key: snapshot[key] for key in (
            "source_run", "track", "organism", "timepoint", "calendar_year",
            "snapshot_id", "arff_path", "arff_sha256", "manifest_path",
        )
    } | {
        "n_genes": counts.get("genes", ""),
        "n_lethal": counts.get("lethal", ""),
        "n_viable": counts.get("viable", ""),
        "n_go_features": counts.get("go_features", ""),
        "lethal_prevalence": safe_div(counts.get("lethal"), counts.get("genes")),
        "phenotype_release": releases.get("phenotype", ""),
        "go_annotation_release": releases.get("go_annotations", ""),
        "go_ontology_release": releases.get("go_ontology", ""),
        "label_source": policies.get("label_source", ""),
        "nonlethal_policy": policies.get("nonlethal", ""),
        "mixed_label_policy": policies.get("mixed_label", ""),
        "ambiguous_label_policy": policies.get("ambiguous", ""),
        "filter_multigenes": policies.get("filter_multigenes", ""),
        "fly_driver_filtering": policies.get("fly_driver_filtering", ""),
        "go_relations": json.dumps(policies.get("go_relations", []), sort_keys=True),
        "go_evidence_exclude": json.dumps(policies.get("go_evidence_exclude", []), sort_keys=True),
        "phenotype_input_sha256": ((inputs.get("phenotype") or {}).get("sha256", "")),
        "gaf_input_sha256": ((inputs.get("go_annotations") or {}).get("sha256", "")),
        "obo_input_sha256": ((inputs.get("go_ontology") or {}).get("sha256", "")),
        "source_tree_sha256": manifest.get("source_tree_sha256", ""),
        "tool_version": manifest.get("tool_version", ""),
        "snapshot_semantics": resource_context.get(
            "snapshot_semantics", "declared_composite_cross_section"
        ),
        "phenotype_availability": resource_context.get(
            "phenotype_availability", "not_recorded_in_base_manifest"
        ),
        "go_annotation_availability": resource_context.get(
            "go_annotation_availability", "not_recorded_in_base_manifest"
        ),
        "go_ontology_availability": resource_context.get(
            "go_ontology_availability", "not_recorded_in_base_manifest"
        ),
        "retrieval_route": resource_context.get(
            "retrieval_route", "not_recorded_in_base_manifest"
        ),
    }


def safe_div(numerator, denominator):
    try:
        return float(numerator) / float(denominator) if float(denominator) else ""
    except (TypeError, ValueError):
        return ""


def set_metrics(left: set, right: set, prefix: str) -> dict:
    union = left | right
    return {
        f"{prefix}_reference": len(left),
        f"{prefix}_alternative": len(right),
        f"{prefix}_shared": len(left & right),
        f"{prefix}_reference_only": len(left - right),
        f"{prefix}_alternative_only": len(right - left),
        f"{prefix}_jaccard": len(left & right) / len(union) if union else 1.0,
    }


def compare_to_primary(base: list[dict], extension: list[dict]) -> list[dict]:
    references = {
        (item["organism"], item["calendar_year"]): item
        for item in base if item["track"] == "primary" and item["calendar_year"]
    }
    reference_cache = {}
    rows = []
    for alternative in extension:
        key = (alternative["organism"], alternative["calendar_year"])
        reference = references.get(key)
        if not reference:
            continue
        if key not in reference_cache:
            reference_cache[key] = arff_contract(reference["arff_path"])
        ref = reference_cache[key]
        alt = arff_contract(alternative["arff_path"])
        common = sorted(set(ref["genes"]) & set(alt["genes"]))
        transitions = Counter((ref["genes"][gene], alt["genes"][gene]) for gene in common)
        rows.append({
            "reference_source_run": reference["source_run"],
            "reference_track": reference["track"],
            "alternative_source_run": alternative["source_run"],
            "alternative_track": alternative["track"],
            "organism": alternative["organism"],
            "timepoint": alternative["timepoint"],
            "calendar_year": alternative["calendar_year"],
            **set_metrics(set(ref["genes"]), set(alt["genes"]), "genes"),
            **set_metrics(ref["features"], alt["features"], "go_features"),
            "shared_label_agreement": (
                sum(ref["genes"][gene] == alt["genes"][gene] for gene in common) / len(common)
                if common else ""
            ),
            "shared_label_churn": (
                sum(ref["genes"][gene] != alt["genes"][gene] for gene in common) / len(common)
                if common else ""
            ),
            "reference_viable_to_alternative_lethal": transitions[(0, 1)],
            "reference_lethal_to_alternative_viable": transitions[(1, 0)],
            "shared_viable_both": transitions[(0, 0)],
            "shared_lethal_both": transitions[(1, 1)],
            "alternative_gene_coverage_vs_primary": safe_div(len(alt["genes"]), len(ref["genes"])),
        })
    return rows


def read_csv_rows(path: Path, delimiter: str = ",") -> list[dict]:
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def collect_analysis_rows(run_root: Path, source_run: str, filename: str) -> list[dict]:
    rows = []
    ml_root = run_root / "04_ml"
    if not ml_root.is_dir():
        return rows
    for path in sorted(ml_root.glob(f"*/*/{filename}")):
        track, organism = path.relative_to(ml_root).parts[:2]
        for row in read_csv_rows(path):
            rows.append({"source_run": source_run, "track": track, "organism": organism, **row})
    return rows


def collect_first_analysis_rows(run_root: Path, source_run: str, filenames: tuple[str, ...]) -> list[dict]:
    for filename in filenames:
        rows = collect_analysis_rows(run_root, source_run, filename)
        if rows:
            return rows
    return []


def collect_publication_table(run_root: Path, source_run: str, filename: str) -> list[dict]:
    path = run_root / "06_publication_tables" / filename
    if not path.is_file():
        return []
    return [{"source_run": source_run, **row} for row in read_csv_rows(path, delimiter="\t")]


def add_dataset_interval_fields(rows: list[dict], left_field: str, right_field: str) -> None:
    for row in rows:
        left = calendar_year(row.get(left_field, ""))
        right = calendar_year(row.get(right_field, ""))
        row[f"{left_field}_calendar_year"] = left
        row[f"{right_field}_calendar_year"] = right
        if left and right:
            signed_gap = int(right) - int(left)
            row["calendar_gap_years"] = signed_gap
            row["absolute_calendar_gap_years"] = abs(signed_gap)
            row["consecutive_calendar_years"] = abs(signed_gap) == 1
        else:
            row["calendar_gap_years"] = ""
            row["absolute_calendar_gap_years"] = ""
            row["consecutive_calendar_years"] = ""


def performance_deltas(performance: list[dict], comparisons: list[dict]) -> list[dict]:
    references = {}
    for row in performance:
        if row["source_run"] == "base" and row["track"] == "primary":
            key = (row["organism"], calendar_year(row.get("dataset", "")), row.get("model"), row.get("panel"))
            references[key] = row
    comparison_lookup = {
        (row["alternative_track"], row["organism"], row["timepoint"]): row
        for row in comparisons
    }
    deltas = []
    for row in performance:
        if row["source_run"] != "extension":
            continue
        year = calendar_year(row.get("dataset", ""))
        key = (row["organism"], year, row.get("model"), row.get("panel"))
        reference = references.get(key)
        if not reference:
            continue
        out = {
            "track": row["track"], "organism": row["organism"],
            "timepoint": row.get("dataset", ""), "calendar_year": year,
            "model": row.get("model", ""), "panel": row.get("panel", ""),
        }
        comparison = comparison_lookup.get((row["track"], row["organism"], row.get("dataset", "")), {})
        for key_name in (
            "alternative_gene_coverage_vs_primary", "genes_jaccard", "go_features_jaccard",
            "shared_label_agreement",
        ):
            out[key_name] = comparison.get(key_name, "")
        for metric in PERFORMANCE_METRICS:
            field = f"{metric}_mean"
            try:
                alternative_value = float(row[field])
                reference_value = float(reference[field])
            except (KeyError, TypeError, ValueError):
                continue
            out[f"alternative_{field}"] = alternative_value
            out[f"primary_{field}"] = reference_value
            out[f"delta_{metric}"] = alternative_value - reference_value
        deltas.append(out)
    return deltas


def collect_gaf_comparisons(extension_root: Path) -> list[dict]:
    rows = []
    for path in sorted((extension_root / "03_qc" / "gaf_comparisons").glob("*/*.json")):
        row = load_json(path)
        relative = path.relative_to(extension_root / "03_qc" / "gaf_comparisons")
        row["organism"] = relative.parts[0]
        row["calendar_year"] = path.stem
        rows.append(row)
    return rows


def vocabulary_drift(repo_root: Path) -> tuple[list[dict], list[dict]]:
    directory = repo_root / "data/mouse/phenotype_data/historical_mgi_VOC_MammalianPhenotype"
    releases = {}
    inventory = []
    for path in sorted(directory.glob("VOC_MammalianPhenotype_*.rpt")):
        year = calendar_year(path.name)
        terms = {}
        with open(path, encoding="utf-8", errors="replace") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) >= 2 and re.fullmatch(r"MP:\d+", row[0]):
                    terms[row[0]] = row[1]
        releases[year] = terms
        inventory.append({
            "year": year, "path": str(path), "sha256": sha256_file(path),
            "n_terms": len(terms),
            "n_obsolete_named_terms": sum("obsolete" in name.lower() for name in terms.values()),
        })
    pairwise = []
    years = sorted(releases)
    for left_year, right_year in zip(years, years[1:]):
        left, right = releases[left_year], releases[right_year]
        shared = set(left) & set(right)
        union = set(left) | set(right)
        pairwise.append({
            "earlier_year": left_year, "later_year": right_year,
            "earlier_terms": len(left), "later_terms": len(right),
            "shared_terms": len(shared), "added_terms": len(set(right) - set(left)),
            "removed_terms": len(set(left) - set(right)),
            "renamed_shared_terms": sum(left[term] != right[term] for term in shared),
            "term_id_jaccard": len(shared) / len(union) if union else 1.0,
        })
    return inventory, pairwise


def float_value(row: dict, key: str):
    try:
        value = float(row.get(key, ""))
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def make_figures(output_dir: Path, inventory: list[dict], comparisons: list[dict], deltas: list[dict]) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    created = []
    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]

    organisms = sorted({row["organism"] for row in inventory})
    fig, axes = plt.subplots(len(organisms), 1, figsize=(11, 2.8 * max(1, len(organisms))), squeeze=False)
    for axis, organism in zip(axes[:, 0], organisms):
        subset = [row for row in inventory if row["organism"] == organism and row["calendar_year"]]
        tracks = sorted({row["track"] for row in subset})
        for index, track in enumerate(tracks):
            points = sorted(
                [(int(row["calendar_year"]), float(row["n_genes"])) for row in subset
                 if row["track"] == track and str(row.get("n_genes", ""))],
            )
            if points:
                axis.plot(*zip(*points), marker="o", ms=3, lw=1.2,
                          color=palette[index % len(palette)], label=track)
        axis.set_title(organism.capitalize())
        axis.set_ylabel("Genes in ARFF")
        axis.grid(axis="y", color="#dddddd", linewidth=0.6)
        axis.legend(fontsize=6, ncol=3, frameon=False, loc="upper left")
    axes[-1, 0].set_xlabel("Release year")
    fig.suptitle("Cohort size changes across resource and policy tracks", y=1.002)
    fig.tight_layout()
    path = figure_dir / "figure_cohort_sizes.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    created.append(str(path))

    if comparisons:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        for index, track in enumerate(sorted({row["alternative_track"] for row in comparisons})):
            subset = [row for row in comparisons if row["alternative_track"] == track]
            x = [float(row["genes_jaccard"]) for row in subset]
            y = [float(row["shared_label_agreement"]) for row in subset if row["shared_label_agreement"] != ""]
            x_for_y = [float(row["genes_jaccard"]) for row in subset if row["shared_label_agreement"] != ""]
            axes[0].scatter(x_for_y, y, s=22, alpha=0.75, color=palette[index % len(palette)], label=track)
            axes[1].scatter(
                [float(row["go_features_jaccard"]) for row in subset],
                [float(row["alternative_gene_coverage_vs_primary"]) for row in subset],
                s=22, alpha=0.75, color=palette[index % len(palette)], label=track,
            )
        axes[0].set(xlabel="Gene-set Jaccard vs primary", ylabel="Shared-gene label agreement")
        axes[1].set(xlabel="GO-feature Jaccard vs primary", ylabel="Gene-count ratio vs primary")
        for axis in axes:
            axis.grid(color="#dddddd", linewidth=0.6)
            axis.legend(fontsize=6, frameon=False)
        fig.suptitle("Dataset-construction sensitivity relative to the primary track")
        fig.tight_layout()
        path = figure_dir / "figure_policy_sensitivity.png"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        created.append(str(path))

    plot_rows = [
        row for row in deltas
        if row.get("panel") == "full"
        and float_value(row, "alternative_gene_coverage_vs_primary") is not None
        and float_value(row, "alternative_average_precision_mean") is not None
    ]
    if plot_rows:
        fig, axis = plt.subplots(figsize=(10, 6))
        tracks = sorted({row["track"] for row in plot_rows})
        for index, track in enumerate(tracks):
            subset = [row for row in plot_rows if row["track"] == track]
            axis.scatter(
                [float(row["alternative_gene_coverage_vs_primary"]) for row in subset],
                [float(row["alternative_average_precision_mean"]) for row in subset],
                s=24, alpha=0.7, color=palette[index % len(palette)], label=track,
            )
        axis.set_xlabel("Alternative cohort size / primary cohort size")
        axis.set_ylabel("Within-snapshot average precision")
        axis.set_title("Predictability changes with dataset construction and coverage")
        axis.grid(color="#dddddd", linewidth=0.6)
        axis.legend(fontsize=6, ncol=2, frameon=False)
        fig.tight_layout()
        path = figure_dir / "figure_predictability_vs_coverage.png"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        created.append(str(path))
    return created


def run_fingerprint(root: Path) -> dict:
    metadata = root / "00_run_metadata"
    result = {"root": str(root)}
    for filename in ("run.complete", "run_fingerprint.sha256", "output_checksums.sha256"):
        path = metadata / filename
        result[filename.replace(".", "_")] = {
            "path": str(path), "exists": path.is_file(),
            "sha256": sha256_file(path) if path.is_file() else "",
        }
    return result


def write_report(
    path: Path, base_root: Path, extension_root: Path, inventory: list[dict],
    comparisons: list[dict], performance: list[dict], deltas: list[dict],
    input_audit: list[dict], figures: list[str], phase_source_hashes: list[dict],
    copy_strategy: str,
) -> None:
    tracks = sorted({row["track"] for row in inventory})
    status_counts = Counter(row.get("status", "") for row in input_audit)
    distinct_phase_hashes = sorted({
        row.get("source_tree_sha256", "") for row in phase_source_hashes
        if row.get("source_tree_sha256", "")
    })
    phase_hash_statement = (
        "All phase manifests declare one source-tree hash."
        if len(distinct_phase_hashes) == 1 else
        "Phase manifests declare multiple source-tree hashes; treat this as an "
        "explicit code transition and inspect `phase_source_hash_inventory.tsv`."
        if distinct_phase_hashes else
        "No phase-level source-tree hashes were available."
    )
    lines = [
        "# PhenGO Consolidated Publication Analysis",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Analysis Contract",
        "",
        f"The completed first run has been copied and verified at `{base_root}`.",
        f"The completed V2 extension has been copied and verified at `{extension_root}`. Every ARFF is linked to its strict PhenGO manifest and SHA-256 digest.",
        "This directory is self-contained: its source-run file inventory records a SHA-256 digest for every copied file.",
        f"Copy strategy: `{copy_strategy}`. APFS clones are independent copy-on-write files; `copy2` denotes ordinary copies.",
        phase_hash_statement,
        "Each nominal year is treated as a declared composite cross-section of independently released phenotype, GAF, and ontology resources, not as a synchronized database release.",
        "The primary question is whether data release and curation choices alter cohort composition, labels, GO features, apparent within-snapshot predictability, cross-snapshot transfer, and feature importance.",
        "",
        "Higher accuracy is not treated as evidence that a label policy is more biologically correct. It may instead reflect easier class separation, leakage-prone curation structure, altered prevalence, or a narrower cohort.",
        "",
        "## Coverage",
        "",
        f"- Snapshots indexed: {len(inventory)}",
        f"- Tracks indexed: {len(tracks)} ({', '.join(tracks)})",
        f"- Primary-versus-alternative matched comparisons: {len(comparisons)}",
        f"- Within-snapshot model summary rows: {len(performance)}",
        f"- Performance delta rows matched to primary: {len(deltas)}",
        f"- Alternative input audit statuses: {dict(sorted(status_counts.items()))}",
        "",
        "## Principal Tables",
        "",
        "- `source_run_file_inventory.tsv`: every physically copied V1 and V2 file with its relative path, size, and SHA-256 digest.",
        "- `phase_source_hash_inventory.tsv`: source-tree hash declared by every ARFF, version-sensitivity, and temporal-analysis manifest.",
        "- `resource_availability.tsv`: provider-confirmed and retrieval-history context for non-random archive gaps.",
        "- `snapshot_inventory.tsv`: one row per ARFF with policies, releases, hashes, cohort size, class balance, and GO-feature count.",
        "- `primary_alternative_comparisons.tsv`: gene and feature overlap plus shared-gene label transitions for each alternative against the same organism/year primary snapshot.",
        "- `within_year_cv_all_tracks.csv`: the complete repeated-CV summaries from both runs.",
        "- `performance_deltas_vs_primary.csv`: matched model/panel performance differences with cohort and label-overlap covariates.",
        "- `cross_year_transfer_all_tracks.csv`: all gene-disjoint cross-year transfer summaries.",
        "- `previous_available_snapshot_label_baseline_all_tracks.csv`: adjacent available-snapshot label persistence with explicit calendar gaps.",
        "- `dataset_drift_all_tracks.csv` and `pairwise_dataset_drift_all_tracks.csv`: cohort, label, and feature drift from every track.",
        "- `evaluation_preflight_all_tracks.csv`: feasibility and class-count checks for all model evaluations.",
        "- `feature_importance_stability_all_tracks.csv`: bootstrap feature-selection stability across all tracks.",
        "- `temporal_summary_all_tracks.tsv` and `temporal_enrichment_statistics_all_tracks.tsv`: combined GO and phenotype temporal outputs.",
        "- `gaf_source_comparisons.tsv`: relation-insensitive stable-ID/GO annotation overlap between GO archive and MOD-provided GAFs.",
        "- `alternative_input_audit.tsv`: all discovered alternative files, including explicit duplicate and truncation exclusions.",
        "- `mouse_vocabulary_inventory.tsv` and `mouse_vocabulary_pairwise.tsv`: available MP vocabulary snapshots and term churn; these reports lack hierarchy and are not used as ontology graphs.",
        "",
        "## Figure Interpretation",
        "",
    ]
    lines.extend(f"- `{Path(item).name}`" for item in figures)
    lines.extend([
        "",
        "The coverage figures establish how apparently minor construction choices change which genes are modelled. The policy-sensitivity panel separates gene turnover, GO-feature turnover, and label disagreement. The predictability-versus-coverage panel should be read as a diagnostic: a performance increase accompanied by strong cohort contraction is not automatically a biological improvement.",
        "",
        "## Required Manuscript Caveats",
        "",
        "The IMPC assertion track and historical MGI association tracks do not contain a complete experimentally verified viable class. Their negative class is operational: genes with other phenotype assertions but no selected lethal MP term. Results from these tracks test source-definition sensitivity and must not be described as definitive viable-versus-lethal prediction.",
        "",
        "Wayback captures are included only when the decompressed file is complete and not an exact duplicate of another selected capture. Distinct captures from the same nominal year retain separate capture IDs. Nominal archive year is not assumed to be the exact resource release date.",
        "",
        "The three available historical MP vocabulary reports contain identifiers, names, and definitions but no parent-child hierarchy. They support vocabulary-churn auditing only; contemporary hierarchy is not back-projected onto historical reports.",
        "",
        "Alternative MOD GAF comparisons use stable gene ID plus GO ID as the principal annotation key because GAF 2.1 and 2.2 encode relation qualifiers differently. Evidence-code-aware overlap is reported separately.",
        "",
        "All model comparisons remain observational sensitivity analyses. They cannot independently distinguish biological change from curation change, annotation growth, assay composition, identifier turnover, or ascertainment bias.",
        "",
        "## Recommended Article Sequence",
        "",
        "1. Establish input and cohort drift with the inventory, GAF comparisons, and MP vocabulary audit.",
        "2. Show primary year-to-year instability using the first run's matched panels and transfer matrices.",
        "3. Introduce alternative construction tracks as controlled perturbations and report gene, label, and feature overlap before ML metrics.",
        "4. Compare within-snapshot predictability only alongside coverage, prevalence, and shared-label agreement.",
        "5. Use cross-year transfer and prediction-instability outputs as the strongest evidence that conclusions depend on resource snapshot.",
        "6. Treat feature-rank instability as instability in the model's apparent biological explanation, not proof that the underlying biology changed.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-run")
    parser.add_argument("--audit-report")
    parser.add_argument("--expected-source-hash")
    parser.add_argument("--base-run")
    parser.add_argument("--extension-run")
    parser.add_argument("--repo-root")
    parser.add_argument("--output-dir")
    parser.add_argument("--input-audit")
    parser.add_argument("--availability-ledger")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.audit_run:
        report = audit_scientific_run(
            Path(args.audit_run), args.expected_source_hash
        )
        if args.audit_report:
            write_json(Path(args.audit_report), report)
        print(json.dumps({
            "run_root": report["run_root"],
            "valid": report["valid"],
            "counts": report["counts"],
            "errors": report["errors"][:8],
        }, indent=2, sort_keys=True))
        raise SystemExit(0 if report["valid"] else 1)
    missing = [
        option for option, value in (
            ("--base-run", args.base_run),
            ("--extension-run", args.extension_run),
            ("--repo-root", args.repo_root),
            ("--output-dir", args.output_dir),
        ) if not value
    ]
    if missing:
        parser.error(f"required arguments are missing: {', '.join(missing)}")

    base_root = Path(args.base_run).resolve()
    extension_root = Path(args.extension_run).resolve()
    output_dir = Path(args.output_dir).resolve()
    try:
        base_integrity = require_valid_scientific_run(base_root)
        extension_integrity = require_valid_scientific_run(extension_root)
    except ValueError as exc:
        parser.error(str(exc))
    if base_integrity["source_tree_sha256"] != extension_integrity["source_tree_sha256"]:
        parser.error(
            "Base and extension runs were generated by different source-tree revisions"
        )
    overlaps_source = (
        output_dir in {base_root, extension_root}
        or output_dir in base_root.parents
        or output_dir in extension_root.parents
        or base_root in output_dir.parents
        or extension_root in output_dir.parents
    )
    if overlaps_source:
        parser.error("Consolidated output must be separate from both source runs")
    if output_dir.exists() and any(output_dir.iterdir()):
        if not args.overwrite:
            parser.error(f"Output is not empty: {output_dir}")
        shutil.rmtree(output_dir)
    logical_copy_bytes = sum(
        path.stat().st_size
        for root in (base_root, extension_root)
        for path in root.rglob("*") if path.is_file()
    )
    copy_strategy = (
        "apfs_clone"
        if supports_apfs_clone((base_root, extension_root), output_dir.parent)
        else "copy2"
    )
    reserve = max(
        CLONE_RESERVE_BYTES if copy_strategy == "apfs_clone" else COPY_RESERVE_BYTES,
        int(logical_copy_bytes * 0.05),
    )
    required_available_bytes = (
        reserve if copy_strategy == "apfs_clone" else logical_copy_bytes + reserve
    )
    if available_bytes(output_dir.parent) < required_available_bytes:
        parser.error(
            "Insufficient free space for a self-contained consolidated copy: "
            f"strategy={copy_strategy}, need at least {required_available_bytes} bytes available"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    copied_root = output_dir / "source_runs"
    copied_root.mkdir()
    copy_inventory = []
    copied_summaries = {}
    for label, source, dirname in (
        ("base", base_root, "v1"),
        ("extension", extension_root, "v2_extension"),
    ):
        records, summary = copy_verified_run(
            source, copied_root / dirname, label, copy_strategy
        )
        copy_inventory.extend(records)
        copied_summaries[label] = summary

    copied_base = copied_root / "v1"
    copied_extension = copied_root / "v2_extension"
    phase_source_hashes = (
        collect_phase_source_hashes(copied_base, "base")
        + collect_phase_source_hashes(copied_extension, "extension")
    )
    base = discover_snapshots(copied_base, "base")
    extension = discover_snapshots(copied_extension, "extension")
    if not base or not extension:
        parser.error("Both source runs must contain strict ARFF snapshots")
    inventory = [inventory_row(item) for item in base + extension]
    comparisons = compare_to_primary(base, extension)
    performance = (
        collect_analysis_rows(copied_base, "base", "within_year_cv_summary.csv")
        + collect_analysis_rows(copied_extension, "extension", "within_year_cv_summary.csv")
    )
    transfer = (
        collect_analysis_rows(copied_base, "base", "cross_year_transfer_summary.csv")
        + collect_analysis_rows(copied_extension, "extension", "cross_year_transfer_summary.csv")
    )
    add_dataset_interval_fields(transfer, "train_dataset", "test_dataset")
    previous_available = (
        collect_first_analysis_rows(
            copied_base, "base",
            ("previous_available_snapshot_label_baseline.csv", "previous_year_label_baseline.csv"),
        )
        + collect_first_analysis_rows(
            copied_extension, "extension",
            ("previous_available_snapshot_label_baseline.csv", "previous_year_label_baseline.csv"),
        )
    )
    add_dataset_interval_fields(previous_available, "train_dataset", "test_dataset")
    dataset_drift = (
        collect_analysis_rows(copied_base, "base", "dataset_drift.csv")
        + collect_analysis_rows(copied_extension, "extension", "dataset_drift.csv")
    )
    pairwise_drift = (
        collect_analysis_rows(copied_base, "base", "pairwise_drift.csv")
        + collect_analysis_rows(copied_extension, "extension", "pairwise_drift.csv")
    )
    evaluation_preflight = (
        collect_analysis_rows(copied_base, "base", "evaluation_preflight.csv")
        + collect_analysis_rows(copied_extension, "extension", "evaluation_preflight.csv")
    )
    instability = (
        collect_analysis_rows(copied_base, "base", "prediction_instability_summary.csv")
        + collect_analysis_rows(copied_extension, "extension", "prediction_instability_summary.csv")
    )
    feature_overlap = (
        collect_analysis_rows(copied_base, "base", "feature_rank_overlap.csv")
        + collect_analysis_rows(copied_extension, "extension", "feature_rank_overlap.csv")
    )
    feature_stability = (
        collect_analysis_rows(copied_base, "base", "feature_importance_stability.csv")
        + collect_analysis_rows(copied_extension, "extension", "feature_importance_stability.csv")
    )
    temporal_summary = (
        collect_publication_table(copied_base, "base", "temporal_summary.tsv")
        + collect_publication_table(copied_extension, "extension", "temporal_summary.tsv")
    )
    temporal_enrichment = (
        collect_publication_table(copied_base, "base", "temporal_enrichment_statistics.tsv")
        + collect_publication_table(copied_extension, "extension", "temporal_enrichment_statistics.tsv")
    )
    deltas = performance_deltas(performance, comparisons)
    gaf_comparisons = collect_gaf_comparisons(copied_extension)
    voc_inventory, voc_pairwise = vocabulary_drift(Path(args.repo_root).resolve())
    if args.input_audit:
        supplied_audit = Path(args.input_audit).resolve()
        try:
            copied_audit = copied_extension / supplied_audit.relative_to(extension_root)
        except ValueError:
            copied_audit = supplied_audit
        input_audit = read_csv_rows(copied_audit, delimiter="\t")
    else:
        input_audit = []
    availability = (
        read_csv_rows(Path(args.availability_ledger), delimiter="\t")
        if args.availability_ledger else []
    )
    availability_context(inventory, availability)

    write_tsv(output_dir / "source_run_file_inventory.tsv", copy_inventory,
              ("source_run", "relative_path", "size_bytes", "sha256", "copied_root"))
    write_tsv(output_dir / "phase_source_hash_inventory.tsv", phase_source_hashes,
              ("source_run", "phase", "manifest_path", "source_tree_sha256"))
    write_json(output_dir / "base_run_integrity.json", base_integrity)
    write_json(output_dir / "extension_run_integrity.json", extension_integrity)
    write_tsv(output_dir / "resource_availability.tsv", availability,
              ("context_id", "organism", "resource", "start_year", "end_year"))
    write_tsv(output_dir / "snapshot_inventory.tsv", inventory,
              ("source_run", "track", "organism", "timepoint", "calendar_year"))
    write_tsv(output_dir / "primary_alternative_comparisons.tsv", comparisons,
              ("alternative_track", "organism", "timepoint", "calendar_year"))
    write_csv(output_dir / "within_year_cv_all_tracks.csv", performance,
              ("source_run", "track", "organism", "dataset", "model", "panel"))
    write_csv(output_dir / "performance_deltas_vs_primary.csv", deltas,
              ("track", "organism", "timepoint", "model", "panel"))
    write_csv(output_dir / "cross_year_transfer_all_tracks.csv", transfer,
              ("source_run", "track", "organism", "train_dataset", "test_dataset", "model", "panel"))
    write_csv(output_dir / "previous_available_snapshot_label_baseline_all_tracks.csv",
              previous_available,
              ("source_run", "track", "organism", "train_dataset", "test_dataset"))
    write_csv(output_dir / "dataset_drift_all_tracks.csv", dataset_drift,
              ("source_run", "track", "organism", "dataset"))
    write_csv(output_dir / "pairwise_dataset_drift_all_tracks.csv", pairwise_drift,
              ("source_run", "track", "organism"))
    write_csv(output_dir / "evaluation_preflight_all_tracks.csv", evaluation_preflight,
              ("source_run", "track", "organism"))
    write_csv(output_dir / "prediction_instability_all_tracks.csv", instability)
    write_csv(output_dir / "feature_rank_overlap_all_tracks.csv", feature_overlap)
    write_csv(output_dir / "feature_importance_stability_all_tracks.csv", feature_stability)
    write_tsv(output_dir / "temporal_summary_all_tracks.tsv", temporal_summary,
              ("source_run", "track", "organism"))
    write_tsv(output_dir / "temporal_enrichment_statistics_all_tracks.tsv", temporal_enrichment,
              ("source_run", "track", "organism"))
    write_tsv(output_dir / "gaf_source_comparisons.tsv", gaf_comparisons,
              ("organism", "calendar_year"))
    write_tsv(output_dir / "alternative_input_audit.tsv", input_audit)
    write_tsv(output_dir / "mouse_vocabulary_inventory.tsv", voc_inventory, ("year",))
    write_tsv(output_dir / "mouse_vocabulary_pairwise.tsv", voc_pairwise,
              ("earlier_year", "later_year"))
    figures = make_figures(output_dir, inventory, comparisons, deltas)
    write_report(
        output_dir / "CONSOLIDATED_REPORT.md", copied_base, copied_extension,
        inventory, comparisons, performance, deltas, input_audit, figures,
        phase_source_hashes, copy_strategy,
    )

    source_runs = {
        "base": {
            "original": run_fingerprint(base_root),
            "copied": run_fingerprint(copied_base),
            "copy_verification": copied_summaries["base"],
        },
        "extension": {
            "original": run_fingerprint(extension_root),
            "copied": run_fingerprint(copied_extension),
            "copy_verification": copied_summaries["extension"],
        },
    }
    output_hashes = {}
    for path in sorted(item for item in output_dir.rglob("*") if item.is_file()):
        relative = path.relative_to(output_dir)
        if path.name == "consolidated_manifest.json" or relative.parts[0] == "source_runs":
            continue
        output_hashes[str(relative)] = sha256_file(path)
    write_json(output_dir / "consolidated_manifest.json", {
        "schema_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "PhenGO publication V2 consolidated multiverse",
        "copy_strategy": copy_strategy,
        "logical_source_bytes": logical_copy_bytes,
        "copy_space_reserve_bytes": reserve,
        "command": list(os.sys.argv),
        "source_runs": source_runs,
        "scientific_run_integrity": {
            "base": base_integrity,
            "extension": extension_integrity,
        },
        "counts": {
            "snapshots": len(inventory), "comparisons": len(comparisons),
            "within_year_summary_rows": len(performance), "transfer_summary_rows": len(transfer),
            "previous_available_snapshot_rows": len(previous_available),
            "dataset_drift_rows": len(dataset_drift),
            "pairwise_drift_rows": len(pairwise_drift),
            "evaluation_preflight_rows": len(evaluation_preflight),
            "performance_delta_rows": len(deltas), "input_audit_rows": len(input_audit),
            "feature_stability_rows": len(feature_stability),
            "temporal_summary_rows": len(temporal_summary),
            "temporal_enrichment_rows": len(temporal_enrichment),
            "availability_context_rows": len(availability),
            "phase_source_hash_rows": len(phase_source_hashes),
            "copied_source_files": len(copy_inventory),
        },
        "outputs": output_hashes,
    })


if __name__ == "__main__":
    main()
