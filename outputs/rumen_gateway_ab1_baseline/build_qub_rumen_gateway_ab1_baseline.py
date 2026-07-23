#!/usr/bin/env python3
"""Build the QUB Rumen Gateway AB1 baseline TSV from the two source workbooks."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_ALL16 = Path("/Users/nicholas/Downloads/all_16S_QUB_formatted_14_12_2025_live_streamlined.xlsx")
DEFAULT_OVERVIEW = Path(
    "/Users/nicholas/Library/CloudStorage/OneDrive-SharedLibraries-Queen'sUniversityBelfast/"
    "RUMEN Gateway QUB-Scientific_Consortium - Documents/03_Partners_sequences/"
    "Current_Sequences_Overview.xlsx"
)
DEFAULT_AB1_DIR = Path(
    "/Users/nicholas/Library/CloudStorage/OneDrive-SharedLibraries-Queen'sUniversityBelfast/"
    "RUMEN Gateway QUB-Scientific_Consortium - Documents/03_Partners_sequences/"
    "QUB/ab1_files"
)

BASELINE_COLUMNS = [
    "baseline_row",
    "match_status",
    "match_score",
    "match_reasons",
    "raw_candidate_count",
    "raw_tie_count",
    "selected_match_note",
    "alternate_raw_candidates",
    "partner_id",
    "partner_isolate_id",
    "partner_isolate_id_status",
    "sh_id",
    "qub_id",
    "sequence_identifier",
    "sequence_primer_from_identifier",
    "sequence_bp",
    "overview_row",
    "overview_organism_type",
    "overview_gtdb_assignment",
    "overview_gtdb_taxon",
    "overview_gtdb_pident",
    "overview_closest_hungate_hit",
    "overview_hungate_hit_pident",
    "overview_jgi_sequenced",
    "overview_mwl",
    "overview_metadata",
    "overview_notes",
    "all16s_row",
    "all16s_db_control_seq",
    "all16s_project",
    "all16s_responsible",
    "all16s_status_info",
    "all16s_bc_status",
    "all16s_sh_id",
    "all16s_project_isolate_id",
    "all16s_short_sequence_code",
    "all16s_original_fasta_id",
    "all16s_order_date",
    "all16s_order_id",
    "all16s_type",
    "all16s_plate_well_id",
    "all16s_used_date",
    "all16s_sequence_order_id",
    "all16s_primer_direction",
    "all16s_sequencing_primer",
    "all16s_base_pairs",
    "all16s_location_of_stock",
    "all16s_method",
    "all16s_isolation_time",
    "all16s_isolation_medium",
    "all16s_year",
    "all16s_growth_medium",
    "all16s_colony_size",
    "all16s_colony_colour",
    "all16s_colony_surface",
    "all16s_maldi_tof",
    "all16s_temperature",
    "all16s_atmosphere",
    "all16s_rf_project_reference",
    "all16s_rf_collection_date",
    "all16s_animal_species",
    "all16s_animal_slaughter_id",
    "all16s_animal_id",
    "all16s_production_system_animal_diet",
    "all16s_animal_week_of_age",
    "ab1_expected_filename",
    "ab1_expected_path",
    "ab1_present",
    "ab1_status",
    "ab1_size_bytes",
    "ab1_modified_time",
    "ab1_directory_listing_status",
    "ab1_directory_listing_error",
]


def text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def clean_sequence_id(value: Any) -> str:
    value_text = text(value)
    if value_text.startswith(">"):
        return value_text[1:]
    return value_text


def emptyish(value: Any) -> bool:
    return text(value).upper() in {"", "NA", "NR", "NAN", "NONE"}


def parse_intish(value: Any) -> int | None:
    value_text = text(value)
    if emptyish(value_text):
        return None
    try:
        return int(float(value_text))
    except ValueError:
        return None


def parse_sequence_primer(sequence_identifier: str) -> str:
    match = re.match(r"^QUB_(.+)_\d+$", text(sequence_identifier), re.IGNORECASE)
    return match.group(1) if match else ""


def canonical_primer(value: Any) -> str:
    return text(value).upper()


def partner_isolate_status(value: str) -> str:
    return "SH_ID" if re.match(r"^SH_\d", text(value), re.IGNORECASE) else "NON_SH_ID"


def read_sources(all16_path: Path, overview_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    all16 = pd.read_excel(all16_path, sheet_name="all_16S_QUB RG", header=1, dtype=str).fillna("")
    all16["_excel_row"] = all16.index + 3

    overview = pd.read_excel(overview_path, sheet_name="Summary", dtype=str).fillna("")
    overview["_excel_row"] = overview.index + 2
    return all16, overview


def inspect_ab1_directory(ab1_dir: Path) -> tuple[str, str, set[str]]:
    try:
        files = {
            path.name
            for path in ab1_dir.iterdir()
            if path.is_file() and path.suffix.lower() == ".ab1"
        }
    except Exception as exc:  # noqa: BLE001 - keep exact OS/TCC error in audit output.
        return "UNAVAILABLE", f"{type(exc).__name__}: {exc}", set()
    return "AVAILABLE", "", files


def project_isolate_contains(project_isolate_id: str, partner_id: str) -> bool:
    if not partner_id:
        return False
    tokens = [part for part in re.split(r"[;\s,]+", project_isolate_id) if part]
    return project_isolate_id == partner_id or partner_id in tokens


def score_raw_candidate(raw_row: dict[str, Any], summary_row: pd.Series) -> tuple[int, list[str]]:
    partner_isolate_id = text(summary_row["Partner Rumen Gateway isolate ID"])
    expected_primer = canonical_primer(parse_sequence_primer(text(summary_row["Sequence_identifier"])))
    expected_bp = parse_intish(summary_row["Sequence_bp"])

    score = 0
    reasons: list[str] = []
    raw_sh_id = text(raw_row.get("SH_ID", ""))
    raw_project_isolate_id = text(raw_row.get("Project_isolate_ID", ""))
    raw_short_code = clean_sequence_id(raw_row.get("Short_sequence_code", ""))
    raw_original_fasta_id = clean_sequence_id(raw_row.get("Original_fasta_ID", ""))

    if partner_isolate_id and raw_sh_id == partner_isolate_id:
        score += 100
        reasons.append("sh_id")
    if partner_isolate_id and raw_short_code == partner_isolate_id:
        score += 80
        reasons.append("short_sequence_code")
    if partner_isolate_id and (
        raw_original_fasta_id == partner_isolate_id or partner_isolate_id in raw_original_fasta_id
    ):
        score += 60
        reasons.append("original_fasta_id_contains")
    if project_isolate_contains(raw_project_isolate_id, partner_isolate_id):
        score += 50
        reasons.append("project_isolate_id")

    if not reasons:
        return 0, []

    raw_primer = canonical_primer(raw_row.get("Sequencing_primer", ""))
    raw_bp = parse_intish(raw_row.get("Base_pairs", ""))
    if expected_primer and raw_primer == expected_primer:
        score += 25
        reasons.append("primer")
    elif expected_primer == "?" and raw_primer in {"", "NR", "?"}:
        score += 15
        reasons.append("unknown_primer")
    if expected_bp is not None and raw_bp == expected_bp:
        score += 25
        reasons.append("base_pairs")
    if text(raw_row.get("BC_Status", "")) == "USED":
        score += 5
        reasons.append("used_bc")

    return score, reasons


def match_status(best: tuple[int, list[str], dict[str, Any]] | None, tie_count: int) -> str:
    if best is None:
        return "NO_RAW_METADATA_MATCH"
    reasons = set(best[1])
    if tie_count > 1:
        return "AMBIGUOUS_RAW_MATCH_SELECTED_EARLIEST_ROW"
    if {"sh_id", "primer", "base_pairs"}.issubset(reasons):
        return "MATCHED_SH_PRIMER_BP"
    if {"sh_id", "primer"}.issubset(reasons):
        return "MATCHED_SH_PRIMER_ONLY"
    if {"sh_id", "base_pairs"}.issubset(reasons):
        return "MATCHED_SH_BP_ONLY"
    if {"short_sequence_code", "base_pairs"}.issubset(reasons) or {
        "original_fasta_id_contains",
        "base_pairs",
    }.issubset(reasons):
        return "MATCHED_NON_SH_RAW_CODE"
    return "MATCHED_ID_ONLY"


def raw_candidate_label(candidate: tuple[int, list[str], dict[str, Any]]) -> str:
    score, reasons, row = candidate
    pieces = [
        f"row={text(row.get('_excel_row', ''))}",
        f"db={text(row.get('Extra DB control seq', ''))}",
        f"short={clean_sequence_id(row.get('Short_sequence_code', ''))}",
        f"fasta={clean_sequence_id(row.get('Original_fasta_ID', ''))}",
        f"primer={text(row.get('Sequencing_primer', ''))}",
        f"bp={text(row.get('Base_pairs', ''))}",
        f"score={score}",
        f"reasons={','.join(reasons)}",
    ]
    return "|".join(pieces)


def selected_match_note(status: str) -> str:
    if status == "AMBIGUOUS_RAW_MATCH_SELECTED_EARLIEST_ROW":
        return "Two or more raw sequencing rows tied on the strongest evidence; earliest all_16S row selected and alternates listed."
    if status == "MATCHED_NON_SH_RAW_CODE":
        return "Overview isolate ID was not recorded as an SH_ID in all_16S; matched by raw sequence code and base pairs."
    if status == "NO_RAW_METADATA_MATCH":
        return "No matching raw sequencing metadata row found in all_16S."
    return ""


def add_raw_fields(output: dict[str, Any], raw_row: dict[str, Any] | None) -> None:
    field_map = {
        "all16s_row": "_excel_row",
        "all16s_db_control_seq": "Extra DB control seq",
        "all16s_project": "Project",
        "all16s_responsible": "Responsible",
        "all16s_status_info": "Status/info",
        "all16s_bc_status": "BC_Status",
        "all16s_sh_id": "SH_ID",
        "all16s_project_isolate_id": "Project_isolate_ID",
        "all16s_short_sequence_code": "Short_sequence_code",
        "all16s_original_fasta_id": "Original_fasta_ID",
        "all16s_order_date": "Order_date",
        "all16s_order_id": "Order_ID",
        "all16s_type": "Type",
        "all16s_plate_well_id": "Plate_well_ID",
        "all16s_used_date": "Used_date",
        "all16s_sequence_order_id": "Sequence_order_ID",
        "all16s_primer_direction": "Primer_direction",
        "all16s_sequencing_primer": "Sequencing_primer",
        "all16s_base_pairs": "Base_pairs",
        "all16s_location_of_stock": "Location of stock",
        "all16s_method": "Method",
        "all16s_isolation_time": "Isolation time",
        "all16s_isolation_medium": "Isolation_medium",
        "all16s_year": "Year",
        "all16s_growth_medium": "Growth_medium",
        "all16s_colony_size": "Colony_size",
        "all16s_colony_colour": "Colony_colour",
        "all16s_colony_surface": "Colony_surface",
        "all16s_maldi_tof": "MALDI TOF",
        "all16s_temperature": "Temperature",
        "all16s_atmosphere": "Atmosphere",
        "all16s_rf_project_reference": "RF_Project_reference",
        "all16s_rf_collection_date": "RF_collection_date",
        "all16s_animal_species": "Animal_species",
        "all16s_animal_slaughter_id": "Animal_slaugther_ID",
        "all16s_animal_id": "Animal_ID",
        "all16s_production_system_animal_diet": "Production_system/Animal diet",
        "all16s_animal_week_of_age": "Animal_week_of_age",
    }
    for output_column, source_column in field_map.items():
        output[output_column] = text(raw_row.get(source_column, "")) if raw_row else ""


def ab1_status_for(raw_row: dict[str, Any] | None, ab1_dir: Path) -> dict[str, Any]:
    if raw_row is None:
        return {
            "ab1_expected_filename": "",
            "ab1_expected_path": "",
            "ab1_present": "FALSE",
            "ab1_status": "NO_RAW_METADATA_MATCH",
            "ab1_size_bytes": "",
            "ab1_modified_time": "",
        }

    original_fasta_id = clean_sequence_id(raw_row.get("Original_fasta_ID", ""))
    if emptyish(original_fasta_id):
        return {
            "ab1_expected_filename": "",
            "ab1_expected_path": "",
            "ab1_present": "FALSE",
            "ab1_status": "NO_EXPECTED_FILENAME_FROM_ALL16S",
            "ab1_size_bytes": "",
            "ab1_modified_time": "",
        }

    expected_filename = f"{original_fasta_id}.ab1"
    expected_path = ab1_dir / expected_filename
    try:
        stat_result = expected_path.stat()
    except FileNotFoundError:
        return {
            "ab1_expected_filename": expected_filename,
            "ab1_expected_path": str(expected_path),
            "ab1_present": "FALSE",
            "ab1_status": "MISSING_EXPECTED_FILE",
            "ab1_size_bytes": "",
            "ab1_modified_time": "",
        }
    except Exception as exc:  # noqa: BLE001 - preserve exact access error.
        return {
            "ab1_expected_filename": expected_filename,
            "ab1_expected_path": str(expected_path),
            "ab1_present": "FALSE",
            "ab1_status": f"FILE_CHECK_ERROR:{type(exc).__name__}:{exc}",
            "ab1_size_bytes": "",
            "ab1_modified_time": "",
        }

    return {
        "ab1_expected_filename": expected_filename,
        "ab1_expected_path": str(expected_path),
        "ab1_present": "TRUE",
        "ab1_status": "PRESENT",
        "ab1_size_bytes": str(stat_result.st_size),
        "ab1_modified_time": datetime.fromtimestamp(stat_result.st_mtime).isoformat(timespec="seconds"),
    }


def build_baseline(
    all16: pd.DataFrame,
    overview: pd.DataFrame,
    ab1_dir: Path,
    directory_status: str,
    directory_error: str,
) -> pd.DataFrame:
    qub_overview = overview[overview["Partner_ID"].eq("QUB")].copy()
    raw_rows = all16[all16["Project"].eq("Rumen_Gateway")].to_dict("records")

    output_rows: list[dict[str, Any]] = []
    for baseline_index, (_, summary_row) in enumerate(qub_overview.iterrows(), start=1):
        candidates: list[tuple[int, list[str], dict[str, Any]]] = []
        for raw_row in raw_rows:
            score, reasons = score_raw_candidate(raw_row, summary_row)
            if score > 0:
                candidates.append((score, reasons, raw_row))
        candidates.sort(key=lambda item: (-item[0], int(item[2]["_excel_row"])))
        best = candidates[0] if candidates else None
        tie_count = sum(1 for candidate in candidates if best and candidate[0] == best[0])
        status = match_status(best, tie_count)
        best_raw = best[2] if best else None

        partner_isolate_id = text(summary_row["Partner Rumen Gateway isolate ID"])
        output: dict[str, Any] = {
            "baseline_row": baseline_index,
            "match_status": status,
            "match_score": str(best[0]) if best else "",
            "match_reasons": ",".join(best[1]) if best else "",
            "raw_candidate_count": len(candidates),
            "raw_tie_count": tie_count,
            "selected_match_note": selected_match_note(status),
            "alternate_raw_candidates": "; ".join(raw_candidate_label(candidate) for candidate in candidates[1:]),
            "partner_id": text(summary_row["Partner_ID"]),
            "partner_isolate_id": partner_isolate_id,
            "partner_isolate_id_status": partner_isolate_status(partner_isolate_id),
            "sh_id": partner_isolate_id if partner_isolate_status(partner_isolate_id) == "SH_ID" else "",
            "qub_id": text(summary_row["RG_Partner Rumen Gateway isolate ID"]),
            "sequence_identifier": text(summary_row["Sequence_identifier"]),
            "sequence_primer_from_identifier": parse_sequence_primer(text(summary_row["Sequence_identifier"])),
            "sequence_bp": text(summary_row["Sequence_bp"]),
            "overview_row": text(summary_row["_excel_row"]),
            "overview_organism_type": text(summary_row["Organism_Type"]),
            "overview_gtdb_assignment": text(summary_row["GTDB_Assignment"]),
            "overview_gtdb_taxon": text(summary_row["GTDB_Taxon"]),
            "overview_gtdb_pident": text(summary_row["GTDB_PIDENT"]),
            "overview_closest_hungate_hit": text(summary_row["Closest_Hungate_Hit"]),
            "overview_hungate_hit_pident": text(summary_row["Hungate_Hit_PIDENT"]),
            "overview_jgi_sequenced": text(summary_row["JGI_Sequenced"]),
            "overview_mwl": text(summary_row["MWL?"]),
            "overview_metadata": text(summary_row["Metadata?"]),
            "overview_notes": text(summary_row["Notes"]),
        }
        add_raw_fields(output, best_raw)
        output.update(ab1_status_for(best_raw, ab1_dir))
        output["ab1_directory_listing_status"] = directory_status
        output["ab1_directory_listing_error"] = directory_error
        output_rows.append(output)

    return pd.DataFrame(output_rows, columns=BASELINE_COLUMNS)


def write_tsv(path: Path, rows: pd.DataFrame) -> None:
    rows.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)


def build_summary(baseline: pd.DataFrame, listed_ab1_files: set[str], directory_status: str, directory_error: str) -> pd.DataFrame:
    expected_files = set(baseline.loc[baseline["ab1_expected_filename"].astype(bool), "ab1_expected_filename"])
    metrics: list[tuple[str, Any]] = [
        ("overview_qub_assigned_sequence_rows", len(baseline)),
        ("unique_partner_isolate_ids", baseline["partner_isolate_id"].nunique()),
        ("unique_sh_like_ids", baseline.loc[baseline["partner_isolate_id_status"].eq("SH_ID"), "partner_isolate_id"].nunique()),
        ("non_sh_partner_isolate_rows", int(baseline["partner_isolate_id_status"].eq("NON_SH_ID").sum())),
        ("unique_qub_ids", baseline["qub_id"].nunique()),
        ("unique_sequence_identifiers", baseline["sequence_identifier"].nunique()),
        ("ab1_expected_files", len(expected_files)),
        ("ab1_present_rows", int(baseline["ab1_present"].eq("TRUE").sum())),
        ("ab1_missing_rows", int(baseline["ab1_status"].eq("MISSING_EXPECTED_FILE").sum())),
        ("ab1_directory_listing_status", directory_status),
        ("ab1_directory_listing_error", directory_error),
        ("listed_ab1_files", len(listed_ab1_files) if directory_status == "AVAILABLE" else "UNAVAILABLE"),
        (
            "listed_ab1_files_not_assigned_to_qub",
            len(listed_ab1_files - expected_files) if directory_status == "AVAILABLE" else "UNAVAILABLE",
        ),
    ]
    for status, count in baseline["match_status"].value_counts().sort_index().items():
        metrics.append((f"match_status__{status}", int(count)))
    for status, count in baseline["ab1_status"].value_counts().sort_index().items():
        metrics.append((f"ab1_status__{status}", int(count)))
    return pd.DataFrame(metrics, columns=["metric", "value"])


def build_extra_files(
    baseline: pd.DataFrame,
    listed_ab1_files: set[str],
    ab1_dir: Path,
    directory_status: str,
    directory_error: str,
) -> pd.DataFrame:
    if directory_status != "AVAILABLE":
        return pd.DataFrame(
            [
                {
                    "filename": "",
                    "path": "",
                    "status": "DIRECTORY_LISTING_UNAVAILABLE",
                    "detail": directory_error,
                }
            ]
        )
    expected_files = set(baseline.loc[baseline["ab1_expected_filename"].astype(bool), "ab1_expected_filename"])
    extra_files = sorted(listed_ab1_files - expected_files)
    return pd.DataFrame(
        [
            {
                "filename": filename,
                "path": str(ab1_dir / filename),
                "status": "NOT_ASSIGNED_TO_QUB_BASELINE",
                "detail": "",
            }
            for filename in extra_files
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all16", type=Path, default=DEFAULT_ALL16)
    parser.add_argument("--overview", type=Path, default=DEFAULT_OVERVIEW)
    parser.add_argument("--ab1-dir", type=Path, default=DEFAULT_AB1_DIR)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    all16, overview = read_sources(args.all16, args.overview)
    directory_status, directory_error, listed_ab1_files = inspect_ab1_directory(args.ab1_dir)
    baseline = build_baseline(all16, overview, args.ab1_dir, directory_status, directory_error)

    baseline_path = args.output_dir / "qub_rumen_gateway_ab1_baseline.tsv"
    missing_path = args.output_dir / "qub_rumen_gateway_missing_ab1_files.tsv"
    ambiguous_path = args.output_dir / "qub_rumen_gateway_ambiguous_raw_matches.tsv"
    extra_path = args.output_dir / "qub_rumen_gateway_extra_ab1_files.tsv"
    summary_path = args.output_dir / "qub_rumen_gateway_ab1_audit_summary.tsv"

    write_tsv(baseline_path, baseline)
    write_tsv(missing_path, baseline[baseline["ab1_status"].eq("MISSING_EXPECTED_FILE")])
    write_tsv(
        ambiguous_path,
        baseline[baseline["match_status"].eq("AMBIGUOUS_RAW_MATCH_SELECTED_EARLIEST_ROW")],
    )
    write_tsv(extra_path, build_extra_files(baseline, listed_ab1_files, args.ab1_dir, directory_status, directory_error))
    write_tsv(summary_path, build_summary(baseline, listed_ab1_files, directory_status, directory_error))

    print(f"Wrote {baseline_path}")
    print(f"Wrote {missing_path}")
    print(f"Wrote {ambiguous_path}")
    print(f"Wrote {extra_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
