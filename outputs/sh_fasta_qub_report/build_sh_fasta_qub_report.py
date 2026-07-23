#!/usr/bin/env python3
"""Report Original_fasta_ID and QUB IDs for each Rumen Gateway all_16S SH_ID."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


DEFAULT_ALL16_PATH = Path("/Users/nicholas/Downloads/all_16S_QUB_formatted_14_12_2025_live_streamlined (2).xlsx")
DEFAULT_OVERVIEW_PATH = Path(
    "/Users/nicholas/Library/CloudStorage/OneDrive-SharedLibraries-Queen'sUniversityBelfast/"
    "RUMEN Gateway QUB-Scientific_Consortium - Documents/03_Partners_sequences/"
    "Current_Sequences_Overview.xlsx"
)
RUMEN_GATEWAY_PROJECT = "Rumen_Gateway"

SUMMARY_COLUMNS = [
    "sh_id",
    "status",
    "all16s_row_count",
    "all16s_rows",
    "original_fasta_id_count",
    "original_fasta_ids",
    "original_fasta_ids_raw",
    "qub_id_count",
    "qub_ids",
    "sequence_identifier_count",
    "sequence_identifiers",
    "overview_row_count",
    "overview_rows",
    "overview_partner_ids",
    "notes",
]

LONG_COLUMNS = [
    "sh_id",
    "all16s_row",
    "all16s_db_control_seq",
    "all16s_project",
    "all16s_bc_status",
    "all16s_short_sequence_code",
    "original_fasta_id",
    "original_fasta_id_raw",
    "all16s_project_isolate_id",
    "all16s_order_id",
    "all16s_sequence_order_id",
    "all16s_primer_direction",
    "all16s_sequencing_primer",
    "all16s_base_pairs",
    "qub_ids_for_sh",
    "sequence_identifiers_for_sh",
    "overview_rows_for_sh",
    "status",
]

UNMATCHED_OVERVIEW_COLUMNS = [
    "sh_id",
    "overview_row_count",
    "overview_rows",
    "qub_id_count",
    "qub_ids",
    "sequence_identifier_count",
    "sequence_identifiers",
    "status",
]


def text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def clean_identifier(value: Any) -> str:
    value_text = text(value)
    return value_text[1:] if value_text.startswith(">") else value_text


def is_emptyish(value: Any) -> bool:
    return text(value).upper() in {"", "NA", "NR", "NAN", "NONE"}


def is_sh_id(value: Any) -> bool:
    return bool(re.match(r"^SH_\d", text(value), flags=re.IGNORECASE))


def unique_preserve_order(values: Iterable[Any], *, clean: bool = False, drop_emptyish: bool = True) -> list[str]:
    seen: set[str] = set()
    output: list[str] = []
    for value in values:
        value_text = clean_identifier(value) if clean else text(value)
        if drop_emptyish and is_emptyish(value_text):
            continue
        if value_text not in seen:
            seen.add(value_text)
            output.append(value_text)
    return output


def join_values(values: Iterable[Any], *, clean: bool = False, drop_emptyish: bool = True) -> str:
    return "; ".join(unique_preserve_order(values, clean=clean, drop_emptyish=drop_emptyish))


def read_workbooks(all16_path: Path, overview_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    all16 = pd.read_excel(all16_path, sheet_name="all_16S_QUB RG", header=1, dtype=str).fillna("")
    all16["_excel_row"] = all16.index + 3

    overview = pd.read_excel(overview_path, sheet_name="Summary", dtype=str).fillna("")
    overview["_excel_row"] = overview.index + 2
    return all16, overview


def validate_columns(all16: pd.DataFrame, overview: pd.DataFrame) -> None:
    required_all16 = {
        "Original_fasta_ID",
        "SH_ID",
        "Extra DB control seq",
        "Project",
        "BC_Status",
        "Short_sequence_code",
    }
    required_overview = {
        "Partner_ID",
        "Partner Rumen Gateway isolate ID",
        "RG_Partner Rumen Gateway isolate ID",
        "Sequence_identifier",
    }
    missing_all16 = sorted(required_all16 - set(all16.columns))
    missing_overview = sorted(required_overview - set(overview.columns))
    if missing_all16 or missing_overview:
        raise ValueError(
            "Missing required columns. "
            f"all_16S missing={missing_all16}; Current_Sequences_Overview missing={missing_overview}"
        )


def overview_by_sh(overview: pd.DataFrame) -> dict[str, pd.DataFrame]:
    overview_sh = overview[
        overview["Partner Rumen Gateway isolate ID"].map(is_sh_id)
        & overview["RG_Partner Rumen Gateway isolate ID"].str.match(r"^QUB_\d+", na=False)
    ].copy()
    return {sh_id: group.copy() for sh_id, group in overview_sh.groupby("Partner Rumen Gateway isolate ID", sort=True)}


def status_for(fasta_ids: list[str], qub_ids: list[str]) -> str:
    if fasta_ids and qub_ids:
        return "OK"
    if fasta_ids and not qub_ids:
        return "MISSING_QUB_ID"
    if qub_ids and not fasta_ids:
        return "MISSING_ORIGINAL_FASTA_ID"
    return "MISSING_FASTA_AND_QUB_ID"


def rumen_gateway_sh_rows(all16: pd.DataFrame) -> pd.DataFrame:
    """Return all_16S rows whose column B project is Rumen_Gateway and whose column H is an SH ID."""
    return all16[all16["Project"].eq(RUMEN_GATEWAY_PROJECT) & all16["SH_ID"].map(is_sh_id)].copy()


def build_summary(all16: pd.DataFrame, overview_lookup: dict[str, pd.DataFrame]) -> pd.DataFrame:
    all16_sh = rumen_gateway_sh_rows(all16)
    rows: list[dict[str, Any]] = []

    for sh_id, group in all16_sh.groupby("SH_ID", sort=True):
        overview_rows = overview_lookup.get(sh_id, pd.DataFrame())
        fasta_ids = unique_preserve_order(group["Original_fasta_ID"], clean=True)
        fasta_ids_raw = unique_preserve_order(group["Original_fasta_ID"], clean=False)
        qub_ids = unique_preserve_order(
            overview_rows.get("RG_Partner Rumen Gateway isolate ID", pd.Series(dtype=str)),
            clean=False,
        )
        sequence_identifiers = unique_preserve_order(
            overview_rows.get("Sequence_identifier", pd.Series(dtype=str)),
            clean=False,
        )
        status = status_for(fasta_ids, qub_ids)
        rows.append(
            {
                "sh_id": sh_id,
                "status": status,
                "all16s_row_count": len(group),
                "all16s_rows": join_values(group["_excel_row"], drop_emptyish=False),
                "original_fasta_id_count": len(fasta_ids),
                "original_fasta_ids": "; ".join(fasta_ids),
                "original_fasta_ids_raw": "; ".join(fasta_ids_raw),
                "qub_id_count": len(qub_ids),
                "qub_ids": "; ".join(qub_ids),
                "sequence_identifier_count": len(sequence_identifiers),
                "sequence_identifiers": "; ".join(sequence_identifiers),
                "overview_row_count": len(overview_rows),
                "overview_rows": join_values(overview_rows.get("_excel_row", pd.Series(dtype=str)), drop_emptyish=False),
                "overview_partner_ids": join_values(overview_rows.get("Partner_ID", pd.Series(dtype=str))),
                "notes": "" if status == "OK" else "Check source workbook coverage for this SH_ID.",
            }
        )
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def build_long(all16: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    all16_sh = rumen_gateway_sh_rows(all16)
    summary_lookup = summary.set_index("sh_id").to_dict("index")
    rows: list[dict[str, Any]] = []

    for _, source_row in all16_sh.iterrows():
        sh_id = text(source_row["SH_ID"])
        summary_row = summary_lookup[sh_id]
        rows.append(
            {
                "sh_id": sh_id,
                "all16s_row": text(source_row["_excel_row"]),
                "all16s_db_control_seq": text(source_row.get("Extra DB control seq", "")),
                "all16s_project": text(source_row.get("Project", "")),
                "all16s_bc_status": text(source_row.get("BC_Status", "")),
                "all16s_short_sequence_code": clean_identifier(source_row.get("Short_sequence_code", "")),
                "original_fasta_id": clean_identifier(source_row.get("Original_fasta_ID", "")),
                "original_fasta_id_raw": text(source_row.get("Original_fasta_ID", "")),
                "all16s_project_isolate_id": text(source_row.get("Project_isolate_ID", "")),
                "all16s_order_id": text(source_row.get("Order_ID", "")),
                "all16s_sequence_order_id": text(source_row.get("Sequence_order_ID", "")),
                "all16s_primer_direction": text(source_row.get("Primer_direction", "")),
                "all16s_sequencing_primer": text(source_row.get("Sequencing_primer", "")),
                "all16s_base_pairs": text(source_row.get("Base_pairs", "")),
                "qub_ids_for_sh": summary_row["qub_ids"],
                "sequence_identifiers_for_sh": summary_row["sequence_identifiers"],
                "overview_rows_for_sh": summary_row["overview_rows"],
                "status": summary_row["status"],
            }
        )
    return pd.DataFrame(rows, columns=LONG_COLUMNS)


def build_unmatched_overview_sh(overview_lookup: dict[str, pd.DataFrame], summary: pd.DataFrame) -> pd.DataFrame:
    all16_sh_ids = set(summary["sh_id"])
    rows: list[dict[str, Any]] = []
    for sh_id, group in overview_lookup.items():
        if sh_id in all16_sh_ids:
            continue
        rows.append(
            {
                "sh_id": sh_id,
                "overview_row_count": len(group),
                "overview_rows": join_values(group["_excel_row"], drop_emptyish=False),
                "qub_id_count": len(unique_preserve_order(group["RG_Partner Rumen Gateway isolate ID"])),
                "qub_ids": join_values(group["RG_Partner Rumen Gateway isolate ID"]),
                "sequence_identifier_count": len(unique_preserve_order(group["Sequence_identifier"])),
                "sequence_identifiers": join_values(group["Sequence_identifier"]),
                "status": "OVERVIEW_QUB_SH_NOT_FOUND_IN_RUMEN_GATEWAY_ALL16S",
            }
        )
    return pd.DataFrame(rows, columns=UNMATCHED_OVERVIEW_COLUMNS)


def build_audit(summary: pd.DataFrame, long_rows: pd.DataFrame, unmatched_overview_sh: pd.DataFrame) -> pd.DataFrame:
    metrics: list[tuple[str, Any]] = [
        ("all16s_project_filter", RUMEN_GATEWAY_PROJECT),
        ("unique_rumen_gateway_sh_ids_in_all16s", len(summary)),
        ("rumen_gateway_all16s_rows_with_valid_sh_id", len(long_rows)),
        ("unique_sh_ids_with_original_fasta_id", int((summary["original_fasta_id_count"].astype(int) > 0).sum())),
        ("unique_sh_ids_without_original_fasta_id", int((summary["original_fasta_id_count"].astype(int) == 0).sum())),
        ("unique_sh_ids_with_qub_id", int((summary["qub_id_count"].astype(int) > 0).sum())),
        ("unique_sh_ids_without_qub_id", int((summary["qub_id_count"].astype(int) == 0).sum())),
        ("unique_sh_ids_missing_qub_id", int(summary["status"].eq("MISSING_QUB_ID").sum())),
        ("unique_sh_ids_missing_original_fasta_id", int(summary["status"].eq("MISSING_ORIGINAL_FASTA_ID").sum())),
        ("overview_qub_sh_ids_not_found_in_rumen_gateway_all16s", len(unmatched_overview_sh)),
        ("total_original_fasta_ids_unique_by_sh", int(summary["original_fasta_id_count"].astype(int).sum())),
        ("total_qub_ids_unique_by_sh", int(summary["qub_id_count"].astype(int).sum())),
    ]
    for status, count in Counter(summary["status"]).items():
        metrics.append((f"summary_status__{status}", count))
    return pd.DataFrame(metrics, columns=["metric", "value"])


def write_tsv(path: Path, rows: pd.DataFrame) -> None:
    rows.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all16", type=Path, default=DEFAULT_ALL16_PATH)
    parser.add_argument("--overview", type=Path, default=DEFAULT_OVERVIEW_PATH)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    all16, overview = read_workbooks(args.all16, args.overview)
    validate_columns(all16, overview)

    overview_lookup = overview_by_sh(overview)
    summary = build_summary(all16, overview_lookup)
    long_rows = build_long(all16, summary)
    unmatched_overview_sh = build_unmatched_overview_sh(overview_lookup, summary)
    audit = build_audit(summary, long_rows, unmatched_overview_sh)

    write_tsv(args.output_dir / "sh_fasta_qub_summary.tsv", summary)
    write_tsv(args.output_dir / "sh_fasta_qub_long.tsv", long_rows)
    write_tsv(args.output_dir / "sh_missing_qub_ids.tsv", summary[summary["status"].eq("MISSING_QUB_ID")])
    write_tsv(args.output_dir / "overview_sh_not_found_in_all16s.tsv", unmatched_overview_sh)
    write_tsv(args.output_dir / "sh_fasta_qub_audit_summary.tsv", audit)

    print(f"Wrote {args.output_dir / 'sh_fasta_qub_summary.tsv'}")
    print(f"Wrote {args.output_dir / 'sh_fasta_qub_long.tsv'}")
    print(f"Wrote {args.output_dir / 'sh_missing_qub_ids.tsv'}")
    print(f"Wrote {args.output_dir / 'overview_sh_not_found_in_all16s.tsv'}")
    print(f"Wrote {args.output_dir / 'sh_fasta_qub_audit_summary.tsv'}")


if __name__ == "__main__":
    main()
