#!/usr/bin/env python3
"""Generate deterministic synthetic previews of PhenGO publication outputs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


ORGANISMS = ["fish", "fly", "mouse", "worm", "yeast"]
MODELS = ["dt", "rf", "gb", "lr", "svm", "nn"]
MODEL_LABELS = {
    "dt": "Decision tree",
    "rf": "Random forest",
    "gb": "Gradient boost",
    "lr": "Logistic regression",
    "svm": "SVM",
    "nn": "MLP",
}
PRIMARY_YEARS = {
    organism: list(range(2015, 2026)) for organism in ORGANISMS[:-1]
}
PRIMARY_YEARS["yeast"] = [2015, 2016, 2017, 2024, 2025]
ORG_COLORS = {
    "fish": "#0072B2",
    "fly": "#D55E00",
    "mouse": "#009E73",
    "worm": "#CC79A7",
    "yeast": "#6B5B95",
}
MODEL_COLORS = {
    "dt": "#6C757D",
    "rf": "#0072B2",
    "gb": "#E69F00",
    "lr": "#009E73",
    "svm": "#D55E00",
    "nn": "#CC79A7",
}
TRACK_COLORS = {
    "primary": "#D55E00",
    "fixed_2025": "#0072B2",
    "is_a_only": "#CC79A7",
    "fixed_2025_is_a": "#009E73",
}
MOCK_STATUS = "synthetic_mock_not_results"


def _primary_years(organism):
    return PRIMARY_YEARS[organism]


def _mock_banner(fig):
    fig.text(
        0.5,
        0.955,
        "SYNTHETIC MOCK - NOT RESULTS",
        ha="center",
        va="top",
        color="#9C0006",
        fontsize=12,
        fontweight="bold",
    )


def _save(fig, path):
    _mock_banner(fig)
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D9E0E5", linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)


def _write_tsv(frame, path):
    frame = frame.copy()
    frame.insert(0, "data_status", MOCK_STATUS)
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
        float_format="%.4f",
    )


def _resource_counts(organism, year):
    oi = ORGANISMS.index(organism)
    elapsed = year - 2010
    gene_base = [5400, 10300, 7200, 13500, 4700][oi]
    feature_base = [4200, 5100, 5900, 4800, 3900][oi]
    genes = int(gene_base + 125 * elapsed + 260 * np.sin(0.61 * elapsed + oi))
    features = int(feature_base + 235 * elapsed + 310 * np.cos(0.47 * elapsed + oi))
    annotations = int(genes * (15.5 + 0.8 * elapsed + oi * 1.7))
    prevalence = np.clip(
        [0.18, 0.24, 0.31, 0.21, 0.16][oi]
        + 0.018 * np.sin(0.78 * elapsed + oi),
        0.08,
        0.45,
    )
    return genes, features, annotations, float(prevalence)


def _performance(organism, year, model, panel="matched_both", track="primary"):
    oi = ORGANISMS.index(organism)
    mi = MODELS.index(model)
    elapsed = year - 2015
    base = [0.62, 0.68, 0.72, 0.60, 0.75][oi]
    model_effect = [-0.055, 0.028, 0.018, 0.0, 0.012, 0.022][mi]
    release_wave = 0.045 * np.sin(0.72 * elapsed + oi * 0.9 + mi * 0.27)
    release_event = 0.035 if year in {2018, 2021, 2024} else 0.0
    panel_effect = 0.018 if panel == "matched_both" else -0.004
    track_effect = {
        "primary": release_wave + release_event,
        "is_a_only": release_wave + release_event - 0.008,
        "fixed_2025": 0.45 * release_wave,
        "fixed_2025_is_a": 0.55 * release_wave - 0.008,
    }[track]
    return float(np.clip(base + model_effect + panel_effect + track_effect, 0.34, 0.93))


def build_tables():
    snapshot_rows = []
    track_years = {
        "primary": lambda org: _primary_years(org),
        "is_a_only": lambda org: _primary_years(org),
        "no_iea": lambda org: _primary_years(org),
        "fixed_2025": lambda org: list(range(2014, 2026)),
        "fixed_2025_is_a": lambda org: list(range(2010, 2026)),
    }
    for track, year_fn in track_years.items():
        for organism in ORGANISMS:
            for year in year_fn(organism):
                genes, features, _, prevalence = _resource_counts(organism, year)
                if track == "no_iea":
                    features = int(features * 0.72)
                elif track in {"is_a_only", "fixed_2025_is_a"}:
                    features = int(features * 0.91)
                snapshot_rows.append({
                    "track": track,
                    "organism": organism,
                    "timepoint": year,
                    "snapshot_id": f"{organism}-{year}-{track}",
                    "strict_snapshot": True,
                    "genes": genes,
                    "lethal": int(genes * prevalence),
                    "viable": genes - int(genes * prevalence),
                    "go_features": features,
                    "label_source": "release_records" if track in {
                        "primary", "is_a_only", "no_iea"
                    } else "fixed_2025_gene_sets",
                    "go_relations": "is_a" if track in {
                        "is_a_only", "fixed_2025_is_a"
                    } else "is_a;part_of",
                    "go_evidence_exclude": "IEA" if track == "no_iea" else "",
                })

    drift_rows = []
    for organism in ORGANISMS:
        previous = None
        for year in range(2010, 2026):
            genes, features, annotations, prevalence = _resource_counts(organism, year)
            gene_jaccard = "" if previous is None else 0.96 - 0.012 * abs(
                np.sin(year + ORGANISMS.index(organism))
            )
            drift_rows.append({
                "track": "fixed_2025_is_a",
                "organism": organism,
                "dataset": f"{organism}-{year}",
                "year": year,
                "n_genes": genes,
                "n_features": features,
                "gaf_annotation_rows": annotations,
                "positive_rate": prevalence,
                "density": min(0.16, annotations / max(1, genes * features)),
                "previous_year_gene_jaccard": gene_jaccard,
            })
            previous = year

    within_rows = []
    for organism in ORGANISMS:
        for year in _primary_years(organism):
            for panel in ["full", "matched_both"]:
                for model in MODELS:
                    ap = _performance(organism, year, model, panel)
                    ci = 0.026 + 0.004 * ((year + MODELS.index(model)) % 3)
                    within_rows.append({
                        "track": "primary",
                        "organism": organism,
                        "panel": panel,
                        "dataset": f"{organism}-{year}",
                        "year": year,
                        "model": model,
                        "n_folds": 5,
                        "n_repeats": 5,
                        "average_precision_mean": ap,
                        "average_precision_sd": ci / 2.2,
                        "average_precision_ci95_low": max(0.0, ap - ci),
                        "average_precision_ci95_high": min(1.0, ap + ci),
                        "mcc_mean": np.clip(ap - 0.285, 0.05, 0.78),
                        "balanced_accuracy_mean": np.clip(0.5 + (ap - 0.35) * 0.53, 0.5, 0.9),
                        "lethal_f1_mean": np.clip(ap - 0.095, 0.1, 0.88),
                        "fit_warnings": 0 if model != "nn" or year % 5 else 1,
                    })

    transfer_rows = []
    for organism in ORGANISMS:
        years = _primary_years(organism)
        for panel in ["full", "matched_both"]:
            for model in MODELS:
                for train_year in years:
                    diagonal = _performance(organism, train_year, model, panel)
                    for test_year in years:
                        distance = abs(test_year - train_year)
                        direction = 0.006 if test_year > train_year else -0.004
                        penalty = distance * (0.009 if panel == "matched_both" else 0.013)
                        ap = np.clip(
                            diagonal - penalty + direction * np.sign(test_year - train_year),
                            0.25,
                            0.93,
                        )
                        transfer_rows.append({
                            "track": "primary",
                            "organism": organism,
                            "panel": panel,
                            "model": model,
                            "train_dataset": f"{organism}-{train_year}",
                            "test_dataset": f"{organism}-{test_year}",
                            "train_year": train_year,
                            "test_year": test_year,
                            "average_precision": ap,
                            "mcc": np.clip(ap - 0.30, 0.01, 0.75),
                            "balanced_accuracy": np.clip(0.5 + (ap - 0.35) * 0.5, 0.5, 0.9),
                        })

    decomposition_rows = []
    for organism in ORGANISMS:
        for year in _primary_years(organism):
            for track in ["primary", "fixed_2025", "is_a_only", "fixed_2025_is_a"]:
                aps = [_performance(organism, year, model, "matched_both", track) for model in MODELS]
                decomposition_rows.append({
                    "organism": organism,
                    "year": year,
                    "track": track,
                    "model_consensus_average_precision": float(np.mean(aps)),
                    "between_model_sd": float(np.std(aps, ddof=1)),
                })

    instability_rows = []
    for organism in ORGANISMS:
        oi = ORGANISMS.index(organism)
        for model in MODELS:
            mi = MODELS.index(model)
            instability_rows.append({
                "track": "primary",
                "organism": organism,
                "panel": "matched_both",
                "model": model,
                "genes_scored_in_multiple_years": 2200 + oi * 510,
                "median_probability_range": 0.16 + 0.021 * mi + 0.012 * oi,
                "genes_with_prediction_flip_fraction": 0.09 + 0.018 * mi + 0.011 * oi,
                "median_prediction_flip_count": 1 + int((mi + oi) / 4),
            })

    term_defs = [
        ("GO:0007049", "cell cycle"),
        ("GO:0006281", "DNA repair"),
        ("GO:0006412", "translation"),
        ("GO:0006915", "apoptotic process"),
        ("GO:0006096", "glycolytic process"),
        ("GO:0006629", "lipid metabolic process"),
        ("GO:0006355", "regulation of transcription"),
        ("GO:0006974", "response to DNA damage"),
    ]
    temporal_rows = []
    for ti, (go_id, name) in enumerate(term_defs):
        for year in range(2015, 2026):
            signal = max(0.0, 1.2 + 1.5 * np.sin((year - 2015) * 0.57 + ti * 0.83))
            fdr = min(1.0, 10 ** (-signal))
            temporal_rows.append({
                "track": "primary",
                "organism": "worm",
                "year": year,
                "go_id": go_id,
                "go_name": name,
                "odds_ratio": 1.0 + 0.7 * signal,
                "fdr": fdr,
                "minus_log10_fdr": signal,
                "lethal_gene_count": int(12 + signal * 7 + ti),
            })

    snapshot = pd.DataFrame(snapshot_rows)
    drift = pd.DataFrame(drift_rows)
    within = pd.DataFrame(within_rows)
    transfer = pd.DataFrame(transfer_rows)
    decomposition = pd.DataFrame(decomposition_rows)
    instability = pd.DataFrame(instability_rows)
    temporal = pd.DataFrame(temporal_rows)

    claim_rows = []
    for organism in ORGANISMS:
        primary = decomposition[
            (decomposition.organism == organism) & (decomposition.track == "primary")
        ]
        fixed = decomposition[
            (decomposition.organism == organism) & (decomposition.track == "fixed_2025")
        ]
        tx = transfer[
            (transfer.organism == organism)
            & (transfer.panel == "matched_both")
            & (transfer.model == "lr")
        ].copy()
        tx["distance"] = (tx.train_year - tx.test_year).abs()
        claim_rows.append({
            "organism": organism,
            "primary_within_year_ap_range": primary.model_consensus_average_precision.max()
            - primary.model_consensus_average_precision.min(),
            "fixed_label_ap_range": fixed.model_consensus_average_precision.max()
            - fixed.model_consensus_average_precision.min(),
            "lr_mean_far_transfer_penalty": tx[tx.distance >= 4].average_precision.mean()
            - tx[tx.distance == 0].average_precision.mean(),
            "minimum_consecutive_gene_jaccard": pd.to_numeric(
                drift[drift.organism == organism].previous_year_gene_jaccard,
                errors="coerce",
            ).min(),
            "interpretive_role": "Illustrates the planned effect-size summary; not evidence.",
        })

    return {
        "snapshot_inventory.tsv": snapshot,
        "dataset_drift.tsv": drift,
        "within_year_model_performance.tsv": within,
        "cross_year_transfer_performance.tsv": transfer,
        "prediction_instability_summary.tsv": instability,
        "track_decomposition.tsv": decomposition,
        "temporal_enrichment_statistics.tsv": temporal,
        "manuscript_claim_summary.tsv": pd.DataFrame(claim_rows),
    }


def figure_design(output):
    years = list(range(2010, 2026))
    tracks = ["primary", "is_a_only", "no_iea", "fixed_2025", "fixed_2025_is_a"]
    fig, axes = plt.subplots(1, 5, figsize=(18, 4.9), sharey=True)
    for ax, organism in zip(axes, ORGANISMS):
        matrix = np.zeros((len(tracks), len(years)))
        for ri, track in enumerate(tracks):
            available = (
                _primary_years(organism)
                if track in {"primary", "is_a_only", "no_iea"}
                else list(range(2014, 2026)) if track == "fixed_2025"
                else years
            )
            for year in available:
                matrix[ri, years.index(year)] = 1
        cmap = LinearSegmentedColormap.from_list("coverage", ["#ECEFF1", ORG_COLORS[organism]])
        ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
        ax.set_title(organism.capitalize(), fontweight="bold")
        ax.set_xticks(range(0, len(years), 3), [years[i] for i in range(0, len(years), 3)], rotation=45)
        ax.set_yticks(range(len(tracks)), tracks)
        ax.set_xlabel("Resource year")
        for spine in ax.spines.values():
            spine.set_color("#B0BEC5")
    fig.suptitle("Figure 1. Planned experimental coverage", y=1.04, fontsize=15, fontweight="bold")
    fig.subplots_adjust(top=0.80, bottom=0.20, wspace=0.20)
    _save(fig, output / "figure_1_experimental_design.png")


def figure_drift(drift, output):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    specs = [
        ("n_genes", "Genes represented"),
        ("n_features", "GO features retained"),
        ("positive_rate", "Lethal-class prevalence"),
    ]
    for ax, (column, label) in zip(axes, specs):
        for organism in ORGANISMS:
            part = drift[drift.organism == organism]
            ax.plot(part.year, part[column], marker="o", markersize=2.8, linewidth=1.5,
                    color=ORG_COLORS[organism], label=organism)
        ax.set_xlabel("GO resource year")
        ax.set_ylabel(label)
        ax.set_xticks([2010, 2013, 2016, 2019, 2022, 2025])
        _style_axis(ax)
    axes[0].legend(frameon=False, ncol=1, fontsize=8)
    fig.suptitle("Figure 2. Resource composition changes across archived releases", y=1.04,
                 fontsize=15, fontweight="bold")
    _save(fig, output / "figure_2_resource_drift.png")


def figure_within_year(within, output):
    data = within[within.panel == "matched_both"]
    fig, axes = plt.subplots(1, 5, figsize=(19, 4.8), sharey=True)
    for ax, organism in zip(axes, ORGANISMS):
        part = data[data.organism == organism]
        for model in MODELS:
            model_part = part[part.model == model]
            ax.plot(model_part.year, model_part.average_precision_mean, marker="o", markersize=2.5,
                    linewidth=1.35, color=MODEL_COLORS[model], label=MODEL_LABELS[model])
            ax.fill_between(
                model_part.year,
                model_part.average_precision_ci95_low,
                model_part.average_precision_ci95_high,
                color=MODEL_COLORS[model],
                alpha=0.045,
                linewidth=0,
            )
        ax.set_title(organism.capitalize(), fontweight="bold")
        ax.set_xlabel("Phenotype/GO year")
        ax.set_ylim(0.42, 0.9)
        ax.tick_params(axis="x", rotation=45)
        _style_axis(ax)
    axes[0].set_ylabel("Average precision (matched genes + terms)")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, -0.09))
    fig.suptitle("Figure 3. Within-year model performance is release-sensitive", y=1.04,
                 fontsize=15, fontweight="bold")
    _save(fig, output / "figure_3_within_year_performance.png")


def _heatmap(ax, matrix, years, title, vmin=0.35, vmax=0.86):
    image = ax.imshow(matrix, cmap="viridis", vmin=vmin, vmax=vmax, aspect="equal", origin="upper")
    ticks = list(range(0, len(years), max(1, len(years) // 5)))
    ax.set_xticks(ticks, [years[i] for i in ticks], rotation=45)
    ax.set_yticks(ticks, [years[i] for i in ticks])
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.set_xlabel("Test year")
    ax.set_ylabel("Train year")
    return image


def figure_transfer(transfer, output):
    organism = "fly"
    years = _primary_years(organism)
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 9.6))
    image = None
    for ax, model in zip(axes.flat, MODELS):
        part = transfer[
            (transfer.organism == organism)
            & (transfer.panel == "matched_both")
            & (transfer.model == model)
        ]
        matrix = part.pivot(index="train_year", columns="test_year", values="average_precision")
        image = _heatmap(ax, matrix.values, years, MODEL_LABELS[model])
    fig.subplots_adjust(top=0.82, right=0.84, wspace=0.24, hspace=0.30)
    colorbar_axis = fig.add_axes([0.875, 0.22, 0.022, 0.56])
    fig.colorbar(image, cax=colorbar_axis, label="Average precision")
    fig.suptitle("Figure 4A. Cross-year transfer for all six models (illustrative fly)",
                 y=1.01, fontsize=15, fontweight="bold")
    _save(fig, output / "figure_4a_transfer_all_models.png")

    fig, axes = plt.subplots(2, 3, figsize=(13, 8.8))
    image = None
    for ax, organism in zip(axes.flat, ORGANISMS):
        years = _primary_years(organism)
        part = transfer[
            (transfer.organism == organism)
            & (transfer.panel == "matched_both")
            & (transfer.model == "lr")
        ]
        matrix = part.pivot(index="train_year", columns="test_year", values="average_precision")
        image = _heatmap(ax, matrix.values, years, organism.capitalize())
    axes.flat[-1].axis("off")
    fig.subplots_adjust(top=0.80, right=0.84, wspace=0.24, hspace=0.34)
    colorbar_axis = fig.add_axes([0.875, 0.22, 0.022, 0.56])
    fig.colorbar(image, cax=colorbar_axis, label="Average precision")
    fig.suptitle("Figure 4B. Directional transfer maps across model organisms",
                 y=1.01, fontsize=15, fontweight="bold")
    fig.text(
        0.5,
        0.895,
        "Illustrative logistic regression; all models are retained in Supplementary Data",
        ha="center",
        fontsize=10,
        color="#455A64",
    )
    _save(fig, output / "figure_4b_transfer_all_organisms.png")


def figure_decomposition(decomposition, output):
    fig, axes = plt.subplots(2, 5, figsize=(18.5, 8.6), sharex="col", sharey=True)
    labels = {
        "primary": "Release labels + GO",
        "fixed_2025": "Fixed labels, is_a + part_of",
        "is_a_only": "Release labels + GO, is_a",
        "fixed_2025_is_a": "Fixed labels, is_a",
    }
    comparisons = [
        ("is_a + part_of", ["primary", "fixed_2025"]),
        ("is_a only", ["is_a_only", "fixed_2025_is_a"]),
    ]
    for row_index, (policy, tracks) in enumerate(comparisons):
        for ax, organism in zip(axes[row_index], ORGANISMS):
            for track in tracks:
                part = decomposition[
                    (decomposition.organism == organism) & (decomposition.track == track)
                ]
                ax.plot(
                    part.year,
                    part.model_consensus_average_precision,
                    marker="o",
                    markersize=2.8,
                    linewidth=1.55,
                    color=TRACK_COLORS[track],
                    label=labels[track],
                )
            if row_index == 0:
                ax.set_title(organism.capitalize(), fontweight="bold")
            if row_index == 1:
                ax.set_xlabel("Year")
            ticks = [
                year for year in [2015, 2017, 2019, 2021, 2023, 2025]
                if year in set(part.year)
            ]
            ax.set_xticks(ticks, [str(year) for year in ticks], rotation=45)
            _style_axis(ax)
        axes[row_index, 0].set_ylabel(
            f"{policy}\nModel-consensus average precision"
        )
    handles = []
    legend_labels = []
    for ax in axes.flat:
        for handle, label in zip(*ax.get_legend_handles_labels()):
            if label not in legend_labels:
                handles.append(handle)
                legend_labels.append(label)
    fig.legend(
        handles,
        legend_labels,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, -0.035),
    )
    fig.suptitle("Figure 5. Decomposing phenotype-label and GO-resource sensitivity", y=1.04,
                 fontsize=15, fontweight="bold")
    _save(fig, output / "figure_5_track_decomposition.png")


def figure_instability(instability, output):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), sharey=True)
    x = np.arange(len(MODELS))
    width = 0.14
    for oi, organism in enumerate(ORGANISMS):
        part = instability[instability.organism == organism].set_index("model").loc[MODELS]
        offset = (oi - 2) * width
        axes[0].bar(x + offset, part.median_probability_range, width, color=ORG_COLORS[organism],
                    label=organism)
        axes[1].bar(x + offset, part.genes_with_prediction_flip_fraction, width,
                    color=ORG_COLORS[organism])
    for ax, ylabel in zip(
        axes,
        ["Median probability range", "Genes with at least one class flip (fraction)"],
    ):
        ax.set_xticks(x, [MODEL_LABELS[model] for model in MODELS], rotation=35, ha="right")
        ax.set_ylabel(ylabel)
        _style_axis(ax)
    axes[0].legend(frameon=False, ncol=1, fontsize=8)
    fig.suptitle("Figure 6. Individual-gene predictions vary across resource years", y=1.04,
                 fontsize=15, fontweight="bold")
    _save(fig, output / "figure_6_prediction_instability.png")


def figure_go_temporal(temporal, output):
    matrix = temporal.pivot(index="go_name", columns="year", values="minus_log10_fdr")
    fig, ax = plt.subplots(figsize=(12.5, 5.8))
    image = ax.imshow(matrix.values, aspect="auto", cmap="magma", vmin=0, vmax=3)
    ax.set_xticks(range(len(matrix.columns)), matrix.columns, rotation=45)
    ax.set_yticks(range(len(matrix.index)), matrix.index)
    ax.set_xlabel("WormBase/GO resource year")
    ax.set_ylabel("GO term")
    fig.colorbar(image, ax=ax, label="-log10(FDR)")
    fig.suptitle("Figure 7. GO associations appear, disappear, or change strength by release",
                 y=1.03, fontsize=15, fontweight="bold")
    _save(fig, output / "figure_7_go_temporal_enrichment.png")


def write_report(output, tables):
    inventory = tables["snapshot_inventory.tsv"]
    table_lines = []
    for track in ["primary", "is_a_only", "no_iea", "fixed_2025", "fixed_2025_is_a"]:
        part = inventory[inventory.track == track]
        table_lines.append(
            f"| `{track}` | {len(part)} | {part.timepoint.min()}-{part.timepoint.max()} | "
            f"{part.go_relations.iloc[0]} |"
        )
    report = f"""# PhenGO paper-output mock

> **SYNTHETIC MOCK - NOT RESULTS.** Every numerical result in this directory was
> generated deterministically for layout and analysis-design review. No value may
> be quoted as a PhenGO finding.

## What this preview is for

The proposed paper should establish resource change first, then show its effects
on ordinary within-year performance, cross-year transfer, individual predictions,
and GO interpretation. The fixed-label analyses separate phenotype-target change
from GO-resource change. All six model families remain visible; no model or year is
selected after inspecting the result.

## Proposed main-paper figures

1. **Experimental design:** the exact organism, year, and track coverage.
2. **Resource drift:** genes, retained GO terms, annotation volume, prevalence,
   Jaccard overlap, and label churn before any ML result is shown.
3. **Within-year performance:** repeated-CV average precision with uncertainty for
   all six comparable models and every organism-year.
4. **Cross-year transfer:** directional train-year/test-year heatmaps. The main text
   can show one pre-declared model plus an organism overview; all 30 organism-model
   matrices belong in the supplement.
5. **Effect decomposition:** release-matched labels versus fixed 2025 labels over
   the same 2015-2025 interval, with the historical is_a series kept separate.
6. **Prediction instability:** per-gene probability ranges and classification flips.
7. **GO interpretation drift:** complete FDR-corrected timelines, including terms
   that disappear rather than only persistent significant terms.

## Planned track coverage

| Track | Mock snapshots | Year span | GO relations |
| --- | ---: | --- | --- |
{chr(10).join(table_lines)}

## Proposed main-paper tables

| Table | Purpose |
| --- | --- |
| `snapshot_inventory.tsv` | Provenance spine: releases, policies, class counts, feature counts, hashes in the real run. |
| `manuscript_claim_summary.tsv` | Organism-level effect sizes: performance range, far-transfer penalty, and resource overlap. |
| `track_decomposition.tsv` | Same-year comparison of release-derived and fixed-label tracks. |
| `prediction_instability_summary.tsv` | Probability spread and decision flips by organism and model. |

The full tables `within_year_model_performance.tsv`,
`cross_year_transfer_performance.tsv`, `dataset_drift.tsv`, and
`temporal_enrichment_statistics.tsv` are better supplied as Supplementary Data,
with selected rows summarized in the manuscript.

## Visual previews

![Experimental design](figure_1_experimental_design.png)

![Resource drift](figure_2_resource_drift.png)

![Within-year performance](figure_3_within_year_performance.png)

![All-model transfer](figure_4a_transfer_all_models.png)

![All-organism transfer](figure_4b_transfer_all_organisms.png)

![Track decomposition](figure_5_track_decomposition.png)

![Prediction instability](figure_6_prediction_instability.png)

![GO temporal enrichment](figure_7_go_temporal_enrichment.png)

## Decisions this mock is meant to settle

- Whether average precision and MCC should be co-primary, or average precision
  primary with MCC as a robustness metric.
- Whether Figure 4A should use logistic regression as a pre-declared transparent
  instrument or present a model-consensus matrix.
- Whether the five organism panels belong together in the main paper or as one
  main example plus organism-specific supplements.
- Whether GO-term timelines should be in the main paper or remain secondary to
  the central performance and transfer results.
- How much of the 2010-2013 `is_a` historical extension belongs in the main result,
  given that the directly policy-matched `is_a + part_of` series begins in 2014.
"""
    (output / "README.md").write_text(report, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="docs/mock_publication_outputs",
        help="Directory for synthetic tables and figure previews",
    )
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)

    tables = build_tables()
    for filename, frame in tables.items():
        _write_tsv(frame, output / filename)

    figure_design(output)
    figure_drift(tables["dataset_drift.tsv"], output)
    figure_within_year(tables["within_year_model_performance.tsv"], output)
    figure_transfer(tables["cross_year_transfer_performance.tsv"], output)
    figure_decomposition(tables["track_decomposition.tsv"], output)
    figure_instability(tables["prediction_instability_summary.tsv"], output)
    figure_go_temporal(tables["temporal_enrichment_statistics.tsv"], output)
    write_report(output, tables)

    manifest = {
        "data_status": MOCK_STATUS,
        "purpose": "layout and analysis-design review only",
        "tables": sorted(tables),
        "figures": sorted(path.name for path in output.glob("figure_*.png")),
    }
    (output / "mock_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output_dir": str(output), **manifest}, indent=2))


if __name__ == "__main__":
    main()
