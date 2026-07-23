#!/usr/bin/env bash
set -Eeuo pipefail

# Unified, one-shot PhenGO publication workflow.
#
# This replaces the two-stage run_publication_pipeline.sh (V1) +
# run_publication_pipeline_v2.sh (V2) system with a single pass:
#   - all 15 publication tracks (V1's 5 GO-relation/label-policy tracks
#     plus V2's 10 alternative-data-provenance tracks) are built against
#     one ledger read, in one output tree
#   - derived per-organism/year label artifacts (fly assignments, worm
#     lethal-term closures) are computed once per unique underlying input
#     and shared across every track that uses that same input, instead of
#     being recomputed once per track namespace
#   - preflight (scientific-import check, test suite, environment/provenance
#     capture, input checksums) runs once, not once per stage
#   - modeling (version_sensitivity, GO_temporal_analysis, reports) runs
#     once across the full 15-track set
#   - one run fingerprint, one checksum manifest, one run.complete marker
#
# Removed relative to V1+V2: the "build/resume/integrity-audit/rebuild V1
# before starting V2" gate, the separate V2 extension output tree, and the
# publication_multiverse.py consolidation step -- there is only one tree now,
# so there is nothing to audit-then-merge.
#
# Architecture: works on both arm64 and x86_64 Python environments. The
# previous hard "must be arm64" checks are gone; the detected architecture
# is recorded in provenance and folded into the run fingerprint instead, so
# a --resume under a different architecture triggers a clean rebuild rather
# than silently mixing numeric results computed under different BLAS/vector
# backends.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_DIR}"

# ---------------------------------------------------------------------------
# Portable helpers (macOS ships `shasum`; most Linux distros ship `sha256sum`)
# ---------------------------------------------------------------------------
sha256() {
    if command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$@"
    elif command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$@"
    else
        printf 'Neither shasum nor sha256sum is available on this system.\n' >&2
        exit 1
    fi
}
sha256_of_stdin() { sha256 | awk '{print $1}'; }

# ---------------------------------------------------------------------------
# Configuration (env-overridable, same PHENGO_* variable names as before
# where a direct equivalent existed)
# ---------------------------------------------------------------------------
PYTHON_DEFAULT="/Users/nicholas/miniconda3/envs/PhenGO/bin/python"
#OUTPUT_DEFAULT="/Users/nicholas/Nextcloud/Current_Work/PhenGO/New_Outputs/publication_run_$(date '+%Y-%m-%d')"
OUTPUT_DEFAULT="/Users/nicholas/Nextcloud/Current_Work/PhenGO/New_Outputs/publication_run_2026-07-22"


PYTHON="${PHENGO_PYTHON:-${PYTHON_DEFAULT}}"
OUTPUT_DIR="${PHENGO_OUTPUT_DIR:-${OUTPUT_DEFAULT}}"
LEDGER="${PHENGO_LEDGER:-${REPO_DIR}/data/publication_snapshots.tsv}"

ALL_TRACKS="primary is_a_only no_iea fixed_2025 fixed_2025_is_a primary_mod_gaf legacy_like_go_archive legacy_like_mod_gaf fail_closed_driver_aware_go_archive fail_closed_driver_aware_mod_gaf worm_fixed_2025_terms_go_archive worm_fixed_2025_terms_mod_gaf mouse_impc_assertions mouse_mgi_genepheno_fixed_terms mouse_mgi_phenogenomp_fixed_terms"

TRACKS="${PHENGO_TRACKS:-${ALL_TRACKS}}"
ORGANISMS="${PHENGO_ORGANISMS:-fish fly mouse worm yeast}"
YEARS="${PHENGO_YEARS:-all}"
MODELS="${PHENGO_MODELS:-all}"
PANELS="${PHENGO_PANELS:-full matched_features matched_genes matched_both}"
SINGLE_SNAPSHOT_MODELS="${PHENGO_SINGLE_SNAPSHOT_MODELS:-all}"
SINGLE_SNAPSHOT_YEARS="${PHENGO_SINGLE_SNAPSHOT_YEARS:-2025}"
CV_FOLDS="${PHENGO_CV_FOLDS:-5}"
CV_REPEATS="${PHENGO_CV_REPEATS:-5}"
CALIBRATION="${PHENGO_CALIBRATION:-sigmoid}"
CALIBRATION_CV="${PHENGO_CALIBRATION_CV:-3}"
N_ESTIMATORS="${PHENGO_N_ESTIMATORS:-300}"
N_JOBS="${PHENGO_N_JOBS:-8}"
MLP_HIDDEN_UNITS="${PHENGO_MLP_HIDDEN_UNITS:-128 64}"
MLP_ALPHA="${PHENGO_MLP_ALPHA:-0.0001}"
MLP_BATCH_SIZE="${PHENGO_MLP_BATCH_SIZE:-32}"
MLP_LEARNING_RATE="${PHENGO_MLP_LEARNING_RATE:-0.001}"
MLP_MAX_ITER="${PHENGO_MLP_MAX_ITER:-300}"
MLP_EARLY_STOPPING="${PHENGO_MLP_EARLY_STOPPING:-1}"
MLP_VALIDATION_FRACTION="${PHENGO_MLP_VALIDATION_FRACTION:-0.1}"
MLP_PATIENCE="${PHENGO_MLP_PATIENCE:-15}"
IMPORTANCE_REPEATS="${PHENGO_IMPORTANCE_REPEATS:-20}"
TOP_K="${PHENGO_TOP_K:-25}"
SEED="${PHENGO_SEED:-42}"
MIN_GO_GENE_COUNT="${PHENGO_MIN_GO_GENE_COUNT:-2}"
MAX_GO_GENE_FRACTION="${PHENGO_MAX_GO_GENE_FRACTION:-0.99}"
RUN_TESTS="${PHENGO_RUN_TESTS:-1}"
RUN_SINGLE_SNAPSHOT_ML="${PHENGO_RUN_SINGLE_SNAPSHOT_ML:-1}"
NN_EPOCHS="${PHENGO_NN_EPOCHS:-100}"
NN_BATCH_SIZE="${PHENGO_NN_BATCH_SIZE:-32}"
NN_HIDDEN_UNITS="${PHENGO_NN_HIDDEN_UNITS:-128 64}"
NN_DROPOUT="${PHENGO_NN_DROPOUT:-0.3}"
SINGLE_TEST_SIZE="${PHENGO_SINGLE_TEST_SIZE:-0.2}"
PERM_REPEATS="${PHENGO_PERM_REPEATS:-10}"
DEEP_PREFLIGHT="${PHENGO_DEEP_PREFLIGHT:-1}"
REQUIRE_ANCHORS="${PHENGO_REQUIRE_ANCHORS:-1}"
RESUME=0
DRY_RUN=0
QUICK=0

INPUT_HELPER="${SCRIPT_DIR}/v2/alternative_inputs.py"
MOUSE_LETHAL_TERMS="${REPO_DIR}/src/PhenGO/data/mouse/mouse_lethal_terms.txt.gz"

usage() {
    cat <<'EOF'
Usage: scripts/run_publication_pipeline_unified.sh [options]

Options:
  --output DIR            Output root (default: New_Outputs/publication_run_<date>)
  --ledger FILE           Snapshot ledger (default: data/publication_snapshots.tsv)
  --tracks "LIST"         Space-separated tracks to run (default: all 15)
  --organisms "LIST"      Space-separated organisms to run
  --years "LIST"          Space-separated years, or all (default: all)
  --resume                Skip completed units in an existing output directory
  --quick                 Fast smoke run across a representative track/organism subset
  --dry-run               Validate and print commands without running computations
  --skip-tests            Do not run the repository test suite during preflight
  --require-anchors       Fail if any annotation-anchor input is absent (default)
  --allow-missing-anchors Record and skip absent annotation anchors
  -h, --help              Show this help

Publication defaults can also be overridden with PHENGO_* environment variables.
Runs on both arm64 and x86_64 Python environments; the detected architecture
is recorded in provenance and folded into the run fingerprint. The full run
is computationally intensive; --resume is safe after interruption, and will
automatically do a clean rebuild if any part of the run fingerprint changed.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --ledger) LEDGER="$2"; shift 2 ;;
        --tracks) TRACKS="$2"; shift 2 ;;
        --organisms) ORGANISMS="$2"; shift 2 ;;
        --years) YEARS="$2"; shift 2 ;;
        --resume) RESUME=1; shift ;;
        --quick) QUICK=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        --skip-tests) RUN_TESTS=0; shift ;;
        --require-anchors) REQUIRE_ANCHORS=1; shift ;;
        --allow-missing-anchors) REQUIRE_ANCHORS=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ ${QUICK} -eq 1 ]]; then
    TRACKS="primary fixed_2025 primary_mod_gaf legacy_like_go_archive mouse_impc_assertions"
    ORGANISMS="fly mouse worm"
    YEARS="2024 2025"
    MODELS="lr"
    SINGLE_SNAPSHOT_MODELS="lr"
    PANELS="full matched_both"
    CV_FOLDS=2
    CV_REPEATS=2
    CALIBRATION="none"
    N_ESTIMATORS=20
    IMPORTANCE_REPEATS=2
    TOP_K=5
    N_JOBS=2
fi

if [[ -z "${TRACKS// }" || -z "${ORGANISMS// }" || -z "${YEARS// }" ]]; then
    printf 'Track, organism, and year selections must not be empty.\n' >&2
    exit 2
fi
for track in ${TRACKS}; do
    case "${track}" in
        primary|is_a_only|no_iea|fixed_2025|fixed_2025_is_a|\
        primary_mod_gaf|legacy_like_go_archive|legacy_like_mod_gaf|\
        fail_closed_driver_aware_go_archive|fail_closed_driver_aware_mod_gaf|\
        worm_fixed_2025_terms_go_archive|worm_fixed_2025_terms_mod_gaf|\
        mouse_impc_assertions|mouse_mgi_genepheno_fixed_terms|\
        mouse_mgi_phenogenomp_fixed_terms) ;;
        *) printf 'Unknown publication track: %s\n' "${track}" >&2; exit 2 ;;
    esac
done
for organism in ${ORGANISMS}; do
    case "${organism}" in
        fish|fly|mouse|worm|yeast) ;;
        *) printf 'Unknown organism: %s\n' "${organism}" >&2; exit 2 ;;
    esac
done
for model in ${MODELS} ${SINGLE_SNAPSHOT_MODELS}; do
    case "${model}" in
        all|dt|rf|gb|lr|svm|nn) ;;
        *) printf 'Unknown model: %s\n' "${model}" >&2; exit 2 ;;
    esac
done
for panel in ${PANELS}; do
    case "${panel}" in
        full|matched_features|matched_genes|matched_both) ;;
        *) printf 'Unknown analysis panel: %s\n' "${panel}" >&2; exit 2 ;;
    esac
done
if [[ "${YEARS}" != "all" ]]; then
    for year in ${YEARS}; do
        if [[ ! "${year}" =~ ^[0-9]{4}$ ]]; then
            printf 'Invalid year selection: %s\n' "${year}" >&2
            exit 2
        fi
    done
fi

if [[ ! -x "${PYTHON}" ]]; then
    printf 'PhenGO Python is not executable: %s\n' "${PYTHON}" >&2
    exit 1
fi
if [[ ! -f "${LEDGER}" ]]; then
    printf 'Snapshot ledger does not exist: %s\n' "${LEDGER}" >&2
    exit 1
fi
if [[ ! -f "${INPUT_HELPER}" ]]; then
    printf 'Required alternative-input helper is missing: %s\n' "${INPUT_HELPER}" >&2
    exit 1
fi
if [[ ! -f "${MOUSE_LETHAL_TERMS}" ]]; then
    printf 'Required mouse lethal-term list is missing: %s\n' "${MOUSE_LETHAL_TERMS}" >&2
    exit 1
fi
if [[ ! "${MLP_EARLY_STOPPING}" =~ ^[01]$ ]]; then
    printf 'PHENGO_MLP_EARLY_STOPPING must be 0 or 1.\n' >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Architecture: record, do not gate. arm64 and x86_64 are both supported;
# anything else gets a warning (some scientific wheels may be unavailable)
# but is not a hard failure.
# ---------------------------------------------------------------------------
ARCHITECTURE="$("${PYTHON}" -c 'import platform; print(platform.machine())')"
case "${ARCHITECTURE}" in
    arm64|aarch64|x86_64|amd64) ;;
    *) printf 'Warning: unrecognised architecture %s -- continuing, but some scientific packages may not have prebuilt wheels.\n' "${ARCHITECTURE}" >&2 ;;
esac

DRY_RUN_MARKER="${OUTPUT_DIR}/00_run_metadata/dry_run.marker"
SAFE_PREFLIGHT_RESTART=0
if [[ -d "${OUTPUT_DIR}" && ! -f "${OUTPUT_DIR}/00_run_metadata/run_fingerprint.sha256" ]]; then
    scientific_output="$({ find \
        "${OUTPUT_DIR}/01_derived_labels" "${OUTPUT_DIR}/02_arff" \
        "${OUTPUT_DIR}/03_qc" "${OUTPUT_DIR}/04_ml" \
        "${OUTPUT_DIR}/04_single_snapshot_ml" "${OUTPUT_DIR}/05_temporal" \
        "${OUTPUT_DIR}/06_publication_tables" "${OUTPUT_DIR}/07_reports" \
        -type f -print -quit 2>/dev/null || true; })"
    [[ -z "${scientific_output}" ]] && SAFE_PREFLIGHT_RESTART=1
fi
if [[ -d "${OUTPUT_DIR}" && -n "$(ls -A "${OUTPUT_DIR}" 2>/dev/null)" && ${RESUME} -ne 1 && ${DRY_RUN} -ne 1 && ! -f "${DRY_RUN_MARKER}" && ${SAFE_PREFLIGHT_RESTART} -ne 1 ]]; then
    printf 'Output directory is not empty: %s\nUse --resume or choose a new directory.\n' "${OUTPUT_DIR}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"/{00_run_metadata,01_derived_labels,02_arff,03_qc,04_ml,04_single_snapshot_ml,05_temporal,06_publication_tables,07_reports,logs}
if [[ ${DRY_RUN} -eq 1 ]]; then
    printf 'planning_only_utc\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" > "${DRY_RUN_MARKER}"
else
    rm -f "${DRY_RUN_MARKER}"
fi

LOCK_FILE="${OUTPUT_DIR}/00_run_metadata/publication_pipeline.lock"
if [[ ${DRY_RUN} -eq 0 ]]; then
    if [[ -s "${LOCK_FILE}" ]]; then
        lock_pid="$(awk -F'\t' '$1 == "pid" {print $2}' "${LOCK_FILE}" 2>/dev/null || true)"
        if [[ "${lock_pid}" =~ ^[0-9]+$ ]] && kill -0 "${lock_pid}" 2>/dev/null; then
            printf 'Publication pipeline is already running with PID %s: %s\n' \
                "${lock_pid}" "${OUTPUT_DIR}" >&2
            exit 1
        fi
    fi
    printf 'pid\t%s\nstarted_utc\t%s\n' "$$" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" > "${LOCK_FILE}"
    cleanup_publication_lock() { rm -f "${LOCK_FILE}"; }
    trap cleanup_publication_lock EXIT
fi

PREVIOUS_RUN_FINGERPRINT=""
if [[ -f "${OUTPUT_DIR}/00_run_metadata/run_fingerprint.sha256" ]]; then
    PREVIOUS_RUN_FINGERPRINT="$(<"${OUTPUT_DIR}/00_run_metadata/run_fingerprint.sha256")"
fi

RUN_LOG="${OUTPUT_DIR}/logs/publication_pipeline.log"
COMMAND_LOG="${OUTPUT_DIR}/00_run_metadata/commands.tsv"
OPTIONAL_REPORT="${OUTPUT_DIR}/00_run_metadata/optional_anchor_status.tsv"
INPUT_LIST="${OUTPUT_DIR}/00_run_metadata/input_paths.txt"
INPUT_AUDIT="${OUTPUT_DIR}/03_qc/alternative_input_audit.tsv"
if [[ ${RESUME} -eq 0 ]]; then
    : > "${RUN_LOG}"
fi
if [[ ${RESUME} -eq 0 || ! -s "${COMMAND_LOG}" ]]; then
    printf 'timestamp_utc\tworking_directory\tcommand\n' > "${COMMAND_LOG}"
fi
if [[ ${RESUME} -eq 0 || ! -s "${OPTIONAL_REPORT}" ]]; then
    printf 'organism\tyear\tstatus\tdetail\n' > "${OPTIONAL_REPORT}"
fi

export PYTHONHASHSEED="${SEED}"
export OMP_NUM_THREADS="${PHENGO_OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${PHENGO_OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${PHENGO_MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${PHENGO_VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${PHENGO_NUMEXPR_NUM_THREADS:-1}"
export MPLBACKEND=Agg
export MPLCONFIGDIR="${OUTPUT_DIR}/00_run_metadata/matplotlib_cache"
mkdir -p "${MPLCONFIGDIR}"

timestamp() { date -u '+%Y-%m-%dT%H:%M:%SZ'; }
quote_command() { local item; for item in "$@"; do printf '%q ' "${item}"; done; }
record_command() {
    local rendered
    rendered="$(quote_command "$@")"
    printf '%s\t%s\t%s\n' "$(timestamp)" "${REPO_DIR}" "${rendered}" >> "${COMMAND_LOG}"
}
run_cmd() {
    record_command "$@"
    printf '\n[%s] ' "$(timestamp)" | tee -a "${RUN_LOG}"
    quote_command "$@" | tee -a "${RUN_LOG}"
    printf '\n' | tee -a "${RUN_LOG}"
    if [[ ${DRY_RUN} -eq 0 ]]; then
        "$@" 2>&1 | tee -a "${RUN_LOG}"
    fi
}
run_capture() {
    local capture="$1"
    shift
    record_command "$@"
    printf '\n[%s] ' "$(timestamp)" | tee -a "${RUN_LOG}"
    quote_command "$@" | tee -a "${RUN_LOG}"
    printf '\n' | tee -a "${RUN_LOG}"
    if [[ ${DRY_RUN} -eq 0 ]]; then
        "$@" 2>&1 | tee -a "${RUN_LOG}" | tee "${capture}"
    fi
}

selected_organism() { [[ " ${ORGANISMS} " == *" $1 "* ]]; }
selected_track() { [[ " ${TRACKS} " == *" $1 "* ]]; }
selected_year() { [[ "${YEARS}" == "all" || " ${YEARS} " == *" $1 "* ]]; }
selected_single_snapshot_year() { [[ "${SINGLE_SNAPSHOT_YEARS}" == "all" || " ${SINGLE_SNAPSHOT_YEARS} " == *" $1 "* ]]; }
absolute_path() {
    if [[ -z "$1" ]]; then printf ''
    elif [[ "$1" = /* ]]; then printf '%s' "$1"
    else printf '%s/%s' "${REPO_DIR}" "$1"
    fi
}
skip_item() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" >> "${OUTPUT_DIR}/00_run_metadata/skipped_inputs.tsv"; }

emit_ledger_rows() {
    local cohort="$1"
    "${PYTHON}" - "${LEDGER}" "${cohort}" <<'PY'
import csv, sys
ledger, cohort = sys.argv[1:]
fields = [
    "include_mode", "cohort", "organism", "year", "snapshot_id",
    "phenotype_file", "viability_file", "gaf_file", "go_obo_file",
    "phenotype_ontology_file", "phenotype_release", "go_annotation_release",
    "go_ontology_release", "phenotype_ontology_release", "label_source",
    "nonlethal_policy", "mixed_label_policy", "ambiguous_label_policy",
    "retrieval_date",
]
with open(ledger, encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["cohort"] == cohort:
            print("\x1f".join(row.get(field, "") for field in fields))
PY
}

# ---------------------------------------------------------------------------
# Preflight: required-input existence check (single pass, was duplicated
# once per stage before)
# ---------------------------------------------------------------------------
required_file_missing=0
: > "${INPUT_LIST}"
: > "${OUTPUT_DIR}/00_run_metadata/skipped_inputs.tsv"
printf 'track\torganism\ttimepoint\treason\tdetail\n' > "${OUTPUT_DIR}/00_run_metadata/skipped_inputs.tsv"
while IFS=$'\x1f' read -r include_mode cohort organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
    selected_organism "${organism}" || continue
    selected_year "${year}" || continue
    for value in "${phenotype_file}" "${viability_file}" "${gaf_file}" "${go_obo_file}" "${phenotype_ontology_file}"; do
        [[ -n "${value}" ]] || continue
        path="$(absolute_path "${value}")"
        printf '%s\n' "${path}" >> "${INPUT_LIST}"
        if [[ ! -f "${path}" ]]; then
            printf 'Missing required input for %s: %s\n' "${snapshot_id}" "${path}" >&2
            required_file_missing=1
        fi
    done
done < <(emit_ledger_rows primary)

while IFS=$'\x1f' read -r include_mode cohort organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
    selected_organism "${organism}" || continue
    missing=()
    for value in "${gaf_file}" "${go_obo_file}"; do
        path="$(absolute_path "${value}")"
        if [[ ! -f "${path}" ]]; then missing+=("${path}"); else printf '%s\n' "${path}" >> "${INPUT_LIST}"; fi
    done
    if [[ ${#missing[@]} -eq 0 ]]; then
        printf '%s\t%s\tavailable\tall annotation-anchor inputs present\n' "${organism}" "${year}" >> "${OPTIONAL_REPORT}"
    else
        printf '%s\t%s\tskipped\t%s\n' "${organism}" "${year}" "${missing[*]}" >> "${OPTIONAL_REPORT}"
        if [[ ${REQUIRE_ANCHORS} -eq 1 ]]; then required_file_missing=1; fi
    fi
done < <(emit_ledger_rows annotation_anchor)

if [[ ${required_file_missing} -ne 0 ]]; then
    printf 'Preflight failed because required inputs are missing.\n' >&2
    exit 1
fi
sort -u "${INPUT_LIST}" -o "${INPUT_LIST}"
cp "${LEDGER}" "${OUTPUT_DIR}/00_run_metadata/publication_snapshots.tsv"

cat > "${OUTPUT_DIR}/00_run_metadata/run_config.env" <<EOF
PHENGO_PYTHON=${PYTHON}
PHENGO_OUTPUT_DIR=${OUTPUT_DIR}
PHENGO_LEDGER=${LEDGER}
PHENGO_TRACKS=${TRACKS}
PHENGO_ORGANISMS=${ORGANISMS}
PHENGO_YEARS=${YEARS}
PHENGO_MODELS=${MODELS}
PHENGO_PANELS=${PANELS}
PHENGO_SINGLE_SNAPSHOT_MODELS=${SINGLE_SNAPSHOT_MODELS}
PHENGO_SINGLE_SNAPSHOT_YEARS=${SINGLE_SNAPSHOT_YEARS}
PHENGO_MLP_HIDDEN_UNITS=${MLP_HIDDEN_UNITS}
PHENGO_MLP_ALPHA=${MLP_ALPHA}
PHENGO_MLP_BATCH_SIZE=${MLP_BATCH_SIZE}
PHENGO_MLP_LEARNING_RATE=${MLP_LEARNING_RATE}
PHENGO_MLP_MAX_ITER=${MLP_MAX_ITER}
PHENGO_MLP_EARLY_STOPPING=${MLP_EARLY_STOPPING}
PHENGO_MLP_VALIDATION_FRACTION=${MLP_VALIDATION_FRACTION}
PHENGO_MLP_PATIENCE=${MLP_PATIENCE}
PHENGO_CV_FOLDS=${CV_FOLDS}
PHENGO_CV_REPEATS=${CV_REPEATS}
PHENGO_CALIBRATION=${CALIBRATION}
PHENGO_CALIBRATION_CV=${CALIBRATION_CV}
PHENGO_N_ESTIMATORS=${N_ESTIMATORS}
PHENGO_N_JOBS=${N_JOBS}
PHENGO_IMPORTANCE_REPEATS=${IMPORTANCE_REPEATS}
PHENGO_TOP_K=${TOP_K}
PHENGO_SEED=${SEED}
PHENGO_MIN_GO_GENE_COUNT=${MIN_GO_GENE_COUNT}
PHENGO_MAX_GO_GENE_FRACTION=${MAX_GO_GENE_FRACTION}
PHENGO_RUN_SINGLE_SNAPSHOT_ML=${RUN_SINGLE_SNAPSHOT_ML}
PHENGO_NN_EPOCHS=${NN_EPOCHS}
PHENGO_NN_BATCH_SIZE=${NN_BATCH_SIZE}
PHENGO_NN_HIDDEN_UNITS=${NN_HIDDEN_UNITS}
PHENGO_NN_DROPOUT=${NN_DROPOUT}
PHENGO_SINGLE_TEST_SIZE=${SINGLE_TEST_SIZE}
PHENGO_PERM_REPEATS=${PERM_REPEATS}
ARCHITECTURE=${ARCHITECTURE}
EOF

{
    printf 'started_utc\t%s\n' "$(timestamp)"
    printf 'uname\t'; uname -a
    printf 'python\t'; "${PYTHON}" -VV
    printf 'python_machine\t%s\n' "${ARCHITECTURE}"
    printf 'phenGO_version\t'; "${PYTHON}" -c 'from PhenGO.constants import PhenGO_VERSION; print(PhenGO_VERSION)'
    if command -v sw_vers >/dev/null 2>&1; then sw_vers; fi
    if command -v sysctl >/dev/null 2>&1; then
        sysctl -n machdep.cpu.brand_string 2>/dev/null || true
        sysctl -n hw.memsize 2>/dev/null || true
        sysctl -n hw.ncpu 2>/dev/null || true
    elif command -v lscpu >/dev/null 2>&1; then
        lscpu 2>/dev/null || true
        nproc 2>/dev/null || true
        awk '/MemTotal/{print}' /proc/meminfo 2>/dev/null || true
    fi
} > "${OUTPUT_DIR}/00_run_metadata/environment.txt"
"${PYTHON}" -m pip freeze > "${OUTPUT_DIR}/00_run_metadata/python_packages.txt"
git -C "${REPO_DIR}" rev-parse HEAD > "${OUTPUT_DIR}/00_run_metadata/git_commit.txt"
git -C "${REPO_DIR}" status --short > "${OUTPUT_DIR}/00_run_metadata/git_status.txt"
git -C "${REPO_DIR}" diff --stat > "${OUTPUT_DIR}/00_run_metadata/git_diff_stat.txt"

if [[ ${DRY_RUN} -eq 0 && ${DEEP_PREFLIGHT} -eq 1 ]]; then
    : > "${OUTPUT_DIR}/00_run_metadata/input_checksums.sha256"
    while IFS= read -r path; do
        if [[ "${path}" == *.gz ]]; then gzip -t "${path}"; fi
        sha256 "${path}" >> "${OUTPUT_DIR}/00_run_metadata/input_checksums.sha256"
    done < "${INPUT_LIST}"
fi

# ---------------------------------------------------------------------------
# Single run fingerprint across every knob that affects scientific output
# ---------------------------------------------------------------------------
SOURCE_TREE_HASH="$("${PYTHON}" -c "from PhenGO.provenance import source_tree_sha256; print(source_tree_sha256('${REPO_DIR}'))")"
DEPENDENCY_STATE="$("${PYTHON}" -c 'import json; from PhenGO.provenance import dependency_versions; print(json.dumps(dependency_versions(), sort_keys=True))')"
PIPELINE_HASH="$(sha256 "${BASH_SOURCE[0]}" "${LEDGER}" "${INPUT_HELPER}" "${MOUSE_LETHAL_TERMS}" | sha256_of_stdin)"
INPUT_HASH="not-computed"
if [[ -s "${OUTPUT_DIR}/00_run_metadata/input_checksums.sha256" ]]; then
    INPUT_HASH="$(sha256 "${OUTPUT_DIR}/00_run_metadata/input_checksums.sha256" | awk '{print $1}')"
fi
RUN_FINGERPRINT="$(printf '%s\n' \
    "${SOURCE_TREE_HASH}" "${DEPENDENCY_STATE}" "${PIPELINE_HASH}" "${INPUT_HASH}" "${PYTHON}" "${ARCHITECTURE}" \
    "${TRACKS}" "${ORGANISMS}" "${YEARS}" "${MODELS}" "${PANELS}" \
    "${SINGLE_SNAPSHOT_MODELS}" "${SINGLE_SNAPSHOT_YEARS}" \
    "${CV_FOLDS}" "${CV_REPEATS}" "${CALIBRATION}" "${CALIBRATION_CV}" \
    "${N_ESTIMATORS}" "${N_JOBS}" "${MLP_HIDDEN_UNITS}" "${MLP_ALPHA}" \
    "${MLP_BATCH_SIZE}" "${MLP_LEARNING_RATE}" "${MLP_MAX_ITER}" \
    "${MLP_EARLY_STOPPING}" "${MLP_VALIDATION_FRACTION}" "${MLP_PATIENCE}" \
    "${IMPORTANCE_REPEATS}" "${TOP_K}" "${SEED}" "${MIN_GO_GENE_COUNT}" \
    "${MAX_GO_GENE_FRACTION}" "${RUN_SINGLE_SNAPSHOT_ML}" "${NN_EPOCHS}" \
    "${NN_BATCH_SIZE}" "${NN_HIDDEN_UNITS}" "${NN_DROPOUT}" \
    "${SINGLE_TEST_SIZE}" "${PERM_REPEATS}" | sha256_of_stdin)"

if [[ ${RESUME} -eq 1 && "${PREVIOUS_RUN_FINGERPRINT}" != "${RUN_FINGERPRINT}" ]]; then
    printf 'Resume fingerprint changed (config, code, inputs, or architecture); clearing and regenerating all scientific outputs.\n' | tee -a "${RUN_LOG}"
    RESUME=0
    if [[ ${DRY_RUN} -eq 0 ]]; then
        for generated_dir in \
            "${OUTPUT_DIR}/01_derived_labels" "${OUTPUT_DIR}/02_arff" \
            "${OUTPUT_DIR}/03_qc" "${OUTPUT_DIR}/04_ml" \
            "${OUTPUT_DIR}/04_single_snapshot_ml" "${OUTPUT_DIR}/05_temporal" \
            "${OUTPUT_DIR}/06_publication_tables" "${OUTPUT_DIR}/07_reports"; do
            find "${generated_dir}" -mindepth 1 -delete
        done
        rm -f "${OUTPUT_DIR}/00_run_metadata/output_checksums.sha256"
    fi
fi
printf '%s\n' "${RUN_FINGERPRINT}" > "${OUTPUT_DIR}/00_run_metadata/run_fingerprint.sha256"
if [[ ${DRY_RUN} -eq 0 ]]; then
    rm -f "${OUTPUT_DIR}/00_run_metadata/run.complete"
fi

run_cmd "${PYTHON}" -c 'import numpy, pandas, scipy, sklearn, networkx, tensorflow, matplotlib; print("scientific Python imports: OK")'
if [[ ${RUN_TESTS} -eq 1 ]]; then
    run_cmd env PYTHONDONTWRITEBYTECODE=1 "${PYTHON}" -m pytest -q
fi
run_cmd "${PYTHON}" "${INPUT_HELPER}" audit --repo-root "${REPO_DIR}" \
    --output-tsv "${INPUT_AUDIT}" --output-json "${OUTPUT_DIR}/03_qc/alternative_input_audit.json"

assert_audit_eligible() {
    local source="$1" relative="${source#${REPO_DIR}/}"
    [[ ${DRY_RUN} -eq 1 ]] && return
    "${PYTHON}" - "${INPUT_AUDIT}" "${relative}" <<'PY'
import csv, sys
audit, target = sys.argv[1:]
with open(audit, encoding="utf-8", newline="") as handle:
    matches = [row for row in csv.DictReader(handle, delimiter="\t") if row["relative_path"] == target]
if len(matches) != 1 or matches[0]["status"] != "eligible":
    status = matches[0]["status"] if matches else "not_in_audit"
    raise SystemExit(f"Selected input is not eligible: {target} ({status})")
PY
}

# ---------------------------------------------------------------------------
# Shared derived-label caches, keyed by the *actual underlying input* rather
# than by track name. Previously, V1 (for primary/is_a_only/no_iea) and V2
# (for legacy_like_go_archive/fail_closed_driver_aware_go_archive) each
# recomputed fly assignments and worm lethal-term closures independently even
# when both were reading the exact same GAF/ontology file for the same year --
# same for worm. These caches are computed once per (year, input-file) and
# reused by every track that needs that combination.
# ---------------------------------------------------------------------------
ensure_mouse_release_sets() {
    local year="$1" source="$2" outdir="${OUTPUT_DIR}/01_derived_labels/mouse/${year}"
    mkdir -p "${outdir}"
    if [[ ${RESUME} -eq 1 && -s "${outdir}/lethal.tsv" && -s "${outdir}/viable.tsv" ]]; then return; fi
    run_capture "${outdir}/conversion_summary.json" \
        "${PYTHON}" -m PhenGO.scripts.publication_inputs mouse-impc \
        --input "${source}" --lethal-output "${outdir}/lethal.tsv" \
        --viable-output "${outdir}/viable.tsv" --excluded-output "${outdir}/excluded.tsv"
}

FLY_ASSIGNMENTS=""
ensure_fly_assignments() {
    # kind is "archive" or "mod" -- the two distinct GAF provenances used
    # across every fly track; computed once per (year, kind), not per track.
    local year="$1" kind="$2" gaf="$3" outdir="${OUTPUT_DIR}/01_derived_labels/fly/${year}/${kind}"
    FLY_ASSIGNMENTS="${outdir}/assignments.tsv"
    mkdir -p "${outdir}"
    if [[ ${RESUME} -eq 1 && -s "${FLY_ASSIGNMENTS}" ]]; then return; fi
    run_capture "${outdir}/conversion_summary.json" \
        "${PYTHON}" -m PhenGO.scripts.publication_inputs fly-assignments \
        --gaf "${gaf}" --output "${FLY_ASSIGNMENTS}" \
        --excluded-output "${outdir}/excluded.tsv" --taxon 7227
}

WORM_TERMS=""
ensure_worm_terms() {
    # key is either a year (current-year ontology) or the literal
    # "fixed_2025" (single 2025-derived closure reused across every year).
    local key="$1" ontology="$2" outdir="${OUTPUT_DIR}/01_derived_labels/worm/${key}"
    WORM_TERMS="${outdir}/lethal_terms.tsv"
    mkdir -p "${outdir}"
    if [[ ${RESUME} -eq 1 && -s "${WORM_TERMS}" ]]; then return; fi
    run_cmd "${PYTHON}" -m PhenGO.scripts.get_phenotype_terms \
        --obo-file "${ontology}" \
        --term-list "${REPO_DIR}/data/worm/lethal_terms/root_lethal_phenotype_terms.txt" \
        --results "${WORM_TERMS}"
}

PREPARED_GAF=""
PREPARED_GAF_LABEL=""
prepare_mod_gaf() {
    local organism="$1" year="$2" archive_gaf="$3" source taxon database outdir
    case "${organism}" in
        fish)
            source="$(find "${REPO_DIR}/data/fish/gene_association/zfin/${year}" -maxdepth 1 -type f -name '*.gz' 2>/dev/null | LC_ALL=C sort | head -1 || true)"
            taxon=7955; database=ZFIN ;;
        fly)
            source="${REPO_DIR}/data/fly/gene_association/fb/${year}/gene_association.fb.gz"
            taxon=7227; database=FB ;;
        worm)
            source="$(find "${REPO_DIR}/data/worm/gene_association/wb/${year}" -maxdepth 1 -type f -name '*.clean.gz' 2>/dev/null | LC_ALL=C sort | head -1 || true)"
            taxon=6239; database=WB ;;
        yeast)
            source="$(find "${REPO_DIR}/data/yeast/gene_association/sgd/${year}" -maxdepth 1 -type f -name '*.gz' 2>/dev/null | LC_ALL=C sort | head -1 || true)"
            taxon=559292; database=SGD ;;
        *) printf 'No alternative MOD GAF mapping for %s\n' "${organism}" >&2; return 1 ;;
    esac
    if [[ ! -f "${source}" ]]; then
        printf 'Alternative MOD GAF is missing: %s/%s\n' "${organism}" "${year}" >&2
        return 1
    fi
    assert_audit_eligible "${source}"
    outdir="${OUTPUT_DIR}/01_derived_inputs/gaf/${organism}/${year}"
    PREPARED_GAF="${outdir}/${organism}_${year}_taxon_filtered.gaf.gz"
    PREPARED_GAF_LABEL="$(basename "${source}")"
    mkdir -p "${outdir}" "${OUTPUT_DIR}/03_qc/gaf_comparisons/${organism}"
    if [[ ! (${RESUME} -eq 1 && -s "${PREPARED_GAF}" && -s "${outdir}/filter_summary.json") ]]; then
        run_cmd "${PYTHON}" "${INPUT_HELPER}" filter-gaf --input "${source}" \
            --output "${PREPARED_GAF}" --taxon "${taxon}" --database "${database}" \
            --summary "${outdir}/filter_summary.json"
    fi
    if [[ ! (${RESUME} -eq 1 && -s "${OUTPUT_DIR}/03_qc/gaf_comparisons/${organism}/${year}.json") ]]; then
        run_cmd "${PYTHON}" "${INPUT_HELPER}" compare-gaf --left "${archive_gaf}" \
            --right "${PREPARED_GAF}" --taxon "${taxon}" \
            --output "${OUTPUT_DIR}/03_qc/gaf_comparisons/${organism}/${year}.json"
    fi
}

FLY_LETHAL=""
FLY_VIABLE=""
FLY_LABEL_AUDIT=""
prepare_fly_fail_closed_labels() {
    local year="$1" phenotype="$2" outdir="${OUTPUT_DIR}/01_derived_labels/fly_fail_closed/${year}"
    FLY_LETHAL="${outdir}/lethal.tsv.gz"
    FLY_VIABLE="${outdir}/viable.tsv.gz"
    FLY_LABEL_AUDIT="${outdir}/row_audit.tsv.gz"
    mkdir -p "${outdir}"
    if [[ ! (${RESUME} -eq 1 && -s "${FLY_LETHAL}" && -s "${FLY_VIABLE}" && -s "${FLY_LABEL_AUDIT}") ]]; then
        run_cmd "${PYTHON}" "${INPUT_HELPER}" fly-labels --input "${phenotype}" \
            --lethal-output "${FLY_LETHAL}" --viable-output "${FLY_VIABLE}" \
            --excluded-output "${outdir}/excluded_genes.tsv.gz" \
            --audit-output "${FLY_LABEL_AUDIT}" --summary "${outdir}/summary.json"
    fi
}

# ---------------------------------------------------------------------------
# build_arff: single unified ARFF constructor for every track. go_relations
# and evidence-code handling depend on the track name; provenance metadata
# (availability/retrieval route) defaults to "ordinary local file" and is
# overridden for tracks known to use alternative-provenance data.
# aux_kind/aux_path attach whichever label mechanism the track needs:
#   fly_assignments   aux_path=<assignments.tsv>
#   fly_gene_sets     aux_path=<lethal>|<viable>|<assignments>
#   worm_terms        aux_path=<lethal_terms.tsv>
#   gene_sets         aux_path=<lethal>|<viable>
#   none              (no additional label flags)
# ---------------------------------------------------------------------------
build_arff() {
    local track="$1" organism="$2" timepoint="$3" snapshot_id="$4" phenotype_file="$5"
    local gaf="$6" go_obo="$7" phenotype_release="$8" go_annotation_release="$9"
    shift 9
    local go_ontology_release="$1" phenotype_ontology_release="$2" retrieval_date="$3"
    local label_source="$4" nonlethal_policy="$5" mixed_policy="$6" ambiguous_policy="$7"
    local aux_kind="$8" aux_path="$9"
    shift 9
    local extra=("$@") output qc arff phenotype_availability go_annotation_availability retrieval_route

    output="${OUTPUT_DIR}/02_arff/${track}/${organism}/${timepoint}"
    qc="${OUTPUT_DIR}/03_qc/${track}/${organism}/${timepoint}"
    arff="${output}/${organism}_PhenGO.arff"

    phenotype_availability="available_local"
    go_annotation_availability="available_local_go_archive"
    retrieval_route="local_archived_or_provider_file"
    [[ "${track}" == *mod_gaf* ]] && go_annotation_availability="available_local_mod_provider_file"
    if [[ "${track}" == mouse_mgi_* ]]; then
        phenotype_availability="web_archive_recovered"
        retrieval_route="wayback_machine_capture"
    elif [[ "${track}" == "mouse_impc_assertions" ]]; then
        phenotype_availability="available_local_incomplete_export"
        retrieval_route="local_IMPC_export"
    fi

    if [[ ${RESUME} -eq 1 && -s "${arff}" && -s "${output}/PhenGO_manifest.json" ]]; then return; fi
    mkdir -p "${output}" "${qc}"

    local command=(
        "${PYTHON}" -m PhenGO.core.PhenGO
        -species "${organism}" -gene_association_file "${gaf}" -go_obo_file "${go_obo}"
        -output_dir "${output}" -snapshot_id "${snapshot_id}"
        -phenotype_release "${phenotype_release}"
        -go_annotation_release "${go_annotation_release}"
        -go_ontology_release "${go_ontology_release}" -retrieval_date "${retrieval_date}"
        -snapshot_semantics declared_composite_cross_section
        -phenotype_availability "${phenotype_availability}"
        -go_annotation_availability "${go_annotation_availability}"
        -go_ontology_availability available_local_annual_GO_file
        -retrieval_route "${retrieval_route}"
        -label_source "${label_source}" -nonlethal_policy "${nonlethal_policy}"
        -mixed_label_policy "${mixed_policy}" -ambiguous_label_policy "${ambiguous_policy}"
        -min_go_gene_count "${MIN_GO_GENE_COUNT}"
        -max_go_gene_fraction "${MAX_GO_GENE_FRACTION}" -strict_snapshot -overwrite
    )
    [[ -n "${phenotype_file}" ]] && command+=( -phenotype_file "${phenotype_file}" )
    [[ -n "${phenotype_ontology_release}" ]] && command+=( -phenotype_ontology_release "${phenotype_ontology_release}" )

    case "${aux_kind}" in
        fly_assignments) command+=( -fly_assignments "${aux_path}" ) ;;
        fly_gene_sets)
            local lethal="${aux_path%%|*}" remainder="${aux_path#*|}"
            local viable="${remainder%%|*}" assignments="${remainder#*|}"
            command+=( -lethal_gene_set "${lethal}" -viable_gene_set "${viable}" \
                -fly_assignments "${assignments}" )
            ;;
        worm_terms) command+=( -worm_phenotypes "${aux_path}" ) ;;
        gene_sets)
            command+=( -lethal_gene_set "${aux_path%%|*}" -viable_gene_set "${aux_path#*|}" ) ;;
        none) ;;
        *) printf 'Unknown ARFF auxiliary kind: %s\n' "${aux_kind}" >&2; return 1 ;;
    esac

    case "${track}" in
        is_a_only|fixed_2025_is_a) command+=( -go_relations is_a ) ;;
        no_iea) command+=( -go_relations is_a part_of -exclude_go_evidence_codes IEA ) ;;
        *) command+=( -go_relations is_a part_of ) ;;
    esac

    if [[ ${#extra[@]} -gt 0 ]]; then command+=( "${extra[@]}" ); fi
    run_cmd "${command[@]}"
    run_cmd "${PYTHON}" -m PhenGO.scripts.arff_validator -i "${arff}" -o "${qc}" --fail-on-error
}

# ---------------------------------------------------------------------------
# Main per-organism/year pass: builds every non-mouse-historical track in
# one loop over the ledger's primary cohort, computing each shared derived
# label artifact at most once per (organism, year).
# ---------------------------------------------------------------------------
while IFS=$'\x1f' read -r include_mode cohort organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
    selected_organism "${organism}" || continue
    selected_year "${year}" || continue
    phenotype="$(absolute_path "${phenotype_file}")"
    archive_gaf="$(absolute_path "${gaf_file}")"
    go_obo="$(absolute_path "${go_obo_file}")"
    phenotype_ontology="$(absolute_path "${phenotype_ontology_file}")"

    need_v1="0"
    for t in primary is_a_only no_iea; do selected_track "${t}" && need_v1=1; done
    need_fixed="0"
    { selected_track fixed_2025 || selected_track fixed_2025_is_a; } && need_fixed=1

    alt_needed=0
    case "${organism}" in
        fish) selected_track primary_mod_gaf && alt_needed=1 ;;
        fly)
            for candidate in primary_mod_gaf legacy_like_mod_gaf fail_closed_driver_aware_mod_gaf; do
                selected_track "${candidate}" && alt_needed=1
            done
            ;;
        worm)
            for candidate in primary_mod_gaf worm_fixed_2025_terms_mod_gaf; do
                selected_track "${candidate}" && alt_needed=1
            done
            ;;
        yeast) selected_track primary_mod_gaf && alt_needed=1 ;;
    esac
    alt_gaf=""
    if [[ ${alt_needed} -eq 1 ]]; then
        if prepare_mod_gaf "${organism}" "${year}" "${archive_gaf}"; then
            alt_gaf="${PREPARED_GAF}"
        elif [[ "${organism}" == "fish" ]]; then
            skip_item primary_mod_gaf fish "${year}" unavailable_provider_gaf \
                "No local ZFIN provider GAF for this year"
        else
            exit 1
        fi
    fi

    # Shared per-(organism,year) derived label caches -- computed once
    # regardless of how many selected tracks need them.
    archive_fly_assignments=""
    mod_fly_assignments=""
    if [[ "${organism}" == "fly" ]]; then
        need_archive_assignments=0
        { [[ ${need_v1} -eq 1 || ${need_fixed} -eq 1 ]] || selected_track legacy_like_go_archive || \
          selected_track fail_closed_driver_aware_go_archive; } && need_archive_assignments=1
        if [[ ${need_archive_assignments} -eq 1 ]]; then
            ensure_fly_assignments "${year}" archive "${archive_gaf}"
            archive_fly_assignments="${FLY_ASSIGNMENTS}"
        fi
        if [[ -n "${alt_gaf}" ]]; then
            ensure_fly_assignments "${year}" mod "${alt_gaf}"
            mod_fly_assignments="${FLY_ASSIGNMENTS}"
        fi
    fi
    current_worm_terms=""
    if [[ "${organism}" == "worm" && ( ${need_v1} -eq 1 || $(selected_track primary_mod_gaf && echo 1 || echo 0) -eq 1 ) ]]; then
        ensure_worm_terms "${year}" "${phenotype_ontology}"
        current_worm_terms="${WORM_TERMS}"
    fi
    if [[ "${organism}" == "mouse" && "${label_source}" == "gene_sets" && ${need_v1} -eq 1 ]]; then
        ensure_mouse_release_sets "${year}" "$(absolute_path "${viability_file}")"
    fi

    # -- V1 tracks: primary / is_a_only / no_iea --
    for track in primary is_a_only no_iea; do
        selected_track "${track}" || continue
        track_snapshot="${snapshot_id}"; [[ "${track}" != "primary" ]] && track_snapshot="${snapshot_id}-${track}"
        aux_kind="none"; aux_path=""
        case "${organism}" in
            fly) aux_kind="fly_assignments"; aux_path="${archive_fly_assignments}" ;;
            worm)
                if [[ "${label_source}" == "release_records" ]]; then
                    aux_kind="worm_terms"; aux_path="${current_worm_terms}"
                fi
                ;;
            mouse)
                if [[ "${label_source}" == "gene_sets" ]]; then
                    aux_kind="gene_sets"
                    aux_path="${OUTPUT_DIR}/01_derived_labels/mouse/${year}/lethal.tsv|${OUTPUT_DIR}/01_derived_labels/mouse/${year}/viable.tsv"
                fi
                ;;
        esac
        build_arff "${track}" "${organism}" "${year}" "${track_snapshot}" \
            "${phenotype}" "${archive_gaf}" "${go_obo}" "${phenotype_release}" \
            "${go_annotation_release}" "${go_ontology_release}" "${phenotype_ontology_release}" \
            "${retrieval_date}" "${label_source}" "${nonlethal_policy}" "${mixed_label_policy}" \
            "${ambiguous_label_policy}" "${aux_kind}" "${aux_path}"
    done

    # -- V2 track: primary_mod_gaf --
    if selected_track primary_mod_gaf && [[ -n "${alt_gaf}" || "${organism}" != "fish" ]]; then
        case "${organism}" in
            fish)
                if [[ -n "${alt_gaf}" ]]; then
                    build_arff primary_mod_gaf fish "${year}" "fish-${year}-primary-mod-gaf" \
                        "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}" \
                        "${go_annotation_release}; ZFIN provider file ${PREPARED_GAF_LABEL}" \
                        "${go_ontology_release}" "" "${retrieval_date}" "${label_source}" \
                        "${nonlethal_policy}" "${mixed_label_policy}" "${ambiguous_label_policy}" none ""
                fi
                ;;
            fly)
                build_arff primary_mod_gaf fly "${year}" "fly-${year}-primary-mod-gaf" \
                    "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}" \
                    "${go_annotation_release}; FlyBase provider file ${PREPARED_GAF_LABEL}" "${go_ontology_release}" \
                    "${phenotype_ontology_release}" "${retrieval_date}" "${label_source}" \
                    "${nonlethal_policy}" "${mixed_label_policy}" "${ambiguous_label_policy}" \
                    fly_assignments "${mod_fly_assignments}"
                ;;
            worm)
                build_arff primary_mod_gaf worm "${year}" "worm-${year}-primary-mod-gaf" \
                    "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}" \
                    "${go_annotation_release}; WormBase provider file ${PREPARED_GAF_LABEL}; taxon-filtered" \
                    "${go_ontology_release}" "${phenotype_ontology_release}" "${retrieval_date}" \
                    release_records "${nonlethal_policy}" "${mixed_label_policy}" \
                    "${ambiguous_label_policy}" worm_terms "${current_worm_terms}"
                ;;
            yeast)
                build_arff primary_mod_gaf yeast "${year}" "yeast-${year}-primary-mod-gaf" \
                    "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}" \
                    "${go_annotation_release}; SGD provider file ${PREPARED_GAF_LABEL}" \
                    "${go_ontology_release}" "" "${retrieval_date}" release_records \
                    "${nonlethal_policy}" "${mixed_label_policy}" "${ambiguous_label_policy}" none ""
                ;;
        esac
    fi

    # -- V2 tracks: legacy_like_go_archive / legacy_like_mod_gaf --
    if selected_track legacy_like_go_archive; then
        case "${organism}" in
            fly)
                build_arff legacy_like_go_archive fly "${year}" "fly-${year}-legacy-like-go-archive" \
                    "${phenotype}" "${archive_gaf}" "${go_obo}" "${phenotype_release}; legacy-like policy" \
                    "${go_annotation_release}" "${go_ontology_release}" "${phenotype_ontology_release}" \
                    "${retrieval_date}" release_records observed_viable lethal_wins exclude \
                    fly_assignments "${archive_fly_assignments}" -allow_multigenes
                ;;
            fish)
                build_arff legacy_like_go_archive fish "${year}" "fish-${year}-legacy-like-go-archive" \
                    "${phenotype}" "${archive_gaf}" "${go_obo}" "${phenotype_release}; legacy-like policy" \
                    "${go_annotation_release}" "${go_ontology_release}" "" "${retrieval_date}" \
                    release_records observed_viable lethal_wins viable none ""
                ;;
            yeast)
                build_arff legacy_like_go_archive yeast "${year}" "yeast-${year}-legacy-like-go-archive" \
                    "${phenotype}" "${archive_gaf}" "${go_obo}" "${phenotype_release}; legacy-like policy" \
                    "${go_annotation_release}" "${go_ontology_release}" "" "${retrieval_date}" \
                    release_records observed_viable lethal_wins exclude none ""
                ;;
        esac
    fi
    if selected_track legacy_like_mod_gaf && [[ "${organism}" == "fly" ]]; then
        build_arff legacy_like_mod_gaf fly "${year}" "fly-${year}-legacy-like-mod-gaf" \
            "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}; legacy-like policy" \
            "${go_annotation_release}; FlyBase provider file ${PREPARED_GAF_LABEL}" "${go_ontology_release}" \
            "${phenotype_ontology_release}" "${retrieval_date}" release_records observed_viable \
            lethal_wins exclude fly_assignments "${mod_fly_assignments}" -allow_multigenes
    fi

    # -- V2 tracks: fail_closed_driver_aware_go_archive / _mod_gaf (fly only) --
    if [[ "${organism}" == "fly" ]] && { selected_track fail_closed_driver_aware_go_archive || selected_track fail_closed_driver_aware_mod_gaf; }; then
        prepare_fly_fail_closed_labels "${year}" "${phenotype}"
        if selected_track fail_closed_driver_aware_go_archive; then
            build_arff fail_closed_driver_aware_go_archive fly "${year}" "fly-${year}-fail-closed-driver-go-archive" \
                "${phenotype}" "${archive_gaf}" "${go_obo}" "${phenotype_release}; fail-closed compound policy" \
                "${go_annotation_release}" "${go_ontology_release}" \
                "${phenotype_ontology_release}; every compound component audited; unresolved excluded" "${retrieval_date}" \
                gene_sets explicit_only exclude exclude fly_gene_sets \
                "${FLY_LETHAL}|${FLY_VIABLE}|${archive_fly_assignments}"
        fi
        if selected_track fail_closed_driver_aware_mod_gaf; then
            build_arff fail_closed_driver_aware_mod_gaf fly "${year}" "fly-${year}-fail-closed-driver-mod-gaf" \
                "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}; fail-closed compound policy" \
                "${go_annotation_release}; FlyBase provider file ${PREPARED_GAF_LABEL}" "${go_ontology_release}" \
                "${phenotype_ontology_release}; every compound component audited; unresolved excluded" "${retrieval_date}" \
                gene_sets explicit_only exclude exclude fly_gene_sets \
                "${FLY_LETHAL}|${FLY_VIABLE}|${mod_fly_assignments}"
        fi
    fi

    # -- V2 tracks: worm_fixed_2025_terms_go_archive / _mod_gaf --
    if [[ "${organism}" == "worm" ]]; then
        fixed_ontology_rel="$("${PYTHON}" - "${LEDGER}" <<'PY'
import csv, sys
with open(sys.argv[1], encoding="utf-8", newline="") as handle:
    matches = [row["phenotype_ontology_file"] for row in csv.DictReader(handle, delimiter="\t")
               if row.get("cohort") == "primary" and row.get("organism") == "worm" and row.get("year") == "2025"]
if len(matches) != 1:
    raise SystemExit("Expected exactly one primary worm 2025 ontology")
print(matches[0])
PY
)"
        if selected_track worm_fixed_2025_terms_go_archive || selected_track worm_fixed_2025_terms_mod_gaf; then
            ensure_worm_terms fixed_2025 "$(absolute_path "${fixed_ontology_rel}")"
            fixed_terms="${WORM_TERMS}"
        fi
        if selected_track worm_fixed_2025_terms_go_archive; then
            build_arff worm_fixed_2025_terms_go_archive worm "${year}" "worm-${year}-fixed-2025-terms-go-archive" \
                "${phenotype}" "${archive_gaf}" "${go_obo}" "${phenotype_release}" \
                "${go_annotation_release}" "${go_ontology_release}" \
                "Worm phenotype lethal closure from 2025 applied to annual records" "${retrieval_date}" \
                release_records "${nonlethal_policy}" "${mixed_label_policy}" "${ambiguous_label_policy}" \
                worm_terms "${fixed_terms}"
        fi
        if selected_track worm_fixed_2025_terms_mod_gaf; then
            build_arff worm_fixed_2025_terms_mod_gaf worm "${year}" "worm-${year}-fixed-2025-terms-mod-gaf" \
                "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}" \
                "${go_annotation_release}; WormBase provider file ${PREPARED_GAF_LABEL}; taxon-filtered" \
                "${go_ontology_release}" "Worm phenotype lethal closure from 2025 applied to annual records" \
                "${retrieval_date}" release_records "${nonlethal_policy}" "${mixed_label_policy}" \
                "${ambiguous_label_policy}" worm_terms "${fixed_terms}"
        fi
    fi
done < <(emit_ledger_rows primary)

# ---------------------------------------------------------------------------
# fixed_2025 / fixed_2025_is_a: derive fixed labels from the primary track's
# 2025 ARFF, then rebuild every annotation-anchor year against them.
# ---------------------------------------------------------------------------
if selected_track fixed_2025 || selected_track fixed_2025_is_a; then
    if ! selected_track primary; then
        printf 'Fixed-label tracks require the primary track in the same run.\n' >&2
        exit 1
    fi
    for organism in ${ORGANISMS}; do
        source_dir="${OUTPUT_DIR}/02_arff/primary/${organism}/2025"
        fixed_dir="${OUTPUT_DIR}/01_derived_labels/fixed_2025/${organism}"
        mkdir -p "${fixed_dir}"
        if [[ ! (${RESUME} -eq 1 && -s "${fixed_dir}/lethal.tsv" && -s "${fixed_dir}/viable.tsv") ]]; then
            run_capture "${fixed_dir}/conversion_summary.json" \
                "${PYTHON}" -m PhenGO.scripts.publication_inputs arff-labels \
                --arff "${source_dir}/${organism}_PhenGO.arff" \
                --identifiers "${source_dir}/gene_identifiers.tsv" \
                --lethal-output "${fixed_dir}/lethal.tsv" \
                --viable-output "${fixed_dir}/viable.tsv"
        fi
    done

    for fixed_track in fixed_2025 fixed_2025_is_a; do
        selected_track "${fixed_track}" || continue
        while IFS=$'\x1f' read -r include_mode cohort organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
            selected_organism "${organism}" || continue
            selected_year "${year}" || continue
            gaf="$(absolute_path "${gaf_file}")"
            aux_kind="gene_sets"
            aux_path="${OUTPUT_DIR}/01_derived_labels/fixed_2025/${organism}/lethal.tsv|${OUTPUT_DIR}/01_derived_labels/fixed_2025/${organism}/viable.tsv"
            if [[ "${organism}" == "fly" ]]; then
                ensure_fly_assignments "${year}" archive "${gaf}"
                aux_kind="fly_gene_sets"
                aux_path="${aux_path}|${FLY_ASSIGNMENTS}"
            fi
            build_arff "${fixed_track}" "${organism}" "${year}" "${organism}-${year}-${fixed_track}" \
                "" "${gaf}" "$(absolute_path "${go_obo_file}")" "Fixed 2025 PhenGO labels for ${organism}" \
                "${go_annotation_release}" "${go_ontology_release}" "${phenotype_ontology_release}" \
                "${retrieval_date}" gene_sets explicit_only exclude exclude "${aux_kind}" "${aux_path}"
        done < <(emit_ledger_rows primary)

        while IFS=$'\x1f' read -r include_mode cohort organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
            selected_organism "${organism}" || continue
            selected_year "${year}" || continue
            if [[ "${fixed_track}" == "fixed_2025" && "${year}" -lt 2014 ]]; then continue; fi
            gaf="$(absolute_path "${gaf_file}")"
            go_obo="$(absolute_path "${go_obo_file}")"
            if [[ ! -f "${gaf}" || ! -f "${go_obo}" ]]; then continue; fi
            aux_kind="gene_sets"
            aux_path="${OUTPUT_DIR}/01_derived_labels/fixed_2025/${organism}/lethal.tsv|${OUTPUT_DIR}/01_derived_labels/fixed_2025/${organism}/viable.tsv"
            if [[ "${organism}" == "fly" ]]; then
                ensure_fly_assignments "${year}" archive "${gaf}"
                aux_kind="fly_gene_sets"
                aux_path="${aux_path}|${FLY_ASSIGNMENTS}"
            fi
            build_arff "${fixed_track}" "${organism}" "${year}" "${snapshot_id}-${fixed_track}" \
                "" "${gaf}" "${go_obo}" "${phenotype_release}" \
                "${go_annotation_release}" "${go_ontology_release}" "${phenotype_ontology_release}" \
                "${retrieval_date}" gene_sets explicit_only exclude exclude "${aux_kind}" "${aux_path}"
        done < <(emit_ledger_rows annotation_anchor)
    done
fi

# ---------------------------------------------------------------------------
# Mouse historical-capture tracks (unchanged logic from V2, ledger-driven)
# ---------------------------------------------------------------------------
emit_mgi_captures() {
    "${PYTHON}" - "${LEDGER}" "${REPO_DIR}" <<'PY'
import csv, os, sys
ledger, root = sys.argv[1:]
captures = [
    ("mouse_mgi_genepheno_fixed_terms", "2016_gzip", "2016", "genepheno", "data/mouse/phenotype_data/historical_mgi_MGI_GenePheno/MGI_GenePheno_2016.rpt.gz"),
    ("mouse_mgi_genepheno_fixed_terms", "2017_plain", "2017", "genepheno", "data/mouse/phenotype_data/historical_mgi_MGI_GenePheno/MGI_GenePheno_2017.rpt"),
    ("mouse_mgi_genepheno_fixed_terms", "2017_gzip", "2017", "genepheno", "data/mouse/phenotype_data/historical_mgi_MGI_GenePheno/MGI_GenePheno_2017.rpt.gz"),
    ("mouse_mgi_genepheno_fixed_terms", "2020_gzip", "2020", "genepheno", "data/mouse/phenotype_data/historical_mgi_MGI_GenePheno/MGI_GenePheno_2020.rpt.gz"),
    ("mouse_mgi_genepheno_fixed_terms", "2021_plain", "2021", "genepheno", "data/mouse/phenotype_data/historical_mgi_MGI_GenePheno/MGI_GenePheno_2021.rpt"),
    ("mouse_mgi_genepheno_fixed_terms", "2021_gzip", "2021", "genepheno", "data/mouse/phenotype_data/historical_mgi_MGI_GenePheno/MGI_GenePheno_2021.rpt.gz"),
    ("mouse_mgi_genepheno_fixed_terms", "2022_gzip", "2022", "genepheno", "data/mouse/phenotype_data/historical_mgi_MGI_GenePheno/MGI_GenePheno_2022.rpt.gz"),
    ("mouse_mgi_genepheno_fixed_terms", "2024_gzip", "2024", "genepheno", "data/mouse/phenotype_data/historical_mgi_MGI_GenePheno/MGI_GenePheno_2024.rpt.gz"),
    ("mouse_mgi_genepheno_fixed_terms", "2025_gzip", "2025", "genepheno", "data/mouse/phenotype_data/historical_mgi_MGI_GenePheno/MGI_GenePheno_2025.rpt.gz"),
    ("mouse_mgi_phenogenomp_fixed_terms", "2017_gzip", "2017", "phenogenomp", "data/mouse/phenotype_data/historical_mgi_MGI_PhenoGenoMP/MGI_PhenoGenoMP_2017.rpt.gz"),
    ("mouse_mgi_phenogenomp_fixed_terms", "2020_gzip", "2020", "phenogenomp", "data/mouse/phenotype_data/historical_mgi_MGI_PhenoGenoMP/MGI_PhenoGenoMP_2020.rpt.gz"),
    ("mouse_mgi_phenogenomp_fixed_terms", "2021_gzip", "2021", "phenogenomp", "data/mouse/phenotype_data/historical_mgi_MGI_PhenoGenoMP/MGI_PhenoGenoMP_2021.rpt.gz"),
    ("mouse_mgi_phenogenomp_fixed_terms", "2024_gzip", "2024", "phenogenomp", "data/mouse/phenotype_data/historical_mgi_MGI_PhenoGenoMP/MGI_PhenoGenoMP_2024.rpt.gz"),
    ("mouse_mgi_phenogenomp_fixed_terms", "2025_gzip", "2025", "phenogenomp", "data/mouse/phenotype_data/historical_mgi_MGI_PhenoGenoMP/MGI_PhenoGenoMP_2025.rpt.gz"),
]
with open(ledger, encoding="utf-8", newline="") as handle:
    mouse = {row["year"]: row for row in csv.DictReader(handle, delimiter="\t") if row["cohort"] == "primary" and row["organism"] == "mouse"}
for track, timepoint, year, fmt, relative in captures:
    row = mouse.get(year)
    if row:
        fields = [track, timepoint, year, fmt, os.path.join(root, relative), row["gaf_file"], row["go_obo_file"], row["go_annotation_release"], row["go_ontology_release"], row["retrieval_date"]]
        print("\x1f".join(fields))
PY
}

if selected_organism mouse && { selected_track mouse_mgi_genepheno_fixed_terms || selected_track mouse_mgi_phenogenomp_fixed_terms; }; then
    while IFS=$'\x1f' read -r track timepoint year report_format source gaf_file go_obo_file go_annotation_release go_ontology_release retrieval_date; do
        selected_track "${track}" || continue
        selected_year "${year}" || continue
        if [[ ! -f "${source}" ]]; then
            skip_item "${track}" mouse "${timepoint}" missing_input "${source}"
            printf 'Required selected historical MGI input is missing for %s/%s: %s\n' \
                "${track}" "${timepoint}" "${source}" >&2
            exit 1
        fi
        assert_audit_eligible "${source}"
        labels="${OUTPUT_DIR}/01_derived_labels/${track}/mouse/${timepoint}"
        mkdir -p "${labels}"
        if [[ ! (${RESUME} -eq 1 && -s "${labels}/lethal.tsv" && -s "${labels}/operational_nonlethal.tsv") ]]; then
            run_cmd "${PYTHON}" "${INPUT_HELPER}" mgi-labels --input "${source}" \
                --format "${report_format}" --lethal-terms "${MOUSE_LETHAL_TERMS}" \
                --lethal-output "${labels}/lethal.tsv" --viable-output "${labels}/operational_nonlethal.tsv" \
                --excluded-output "${labels}/excluded.tsv" --summary "${labels}/summary.json"
        fi
        build_arff "${track}" mouse "${timepoint}" "mouse-${timepoint}-${report_format}-fixed-terms" \
            "${source}" "$(absolute_path "${gaf_file}")" "$(absolute_path "${go_obo_file}")" \
            "Historical MGI ${report_format} capture ${timepoint}; fixed MP lethal terms; operational other-phenotype negative class" \
            "${go_annotation_release}" "${go_ontology_release}" "" "${retrieval_date}" \
            gene_sets explicit_only exclude exclude gene_sets \
            "${labels}/lethal.tsv|${labels}/operational_nonlethal.tsv"
    done < <(emit_mgi_captures)
fi

if selected_track mouse_impc_assertions && selected_organism mouse; then
    while IFS=$'\x1f' read -r organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
        [[ "${organism}" == "mouse" ]] || continue
        selected_year "${year}" || continue
        source="$(find "${REPO_DIR}/data/mouse/phenotype_data/ebi/${year}" -maxdepth 1 -type f -name '*.csv.gz' | LC_ALL=C sort | head -1)"
        if [[ ! -f "${source}" ]]; then
            skip_item mouse_impc_assertions mouse "${year}" missing_input "${source}"
            printf 'Required selected IMPC input is missing for %s: %s\n' "${year}" "${source}" >&2
            exit 1
        fi
        assert_audit_eligible "${source}"
        labels="${OUTPUT_DIR}/01_derived_labels/mouse_impc_assertions/mouse/${year}"
        mkdir -p "${labels}"
        if [[ ! (${RESUME} -eq 1 && -s "${labels}/lethal.tsv" && -s "${labels}/operational_nonlethal.tsv") ]]; then
            run_cmd "${PYTHON}" "${INPUT_HELPER}" impc-labels --input "${source}" \
                --lethal-terms "${MOUSE_LETHAL_TERMS}" --lethal-output "${labels}/lethal.tsv" \
                --viable-output "${labels}/operational_nonlethal.tsv" \
                --excluded-output "${labels}/excluded.tsv" --summary "${labels}/summary.json"
        fi
        build_arff mouse_impc_assertions mouse "${year}" "mouse-${year}-impc-assertions" \
            "${source}" "$(absolute_path "${gaf_file}")" "$(absolute_path "${go_obo_file}")" \
            "IMPC assertions nominal ${year}; operational other-phenotype negative class" \
            "${go_annotation_release}" "${go_ontology_release}" "" "${retrieval_date}" \
            gene_sets explicit_only exclude exclude gene_sets \
            "${labels}/lethal.tsv|${labels}/operational_nonlethal.tsv"
    done < <(emit_ledger_rows primary)
fi

# ---------------------------------------------------------------------------
# Single-snapshot ML (unchanged from V1, uses the primary track's ARFFs)
# ---------------------------------------------------------------------------
if [[ ${RUN_SINGLE_SNAPSHOT_ML} -eq 1 ]]; then
    read -r -a single_model_args <<< "${SINGLE_SNAPSHOT_MODELS}"
    read -r -a hidden_unit_args <<< "${NN_HIDDEN_UNITS}"
    while IFS=$'\x1f' read -r include_mode cohort organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
        selected_organism "${organism}" || continue
        selected_single_snapshot_year "${year}" || continue
        arff="${OUTPUT_DIR}/02_arff/primary/${organism}/${year}/${organism}_PhenGO.arff"
        [[ -f "${arff}" ]] || continue
        predict_output="${OUTPUT_DIR}/04_single_snapshot_ml/${organism}/${year}/Predict"
        mkdir -p "${predict_output}"
        if [[ ${RESUME} -eq 1 && -s "${predict_output}/PhenGO_Predict_params.txt" ]]; then continue; fi
        run_cmd "${PYTHON}" -m PhenGO.predict.cli \
            -arff_file "${arff}" -output_dir "${predict_output}" \
            -model "${single_model_args[@]}" -epochs "${NN_EPOCHS}" \
            -batch_size "${NN_BATCH_SIZE}" -hidden_units "${hidden_unit_args[@]}" \
            -dropout "${NN_DROPOUT}" -n_estimators "${N_ESTIMATORS}" \
            -calibration "${CALIBRATION}" -calibration_cv "${CALIBRATION_CV}" \
            -n_jobs "${N_JOBS}" -test_size "${SINGLE_TEST_SIZE}" \
            -perm_repeats "${PERM_REPEATS}" -seed "${SEED}" -overwrite
    done < <(emit_ledger_rows primary)

    for organism in ${ORGANISMS}; do
        predict_parent="${OUTPUT_DIR}/04_single_snapshot_ml/${organism}"
        [[ -d "${predict_parent}" ]] || continue
        report_output="${OUTPUT_DIR}/07_reports/single_snapshot/${organism}"
        mkdir -p "${report_output}"
        run_cmd "${PYTHON}" -m PhenGO.scripts.predict_report \
            -input_dir "${predict_parent}" -output_dir "${report_output}" --top-n "${TOP_K}"
    done
fi

# ---------------------------------------------------------------------------
# Single modeling pass across every selected track (was two separate loops,
# one per script, previously)
# ---------------------------------------------------------------------------
read -r -a model_args <<< "${MODELS}"
read -r -a panel_args <<< "${PANELS}"
read -r -a mlp_hidden_unit_args <<< "${MLP_HIDDEN_UNITS}"
mlp_early_stopping_args=( -nn_early_stopping )
[[ "${MLP_EARLY_STOPPING}" == "0" ]] && mlp_early_stopping_args=( -no_nn_early_stopping )

for track in ${TRACKS}; do
    for organism in ${ORGANISMS}; do
        arff_parent="${OUTPUT_DIR}/02_arff/${track}/${organism}"
        [[ -d "${arff_parent}" ]] || continue
        dataset_count="$(find "${arff_parent}" -type f -name '*_PhenGO.arff' | wc -l | tr -d ' ')"
        if [[ ${DRY_RUN} -eq 0 && "${dataset_count}" -lt 2 ]]; then
            skip_item "${track}" "${organism}" all insufficient_snapshots "${dataset_count} ARFF(s)"
            continue
        fi
        ml_output="${OUTPUT_DIR}/04_ml/${track}/${organism}"
        temporal_output="${OUTPUT_DIR}/05_temporal/${track}/${organism}"
        mkdir -p "${ml_output}" "${temporal_output}" "${OUTPUT_DIR}/07_reports/${track}"
        if [[ ! (${RESUME} -eq 1 && -s "${ml_output}/version_sensitivity_manifest.json") ]]; then
            run_cmd "${PYTHON}" -m PhenGO.predict.version_sensitivity \
                -input_dir "${arff_parent}" -output_dir "${ml_output}" \
                -models "${model_args[@]}" -panels "${panel_args[@]}" \
                -cv_folds "${CV_FOLDS}" -cv_repeats "${CV_REPEATS}" \
                -calibration "${CALIBRATION}" -calibration_cv "${CALIBRATION_CV}" \
                -n_estimators "${N_ESTIMATORS}" -n_jobs "${N_JOBS}" \
                -nn_hidden_units "${mlp_hidden_unit_args[@]}" -nn_alpha "${MLP_ALPHA}" \
                -nn_batch_size "${MLP_BATCH_SIZE}" -nn_learning_rate_init "${MLP_LEARNING_RATE}" \
                -nn_max_iter "${MLP_MAX_ITER}" "${mlp_early_stopping_args[@]}" \
                -nn_validation_fraction "${MLP_VALIDATION_FRACTION}" \
                -nn_n_iter_no_change "${MLP_PATIENCE}" \
                -importance_repeats "${IMPORTANCE_REPEATS}" -top_k "${TOP_K}" \
                -seed "${SEED}" -overwrite
        fi
        if [[ ! (${RESUME} -eq 1 && -s "${temporal_output}/${organism}_${track}_summary_timeline.csv") ]]; then
            run_cmd "${PYTHON}" -m PhenGO.scripts.GO_temporal_analysis \
                -input_dir "${arff_parent}" -output_dir "${temporal_output}" \
                -o "${organism}_${track}" --top-n "${TOP_K}" --max-fdr 0.05 -overwrite
        fi
        run_cmd "${PYTHON}" -m PhenGO.scripts.report_generator \
            -i "${ml_output}" -o "${OUTPUT_DIR}/07_reports/${track}/${organism}_ml_report.html"
    done
done

run_cmd "${PYTHON}" -m PhenGO.scripts.publication_summary \
    --run-root "${OUTPUT_DIR}" --output-dir "${OUTPUT_DIR}/06_publication_tables"
run_cmd "${PYTHON}" -m PhenGO.scripts.report_generator \
    -i "${OUTPUT_DIR}/06_publication_tables" \
    -o "${OUTPUT_DIR}/07_reports/publication_tables_report.html"

if [[ ${DRY_RUN} -eq 0 ]]; then
    : > "${OUTPUT_DIR}/00_run_metadata/output_checksums.sha256"
    find "${OUTPUT_DIR}/01_derived_inputs" "${OUTPUT_DIR}/01_derived_labels" "${OUTPUT_DIR}/02_arff" \
        "${OUTPUT_DIR}/03_qc" "${OUTPUT_DIR}/04_ml" "${OUTPUT_DIR}/04_single_snapshot_ml" \
        "${OUTPUT_DIR}/05_temporal" "${OUTPUT_DIR}/06_publication_tables" \
        -type f -print 2>/dev/null | LC_ALL=C sort | while IFS= read -r path; do
            sha256 "${path}"
        done >> "${OUTPUT_DIR}/00_run_metadata/output_checksums.sha256"
    printf 'completed_utc\t%s\narchitecture\t%s\n' "$(timestamp)" "${ARCHITECTURE}" \
        > "${OUTPUT_DIR}/00_run_metadata/run.complete"
fi

if [[ ${DRY_RUN} -eq 1 ]]; then
    printf '\nPhenGO unified publication dry-run plan complete.\nPlanned output: %s\n' "${OUTPUT_DIR}"
else
    printf '\nPhenGO unified publication workflow complete.\nOutput: %s\n' "${OUTPUT_DIR}"
fi
