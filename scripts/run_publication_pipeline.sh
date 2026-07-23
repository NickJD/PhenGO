#!/usr/bin/env bash
set -Eeuo pipefail

# One-shot, publication-oriented PhenGO workflow for Apple Silicon.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_DIR}"
PYTHON_DEFAULT="/Users/nicholas/miniconda3/envs/PhenGO/bin/python"
OUTPUT_DEFAULT="/Users/nicholas/Nextcloud/Current_Work/PhenGO/New_Outputs"

PYTHON="${PHENGO_PYTHON:-${PYTHON_DEFAULT}}"
OUTPUT_DIR="${PHENGO_OUTPUT_DIR:-${OUTPUT_DEFAULT}}"
LEDGER="${PHENGO_LEDGER:-${REPO_DIR}/data/publication_snapshots.tsv}"
TRACKS="${PHENGO_TRACKS:-primary is_a_only no_iea fixed_2025 fixed_2025_is_a}"
ORGANISMS="${PHENGO_ORGANISMS:-fish fly mouse worm yeast}"
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

usage() {
    cat <<'EOF'
Usage: scripts/run_publication_pipeline.sh [options]

Options:
  --output DIR          Output root (default: New_Outputs in Nextcloud)
  --ledger FILE         Snapshot ledger (default: data/publication_snapshots.tsv)
  --tracks "LIST"       Space-separated tracks to run
  --organisms "LIST"    Space-separated organisms to run
  --resume              Skip completed units in an existing output directory
  --quick               Two-fold smoke analysis using LR and primary/fixed tracks
  --dry-run             Validate and print commands without running computations
  --skip-tests          Do not run the repository test suite during preflight
  --require-anchors     Fail if any annotation-anchor input is absent (default)
  --allow-missing-anchors
                        Record and skip absent annotation anchors
  -h, --help            Show this help

Publication defaults can also be overridden with PHENGO_* environment variables.
The full run is computationally intensive; --resume is safe after interruption.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --ledger) LEDGER="$2"; shift 2 ;;
        --tracks) TRACKS="$2"; shift 2 ;;
        --organisms) ORGANISMS="$2"; shift 2 ;;
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
    TRACKS="primary fixed_2025"
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

if [[ -z "${TRACKS// }" || -z "${ORGANISMS// }" ]]; then
    printf 'Track and organism selections must not be empty.\n' >&2
    exit 2
fi
for track in ${TRACKS}; do
    case "${track}" in
        primary|is_a_only|no_iea|fixed_2025|fixed_2025_is_a) ;;
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

if [[ ! -x "${PYTHON}" ]]; then
    printf 'PhenGO Python is not executable: %s\n' "${PYTHON}" >&2
    exit 1
fi
if [[ ! -f "${LEDGER}" ]]; then
    printf 'Snapshot ledger does not exist: %s\n' "${LEDGER}" >&2
    exit 1
fi
if [[ ! "${MLP_EARLY_STOPPING}" =~ ^[01]$ ]]; then
    printf 'PHENGO_MLP_EARLY_STOPPING must be 0 or 1.\n' >&2
    exit 1
fi

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
if [[ -d "${OUTPUT_DIR}" ]] && [[ -n "$(ls -A "${OUTPUT_DIR}" 2>/dev/null)" ]] && [[ ${RESUME} -ne 1 ]] && [[ ${DRY_RUN} -ne 1 ]] && [[ ! -f "${DRY_RUN_MARKER}" ]] && [[ ${SAFE_PREFLIGHT_RESTART} -ne 1 ]]; then
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

quote_command() {
    local item
    for item in "$@"; do printf '%q ' "${item}"; done
}

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

selected_organism() {
    [[ " ${ORGANISMS} " == *" $1 "* ]]
}

selected_track() {
    [[ " ${TRACKS} " == *" $1 "* ]]
}

selected_single_snapshot_year() {
    [[ "${SINGLE_SNAPSHOT_YEARS}" == "all" ]] ||
        [[ " ${SINGLE_SNAPSHOT_YEARS} " == *" $1 "* ]]
}

absolute_path() {
    if [[ -z "$1" ]]; then
        printf ''
    elif [[ "$1" = /* ]]; then
        printf '%s' "$1"
    else
        printf '%s/%s' "${REPO_DIR}" "$1"
    fi
}

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

required_file_missing=0
: > "${INPUT_LIST}"
while IFS=$'\x1f' read -r include_mode cohort organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
    selected_organism "${organism}" || continue
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

architecture="$("${PYTHON}" -c 'import platform; print(platform.machine())')"
if [[ "${architecture}" != "arm64" ]]; then
    printf 'Expected an ARM64 Python environment, observed %s from %s\n' "${architecture}" "${PYTHON}" >&2
    exit 1
fi

cat > "${OUTPUT_DIR}/00_run_metadata/run_config.env" <<EOF
PHENGO_PYTHON=${PYTHON}
PHENGO_OUTPUT_DIR=${OUTPUT_DIR}
PHENGO_LEDGER=${LEDGER}
PHENGO_TRACKS=${TRACKS}
PHENGO_ORGANISMS=${ORGANISMS}
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
EOF

{
    printf 'started_utc\t%s\n' "$(timestamp)"
    printf 'uname\t'; uname -a
    printf 'python\t'; "${PYTHON}" -VV
    printf 'python_machine\t%s\n' "${architecture}"
    printf 'phenGO_version\t'; "${PYTHON}" -c 'from PhenGO.constants import PhenGO_VERSION; print(PhenGO_VERSION)'
    if command -v sw_vers >/dev/null 2>&1; then sw_vers; fi
    if command -v sysctl >/dev/null 2>&1; then
        sysctl -n machdep.cpu.brand_string 2>/dev/null || true
        sysctl -n hw.memsize 2>/dev/null || true
        sysctl -n hw.ncpu 2>/dev/null || true
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
        shasum -a 256 "${path}" >> "${OUTPUT_DIR}/00_run_metadata/input_checksums.sha256"
    done < "${INPUT_LIST}"
fi

SOURCE_TREE_HASH="$("${PYTHON}" -c "from PhenGO.provenance import source_tree_sha256; print(source_tree_sha256('${REPO_DIR}'))")"
DEPENDENCY_STATE="$("${PYTHON}" -c 'import json; from PhenGO.provenance import dependency_versions; print(json.dumps(dependency_versions(), sort_keys=True))')"
PIPELINE_HASH="$(shasum -a 256 "${BASH_SOURCE[0]}" "${LEDGER}" | shasum -a 256 | awk '{print $1}')"
INPUT_HASH="not-computed"
if [[ -s "${OUTPUT_DIR}/00_run_metadata/input_checksums.sha256" ]]; then
    INPUT_HASH="$(shasum -a 256 "${OUTPUT_DIR}/00_run_metadata/input_checksums.sha256" | awk '{print $1}')"
fi
RUN_FINGERPRINT="$(printf '%s\n' \
    "${SOURCE_TREE_HASH}" "${DEPENDENCY_STATE}" "${PIPELINE_HASH}" "${INPUT_HASH}" "${PYTHON}" \
    "${TRACKS}" "${ORGANISMS}" "${MODELS}" "${PANELS}" \
    "${SINGLE_SNAPSHOT_MODELS}" "${SINGLE_SNAPSHOT_YEARS}" \
    "${CV_FOLDS}" "${CV_REPEATS}" "${CALIBRATION}" "${CALIBRATION_CV}" \
    "${N_ESTIMATORS}" "${N_JOBS}" "${MLP_HIDDEN_UNITS}" "${MLP_ALPHA}" \
    "${MLP_BATCH_SIZE}" "${MLP_LEARNING_RATE}" "${MLP_MAX_ITER}" \
    "${MLP_EARLY_STOPPING}" "${MLP_VALIDATION_FRACTION}" "${MLP_PATIENCE}" \
    "${IMPORTANCE_REPEATS}" "${TOP_K}" "${SEED}" "${MIN_GO_GENE_COUNT}" \
    "${MAX_GO_GENE_FRACTION}" "${RUN_SINGLE_SNAPSHOT_ML}" "${NN_EPOCHS}" \
    "${NN_BATCH_SIZE}" "${NN_HIDDEN_UNITS}" "${NN_DROPOUT}" \
    "${SINGLE_TEST_SIZE}" "${PERM_REPEATS}" | shasum -a 256 | awk '{print $1}')"
if [[ ${RESUME} -eq 1 && "${PREVIOUS_RUN_FINGERPRINT}" != "${RUN_FINGERPRINT}" ]]; then
    printf 'Resume fingerprint changed; clearing and regenerating all scientific outputs.\n' | tee -a "${RUN_LOG}"
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

run_cmd "${PYTHON}" -c 'import numpy, pandas, scipy, sklearn, networkx, tensorflow; print("scientific Python imports: OK")'
if [[ ${RUN_TESTS} -eq 1 ]]; then
    run_cmd env PYTHONDONTWRITEBYTECODE=1 "${PYTHON}" -m pytest -q
fi

ensure_mouse_sets() {
    local year="$1" source="$2" outdir="${OUTPUT_DIR}/01_derived_labels/mouse/${year}"
    mkdir -p "${outdir}"
    if [[ ${RESUME} -eq 1 && -s "${outdir}/lethal.tsv" && -s "${outdir}/viable.tsv" ]]; then return; fi
    run_capture "${outdir}/conversion_summary.json" \
        "${PYTHON}" -m PhenGO.scripts.publication_inputs mouse-impc \
        --input "${source}" --lethal-output "${outdir}/lethal.tsv" \
        --viable-output "${outdir}/viable.tsv" --excluded-output "${outdir}/excluded.tsv"
}

ensure_fly_assignments() {
    local year="$1" gaf="$2" outdir="${OUTPUT_DIR}/01_derived_labels/fly/${year}"
    mkdir -p "${outdir}"
    if [[ ${RESUME} -eq 1 && -s "${outdir}/assignments.tsv" ]]; then return; fi
    run_capture "${outdir}/conversion_summary.json" \
        "${PYTHON}" -m PhenGO.scripts.publication_inputs fly-assignments \
        --gaf "${gaf}" --output "${outdir}/assignments.tsv" \
        --excluded-output "${outdir}/excluded.tsv" --taxon 7227
}

ensure_worm_terms() {
    local year="$1" obo="$2" outdir="${OUTPUT_DIR}/01_derived_labels/worm/${year}"
    mkdir -p "${outdir}"
    if [[ ${RESUME} -eq 1 && -s "${outdir}/lethal_terms.tsv" ]]; then return; fi
    run_cmd "${PYTHON}" -m PhenGO.scripts.get_phenotype_terms \
        --obo-file "${obo}" \
        --term-list "${REPO_DIR}/data/worm/lethal_terms/root_lethal_phenotype_terms.txt" \
        --results "${outdir}/lethal_terms.tsv"
}

while IFS=$'\x1f' read -r include_mode cohort organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
    selected_organism "${organism}" || continue
    gaf="$(absolute_path "${gaf_file}")"
    case "${organism}" in
        mouse) ensure_mouse_sets "${year}" "$(absolute_path "${viability_file}")" ;;
        fly) ensure_fly_assignments "${year}" "${gaf}" ;;
        worm) ensure_worm_terms "${year}" "$(absolute_path "${phenotype_ontology_file}")" ;;
    esac
done < <(emit_ledger_rows primary)

build_arff() {
    local track="$1" organism="$2" year="$3" snapshot_id="$4"
    local phenotype_file="$5" viability_file="$6" gaf="$7" go_obo="$8" phenotype_ontology_file="$9"
    shift 9
    local phenotype_release="$1" go_annotation_release="$2" go_ontology_release="$3" phenotype_ontology_release="$4"
    local label_source="$5" nonlethal_policy="$6" mixed_label_policy="$7" ambiguous_label_policy="$8" retrieval_date="$9"
    local output="${OUTPUT_DIR}/02_arff/${track}/${organism}/${year}"
    local qc="${OUTPUT_DIR}/03_qc/${track}/${organism}/${year}"
    local arff="${output}/${organism}_PhenGO.arff"
    if [[ ${RESUME} -eq 1 && -s "${arff}" && -s "${output}/PhenGO_manifest.json" ]]; then return; fi
    mkdir -p "${output}" "${qc}"

    local command=(
        "${PYTHON}" -m PhenGO.core.PhenGO
        -species "${organism}"
        -gene_association_file "${gaf}"
        -go_obo_file "${go_obo}"
        -output_dir "${output}"
        -snapshot_id "${snapshot_id}"
        -phenotype_release "${phenotype_release}"
        -go_annotation_release "${go_annotation_release}"
        -go_ontology_release "${go_ontology_release}"
        -retrieval_date "${retrieval_date}"
        -snapshot_semantics declared_composite_cross_section
        -phenotype_availability available_local
        -go_annotation_availability available_local_go_archive
        -go_ontology_availability available_local_annual_GO_file
        -retrieval_route local_archived_or_provider_file
        -label_source "${label_source}"
        -nonlethal_policy "${nonlethal_policy}"
        -mixed_label_policy "${mixed_label_policy}"
        -ambiguous_label_policy "${ambiguous_label_policy}"
        -min_go_gene_count "${MIN_GO_GENE_COUNT}"
        -max_go_gene_fraction "${MAX_GO_GENE_FRACTION}"
        -strict_snapshot -overwrite
    )
    if [[ -n "${phenotype_file}" ]]; then command+=( -phenotype_file "${phenotype_file}" ); fi
    if [[ -n "${phenotype_ontology_release}" ]]; then
        command+=( -phenotype_ontology_release "${phenotype_ontology_release}" )
    fi

    if [[ "${label_source}" == "gene_sets" ]]; then
        if [[ "${track}" == fixed_2025* ]]; then
            command+=(
                -lethal_gene_set "${OUTPUT_DIR}/01_derived_labels/fixed_2025/${organism}/lethal.tsv"
                -viable_gene_set "${OUTPUT_DIR}/01_derived_labels/fixed_2025/${organism}/viable.tsv"
            )
        elif [[ "${organism}" == "mouse" ]]; then
            command+=(
                -lethal_gene_set "${OUTPUT_DIR}/01_derived_labels/mouse/${year}/lethal.tsv"
                -viable_gene_set "${OUTPUT_DIR}/01_derived_labels/mouse/${year}/viable.tsv"
            )
        fi
    fi
    case "${organism}" in
        fly)
            ensure_fly_assignments "${year}" "${gaf}"
            command+=( -fly_assignments "${OUTPUT_DIR}/01_derived_labels/fly/${year}/assignments.tsv" )
            ;;
        worm)
            if [[ "${label_source}" == "release_records" ]]; then
                command+=( -worm_phenotypes "${OUTPUT_DIR}/01_derived_labels/worm/${year}/lethal_terms.tsv" )
            fi
            ;;
    esac

    case "${track}" in
        is_a_only|fixed_2025_is_a) command+=( -go_relations is_a ) ;;
        no_iea) command+=( -go_relations is_a part_of -exclude_go_evidence_codes IEA ) ;;
        *) command+=( -go_relations is_a part_of ) ;;
    esac
    run_cmd "${command[@]}"
    run_cmd "${PYTHON}" -m PhenGO.scripts.arff_validator \
        -i "${arff}" -o "${qc}" --fail-on-error
}

for track in primary is_a_only no_iea; do
    selected_track "${track}" || continue
    while IFS=$'\x1f' read -r include_mode cohort organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
        selected_organism "${organism}" || continue
        track_snapshot="${snapshot_id}"
        if [[ "${track}" != "primary" ]]; then track_snapshot="${snapshot_id}-${track}"; fi
        build_arff "${track}" "${organism}" "${year}" "${track_snapshot}" \
            "$(absolute_path "${phenotype_file}")" "$(absolute_path "${viability_file}")" \
            "$(absolute_path "${gaf_file}")" "$(absolute_path "${go_obo_file}")" \
            "$(absolute_path "${phenotype_ontology_file}")" "${phenotype_release}" \
            "${go_annotation_release}" "${go_ontology_release}" "${phenotype_ontology_release}" \
            "${label_source}" "${nonlethal_policy}" "${mixed_label_policy}" \
            "${ambiguous_label_policy}" "${retrieval_date}"
    done < <(emit_ledger_rows primary)
done

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
            build_arff "${fixed_track}" "${organism}" "${year}" "${organism}-${year}-${fixed_track}" \
                "" "" "$(absolute_path "${gaf_file}")" "$(absolute_path "${go_obo_file}")" "" \
                "Fixed 2025 PhenGO labels for ${organism}" "${go_annotation_release}" \
                "${go_ontology_release}" "${phenotype_ontology_release}" gene_sets \
                explicit_only exclude exclude "${retrieval_date}"
        done < <(emit_ledger_rows primary)

        while IFS=$'\x1f' read -r include_mode cohort organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
            selected_organism "${organism}" || continue
            if [[ "${fixed_track}" == "fixed_2025" && "${year}" -lt 2014 ]]; then continue; fi
            gaf="$(absolute_path "${gaf_file}")"
            go_obo="$(absolute_path "${go_obo_file}")"
            if [[ ! -f "${gaf}" || ! -f "${go_obo}" ]]; then continue; fi
            build_arff "${fixed_track}" "${organism}" "${year}" "${snapshot_id}-${fixed_track}" \
                "" "" "${gaf}" "${go_obo}" "" "${phenotype_release}" \
                "${go_annotation_release}" "${go_ontology_release}" "${phenotype_ontology_release}" \
                gene_sets explicit_only exclude exclude "${retrieval_date}"
        done < <(emit_ledger_rows annotation_anchor)
    done
fi

if [[ ${RUN_SINGLE_SNAPSHOT_ML} -eq 1 ]]; then
    read -r -a single_model_args <<< "${SINGLE_SNAPSHOT_MODELS}"
    read -r -a hidden_unit_args <<< "${NN_HIDDEN_UNITS}"
    while IFS=$'\x1f' read -r include_mode cohort organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
        selected_organism "${organism}" || continue
        selected_single_snapshot_year "${year}" || continue
        arff="${OUTPUT_DIR}/02_arff/primary/${organism}/${year}/${organism}_PhenGO.arff"
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

read -r -a model_args <<< "${MODELS}"
read -r -a panel_args <<< "${PANELS}"
read -r -a mlp_hidden_unit_args <<< "${MLP_HIDDEN_UNITS}"
mlp_early_stopping_args=( -nn_early_stopping )
if [[ "${MLP_EARLY_STOPPING}" == "0" ]]; then
    mlp_early_stopping_args=( -no_nn_early_stopping )
fi
for track in ${TRACKS}; do
    for organism in ${ORGANISMS}; do
        arff_parent="${OUTPUT_DIR}/02_arff/${track}/${organism}"
        [[ -d "${arff_parent}" ]] || continue
        dataset_count="$(find "${arff_parent}" -type f -name '*_PhenGO.arff' | wc -l | tr -d ' ')"
        if [[ ${DRY_RUN} -eq 0 && "${dataset_count}" -lt 2 ]]; then
            printf 'Skipping analysis for %s/%s: fewer than two ARFFs.\n' "${track}" "${organism}" | tee -a "${RUN_LOG}"
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
                -nn_hidden_units "${mlp_hidden_unit_args[@]}" \
                -nn_alpha "${MLP_ALPHA}" -nn_batch_size "${MLP_BATCH_SIZE}" \
                -nn_learning_rate_init "${MLP_LEARNING_RATE}" \
                -nn_max_iter "${MLP_MAX_ITER}" \
                "${mlp_early_stopping_args[@]}" \
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
    find "${OUTPUT_DIR}/01_derived_labels" "${OUTPUT_DIR}/02_arff" "${OUTPUT_DIR}/04_ml" "${OUTPUT_DIR}/04_single_snapshot_ml" "${OUTPUT_DIR}/05_temporal" "${OUTPUT_DIR}/06_publication_tables" \
        -type f -print | LC_ALL=C sort | while IFS= read -r path; do
            shasum -a 256 "${path}"
        done >> "${OUTPUT_DIR}/00_run_metadata/output_checksums.sha256"
    printf 'completed_utc\t%s\n' "$(timestamp)" > "${OUTPUT_DIR}/00_run_metadata/run.complete"
fi

if [[ ${DRY_RUN} -eq 1 ]]; then
    printf '\nPhenGO publication dry-run plan complete.\nPlanned output: %s\n' "${OUTPUT_DIR}"
else
    printf '\nPhenGO publication workflow complete.\nOutput: %s\n' "${OUTPUT_DIR}"
fi
