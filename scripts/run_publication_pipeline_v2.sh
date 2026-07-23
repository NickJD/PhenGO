#!/usr/bin/env bash
set -Eeuo pipefail

# Master publication workflow. It builds, resumes, or integrity-rebuilds V1,
# adds V2 experiments, and creates a verified unified publication artifact
# rather than a set of external references.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_DIR}"
PYTHON_DEFAULT="/Users/nicholas/miniconda3/envs/PhenGO/bin/python"
#RUN_ROOT_DEFAULT="/Users/nicholas/Nextcloud/Current_Work/PhenGO/New_Outputs/publication_run_$(date '+%Y-%m-%d')"
RUN_ROOT_DEFAULT="/Users/nicholas/Nextcloud/Current_Work/PhenGO/New_Outputs/publication_run_2026-07-22"
BASE_RUN_DEFAULT="${RUN_ROOT_DEFAULT}"
OUTPUT_DEFAULT="${RUN_ROOT_DEFAULT}_v2_extension"
CONSOLIDATED_DEFAULT="${RUN_ROOT_DEFAULT}_consolidated"

PYTHON="${PHENGO_PYTHON:-${PYTHON_DEFAULT}}"
BASE_RUN="${PHENGO_BASE_RUN:-${BASE_RUN_DEFAULT}}"
OUTPUT_DIR="${PHENGO_V2_OUTPUT_DIR:-${OUTPUT_DEFAULT}}"
CONSOLIDATED_DIR="${PHENGO_CONSOLIDATED_OUTPUT_DIR:-${CONSOLIDATED_DEFAULT}}"
TRACKS="${PHENGO_V2_TRACKS:-primary_mod_gaf legacy_like_go_archive legacy_like_mod_gaf fail_closed_driver_aware_go_archive fail_closed_driver_aware_mod_gaf worm_fixed_2025_terms_go_archive worm_fixed_2025_terms_mod_gaf mouse_impc_assertions mouse_mgi_genepheno_fixed_terms mouse_mgi_phenogenomp_fixed_terms}"
ORGANISMS="${PHENGO_V2_ORGANISMS:-fish fly mouse worm yeast}"
YEARS="${PHENGO_V2_YEARS:-all}"
MODELS="${PHENGO_MODELS:-all}"
PANELS="${PHENGO_PANELS:-full matched_features matched_genes matched_both}"
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
VERIFY_BASE="${PHENGO_VERIFY_BASE:-1}"
RESUME=0
DRY_RUN=0
QUICK=0
EXTENSION_ONLY=0

BASE_PIPELINE="${SCRIPT_DIR}/run_publication_pipeline.sh"
INPUT_HELPER="${SCRIPT_DIR}/v2/alternative_inputs.py"
MULTIVERSE_HELPER="${SCRIPT_DIR}/v2/publication_multiverse.py"
MOUSE_LETHAL_TERMS="${REPO_DIR}/src/PhenGO/data/mouse/mouse_lethal_terms.txt.gz"
AVAILABILITY_LEDGER="${REPO_DIR}/data/resource_availability.tsv"

usage() {
    cat <<'EOF'
Usage: scripts/run_publication_pipeline_v2.sh [options]

Options:
  --run-root PREFIX       Set V1 to PREFIX and derive V2/consolidated suffixes
  --base-run DIR          Completed V1 publication run
  --output DIR            New V2 extension output root
  --consolidated-output DIR
                          Self-contained V1+V2 copy, tables, figures, and reports
  --tracks "LIST"         Space-separated V2 tracks
  --organisms "LIST"      Space-separated organism names
  --years "LIST"          Space-separated years, or all
  --resume                Resume an interrupted run with the same fingerprint
  --quick                 Two-year LR smoke run on representative tracks
  --dry-run               Print commands; permits an incomplete V1 base
  --extension-only        Require a completed V1 instead of building/resuming it
  --skip-tests            Skip focused V2 tests
  --skip-base-checksum    Do not verify every V1 output against its checksum ledger
  -h, --help              Show this help

Without --extension-only, an absent or incomplete V1 is built or resumed first.
If a completed V1 fails the scientific integrity audit, it is rebuilt cleanly
with the current code unless --extension-only was selected.
The full default runs every model family, including the comparable sklearn MLP,
on full, matched-feature, matched-gene, and matched-both panels. It is intended
to run as the sole publication entry point. The consolidated directory creates
independent copies of both completed runs and verifies every copied file.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root)
            BASE_RUN="$2"
            OUTPUT_DIR="${2}_v2_extension"
            CONSOLIDATED_DIR="${2}_consolidated"
            shift 2
            ;;
        --base-run) BASE_RUN="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --consolidated-output) CONSOLIDATED_DIR="$2"; shift 2 ;;
        --tracks) TRACKS="$2"; shift 2 ;;
        --organisms) ORGANISMS="$2"; shift 2 ;;
        --years) YEARS="$2"; shift 2 ;;
        --resume) RESUME=1; shift ;;
        --quick) QUICK=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        --extension-only) EXTENSION_ONLY=1; shift ;;
        --skip-tests) RUN_TESTS=0; shift ;;
        --skip-base-checksum) VERIFY_BASE=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ ${QUICK} -eq 1 ]]; then
    TRACKS="primary_mod_gaf legacy_like_go_archive mouse_impc_assertions"
    ORGANISMS="fly mouse worm"
    YEARS="2024 2025"
    MODELS="lr"
    PANELS="full matched_both"
    CV_FOLDS=2
    CV_REPEATS=2
    CALIBRATION=none
    N_ESTIMATORS=20
    IMPORTANCE_REPEATS=2
    TOP_K=5
    N_JOBS=2
fi

if [[ -z "${TRACKS// }" || -z "${ORGANISMS// }" || -z "${YEARS// }" || -z "${MODELS// }" || -z "${PANELS// }" ]]; then
    printf 'Track, organism, year, model, and panel selections must not be empty.\n' >&2
    exit 2
fi
for track in ${TRACKS}; do
    case "${track}" in
        primary_mod_gaf|legacy_like_go_archive|legacy_like_mod_gaf|\
        fail_closed_driver_aware_go_archive|fail_closed_driver_aware_mod_gaf|\
        worm_fixed_2025_terms_go_archive|worm_fixed_2025_terms_mod_gaf|\
        mouse_impc_assertions|mouse_mgi_genepheno_fixed_terms|\
        mouse_mgi_phenogenomp_fixed_terms) ;;
        *) printf 'Unknown V2 track: %s\n' "${track}" >&2; exit 2 ;;
    esac
done
for organism in ${ORGANISMS}; do
    case "${organism}" in
        fish|fly|mouse|worm|yeast) ;;
        *) printf 'Unknown organism: %s\n' "${organism}" >&2; exit 2 ;;
    esac
done
for model in ${MODELS}; do
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

for required in "${PYTHON}" "${BASE_PIPELINE}" "${INPUT_HELPER}" "${MULTIVERSE_HELPER}" "${MOUSE_LETHAL_TERMS}" "${AVAILABILITY_LEDGER}"; do
    if [[ ! -e "${required}" ]]; then
        printf 'Required V2 dependency is missing: %s\n' "${required}" >&2
        exit 1
    fi
done
if [[ ! -x "${PYTHON}" ]]; then
    printf 'PhenGO Python is not executable: %s\n' "${PYTHON}" >&2
    exit 1
fi
if [[ "$("${PYTHON}" -c 'import platform; print(platform.machine())')" != "arm64" ]]; then
    printf 'V2 requires the ARM64 PhenGO environment: %s\n' "${PYTHON}" >&2
    exit 1
fi
if [[ "${OUTPUT_DIR}" == "${BASE_RUN}" || "${CONSOLIDATED_DIR}" == "${BASE_RUN}" || "${OUTPUT_DIR}" == "${CONSOLIDATED_DIR}" ]]; then
    printf 'Base, extension, and consolidated directories must be distinct.\n' >&2
    exit 1
fi

DRY_RUN_MARKER="${OUTPUT_DIR}/00_run_metadata/dry_run.marker"
if [[ -d "${OUTPUT_DIR}" && -n "$(ls -A "${OUTPUT_DIR}" 2>/dev/null)" && ${RESUME} -ne 1 && ${DRY_RUN} -ne 1 && ! -f "${DRY_RUN_MARKER}" ]]; then
    printf 'Extension output is not empty: %s\nUse --resume or a new directory.\n' "${OUTPUT_DIR}" >&2
    exit 1
fi
if [[ -d "${CONSOLIDATED_DIR}" && -n "$(ls -A "${CONSOLIDATED_DIR}" 2>/dev/null)" && ${RESUME} -ne 1 && ${DRY_RUN} -ne 1 ]]; then
    printf 'Consolidated output is not empty: %s\nUse --resume or a new directory.\n' "${CONSOLIDATED_DIR}" >&2
    exit 1
fi

BASE_COMPLETE="${BASE_RUN}/00_run_metadata/run.complete"
CURRENT_SOURCE_TREE_HASH="$("${PYTHON}" -c "from PhenGO.provenance import source_tree_sha256; print(source_tree_sha256('${REPO_DIR}'))")"
prepare_base_run() {
    local base_lock base_pid base_log base_log_mtime now_epoch process_command
    local process_discovery_available=0
    local -a base_args
    base_lock="${BASE_RUN}/00_run_metadata/publication_pipeline.lock"
    if [[ -s "${base_lock}" ]]; then
        base_pid="$(awk -F'\t' '$1 == "pid" {print $2}' "${base_lock}" 2>/dev/null || true)"
        if [[ "${base_pid}" =~ ^[0-9]+$ ]] && kill -0 "${base_pid}" 2>/dev/null; then
            printf 'V1 is already running with PID %s: %s\n' "${base_pid}" "${BASE_RUN}" >&2
            return 1
        fi
    fi
    if command -v pgrep >/dev/null 2>&1 && command -v ps >/dev/null 2>&1; then
        process_discovery_available=1
        while IFS= read -r base_pid; do
            [[ "${base_pid}" =~ ^[0-9]+$ ]] || continue
            [[ "${base_pid}" == "$$" ]] && continue
            process_command="$(ps -p "${base_pid}" -o command= 2>/dev/null || true)"
            if [[ "${process_command}" == *"run_publication_pipeline.sh"* ]] && \
                    { [[ "${process_command}" == *"--output ${BASE_RUN}"* ]] || \
                      [[ "${process_command}" == *"--output=${BASE_RUN}"* ]]; }; then
                printf 'V1 has an existing running or suspended process with PID %s: %s\n' \
                    "${base_pid}" "${BASE_RUN}" >&2
                return 1
            fi
        done < <(pgrep -f 'run_publication_pipeline\.sh' 2>/dev/null || true)
    fi
    if [[ ! -f "${BASE_COMPLETE}" && ${process_discovery_available} -eq 0 ]]; then
        base_log="${BASE_RUN}/logs/publication_pipeline.log"
        if [[ -f "${base_log}" ]]; then
            base_log_mtime="$(stat -f '%m' "${base_log}" 2>/dev/null || stat -c '%Y' "${base_log}" 2>/dev/null || printf '0')"
            now_epoch="$(date +%s)"
            if [[ "${base_log_mtime}" =~ ^[0-9]+$ ]] && (( now_epoch - base_log_mtime < 300 )); then
                printf 'V1 appears active because its log changed within five minutes: %s\n' "${base_log}" >&2
                return 1
            fi
        fi
    fi

    base_args=( "${BASE_PIPELINE}" --output "${BASE_RUN}" )
    if [[ -f "${BASE_RUN}/00_run_metadata/run_fingerprint.sha256" && ! -f "${BASE_RUN}/00_run_metadata/dry_run.marker" ]]; then
        base_args+=( --resume )
    fi
    [[ ${QUICK} -eq 1 ]] && base_args+=( --quick )
    [[ ${RUN_TESTS} -eq 0 ]] && base_args+=( --skip-tests )
    printf 'Preparing V1 before V2: '
    printf '%q ' "${base_args[@]}"
    printf '\n'
    "${base_args[@]}"
    if [[ ! -f "${BASE_COMPLETE}" ]]; then
        printf 'V1 returned without its completion marker: %s\n' "${BASE_COMPLETE}" >&2
        return 1
    fi
}

if [[ ! -f "${BASE_COMPLETE}" ]]; then
    if [[ ${DRY_RUN} -eq 1 ]]; then
        printf 'DRY RUN: V1 completion marker is absent; V1 would be prepared first: %s\n' "${BASE_COMPLETE}" >&2
        base_dry_args=( "${BASE_PIPELINE}" --output "${BASE_RUN}" --dry-run )
        if [[ -f "${BASE_RUN}/00_run_metadata/run_fingerprint.sha256" && ! -f "${BASE_RUN}/00_run_metadata/dry_run.marker" ]]; then
            base_dry_args+=( --resume )
        fi
        [[ ${QUICK} -eq 1 ]] && base_dry_args+=( --quick )
        [[ ${RUN_TESTS} -eq 0 ]] && base_dry_args+=( --skip-tests )
        printf 'DRY RUN V1 command: '
        printf '%q ' "${base_dry_args[@]}"
        printf '\n'
        base_scientific_output="$({ find \
            "${BASE_RUN}/01_derived_labels" "${BASE_RUN}/02_arff" \
            "${BASE_RUN}/03_qc" "${BASE_RUN}/04_ml" \
            "${BASE_RUN}/04_single_snapshot_ml" "${BASE_RUN}/05_temporal" \
            "${BASE_RUN}/06_publication_tables" "${BASE_RUN}/07_reports" \
            -type f -print -quit 2>/dev/null || true; })"
        if [[ -z "${base_scientific_output}" ]]; then
            "${base_dry_args[@]}"
        else
            printf 'DRY RUN: existing incomplete V1 has scientific outputs; its live tree is not modified.\n' >&2
        fi
    elif [[ ${EXTENSION_ONLY} -eq 1 ]]; then
        printf 'V1 is not complete and --extension-only was selected: %s\n' "${BASE_COMPLETE}" >&2
        exit 1
    else
        prepare_base_run
    fi
fi

BASE_INTEGRITY_REPORT=""
if [[ ${DRY_RUN} -eq 0 ]]; then
    BASE_INTEGRITY_REPORT="$(mktemp "${TMPDIR:-/tmp}/phengo_base_integrity.XXXXXX")"
    if ! "${PYTHON}" "${MULTIVERSE_HELPER}" --audit-run "${BASE_RUN}" \
            --audit-report "${BASE_INTEGRITY_REPORT}" \
            --expected-source-hash "${CURRENT_SOURCE_TREE_HASH}"; then
        if [[ ${EXTENSION_ONLY} -eq 1 ]]; then
            printf 'Completed V1 failed scientific integrity and --extension-only forbids rebuilding it.\n' >&2
            rm -f "${BASE_INTEGRITY_REPORT}"
            exit 1
        fi
        printf 'Completed V1 failed scientific integrity; rebuilding it cleanly with the current code.\n' >&2
        prepare_base_run
        if ! "${PYTHON}" "${MULTIVERSE_HELPER}" --audit-run "${BASE_RUN}" \
                --audit-report "${BASE_INTEGRITY_REPORT}" \
                --expected-source-hash "${CURRENT_SOURCE_TREE_HASH}"; then
            printf 'Rebuilt V1 still failed scientific integrity; V2 will not continue.\n' >&2
            rm -f "${BASE_INTEGRITY_REPORT}"
            exit 1
        fi
    fi
fi
LEDGER="${BASE_RUN}/00_run_metadata/publication_snapshots.tsv"
if [[ ! -f "${LEDGER}" ]]; then
    LEDGER="${REPO_DIR}/data/publication_snapshots.tsv"
fi
if [[ ! -f "${LEDGER}" ]]; then
    printf 'Publication snapshot ledger is missing.\n' >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"/{00_run_metadata,01_derived_inputs,01_derived_labels,02_arff,03_qc,04_ml,05_temporal,06_publication_tables,07_reports,logs}
if [[ -n "${BASE_INTEGRITY_REPORT}" ]]; then
    cp "${BASE_INTEGRITY_REPORT}" "${OUTPUT_DIR}/03_qc/base_run_integrity.json"
    rm -f "${BASE_INTEGRITY_REPORT}"
fi
if [[ ${DRY_RUN} -eq 1 ]]; then
    printf 'planning_only_utc\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" > "${DRY_RUN_MARKER}"
else
    rm -f "${DRY_RUN_MARKER}"
fi
# Keep the live-process lock outside the run tree so consolidation never copies
# transient state into the self-contained publication artifact.
LOCK_FILE="${OUTPUT_DIR}.publication_pipeline_v2.lock"
if [[ ${DRY_RUN} -eq 0 ]]; then
    if [[ -s "${LOCK_FILE}" ]]; then
        lock_pid="$(awk -F'\t' '$1 == "pid" {print $2}' "${LOCK_FILE}" 2>/dev/null || true)"
        if [[ "${lock_pid}" =~ ^[0-9]+$ ]] && kill -0 "${lock_pid}" 2>/dev/null; then
            printf 'V2 pipeline is already running with PID %s: %s\n' \
                "${lock_pid}" "${OUTPUT_DIR}" >&2
            exit 1
        fi
    fi
    printf 'pid\t%s\nstarted_utc\t%s\n' "$$" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" > "${LOCK_FILE}"
fi
cleanup_v2_pipeline() {
    status=$?
    rm -f "${LOCK_FILE}"
    if [[ ${status} -ne 0 ]]; then
        rm -f "${OUTPUT_DIR}/00_run_metadata/run.complete"
    fi
}
trap cleanup_v2_pipeline EXIT
RUN_LOG="${OUTPUT_DIR}/logs/publication_pipeline_v2.log"
COMMAND_LOG="${OUTPUT_DIR}/00_run_metadata/commands.tsv"
SKIP_LOG="${OUTPUT_DIR}/00_run_metadata/skipped_inputs.tsv"
INPUT_AUDIT="${OUTPUT_DIR}/03_qc/alternative_input_audit.tsv"
if [[ ${RESUME} -eq 0 ]]; then
    : > "${RUN_LOG}"
fi
if [[ ${RESUME} -eq 0 || ! -s "${COMMAND_LOG}" ]]; then
    printf 'timestamp_utc\tworking_directory\tcommand\n' > "${COMMAND_LOG}"
fi
if [[ ${RESUME} -eq 0 || ! -s "${SKIP_LOG}" ]]; then
    printf 'track\torganism\ttimepoint\treason\tdetail\n' > "${SKIP_LOG}"
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
selected_track() { [[ " ${TRACKS} " == *" $1 "* ]]; }
selected_organism() { [[ " ${ORGANISMS} " == *" $1 "* ]]; }
selected_year() { [[ "${YEARS}" == "all" || " ${YEARS} " == *" $1 "* ]]; }
absolute_path() { [[ "$1" = /* ]] && printf '%s' "$1" || printf '%s/%s' "${REPO_DIR}" "$1"; }
skip_item() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" >> "${SKIP_LOG}"; }

emit_ledger_rows() {
    "${PYTHON}" - "${LEDGER}" <<'PY'
import csv, sys
fields = [
    "organism", "year", "snapshot_id", "phenotype_file", "viability_file",
    "gaf_file", "go_obo_file", "phenotype_ontology_file", "phenotype_release",
    "go_annotation_release", "go_ontology_release", "phenotype_ontology_release",
    "label_source", "nonlethal_policy", "mixed_label_policy",
    "ambiguous_label_policy", "retrieval_date",
]
with open(sys.argv[1], encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row.get("cohort") == "primary":
            print("\x1f".join(row.get(field, "") for field in fields))
PY
}

if [[ ${VERIFY_BASE} -eq 1 && ${DRY_RUN} -eq 0 ]]; then
    BASE_CHECKSUMS="${BASE_RUN}/00_run_metadata/output_checksums.sha256"
    if [[ ! -s "${BASE_CHECKSUMS}" ]]; then
        printf 'V1 output checksum ledger is missing: %s\n' "${BASE_CHECKSUMS}" >&2
        exit 1
    fi
    run_cmd shasum -a 256 -c "${BASE_CHECKSUMS}"
fi

SOURCE_TREE_HASH="${CURRENT_SOURCE_TREE_HASH}"
BASE_FINGERPRINT="missing"
if [[ -f "${BASE_RUN}/00_run_metadata/run_fingerprint.sha256" ]]; then
    BASE_FINGERPRINT="$(<"${BASE_RUN}/00_run_metadata/run_fingerprint.sha256")"
fi
PIPELINE_HASH="$(shasum -a 256 "${BASH_SOURCE[0]}" "${INPUT_HELPER}" "${MULTIVERSE_HELPER}" "${LEDGER}" "${AVAILABILITY_LEDGER}" | shasum -a 256 | awk '{print $1}')"
RUN_FINGERPRINT="$(printf '%s\n' "${SOURCE_TREE_HASH}" "${BASE_FINGERPRINT}" "${PIPELINE_HASH}" "${TRACKS}" "${ORGANISMS}" "${YEARS}" "${MODELS}" "${PANELS}" "${CV_FOLDS}" "${CV_REPEATS}" "${CALIBRATION}" "${CALIBRATION_CV}" "${N_ESTIMATORS}" "${N_JOBS}" "${SEED}" | shasum -a 256 | awk '{print $1}')"
if [[ ${RESUME} -eq 1 && -f "${OUTPUT_DIR}/00_run_metadata/run_fingerprint.sha256" ]]; then
    PREVIOUS_FINGERPRINT="$(<"${OUTPUT_DIR}/00_run_metadata/run_fingerprint.sha256")"
    if [[ "${PREVIOUS_FINGERPRINT}" != "${RUN_FINGERPRINT}" ]]; then
        printf 'Resume fingerprint differs. Use a new extension output directory.\n' >&2
        exit 1
    fi
fi
printf '%s\n' "${RUN_FINGERPRINT}" > "${OUTPUT_DIR}/00_run_metadata/run_fingerprint.sha256"
cp "${LEDGER}" "${OUTPUT_DIR}/00_run_metadata/publication_snapshots.tsv"
cat > "${OUTPUT_DIR}/00_run_metadata/run_config.env" <<EOF
PHENGO_BASE_RUN=${BASE_RUN}
PHENGO_V2_OUTPUT_DIR=${OUTPUT_DIR}
PHENGO_CONSOLIDATED_OUTPUT_DIR=${CONSOLIDATED_DIR}
PHENGO_V2_TRACKS=${TRACKS}
PHENGO_V2_ORGANISMS=${ORGANISMS}
PHENGO_V2_YEARS=${YEARS}
PHENGO_MODELS=${MODELS}
PHENGO_PANELS=${PANELS}
PHENGO_CV_FOLDS=${CV_FOLDS}
PHENGO_CV_REPEATS=${CV_REPEATS}
PHENGO_CALIBRATION=${CALIBRATION}
PHENGO_CALIBRATION_CV=${CALIBRATION_CV}
PHENGO_N_ESTIMATORS=${N_ESTIMATORS}
PHENGO_N_JOBS=${N_JOBS}
PHENGO_SEED=${SEED}
EOF
{
    printf 'started_utc\t%s\n' "$(timestamp)"
    printf 'base_run\t%s\n' "${BASE_RUN}"
    printf 'base_fingerprint\t%s\n' "${BASE_FINGERPRINT}"
    printf 'base_completion_sha256\t'; [[ -f "${BASE_COMPLETE}" ]] && shasum -a 256 "${BASE_COMPLETE}" | awk '{print $1}' || printf 'not_complete\n'
    printf 'source_tree_sha256\t%s\n' "${SOURCE_TREE_HASH}"
    printf 'python\t'; "${PYTHON}" -VV
    printf 'python_machine\t'; "${PYTHON}" -c 'import platform; print(platform.machine())'
    printf 'git_commit\t'; git -C "${REPO_DIR}" rev-parse HEAD
} > "${OUTPUT_DIR}/00_run_metadata/environment.txt"
"${PYTHON}" -m pip freeze > "${OUTPUT_DIR}/00_run_metadata/python_packages.txt"
git -C "${REPO_DIR}" status --short > "${OUTPUT_DIR}/00_run_metadata/git_status.txt"
if [[ ${DRY_RUN} -eq 0 ]]; then rm -f "${OUTPUT_DIR}/00_run_metadata/run.complete"; fi

run_cmd "${PYTHON}" -c 'import numpy, pandas, scipy, sklearn, networkx, tensorflow, matplotlib; print("V2 scientific imports: OK")'
if [[ ${RUN_TESTS} -eq 1 ]]; then
    run_cmd env PYTHONDONTWRITEBYTECODE=1 "${PYTHON}" -m pytest -q \
        tests/test_v2_alternative_inputs.py tests/test_v2_publication_multiverse.py
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

FLY_ASSIGNMENTS=""
prepare_fly_assignments() {
    local namespace="$1" year="$2" gaf="$3" outdir
    outdir="${OUTPUT_DIR}/01_derived_labels/${namespace}/fly/${year}"
    FLY_ASSIGNMENTS="${outdir}/assignments.tsv"
    mkdir -p "${outdir}"
    if [[ ! (${RESUME} -eq 1 && -s "${FLY_ASSIGNMENTS}") ]]; then
        run_cmd "${PYTHON}" -m PhenGO.scripts.publication_inputs fly-assignments \
            --gaf "${gaf}" --output "${FLY_ASSIGNMENTS}" \
            --excluded-output "${outdir}/excluded.tsv" --taxon 7227
    fi
}

FLY_LETHAL=""
FLY_VIABLE=""
FLY_LABEL_AUDIT=""
prepare_fly_fail_closed_labels() {
    local year="$1" phenotype="$2" outdir
    outdir="${OUTPUT_DIR}/01_derived_labels/fly_fail_closed/${year}"
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

WORM_TERMS=""
prepare_worm_terms() {
    local key="$1" ontology="$2" outdir
    outdir="${OUTPUT_DIR}/01_derived_labels/worm_corrected/${key}"
    WORM_TERMS="${outdir}/lethal_terms.tsv"
    mkdir -p "${outdir}"
    if [[ ! (${RESUME} -eq 1 && -s "${WORM_TERMS}") ]]; then
        run_cmd "${PYTHON}" -m PhenGO.scripts.get_phenotype_terms \
            --obo-file "${ontology}" \
            --term-list "${REPO_DIR}/data/worm/lethal_terms/root_lethal_phenotype_terms.txt" \
            --results "${WORM_TERMS}"
    fi
}

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
        -go_relations is_a part_of -min_go_gene_count "${MIN_GO_GENE_COUNT}"
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
    if [[ ${#extra[@]} -gt 0 ]]; then command+=( "${extra[@]}" ); fi
    run_cmd "${command[@]}"
    run_cmd "${PYTHON}" -m PhenGO.scripts.arff_validator -i "${arff}" -o "${qc}" --fail-on-error
}

while IFS=$'\x1f' read -r organism year snapshot_id phenotype_file viability_file gaf_file go_obo_file phenotype_ontology_file phenotype_release go_annotation_release go_ontology_release phenotype_ontology_release label_source nonlethal_policy mixed_label_policy ambiguous_label_policy retrieval_date; do
    selected_organism "${organism}" || continue
    selected_year "${year}" || continue
    phenotype="$(absolute_path "${phenotype_file}")"
    archive_gaf="$(absolute_path "${gaf_file}")"
    go_obo="$(absolute_path "${go_obo_file}")"

    alt_needed=0
    case "${organism}" in
        fish)
            if selected_track primary_mod_gaf; then alt_needed=1; fi
            ;;
        fly)
            for candidate in primary_mod_gaf legacy_like_mod_gaf fail_closed_driver_aware_mod_gaf; do
                if selected_track "${candidate}"; then alt_needed=1; fi
            done
            ;;
        worm)
            for candidate in primary_mod_gaf worm_fixed_2025_terms_mod_gaf; do
                if selected_track "${candidate}"; then alt_needed=1; fi
            done
            ;;
        yeast)
            if selected_track primary_mod_gaf; then alt_needed=1; fi
            ;;
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

    if selected_track primary_mod_gaf; then
        case "${organism}" in
            fish)
                if [[ -n "${alt_gaf}" ]]; then
                    build_arff primary_mod_gaf fish "${year}" "fish-${year}-primary-mod-gaf" \
                        "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}" \
                        "${go_annotation_release}; ZFIN provider file ${PREPARED_GAF_LABEL}" \
                        "${go_ontology_release}" "" "${retrieval_date}" release_records \
                        "${nonlethal_policy}" "${mixed_label_policy}" "${ambiguous_label_policy}" \
                        none ""
                fi
                ;;
            fly)
                prepare_fly_assignments primary_mod_gaf "${year}" "${alt_gaf}"
                build_arff primary_mod_gaf fly "${year}" "fly-${year}-primary-mod-gaf" \
                    "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}" \
                    "${go_annotation_release}; FlyBase provider file ${PREPARED_GAF_LABEL}" "${go_ontology_release}" \
                    "${phenotype_ontology_release}" "${retrieval_date}" release_records \
                    "${nonlethal_policy}" "${mixed_label_policy}" "${ambiguous_label_policy}" \
                    fly_assignments "${FLY_ASSIGNMENTS}"
                ;;
            worm)
                prepare_worm_terms "${year}" "$(absolute_path "${phenotype_ontology_file}")"
                worm_terms="${WORM_TERMS}"
                build_arff primary_mod_gaf worm "${year}" "worm-${year}-primary-mod-gaf" \
                    "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}" \
                    "${go_annotation_release}; WormBase provider file ${PREPARED_GAF_LABEL}; taxon-filtered" \
                    "${go_ontology_release}" "${phenotype_ontology_release}" "${retrieval_date}" \
                    release_records "${nonlethal_policy}" "${mixed_label_policy}" \
                    "${ambiguous_label_policy}" worm_terms "${worm_terms}"
                ;;
            yeast)
                build_arff primary_mod_gaf yeast "${year}" "yeast-${year}-primary-mod-gaf" \
                    "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}" \
                    "${go_annotation_release}; SGD provider file ${PREPARED_GAF_LABEL}" \
                    "${go_ontology_release}" "" "${retrieval_date}" release_records \
                    "${nonlethal_policy}" "${mixed_label_policy}" "${ambiguous_label_policy}" \
                    none ""
                ;;
        esac
    fi

    if selected_track legacy_like_go_archive; then
        case "${organism}" in
            fly)
                prepare_fly_assignments legacy_like_go_archive "${year}" "${archive_gaf}"
                build_arff legacy_like_go_archive fly "${year}" "fly-${year}-legacy-like-go-archive" \
                    "${phenotype}" "${archive_gaf}" "${go_obo}" "${phenotype_release}; legacy-like policy" \
                    "${go_annotation_release}" "${go_ontology_release}" "${phenotype_ontology_release}" \
                    "${retrieval_date}" release_records observed_viable lethal_wins exclude \
                    fly_assignments "${FLY_ASSIGNMENTS}" -allow_multigenes
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
        prepare_fly_assignments legacy_like_mod_gaf "${year}" "${alt_gaf}"
        build_arff legacy_like_mod_gaf fly "${year}" "fly-${year}-legacy-like-mod-gaf" \
            "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}; legacy-like policy" \
            "${go_annotation_release}; FlyBase provider file ${PREPARED_GAF_LABEL}" "${go_ontology_release}" \
            "${phenotype_ontology_release}" "${retrieval_date}" release_records observed_viable \
            lethal_wins exclude fly_assignments "${FLY_ASSIGNMENTS}" -allow_multigenes
    fi

    if [[ "${organism}" == "fly" ]] && { selected_track fail_closed_driver_aware_go_archive || selected_track fail_closed_driver_aware_mod_gaf; }; then
        prepare_fly_fail_closed_labels "${year}" "${phenotype}"
        if selected_track fail_closed_driver_aware_go_archive; then
            prepare_fly_assignments fail_closed_driver_aware_go_archive "${year}" "${archive_gaf}"
            build_arff fail_closed_driver_aware_go_archive fly "${year}" "fly-${year}-fail-closed-driver-go-archive" \
                "${phenotype}" "${archive_gaf}" "${go_obo}" "${phenotype_release}; fail-closed compound policy" \
                "${go_annotation_release}" "${go_ontology_release}" \
                "${phenotype_ontology_release}; every compound component audited; unresolved excluded" "${retrieval_date}" \
                gene_sets explicit_only exclude exclude fly_gene_sets \
                "${FLY_LETHAL}|${FLY_VIABLE}|${FLY_ASSIGNMENTS}"
        fi
        if selected_track fail_closed_driver_aware_mod_gaf; then
            prepare_fly_assignments fail_closed_driver_aware_mod_gaf "${year}" "${alt_gaf}"
            build_arff fail_closed_driver_aware_mod_gaf fly "${year}" "fly-${year}-fail-closed-driver-mod-gaf" \
                "${phenotype}" "${alt_gaf}" "${go_obo}" "${phenotype_release}; fail-closed compound policy" \
                "${go_annotation_release}; FlyBase provider file ${PREPARED_GAF_LABEL}" "${go_ontology_release}" \
                "${phenotype_ontology_release}; every compound component audited; unresolved excluded" "${retrieval_date}" \
                gene_sets explicit_only exclude exclude fly_gene_sets \
                "${FLY_LETHAL}|${FLY_VIABLE}|${FLY_ASSIGNMENTS}"
        fi
    fi

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
        prepare_worm_terms fixed_2025 "$(absolute_path "${fixed_ontology_rel}")"
        fixed_terms="${WORM_TERMS}"
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
done < <(emit_ledger_rows)

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
    done < <(emit_ledger_rows)
fi

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

while IFS=$'\x1f' read -r track timepoint year report_format source gaf_file go_obo_file go_annotation_release go_ontology_release retrieval_date; do
    selected_track "${track}" || continue
    selected_organism mouse || continue
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

read -r -a model_args <<< "${MODELS}"
read -r -a panel_args <<< "${PANELS}"
read -r -a hidden_unit_args <<< "${MLP_HIDDEN_UNITS}"
early_stopping_args=( -nn_early_stopping )
[[ "${MLP_EARLY_STOPPING}" == "0" ]] && early_stopping_args=( -no_nn_early_stopping )

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
                -nn_hidden_units "${hidden_unit_args[@]}" -nn_alpha "${MLP_ALPHA}" \
                -nn_batch_size "${MLP_BATCH_SIZE}" -nn_learning_rate_init "${MLP_LEARNING_RATE}" \
                -nn_max_iter "${MLP_MAX_ITER}" "${early_stopping_args[@]}" \
                -nn_validation_fraction "${MLP_VALIDATION_FRACTION}" \
                -nn_n_iter_no_change "${MLP_PATIENCE}" \
                -importance_repeats "${IMPORTANCE_REPEATS}" -top_k "${TOP_K}" \
                -seed "${SEED}" -overwrite
        fi
        if [[ ! (${RESUME} -eq 1 && -s "${temporal_output}/temporal_analysis_manifest.json") ]]; then
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
    -i "${OUTPUT_DIR}/06_publication_tables" -o "${OUTPUT_DIR}/07_reports/publication_tables_report.html"

if [[ ${DRY_RUN} -eq 0 ]]; then
    multiverse_args=(
        "${PYTHON}" "${MULTIVERSE_HELPER}" --base-run "${BASE_RUN}"
        --extension-run "${OUTPUT_DIR}" --repo-root "${REPO_DIR}"
        --output-dir "${CONSOLIDATED_DIR}" --input-audit "${INPUT_AUDIT}"
        --availability-ledger "${AVAILABILITY_LEDGER}"
    )
    [[ ${RESUME} -eq 1 ]] && multiverse_args+=( --overwrite )
    record_command "${multiverse_args[@]}"

    : > "${OUTPUT_DIR}/00_run_metadata/output_checksums.sha256"
    find "${OUTPUT_DIR}/01_derived_inputs" "${OUTPUT_DIR}/01_derived_labels" \
        "${OUTPUT_DIR}/02_arff" "${OUTPUT_DIR}/03_qc" "${OUTPUT_DIR}/04_ml" \
        "${OUTPUT_DIR}/05_temporal" "${OUTPUT_DIR}/06_publication_tables" \
        -type f -print | LC_ALL=C sort | while IFS= read -r path; do
            shasum -a 256 "${path}"
        done >> "${OUTPUT_DIR}/00_run_metadata/output_checksums.sha256"
    printf 'completed_utc\t%s\nbase_run\t%s\nbase_fingerprint\t%s\n' \
        "$(timestamp)" "${BASE_RUN}" "${BASE_FINGERPRINT}" > "${OUTPUT_DIR}/00_run_metadata/run.complete"

    # Freeze the extension tree before consolidation. Logging the copy command
    # through run_cmd would append to RUN_LOG while that same tree was copied.
    printf '\n[%s] ' "$(timestamp)"
    quote_command "${multiverse_args[@]}"
    printf '\n'
    "${multiverse_args[@]}"
    consolidated_manifest_sha="$(shasum -a 256 "${CONSOLIDATED_DIR}/consolidated_manifest.json" | awk '{print $1}')"
    printf 'completed_utc\t%s\nbase_run\t%s\nextension_run\t%s\nconsolidated_manifest_sha256\t%s\n' \
        "$(timestamp)" "${BASE_RUN}" "${OUTPUT_DIR}" "${consolidated_manifest_sha}" \
        > "${CONSOLIDATED_DIR}/consolidated.complete"
fi

if [[ ${DRY_RUN} -eq 1 ]]; then
    printf '\nPhenGO publication dry-run plan complete.\nPlanned extension: %s\nPlanned consolidated output: %s\n' \
        "${OUTPUT_DIR}" "${CONSOLIDATED_DIR}"
else
    printf '\nPhenGO publication V2 complete.\nExtension: %s\nConsolidated: %s\n' \
        "${OUTPUT_DIR}" "${CONSOLIDATED_DIR}"
fi
