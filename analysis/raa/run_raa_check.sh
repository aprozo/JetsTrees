#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Preliminary inclusive-jet R_AA using ROOT 6 in Apptainer
# Run from anywhere: bash analysis/raa/run_raa_check.sh
# ------------------------------------------------------------

# Project locations on the host.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
SIF="${BASE}/analysis/unfolding/roounfold.sif"
MACRO="${BASE}/analysis/raa/raa_check.cxx"
INPUT="${BASE}/analysis/unfolding/out_data_HSCalib/unfolded_data.root"
PP_DIR="${BASE}/analysis/pp_reference"
OUT_DIR="${BASE}/analysis/raa/out_raa"

# Plot settings.
MIN_PLOT_PT="5.0"

# The same project directory, as seen inside the container.
CONTAINER_BASE="/work/JetsTrees"
C_MACRO="${CONTAINER_BASE}/analysis/raa/raa_check.cxx"
C_INPUT="${CONTAINER_BASE}/analysis/unfolding/out_data_HSCalib/unfolded_data.root"
C_PP_DIR="${CONTAINER_BASE}/analysis/pp_reference"
C_OUT_DIR="${CONTAINER_BASE}/analysis/raa/out_raa"

# Check required inputs before starting ROOT.
for file in "${SIF}" "${MACRO}" "${INPUT}"; do
    if [[ ! -f "${file}" ]]; then
        echo "[ERROR] Missing file: ${file}" >&2
        exit 1
    fi
done

if [[ ! -d "${PP_DIR}" ]]; then
    echo "[ERROR] Missing pp reference directory: ${PP_DIR}" >&2
    exit 1
fi

mkdir -p "${OUT_DIR}/png"

echo "------------------------------------------------------------"
echo "Preliminary inclusive-jet R_AA"
echo "Au+Au input   : ${INPUT}"
echo "pp references : ${PP_DIR}"
echo "Output        : ${OUT_DIR}"
echo "Minimum pT    : ${MIN_PLOT_PT} GeV/c"
echo "------------------------------------------------------------"

# ROOT macro arguments are paths *inside* the container.
ROOT_CALL="${C_MACRO}(\"${C_INPUT}\",\"${C_PP_DIR}\",\"${C_OUT_DIR}\",${MIN_PLOT_PT})"

# Bind the host project to /work/JetsTrees.
apptainer exec \
    --bind "${BASE}:${CONTAINER_BASE}" \
    "${SIF}" \
    root -l -b -q "${ROOT_CALL}"

echo "[DONE] Output directory: ${OUT_DIR}"