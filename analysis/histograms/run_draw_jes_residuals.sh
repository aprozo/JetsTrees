#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT="${1:-${SCRIPT_DIR}/jes_jer_hists.root}"
PLOT_DIR="${2:-${SCRIPT_DIR}/jes_jer_plots}"
MACRO="${SCRIPT_DIR}/draw_jes_residuals.C"

echo "----------------------------------------"
echo "Drawing JES residual distributions"
echo "Script dir : ${SCRIPT_DIR}"
echo "Input      : ${INPUT}"
echo "Plot dir   : ${PLOT_DIR}"
echo "Macro      : ${MACRO}"
echo "----------------------------------------"

root -l -b -q "${MACRO}(\"${INPUT}\",\"${PLOT_DIR}\")"