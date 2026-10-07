#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

IN_BASENAME="${1:-embedding_merged_MCReco1p5.root}"
OUT_BASENAME="${2:-jes_jer_hists.root}"
PLOT_DIR="${3:-${SCRIPT_DIR}/jes_jer_plots}"

INPUT="${SCRIPT_DIR}/../../trees/${IN_BASENAME}"
OUTPUT="${SCRIPT_DIR}/${OUT_BASENAME}"
MACRO="${SCRIPT_DIR}/make_jes_jer.C"

echo "----------------------------------------"
echo "Running MC/reco and JES/JER diagnostics"
echo "Script dir : $SCRIPT_DIR"
echo "Input      : $INPUT"
echo "Output     : $OUTPUT"
echo "Plot dir   : $PLOT_DIR"
echo "Macro      : $MACRO"
echo "----------------------------------------"

root -l -b -q "${MACRO}(\"${INPUT}\",\"${OUTPUT}\",\"${PLOT_DIR}\")"

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"