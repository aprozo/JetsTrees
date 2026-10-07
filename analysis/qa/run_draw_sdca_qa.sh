#!/bin/bash

set -e

# Directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input ROOT file in JetsTrees/trees/
IN_BASENAME="${1:-data_merged.root}"

# Output directory inside analysis/qa/
OUT_BASENAME="${2:-sdca_qa_plots}"

INPUT="${SCRIPT_DIR}/../../trees/${IN_BASENAME}"
OUT_DIR="${SCRIPT_DIR}/${OUT_BASENAME}"
MACRO="${SCRIPT_DIR}/draw_sdca_qa.C"

mkdir -p "$OUT_DIR"

echo "----------------------------------------"
echo "Drawing signed-DCA QA"
echo "Script dir : $SCRIPT_DIR"
echo "Input      : $INPUT"
echo "Output dir : $OUT_DIR"
echo "Macro      : $MACRO"
echo "----------------------------------------"

root -l -b -q "$MACRO(\"$INPUT\",\"$OUT_DIR\")"

echo "----------------------------------------"
echo "Done."
echo "Plots saved in:"
echo "$OUT_DIR"
echo "----------------------------------------"