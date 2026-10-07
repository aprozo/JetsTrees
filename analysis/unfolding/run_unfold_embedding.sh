#!/usr/bin/env bash
set -euo pipefail

########################
# Paths relative to this script
########################

# Directory where *this* script lives (analysis/unfolding)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# JetsTrees base (two levels up)
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Files relative to this structure
SIF="${SCRIPT_DIR}/roounfold.sif"              # analysis/unfolding/roounfold.sif
MACRO="${SCRIPT_DIR}/unfold_embedding.cxx"    # analysis/unfolding/unfold_embedding.cxx

########################
# Arguments
########################

# 1st arg: method (Bayes or SVD)
METHOD="${1:-BAYES}"

# 2nd arg: input (either basename in trees/ or an absolute path)
if [[ $# -ge 2 ]]; then
  if [[ "$2" = /* ]]; then
    INPUT="$2"  
  else
    INPUT="${BASE}/trees/$2"
  fi
else
  INPUT="${BASE}/trees/embedding_merged_HSCalib.root"
fi

# 3rd arg: output directory (default: analysis/unfolding/out_embedding under BASE)
OUT_DIR="${3:-${SCRIPT_DIR}/out_embedding_${METHOD}}"

########################
# Checks
########################

echo "----------------------------------------"
echo "Running unfolding"
echo "Method     : $METHOD"
echo "SCRIPT_DIR : $SCRIPT_DIR"
echo "BASE       : $BASE"
echo "SIF        : $SIF"
echo "Macro      : $MACRO"
echo "Input      : $INPUT"
echo "Output dir : $OUT_DIR"
echo "----------------------------------------"

[[ -f "$SIF"   ]] || { echo "ERROR: SIF not found:   $SIF";   exit 1; }
[[ -f "$MACRO" ]] || { echo "ERROR: MACRO not found: $MACRO"; exit 1; }
[[ -f "$INPUT" ]] || { echo "ERROR: Input not found: $INPUT"; exit 1; }

mkdir -p "$OUT_DIR"

########################
# Run inside container
########################

apptainer exec -e -B /gpfs01 \
  "$SIF" \
  root -l -b <<EOF
gSystem->Load("libRooUnfold");
.x ${MACRO}+("${INPUT}","${OUT_DIR}","${METHOD}");
.q
EOF

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"