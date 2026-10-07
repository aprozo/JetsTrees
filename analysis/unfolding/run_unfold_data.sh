#!/usr/bin/env bash
set -euo pipefail

########################
# Configuration
########################

# Directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Base dir = two levels up from unfolding/
# (i.e. .../JetsTrees, assuming the same structure)
BASE="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Singularity image and macro, both in the same dir as this script
SIF="${SCRIPT_DIR}/roounfold.sif"
MACRO="${SCRIPT_DIR}/unfold_data_wEff.cxx"

########################
# Arguments
########################

# 1st arg: HT data input
if [[ $# -ge 1 ]]; then
  if [[ "$1" = /* ]]; then
    INPUT_HT="$1"
  else
    INPUT_HT="${BASE}/trees/$1"
  fi
else
  INPUT_HT="${BASE}/trees/data_merged_OrigCalib.root"
fi

# 2nd arg: MB data input
if [[ $# -ge 2 ]]; then
  if [[ "$2" = /* ]]; then
    INPUT_MB="$2"
  else
    INPUT_MB="${BASE}/trees/$2"
  fi
else
  INPUT_MB="${BASE}/trees/dataMB_merged.root"
fi

# 3rd arg: RESPONSE ROOT FILE
RESP_FILE="${3:-${SCRIPT_DIR}/out_embedding_BAYES_OrigCalib/responses_embedding.root}"

# 4th arg: EFFICIENCIES ROOT FILE
EFF_FILE="${4:-${BASE}/analysis/efficiencies/efficiencies.root}"

# 5th arg: output directory
OUT_DIR="${5:-${SCRIPT_DIR}/out_data}"

# 6th arg: number of Bayes iterations
NITER="${6:-4}"

########################
# Checks
########################

echo "----------------------------------------"
echo "Running unfolding on REAL DATA"
echo "SCRIPT_DIR    : $SCRIPT_DIR"
echo "BASE          : $BASE"
echo "SIF           : $SIF"
echo "Macro         : $MACRO"
echo "Input HT data : $INPUT_HT"
echo "Input MB data : $INPUT_MB"
echo "Resp. file    : $RESP_FILE"
echo "Eff. file     : $EFF_FILE"
echo "Output dir    : $OUT_DIR"
echo "Iterations    : $NITER"
echo "----------------------------------------"

[[ -f "$SIF"       ]] || { echo "ERROR: SIF not found:       $SIF";       exit 1; }
[[ -f "$MACRO"     ]] || { echo "ERROR: MACRO not found:     $MACRO";     exit 1; }
[[ -f "$INPUT_HT"  ]] || { echo "ERROR: HT input not found:  $INPUT_HT";  exit 1; }
[[ -f "$INPUT_MB"  ]] || { echo "ERROR: MB input not found:  $INPUT_MB";  exit 1; }
[[ -f "$RESP_FILE" ]] || { echo "ERROR: Resp. file not found: $RESP_FILE"; exit 1; }
[[ -f "$EFF_FILE"  ]] || { echo "ERROR: Eff. file not found: $EFF_FILE"; exit 1; }

mkdir -p "$OUT_DIR"

########################
# Run inside container
########################

# IMPORTANT: bind /gpfs01 so the container sees the same paths
apptainer exec -e -B /gpfs01 \
  "$SIF" \
  root -l -b <<EOF
gSystem->Load("libRooUnfold");
.x ${MACRO}+("${INPUT_HT}","${INPUT_MB}","${RESP_FILE}","${EFF_FILE}","${OUT_DIR}",${NITER});
.q
EOF

echo "----------------------------------------"
echo "Done."
echo "----------------------------------------"
