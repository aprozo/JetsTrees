#!/usr/bin/env bash
set -euo pipefail

# Place beside compare_calibrations.C and roounfold.sif.
# Usage: bash run_compare_calibrations.sh ORIG_FOLDER HS_FOLDER [OUTPUT_FOLDER]
# Relative folder arguments are resolved relative to this wrapper.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIF="${SCRIPT_DIR}/roounfold.sif"
MACRO="${SCRIPT_DIR}/compare_calibrations.C"
if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 ORIG_FOLDER HS_FOLDER [OUTPUT_FOLDER]" >&2
  exit 2
fi
command -v apptainer >/dev/null || { echo 'ERROR: apptainer is unavailable.' >&2; exit 1; }
[[ -f "$SIF" ]] || { echo "ERROR: Missing $SIF" >&2; exit 1; }
[[ -f "$MACRO" ]] || { echo "ERROR: Missing $MACRO" >&2; exit 1; }

absolute_folder() {
  local folder="$1"
  [[ "$folder" = /* ]] || folder="${SCRIPT_DIR}/${folder}"
  (cd "$folder" && pwd -P)
}
ORIG_DIR="$(absolute_folder "$1")"
HS_DIR="$(absolute_folder "$2")"
OUT_DIR="${3:-calibration_comparison}"
[[ "$OUT_DIR" = /* ]] || OUT_DIR="${SCRIPT_DIR}/${OUT_DIR}"
mkdir -p "$OUT_DIR"
OUT_DIR="$(absolute_folder "$OUT_DIR")"
for folder in "$ORIG_DIR" "$HS_DIR"; do
  [[ -f "${folder}/unfolded_data.root" ]] || { echo "ERROR: Missing ${folder}/unfolded_data.root" >&2; exit 1; }
done
# Paths enter a C++ string in the ROOT invocation; reject ambiguous characters.
for value in "$SCRIPT_DIR" "$ORIG_DIR" "$HS_DIR" "$OUT_DIR"; do
  case "$value" in
    *\"*|*\\*|*$'\n'*|*$'\r'*) echo 'ERROR: Paths must not contain quotes, backslashes or newlines.' >&2; exit 1 ;;
  esac
done
# Use a fresh ACLiC directory, avoiding libraries compiled by the host ROOT.
BUILD_DIR="$(mktemp -d "${OUT_DIR}/aclic.XXXXXX")"
trap 'rm -rf -- "$BUILD_DIR"' EXIT
LOG="${OUT_DIR}/compare_calibrations.log"
echo "Original: $ORIG_DIR"
echo "Hanseul baseline: $HS_DIR"
echo "Output: $OUT_DIR"
echo "Container: $SIF"

# Bind both /gpfs01 and resolved directories (also covers paths outside GPFS).
# --env explicitly forwards only the comparison parameters into the clean environment.
apptainer exec -e \
  -B /gpfs01 \
  -B "$SCRIPT_DIR" -B "$ORIG_DIR" -B "$HS_DIR" -B "$OUT_DIR" \
  --env "CALIB_MACRO=$MACRO" \
  --env "CALIB_ORIG=$ORIG_DIR" \
  --env "CALIB_HS=$HS_DIR" \
  --env "CALIB_OUT=$OUT_DIR" \
  --env "CALIB_BUILD=$BUILD_DIR" \
  "$SIF" root -l -b 2>&1 <<'ROOT_COMMANDS' | tee "$LOG"
std::cout << "Container ROOT version: " << gROOT->GetVersion() << std::endl;
gSystem->SetBuildDir(gSystem->Getenv("CALIB_BUILD"), kTRUE);
if (gROOT->LoadMacro(Form("%s+",gSystem->Getenv("CALIB_MACRO"))) < 0) gSystem->Exit(1);
Int_t comparisonError = 0;
gROOT->ProcessLine(Form("compare_calibrations(\"%s\",\"%s\",\"%s\")",gSystem->Getenv("CALIB_ORIG"),gSystem->Getenv("CALIB_HS"),gSystem->Getenv("CALIB_OUT")), &comparisonError);
if (comparisonError != 0) gSystem->Exit(1);
.q
ROOT_COMMANDS

echo "Finished. Log: $LOG"
