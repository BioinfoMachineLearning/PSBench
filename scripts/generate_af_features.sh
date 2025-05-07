# -------------------------------------------------------------------------------------------------------------------------------------
# Following code curated for PoseBench: (https://github.com/BioinfoMachineLearning/PSBench)
# -------------------------------------------------------------------------------------------------------------------------------------


#!/bin/bash

# Default values
FASTA_DIR=""
PREDICTED_DIR=""
PKL_DIR=""
OUTDIR=""
TARGETS=()

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fasta_dir) FASTA_DIR="$2"; shift ;;
        --predicted_dir) PREDICTED_DIR="$2"; shift ;;
        --pkl_dir) PKL_DIR="$2"; shift ;;
        --outdir) OUTDIR="$2"; shift ;;
        --targets) shift; while [[ "$1" != --* && "$1" != "" ]]; do TARGETS+=("$1"); shift; done ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check for required arguments
if [[ -z "$FASTA_DIR" || -z "$PREDICTED_DIR" || -z "$PKL_DIR" || -z "$OUTDIR" || ${#TARGETS[@]} -eq 0 ]]; then
    echo "Usage: $0 --fasta_dir PATH --PREDICTED_DIR PATH --pkl_dir PATH --outdir PATH --targets TARGET1 [TARGET2 ...]"
    exit 1
fi

# Main loop
for target in "${TARGETS[@]}"; do
    echo "Processing $target"

    FASTA="${FASTA_DIR}/${target}.fasta"
    PDB_PATH="${PREDICTED_DIR}/${target}/"
    PKL_PATH="${PKL_DIR}/${target}/"
    OUTDIR="${OUTDIR}"

    python generate_af_features.py \
        --fasta "$FASTA" \
        --predicted_dir "$PDB_PATH" \
        --pkl_dir "$PKL_PATH" \
        --outdir "$OUTDIR"

    echo "Done with $target"
    echo "----------------------------------------"
done
