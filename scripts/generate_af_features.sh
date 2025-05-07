: <<'LICENSE' 
MIT License

Copyright (c) 2025 BioinfoMachineLearning

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. 
LICENSE

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
