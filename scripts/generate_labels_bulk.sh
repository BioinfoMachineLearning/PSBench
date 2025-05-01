#!/bin/bash

# Print help for usage
if [ "$#" -lt 7 ]; then
  echo "Usage: $0 <fasta_dir> <indir_base> <nativedir_base> <outdir_base> <usalign_path> <clustalw_path> <target1> [<target2> ...]"
  echo ""
  echo "Example:"
  echo "$0 /path/to/fasta /path/to/predicted /path/to/native /path/to/output /path/to/USalign /path/to/clustalw T1001 T1002"
  exit 1
fi

# Assign inputs
FASTA_DIR=$1
INDIR_BASE=$2
NATIVEDIR_BASE=$3
OUTDIR_BASE=$4
USALIGN=$5
CLUSTALW=$6
shift 6

TARGETS=("$@")

for target in "${TARGETS[@]}"; do
    echo "Processing $target"
    
    FASTA="${FASTA_DIR}/${target}.fasta"
    INDIR="${INDIR_BASE}/${target}/"
    NATIVE="${NATIVEDIR_BASE}/${target}.pdb"
    OUTDIR="${OUTDIR_BASE}/"

    python generate_labels.py \
        --fasta "$FASTA" \
        --indir "$INDIR" \
        --nativedir "$NATIVE" \
        --outdir "$OUTDIR" \
        --usalign_program "$USALIGN" \
        --clustalw_program "$CLUSTALW"
    
    echo "Done with $target"
    echo "----------------------------------------"
done
