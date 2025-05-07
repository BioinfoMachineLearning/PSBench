# -------------------------------------------------------------------------------------------------------------------------------------
# Following code curated for PoseBench: (https://github.com/BioinfoMachineLearning/PSBench)
# -------------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash

# Default values
TARGETS=()

# Help message
print_help() {
  echo "Usage: $0 --fasta_dir PATH --predicted_dir PATH --native_dir PATH --outdir PATH --usalign PATH --clustalw PATH --targets TARGET1 [TARGET2 ...]"
  echo ""
  echo "Example:"
  echo "$0 --fasta_dir /path/to/fasta_folder --predicted_dir /path/to/predicted_folder --native_dir /path/to/native_folder \\"
  echo "   --outdir /path/to/output_folder --usalign /path/to/USalign_binary --clustalw /path/to/clustalw_binary --targets T1001 T1002"
  exit 1
}

# Parse named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --fasta_dir)
      FASTA_DIR="$2"
      shift 2
      ;;
    --predicted_dir)
      PREDICTED_DIR_BASE="$2"
      shift 2
      ;;
    --native_dir)
      NATIVE_DIR_BASE="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR_BASE="$2"
      shift 2
      ;;
    --usalign)
      USALIGN="$2"
      shift 2
      ;;
    --clustalw)
      CLUSTALW="$2"
      shift 2
      ;;
    --targets)
      shift
      while [[ $# -gt 0 && $1 != --* ]]; do
        TARGETS+=("$1")
        shift
      done
      ;;
    --help|-h)
      print_help
      ;;
    *)
      echo "Unknown argument: $1"
      print_help
      ;;
  esac
done

# Check for required args
if [[ -z "$FASTA_DIR" || -z "$PREDICTED_DIR_BASE" || -z "$NATIVE_DIR_BASE" || -z "$OUTDIR_BASE" || -z "$USALIGN" || -z "$CLUSTALW" || ${#TARGETS[@]} -eq 0 ]]; then
  echo "Error: Missing required arguments."
  print_help
fi

# Run for each target
for target in "${TARGETS[@]}"; do
    echo "Processing $target"

    FASTA="${FASTA_DIR}/${target}.fasta"
    PREDICTED="${PREDICTED_DIR_BASE}/${target}/"
    NATIVE="${NATIVE_DIR_BASE}/${target}.pdb"
    OUTDIR="${OUTDIR_BASE}/"

    python generate_quality_scores.py \
        --fasta "$FASTA" \
        --predicted_dir "$PREDICTED" \
        --native_dir "$NATIVE" \
        --outdir "$OUTDIR" \
        --usalign_program "$USALIGN" \
        --clustalw_program "$CLUSTALW"

    echo "Done with $target"
    echo "----------------------------------------"
done
