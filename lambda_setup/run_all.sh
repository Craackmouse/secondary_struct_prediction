#!/bin/bash
# Master pipeline: prepare → infer → convert → submission.csv
# Run from the lambda_setup/ directory on your Lambda Labs instance.
#
# Usage:
#   bash run_all.sh           # full pipeline with MSA server (recommended)
#   bash run_all.sh --no-msa  # skip MSA server (faster, lower quality)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=============================================="
echo "  RNA 3D Structure Prediction Pipeline"
echo "  Working dir: $SCRIPT_DIR"
echo "=============================================="
echo

# --- Activate conda ---
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rna_pred

# --- Check test_sequences.csv ---
if [[ ! -f "test_sequences.csv" ]]; then
    echo "ERROR: test_sequences.csv not found in $SCRIPT_DIR"
    echo "Upload it with:"
    echo "  scp test_sequences.csv ubuntu@<LAMBDA_IP>:~/lambda_setup/"
    exit 1
fi

# --- Step 1: Prepare Boltz YAML inputs ---
echo "[STEP 1/3] Preparing YAML inputs..."
python prepare_inputs.py
echo

# --- Step 2: Run Boltz inference ---
echo "[STEP 2/3] Running Boltz inference (this will take a while)..."
bash run_boltz.sh "$1"
echo

# --- Step 3: Convert to submission CSV ---
echo "[STEP 3/3] Converting outputs to submission.csv..."
python convert_submission.py
echo

echo "=============================================="
echo "  DONE. Submission file: $SCRIPT_DIR/submission.csv"
echo "  Download with:"
echo "    scp ubuntu@<LAMBDA_IP>:~/lambda_setup/submission.csv ."
echo "=============================================="
