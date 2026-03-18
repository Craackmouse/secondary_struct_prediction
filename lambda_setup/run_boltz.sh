#!/bin/bash
# Run Boltz inference on all test targets.
# Generates 5 diverse structure predictions per target using the MSA server.
#
# Usage: bash run_boltz.sh [--no-msa]   (--no-msa skips MSA server, faster but worse)

set -e

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rna_pred

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUTS_DIR="$SCRIPT_DIR/inputs"
OUTPUTS_DIR="$SCRIPT_DIR/outputs"
LOGS_DIR="$SCRIPT_DIR/logs"
FAILED_LOG="$SCRIPT_DIR/failed_targets.txt"

mkdir -p "$OUTPUTS_DIR" "$LOGS_DIR"
> "$FAILED_LOG"

# MSA flag
USE_MSA="--use_msa_server"
if [[ "$1" == "--no-msa" ]]; then
    USE_MSA=""
    echo "[INFO] MSA server disabled — predictions will be lower quality."
fi

YAMLS=("$INPUTS_DIR"/*.yaml)
TOTAL=${#YAMLS[@]}
echo "=== Boltz Inference: $TOTAL targets, 5 samples each ==="
echo "Outputs → $OUTPUTS_DIR"
echo

START_ALL=$(date +%s)
COUNT=0

for yaml_file in "${YAMLS[@]}"; do
    TARGET=$(basename "$yaml_file" .yaml)
    COUNT=$((COUNT + 1))
    OUT_DIR="$OUTPUTS_DIR/$TARGET"

    # Skip if already done (all 5 CIFs exist)
    N_DONE=$(find "$OUT_DIR" -name "*.cif" 2>/dev/null | wc -l)
    if [[ "$N_DONE" -ge 5 ]]; then
        echo "[$COUNT/$TOTAL] $TARGET — already done ($N_DONE CIFs), skipping."
        continue
    fi

    echo "[$COUNT/$TOTAL] $TARGET — running..."
    START=$(date +%s)

    boltz predict "$yaml_file" \
        $USE_MSA \
        --diffusion_samples 5 \
        --no_kernels \
        --out_dir "$OUT_DIR" \
        2>&1 | tee "$LOGS_DIR/${TARGET}.log" || {
            echo "[FAIL] $TARGET — check logs/${TARGET}.log"
            echo "$TARGET" >> "$FAILED_LOG"
            continue
        }

    END=$(date +%s)
    echo "[$COUNT/$TOTAL] $TARGET — done in $((END - START))s"
    echo
done

END_ALL=$(date +%s)
echo "=== All done in $((END_ALL - START_ALL))s ==="

FAILED=$(wc -l < "$FAILED_LOG")
if [[ "$FAILED" -gt 0 ]]; then
    echo "[WARN] $FAILED target(s) failed — see failed_targets.txt"
    cat "$FAILED_LOG"
fi

echo
echo "Next: python convert_submission.py"
