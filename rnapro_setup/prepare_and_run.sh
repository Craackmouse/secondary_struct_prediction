#!/bin/bash
# Prepare inputs and run RNAPro inference on all test targets.
# Run from home directory: bash ~/repo/rnapro_setup/prepare_and_run.sh
set -e

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rnapro

REPO_DIR=~/repo
RNAPRO_DIR=~/RNAPro
LAMBDA_DIR=~/repo/lambda_setup
OUTPUT_DIR=$RNAPRO_DIR/output
MSA_DIR=$RNAPRO_DIR/msa

mkdir -p "$OUTPUT_DIR" "$MSA_DIR"

# --- Step 1: Prepare sequences CSV for RNAPro ---
echo "[STEP 1] Preparing sequences CSV..."
python3 "$REPO_DIR/rnapro_setup/prepare_rnapro_inputs.py" \
    --test_csv "$LAMBDA_DIR/test_sequences.csv" \
    --out_csv "$RNAPRO_DIR/sequences.csv" \
    --msa_dir "$MSA_DIR"

# --- Step 2: Run RNAPro inference ---
echo "[STEP 2] Running RNAPro inference..."
cd "$RNAPRO_DIR"

export LAYERNORM_TYPE=torch

# Run without templates first (still strong — RibonanzaNet2 + MSA carry most of the weight)
# Templates can be added later for a second run
python3 runner/inference.py \
    --model_name rnapro_base \
    --seeds 42 \
    --dump_dir "$OUTPUT_DIR" \
    --load_checkpoint_path ./rnapro_base.pt \
    --use_msa true \
    --use_template none \
    --model.use_template none \
    --model.use_RibonanzaNet2 true \
    --model.template_embedder.n_blocks 2 \
    --model.ribonanza_net_path ./release_data/ribonanzanet2_checkpoint \
    --model.N_cycle 10 \
    --sample_diffusion.N_sample 5 \
    --sample_diffusion.N_step 200 \
    --load_strict false \
    --num_workers 0 \
    --triangle_attention torch \
    --triangle_multiplicative torch \
    --sequences_csv ./sequences.csv \
    --max_len 5000 \
    --logger logging \
    2>&1 | tee "$REPO_DIR/rnapro_setup/rnapro_run.log"

echo ""
echo "=== RNAPro inference complete ==="
echo "Output: $OUTPUT_DIR"
echo "Next: python3 ~/repo/rnapro_setup/convert_rnapro_submission.py"
