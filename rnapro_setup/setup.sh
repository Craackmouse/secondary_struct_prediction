#!/bin/bash
# RNAPro setup on Lambda Labs A100 instance.
# Run from home directory: bash ~/repo/rnapro_setup/setup.sh
set -e

echo "=== RNAPro Setup ==="

# --- Conda env ---
source "$(conda info --base)/etc/profile.d/conda.sh"
conda create -n rnapro python=3.12 -y
conda activate rnapro

# --- Clone RNAPro ---
cd ~
if [ ! -d "RNAPro" ]; then
    git clone https://github.com/NVIDIA-Digital-Bio/RNAPro.git
fi
cd RNAPro

# --- Install dependencies ---
pip install -r requirements.txt
pip install -e .

# --- Download model weights from HuggingFace ---
echo "=== Downloading model weights ==="
pip install huggingface_hub
python3 -c "
from huggingface_hub import hf_hub_download
import shutil, os
# Download Public-Best checkpoint
path = hf_hub_download('nvidia/RNAPro-Public-Best-500M', filename='rnapro_base.pt')
shutil.copy(path, './rnapro_base.pt')
print(f'Weights saved to ./rnapro_base.pt ({os.path.getsize(\"./rnapro_base.pt\") / 1e9:.1f} GB)')
"

# --- Download RibonanzaNet2 checkpoint ---
echo "=== Downloading RibonanzaNet2 encoder ==="
mkdir -p release_data/ribonanzanet2_checkpoint
cd release_data/ribonanzanet2_checkpoint
if [ ! -f "pytorch_model_fsdp.bin" ]; then
    pip install kaggle
    # Try kaggle API first, fall back to direct download
    python3 -c "
from huggingface_hub import hf_hub_download
import os
# RibonanzaNet2 is also mirrored; try kaggle models API
os.system('kaggle models instances versions download shujun717/ribonanzanet2/pyTorch/alpha/1 -p . --untar')
" || {
    echo "[WARN] Kaggle download failed. Trying alternative..."
    curl -L -o ribonanzanet2.tar.gz \
      "https://www.kaggle.com/api/v1/models/shujun717/ribonanzanet2/pyTorch/alpha/1/download"
    tar -xzvf ribonanzanet2.tar.gz
    rm -f ribonanzanet2.tar.gz
}
fi
cd ~/RNAPro

# --- Generate CCD cache ---
echo "=== Generating CCD cache ==="
python3 preprocess/gen_ccd_cache.py -n 8

echo ""
echo "=== Setup complete ==="
echo "Next: bash ~/repo/rnapro_setup/prepare_and_run.sh"
