#!/bin/bash
set -e

echo "=== RNA Structure Prediction — Lambda Labs H100 Setup ==="
echo "This script sets up the Boltz environment for RNA 3D structure prediction."
echo

# --- Conda env ---
if conda env list | grep -q "rna_pred"; then
    echo "[INFO] rna_pred env already exists, skipping create."
else
    conda create -n rna_pred python=3.10 -y
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rna_pred

# --- Core packages ---
pip install --upgrade pip
pip install boltz pandas biopython gemmi tqdm requests PyYAML

echo
echo "=== Verifying installs ==="
python -c "import boltz; print('boltz OK')"
python -c "import gemmi; print('gemmi OK')"
python -c "import pandas; print('pandas OK')"

echo
echo "=== Setup complete. Next steps: ==="
echo "  1. Upload test_sequences.csv to this directory"
echo "  2. Run: conda activate rna_pred && bash run_all.sh"
