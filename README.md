# RNA 3D Structure Prediction

Entry for the [Stanford RNA 3D Folding Part 2](https://www.kaggle.com/competitions/stanford-rna-3d-folding-2) Kaggle competition.

## Problem Statement

Given an RNA nucleotide sequence, predict the 3D coordinates of the C1' atom for each residue. Each target requires multiple structure predictions (up to 10 samples), scored by TM-align against the ground truth. The competition dataset is designed to be template-free — test targets have no close structural relatives in the PDB.

Key challenges:
- RNA 3D structure prediction is far less mature than protein folding
- Test sequences are deliberately dissimilar to training data
- Long sequences (>512 nt) must be handled via chunking strategies
- Predictions must be physically plausible (valid bond lengths, no steric clashes)

## Approaches

### 1. Protenix + TBM Pipeline (`stanford-rna-3d-folding.ipynb`)

The main submission notebook. Uses a two-phase strategy:

**Phase 1 — Template-Based Modeling (TBM):**
- Aligns each test sequence against a pool of ~3,700 training structures using a global pairwise aligner
- For high-similarity hits (configurable thresholds), maps template coordinates onto the query via sequence alignment
- Generates diversity across the 10 sample slots using geometric perturbations: small noise injection, hinge bending at chain pivots, inter-chain jittering, and spline-based smooth wiggle
- All TBM predictions are refined with `adaptive_rna_constraints`

**Phase 2 — Protenix (deep learning):**
- Targets with insufficient template coverage are routed to [Protenix](https://github.com/bytedance/protenix), ByteDance's open-source AlphaFold3 reproduction
- Long sequences are split into overlapping chunks (512 nt, 128 overlap), predicted independently, then stitched back via Kabsch alignment on overlapping residues with linear blending
- Protenix predictions also receive a gentle physics constraint pass (confidence=0.9)

**Phase 3 — Post-prediction refinements:**
- **Self-consistency reordering**: computes pairwise RMSD across all 10 predictions per target and places the most "consensus" structure (lowest average RMSD to all others) in slot 1
- **ViennaRNA secondary structure restraints** (optional): predicts base pairs via thermodynamic folding and applies soft distance restraints pushing paired C1' atoms toward ~10.4 A (canonical A-form geometry)
- **`adaptive_rna_constraints`**: iterative physics correction that normalizes C1'–C1' bond lengths (~5.95 A), enforces i→i+2 distances (~10.2 A), applies Laplacian smoothing, and resolves steric clashes via self-avoidance

### 2. Tri-Model Hybrid (`rnapro-inference-with-tbm.ipynb`)

Combines three independent prediction sources per target:

| Slot | Source | Description |
|------|--------|-------------|
| 1 | Best TBM | Template-based modeling (enhanced) |
| 2 | Best Boltz-2 | Deep learning ([Boltz](https://github.com/jost-tech/boltz)) |
| 3 | Best RNAPro | Deep learning ([RNAPro](https://github.com/NVIDIA-Digital-Bio/RNAPro), NVIDIA 500M model) |
| 4 | Boltz + TBM blend | 50/50 centroid-aligned average |
| 5 | RNAPro + Boltz blend | 50/50 centroid-aligned average |

Fallback tiers handle cases where individual models fail:
- **Tier A**: all three models available — use full strategy above
- **Tier B**: RNAPro unavailable — TBM, Boltz, blend, geometric variants
- **Tier C**: Boltz unavailable — TBM, RNAPro, blend, geometric variants
- **Tier D**: only TBM available — 5 geometric variants of TBM

### 3. Protenix with Biopython Fix (`rna-3d-protenix-biopython-fix.ipynb`)

Variant of the Protenix+TBM pipeline with a patched Biopython integration for improved atom mask extraction (`get_c1_mask`). Handles edge cases in the C1' coordinate extraction where `centre_atom_mask`, `atom_name`, and `atom_to_tokatom_idx` fallbacks are needed depending on the Protenix version.

### 4. RNAPro Standalone (`rnapro_setup/`)

Setup and inference scripts for running NVIDIA's RNAPro model independently:
- `setup.sh` — environment setup and dependency installation
- `prepare_rnapro_inputs.py` — converts competition CSV to RNAPro input format
- `prepare_and_run.sh` — end-to-end inference pipeline
- `convert_rnapro_submission.py` — converts RNAPro output back to competition submission format

## Physics Corrections

All pipelines share a common `adaptive_rna_constraints` function that enforces RNA geometry:

| Constraint | Target | Description |
|-----------|--------|-------------|
| Bond length | C1'(i)–C1'(i+1) ~5.95 A | Consecutive backbone distance |
| Next-nearest | C1'(i)–C1'(i+2) ~10.2 A | Soft angle-like restraint |
| Smoothing | Laplacian | Reduces local coordinate noise |
| Clash removal | min 3.2 A separation | Self-avoidance for residues >2 apart |
| Base-pair distance | C1'–C1' ~10.4 A | For ViennaRNA-predicted pairs (optional) |

Constraint strength adapts to prediction confidence — strong for low-confidence de-novo structures, gentle for high-confidence deep learning predictions.

## Models Used

| Model | Type | Source |
|-------|------|--------|
| Protenix v1 | AF3 reproduction | [ByteDance](https://github.com/bytedance/protenix) |
| RNAPro 500M | RNA-specific | [NVIDIA](https://github.com/NVIDIA-Digital-Bio/RNAPro) |
| Boltz-2 | AF3-like | [Boltz](https://github.com/jost-tech/boltz) |
| ViennaRNA | 2D structure (thermodynamic) | [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) |

## Requirements

Designed to run on Kaggle (P100 GPU, 30GB RAM). Dependencies are loaded from Kaggle datasets:
- Biopython 1.86
- Biotite
- RDKit
- US-align
- PyTorch (Kaggle pre-installed)

## Repository Structure

```
.
├── stanford-rna-3d-folding.ipynb          # Main submission: Protenix + TBM
├── rnapro-inference-with-tbm.ipynb        # Tri-model hybrid: RNAPro + Boltz + TBM
├── rna-3d-protenix-biopython-fix.ipynb    # Protenix variant with Biopython patch
└── rnapro_setup/
    ├── setup.sh                           # RNAPro environment setup
    ├── prepare_rnapro_inputs.py           # CSV → RNAPro input conversion
    ├── prepare_and_run.sh                 # RNAPro inference pipeline
    └── convert_rnapro_submission.py       # RNAPro output → submission CSV
```
