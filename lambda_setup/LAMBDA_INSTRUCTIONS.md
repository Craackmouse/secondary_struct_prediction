# Lambda Labs H100 — RNA Structure Prediction Setup

## What this does
Runs **Boltz** (open-source AlphaFold3-like model) on 28 RNA test sequences
and generates a Kaggle-ready submission.csv with 5 predicted structures per target.

Expected runtime on H100: ~2–4 hours total (varies by MSA server latency).

---

## Step 1 — Launch a Lambda Labs instance
- GPU: H100 (80 GB) recommended; A100 also works
- Image: Ubuntu 22.04 LTS (Lambda base image)
- Storage: at least 20 GB (for model weights + outputs)

---

## Step 2 — Upload files from your Mac

```bash
# From your local machine (replace <LAMBDA_IP>):
scp -r lambda_setup/ ubuntu@<LAMBDA_IP>:~/lambda_setup/
scp data_extracted/test_sequences.csv ubuntu@<LAMBDA_IP>:~/lambda_setup/
```

---

## Step 3 — SSH in and run setup

```bash
ssh ubuntu@<LAMBDA_IP>

cd ~/lambda_setup
bash setup.sh
```

---

## Step 4 — Run the full pipeline

```bash
cd ~/lambda_setup
conda activate rna_pred
bash run_all.sh
```

This runs three steps automatically:
1. `prepare_inputs.py` — parses test_sequences.csv → Boltz YAML inputs (handles multi-chain, ligands)
2. `run_boltz.sh` — runs Boltz with 5 diffusion samples per target, using ColabFold MSA server
3. `convert_submission.py` — extracts C1' atom coordinates → submission.csv

Progress and errors are logged to `logs/{target_id}.log`.
Failed targets are listed in `failed_targets.txt` and zero-filled in submission.csv.

---

## Step 5 — Download submission

```bash
# From your local machine:
scp ubuntu@<LAMBDA_IP>:~/lambda_setup/submission.csv .
```

Then submit on Kaggle.

---

## Tips

**Resume interrupted run:** `run_boltz.sh` skips targets that already have 5 CIF outputs,
so you can safely re-run after interruption.

**Skip MSA server (faster, lower quality):**
```bash
bash run_all.sh --no-msa
```

**Check progress:**
```bash
ls outputs/          # completed targets
tail -f logs/9MME.log  # live log for long-running target
```

**Monitor GPU:**
```bash
watch -n2 nvidia-smi
```

---

## Troubleshooting

**`ModuleNotFoundError`**: Run `conda activate rna_pred` first.

**MSA server timeout**: Re-run `run_boltz.sh`; it will skip already-completed targets.

**9MME (4640 nt, 8 chains)**: This is the largest target and may take 30–60 min alone
or OOM on A100. On H100 it should fit. If it fails, it'll be zero-filled — acceptable
for a first submission.

**Submission row count**: Should be 9762 rows (matching sample_submission.csv).

---

## File structure after pipeline completes

```
lambda_setup/
├── test_sequences.csv       # competition test data (you upload this)
├── inputs/                  # Boltz YAML inputs (one per target)
├── outputs/                 # Boltz prediction outputs (CIF files)
│   ├── 8ZNQ/
│   │   └── predictions/
│   │       ├── 8ZNQ_model_0.cif
│   │       └── ...
│   └── ...
├── chain_map.json           # chain ordering metadata
├── logs/                    # per-target inference logs
├── failed_targets.txt       # targets that failed (if any)
└── submission.csv           # FINAL SUBMISSION FILE ← this is what you need
```
