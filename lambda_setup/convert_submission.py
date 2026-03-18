#!/usr/bin/env python3
"""
Convert Boltz CIF output files → Kaggle submission CSV.

Submission format:
    ID,resname,resid,x_1,y_1,z_1,...,x_5,y_5,z_5

Where:
    ID      = {target_id}_{resid}  (e.g. 8ZNQ_1)
    resname = nucleotide (A/U/G/C)
    resid   = 1-indexed sequential position across ALL chains
    x_i,y_i,z_i = C1' atom coordinates from the i-th Boltz model
"""

import csv
import json
import sys
from pathlib import Path

import gemmi
import pandas as pd

N_MODELS = 5  # required predictions per target


# --------------------------------------------------------------------------- #
# CIF parsing
# --------------------------------------------------------------------------- #

def extract_c1prime_coords(cif_path: Path) -> list[dict]:
    """
    Parse a mmCIF file and return C1' atom positions per residue, in chain order.

    Returns:
        List of {'chain': str, 'resnum': int, 'resname': str, 'x': float, 'y': float, 'z': float}
        sorted by (chain_label, resnum).  Only standard RNA residues included.
    """
    try:
        structure = gemmi.read_structure(str(cif_path))
    except Exception as e:
        print(f"    [ERROR] Could not read {cif_path.name}: {e}")
        return []

    RNA_RESIDUES = {"A", "U", "G", "C", "ADE", "URA", "GUA", "CYT",
                    "DA", "DU", "DG", "DC"}  # some models use full names

    records = []
    model = structure[0]  # first model
    for chain in model:
        for residue in chain:
            rname = residue.name.strip().upper()
            # Normalise 3-letter codes to 1-letter
            name_map = {"ADE": "A", "URA": "U", "GUA": "G", "CYT": "C",
                        "DA": "A", "DU": "U", "DG": "G", "DC": "C"}
            rname_1 = name_map.get(rname, rname if len(rname) == 1 else None)
            if rname_1 not in ("A", "U", "G", "C"):
                continue  # skip ligands / non-RNA

            c1p = None
            for atom in residue:
                if atom.name == "C1'":
                    c1p = atom
                    break

            if c1p is None:
                # Fill with sentinel; will be zeroed in output
                records.append({
                    "chain": chain.name,
                    "resnum": residue.seqid.num,
                    "resname": rname_1,
                    "x": 0.0, "y": 0.0, "z": 0.0,
                    "missing": True,
                })
            else:
                records.append({
                    "chain": chain.name,
                    "resnum": residue.seqid.num,
                    "resname": rname_1,
                    "x": round(c1p.pos.x, 3),
                    "y": round(c1p.pos.y, 3),
                    "z": round(c1p.pos.z, 3),
                    "missing": False,
                })

    return records


# --------------------------------------------------------------------------- #
# Chain-map to sequential resid
# --------------------------------------------------------------------------- #

def records_to_sequential(records: list[dict],
                           chain_order: list[dict]) -> list[dict]:
    """
    Map chain + resnum → sequential 1-based resid using the chain_order from
    chain_map.json.

    chain_order: [{'chain_id': 'A', 'seq_len': 52}, {'chain_id': 'B', 'seq_len': 52}, ...]

    Returns list of {'resid': int, 'resname': str, 'x': float, 'y': float, 'z': float}.
    """
    # Build offset map: for each Boltz chain ID (A, B, …) → offset in sequential resid
    offset = 0
    chain_offset: dict[str, int] = {}
    for entry in chain_order:
        chain_offset[entry["chain_id"]] = offset
        offset += entry["seq_len"]

    # Sort records by (chain position in chain_order, resnum within chain)
    chain_position = {entry["chain_id"]: i for i, entry in enumerate(chain_order)}

    records_sorted = sorted(
        records,
        key=lambda r: (chain_position.get(r["chain"], 999), r["resnum"])
    )

    results = []
    for rec in records_sorted:
        cid = rec["chain"]
        base_offset = chain_offset.get(cid, 0)
        seq_resid = base_offset + rec["resnum"]
        results.append({
            "resid": seq_resid,
            "resname": rec["resname"],
            "x": rec["x"],
            "y": rec["y"],
            "z": rec["z"],
        })

    return results


# --------------------------------------------------------------------------- #
# Find CIF files for a target
# --------------------------------------------------------------------------- #

def find_cif_files(target_output_dir: Path, target_id: str) -> list[Path]:
    """
    Locate up to N_MODELS CIF files in the Boltz output directory.

    Boltz output layout (varies by version):
        outputs/{target_id}/predictions/{target_id}_model_{i}.cif
        or
        outputs/{target_id}/boltz_results_{target_id}/predictions/...
    """
    cifs: list[Path] = []

    # Search recursively for all .cif files under the target output dir
    all_cifs = sorted(target_output_dir.rglob("*.cif"))

    # Prefer files with 'model_' in the name
    model_cifs = [c for c in all_cifs if "model_" in c.name]
    if model_cifs:
        cifs = model_cifs[:N_MODELS]
    elif all_cifs:
        cifs = all_cifs[:N_MODELS]

    return cifs


# --------------------------------------------------------------------------- #
# Build submission rows for one target
# --------------------------------------------------------------------------- #

def process_target(target_id: str,
                   chain_order: list[dict],
                   outputs_dir: Path,
                   test_rows_map: dict) -> list[dict] | None:
    """
    Returns list of submission row dicts, or None on failure.
    """
    target_out = outputs_dir / target_id
    if not target_out.exists():
        print(f"  [MISSING] {target_id}: output directory not found.")
        return None

    cif_files = find_cif_files(target_out, target_id)

    if not cif_files:
        print(f"  [MISSING] {target_id}: no CIF files found in {target_out}")
        return None

    print(f"  {target_id}: found {len(cif_files)} CIF(s)")

    # Extract coords from each model
    all_model_coords: list[list[dict]] = []
    for cif in cif_files:
        records = extract_c1prime_coords(cif)
        if not records:
            print(f"    [WARN] {cif.name}: no residues extracted")
            continue
        sequential = records_to_sequential(records, chain_order)
        all_model_coords.append(sequential)

    if not all_model_coords:
        print(f"  [ERROR] {target_id}: all CIF extractions failed")
        return None

    # Pad to N_MODELS by repeating last model
    while len(all_model_coords) < N_MODELS:
        all_model_coords.append(all_model_coords[-1])

    # Get expected sequence from test data
    full_seq = test_rows_map.get(target_id, {}).get("sequence", "")
    n_residues = len(full_seq)
    if n_residues == 0:
        # Fall back to max resid seen
        n_residues = max(r["resid"] for model in all_model_coords for r in model)

    # Build one row per residue
    submission_rows = []
    for res_idx in range(1, n_residues + 1):
        resname = full_seq[res_idx - 1] if full_seq else "N"
        row = {
            "ID": f"{target_id}_{res_idx}",
            "resname": resname,
            "resid": res_idx,
        }
        # Fill coords from each model
        for model_num, model_coords in enumerate(all_model_coords[:N_MODELS], start=1):
            # Look up this resid in model
            coord = next((r for r in model_coords if r["resid"] == res_idx), None)
            if coord:
                row[f"x_{model_num}"] = coord["x"]
                row[f"y_{model_num}"] = coord["y"]
                row[f"z_{model_num}"] = coord["z"]
            else:
                row[f"x_{model_num}"] = 0.0
                row[f"y_{model_num}"] = 0.0
                row[f"z_{model_num}"] = 0.0

        submission_rows.append(row)

    return submission_rows


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    base_dir = Path(__file__).parent
    outputs_dir = base_dir / "outputs"
    chain_map_path = base_dir / "chain_map.json"
    test_csv_path = base_dir / "test_sequences.csv"
    sample_csv_path = base_dir / "sample_submission.csv"
    out_csv = base_dir / "submission.csv"

    # Load chain map
    if not chain_map_path.exists():
        print("ERROR: chain_map.json not found. Run prepare_inputs.py first.", file=sys.stderr)
        sys.exit(1)
    with open(chain_map_path) as f:
        chain_map: dict[str, list[dict]] = json.load(f)

    # Load test sequences for resname lookup
    test_rows_map: dict[str, dict] = {}
    if test_csv_path.exists():
        with open(test_csv_path) as f:
            for row in csv.DictReader(f):
                test_rows_map[row["target_id"]] = row

    # Build predictions keyed by ID for all targets
    pred_by_id: dict[str, dict] = {}
    zero_filled_targets: list[str] = []

    for target_id, chain_order in chain_map.items():
        rows = process_target(target_id, chain_order, outputs_dir, test_rows_map)
        if rows is None:
            zero_filled_targets.append(target_id)
        else:
            for row in rows:
                pred_by_id[row["ID"]] = row

    # Use sample_submission.csv as template for exact row ordering
    if not sample_csv_path.exists():
        print("ERROR: sample_submission.csv not found.", file=sys.stderr)
        sys.exit(1)

    COLUMNS = (["ID", "resname", "resid"] +
               [f"{x}_{i}" for i in range(1, N_MODELS + 1) for x in ("x", "y", "z")])

    all_rows: list[dict] = []
    with open(sample_csv_path) as f:
        for sample_row in csv.DictReader(f):
            sid = sample_row["ID"]
            if sid in pred_by_id:
                # Use our coords, keep sample's ID/resname/resid
                row = {"ID": sid, "resname": sample_row["resname"],
                       "resid": sample_row["resid"]}
                for col in COLUMNS:
                    if col not in ("ID", "resname", "resid"):
                        row[col] = pred_by_id[sid].get(col, 0.0)
                all_rows.append(row)
            else:
                # Zero-fill missing predictions
                row = {"ID": sid, "resname": sample_row["resname"],
                       "resid": sample_row["resid"]}
                for m in range(1, N_MODELS + 1):
                    row[f"x_{m}"] = 0.0
                    row[f"y_{m}"] = 0.0
                    row[f"z_{m}"] = 0.0
                all_rows.append(row)

    # Write CSV via pandas (LF line endings, matching Kaggle scorer expectations)
    df = pd.DataFrame(all_rows, columns=COLUMNS)
    df.to_csv(out_csv, index=False)

    predicted_count = len(chain_map) - len(zero_filled_targets)
    print(f"\n=== Submission written: {out_csv} ===")
    print(f"  Total rows: {len(all_rows)}")
    print(f"  Targets with predictions: {predicted_count}")
    if zero_filled_targets:
        print(f"  Zero-filled (no output): {zero_filled_targets}")


if __name__ == "__main__":
    main()
