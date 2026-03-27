#!/usr/bin/env python3
"""
Convert RNAPro CIF output → Kaggle submission CSV.
Uses sample_submission.csv as template for exact row ordering.
"""

import argparse
import csv
import sys
from pathlib import Path

import gemmi

N_MODELS = 5


def extract_c1prime_coords(cif_path: Path) -> dict[int, dict]:
    """Parse CIF, return {resid: {resname, x, y, z}} for RNA C1' atoms."""
    try:
        structure = gemmi.read_structure(str(cif_path))
    except Exception as e:
        print(f"    [ERROR] {cif_path.name}: {e}")
        return {}

    name_map = {"ADE": "A", "URA": "U", "GUA": "G", "CYT": "C",
                "DA": "A", "DU": "U", "DG": "G", "DC": "C"}
    coords = {}
    resid_counter = 0

    model = structure[0]
    for chain in model:
        for residue in chain:
            rname = residue.name.strip().upper()
            rname_1 = name_map.get(rname, rname if len(rname) == 1 else None)
            if rname_1 not in ("A", "U", "G", "C"):
                continue
            resid_counter += 1
            c1p = next((a for a in residue if a.name == "C1'"), None)
            if c1p:
                coords[resid_counter] = {
                    "resname": rname_1,
                    "x": round(c1p.pos.x, 3),
                    "y": round(c1p.pos.y, 3),
                    "z": round(c1p.pos.z, 3),
                }
            else:
                coords[resid_counter] = {
                    "resname": rname_1, "x": 0.0, "y": 0.0, "z": 0.0,
                }
    return coords


def collect_rnapro_predictions(rnapro_output: Path) -> dict[str, dict]:
    """Map submission row ID -> {x_1, y_1, ...} from RNAPro CIF tree."""
    target_cifs: dict[str, list[Path]] = {}
    if not rnapro_output.exists():
        return {}

    for cif in sorted(rnapro_output.rglob("*.cif")):
        tid = cif.stem.split("_")[0] if "_" in cif.stem else cif.stem
        parent_tid = cif.parent.name
        for candidate in (parent_tid, tid):
            if len(candidate) == 4 and candidate[0].isdigit():
                tid = candidate
                break
        target_cifs.setdefault(tid, []).append(cif)

    pred_by_id: dict[str, dict] = {}
    for tid, cifs in target_cifs.items():
        cifs = cifs[:N_MODELS]
        all_model_coords = []
        for cif in cifs:
            coords = extract_c1prime_coords(cif)
            if coords:
                all_model_coords.append(coords)

        while len(all_model_coords) < N_MODELS:
            if all_model_coords:
                all_model_coords.append(all_model_coords[-1])
            else:
                break

        if not all_model_coords:
            print(f"  [SKIP] {tid}: no valid CIFs")
            continue

        print(f"  {tid}: {len(cifs)} CIFs")

        max_resid = max(max(m.keys()) for m in all_model_coords)
        for resid in range(1, max_resid + 1):
            row_id = f"{tid}_{resid}"
            for mi, model in enumerate(all_model_coords, 1):
                c = model.get(resid, {"x": 0.0, "y": 0.0, "z": 0.0})
                pred_by_id.setdefault(row_id, {})[f"x_{mi}"] = c["x"]
                pred_by_id[row_id][f"y_{mi}"] = c["y"]
                pred_by_id[row_id][f"z_{mi}"] = c["z"]

    return pred_by_id


def write_submission_from_template(
    sample_csv: Path,
    out_csv: Path,
    pred_by_id: dict[str, dict],
) -> tuple[int, int]:
    """Fill coordinate columns from pred_by_id; zeros for missing rows."""
    COLUMNS = (["ID", "resname", "resid"] +
               [f"{x}_{i}" for i in range(1, N_MODELS + 1) for x in ("x", "y", "z")])

    all_rows = []
    with open(sample_csv) as f:
        for sample_row in csv.DictReader(f):
            sid = sample_row["ID"]
            row = {"ID": sid, "resname": sample_row["resname"],
                   "resid": sample_row["resid"]}
            if sid in pred_by_id:
                for col in COLUMNS:
                    if col not in ("ID", "resname", "resid"):
                        row[col] = pred_by_id[sid].get(col, 0.0)
            else:
                for m in range(1, N_MODELS + 1):
                    row[f"x_{m}"] = 0.0
                    row[f"y_{m}"] = 0.0
                    row[f"z_{m}"] = 0.0
            all_rows.append(row)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(all_rows)

    predicted = sum(
        1 for r in all_rows
        if any(float(r.get(f"x_{i}", 0)) != 0 for i in range(1, N_MODELS + 1))
    )
    return len(all_rows), predicted


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    default_sample = repo_root / "lambda_setup" / "sample_submission.csv"
    default_out = Path(__file__).resolve().parent / "submission.csv"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample-submission",
        type=Path,
        default=default_sample,
        help="sample_submission.csv from the competition (row order template)",
    )
    parser.add_argument(
        "--rnapro-output",
        type=Path,
        default=Path.home() / "RNAPro" / "output",
        help="Directory containing RNAPro *.cif outputs",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=default_out,
        help="Output submission path",
    )
    args = parser.parse_args()

    if not args.sample_submission.exists():
        print("ERROR: sample_submission.csv not found", file=sys.stderr)
        sys.exit(1)

    pred_by_id = collect_rnapro_predictions(args.rnapro_output)
    n_rows, predicted = write_submission_from_template(
        args.sample_submission, args.out_csv, pred_by_id
    )
    print(f"\n=== Submission: {args.out_csv} ===")
    print(f"  Total rows: {n_rows}")
    print(f"  Rows with predictions: {predicted}")


if __name__ == "__main__":
    main()
