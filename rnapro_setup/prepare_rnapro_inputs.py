#!/usr/bin/env python3
"""
Convert test_sequences.csv → RNAPro sequences.csv format.
Also generates placeholder MSA files (single-sequence MSAs) for targets
without precomputed MSAs.

RNAPro sequences.csv format: target_id,sequence
RNAPro MSA format: {msa_dir}/{target_id}.MSA.fasta
"""

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--test_csv", required=True)
    parser.add_argument("--out_csv", required=True)
    parser.add_argument("--msa_dir", required=True)
    args = parser.parse_args()

    msa_dir = Path(args.msa_dir)
    msa_dir.mkdir(exist_ok=True)

    with open(args.test_csv) as f:
        rows = list(csv.DictReader(f))

    with open(args.out_csv, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["target_id", "sequence"])

        for row in rows:
            tid = row["target_id"]
            seq = row["sequence"].upper().replace("T", "U")
            writer.writerow([tid, seq])

            # Create single-sequence MSA if none exists
            msa_file = msa_dir / f"{tid}.MSA.fasta"
            if not msa_file.exists():
                with open(msa_file, "w") as mf:
                    mf.write(f">{tid}\n{seq}\n")

    print(f"Written {len(rows)} targets to {args.out_csv}")
    print(f"MSA files in {msa_dir}")


if __name__ == "__main__":
    main()
