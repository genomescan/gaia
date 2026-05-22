#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser(description="Convert a directory of bin FASTA files to a DAS Tool-compatible contig2bin TSV.")
    p.add_argument("--bins", required=True, help="Directory containing bin FASTA files.")
    p.add_argument("--output", required=True, help="Output contig2bin TSV.")
    p.add_argument("--tool", default="", help="Optional tool name.")
    return p.parse_args()

def is_fasta(path: Path) -> bool:
    suffixes = "".join(path.suffixes).lower()
    return path.is_file() and (
        path.suffix.lower() in {".fa", ".fna", ".fasta", ".fas"} or
        suffixes.endswith(".fa.gz") or suffixes.endswith(".fna.gz") or suffixes.endswith(".fasta.gz")
    )

def iter_headers(path: Path):
    suffixes = "".join(path.suffixes).lower()
    opener = gzip.open if suffixes.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip().split()[0]

def main():
    args = parse_args()
    bins_dir = Path(args.bins)
    out_path = Path(args.output)

    if not bins_dir.exists():
        raise SystemExit(f"ERROR: bin directory does not exist: {bins_dir}")

    fasta_files = sorted([p for p in bins_dir.iterdir() if is_fasta(p)])
    if not fasta_files:
        raise SystemExit(f"ERROR: no FASTA bin files found in: {bins_dir}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    n_rows = 0
    with out_path.open("w") as out:
        for fasta in fasta_files:
            bin_id = fasta.name
            for ext in (".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta", ".fas"):
                if bin_id.endswith(ext):
                    bin_id = bin_id[:-len(ext)]
                    break
            for contig_id in iter_headers(fasta):
                key = (contig_id, bin_id)
                if key in seen:
                    continue
                seen.add(key)
                out.write(f"{contig_id}\t{bin_id}\n")
                n_rows += 1

    if n_rows == 0:
        raise SystemExit(f"ERROR: no contig headers were parsed from FASTA files in: {bins_dir}")

if __name__ == "__main__":
    main()
