#!/usr/bin/env python3
"""
parse_taxonomy.py – extract top-N species from Kraken2 and/or Centrifuger reports.

Outputs a JSON file with the following structure:
{
  "kraken2": [{"name": ..., "reads": ..., "percent": ...}, ...],
  "centrifuger": [{"name": ..., "reads": ..., "percent": ...}, ...]
}
Each list is sorted descending by reads and capped at --top-n entries.
"""
import argparse
import csv
import json
import os


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _exists_nonempty(path):
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


# ---------------------------------------------------------------------------
# Kraken2 report parser
#
# Format (tab-separated, no header):
#   %reads  clade_reads  direct_reads  rank  taxid  name
# ---------------------------------------------------------------------------

def parse_kraken2(path, top_n=10):
    """Return top-n species-level rows sorted by clade reads (descending)."""
    if not _exists_nonempty(path):
        return []
    species = []
    with open(path, newline="") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rank = parts[3].strip()
            if rank != "S":
                continue
            try:
                pct = round(float(parts[0].strip()), 1)
                reads = int(parts[1].strip())
            except ValueError:
                continue
            name = parts[5].strip()
            species.append({"name": name, "reads": reads, "percent": pct})
    species.sort(key=lambda x: x["reads"], reverse=True)
    return species[:top_n]


# ---------------------------------------------------------------------------
# Centrifuger report parser
#
# centrifuger-quant output (tab-separated, with header):
#   name  taxID  taxRank  genomeSize  numReads  numUniqueReads  abundance
# ---------------------------------------------------------------------------

def parse_centrifuger(path, top_n=10):
    """Return top-n species-level rows sorted by numReads (descending)."""
    if not _exists_nonempty(path):
        return []
    species = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return []
        for row in reader:
            rank = row.get("taxRank", "").strip()
            if rank not in ("species", "S"):
                continue
            try:
                reads = int(row.get("numReads", 0))
                abundance = float(row.get("abundance", 0.0))
            except (ValueError, TypeError):
                continue
            name = row.get("name", "").strip()
            if not name or name in ("unclassified", "root"):
                continue
            species.append({"name": name, "reads": reads, "percent": round(abundance * 100, 1)})
    species.sort(key=lambda x: x["reads"], reverse=True)
    return species[:top_n]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Parse taxonomy reports into top-N JSON.")
    ap.add_argument("--kraken2-report", default="", help="Kraken2 report file path")
    ap.add_argument("--centrifuger-report", default="", help="Centrifuger quant report path")
    ap.add_argument("--top-n", type=int, default=10, help="Number of top species to include")
    ap.add_argument("--output", required=True, help="Output JSON file path")
    args = ap.parse_args()

    result = {
        "kraken2": parse_kraken2(args.kraken2_report, args.top_n),
        "centrifuger": parse_centrifuger(args.centrifuger_report, args.top_n),
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2)


if __name__ == "__main__":
    main()
