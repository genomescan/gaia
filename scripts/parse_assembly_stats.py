#!/usr/bin/env python3
"""
parse_assembly_stats.py – extract assembly quality metrics (N50, contig
count, total length, largest contig) for a single sample.

Reuses existing pipeline outputs instead of re-running any tools:
  - If MetaQUAST was run for the sample (assembly_qc.metaquast: true), its
    ``report.tsv`` already contains N50 and other QC metrics, so those are
    used directly.
  - Otherwise, falls back to computing the same metrics from Flye's
    ``assembly_info.txt`` (always produced by the metaflye_assemble rule),
    using the contig lengths listed there.

Outputs a single JSON file:
{
  "sample": ...,
  "source": "metaquast" | "assembly_info" | "none",
  "contig_count": ...,
  "total_length": ...,
  "n50": ...,
  "largest_contig": ...
}
"""
import argparse
import csv
import json
import os


def _exists_nonempty(path):
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def _n50(lengths):
    """Compute N50 from a list of contig/sequence lengths."""
    if not lengths:
        return 0
    lengths_sorted = sorted(lengths, reverse=True)
    total = sum(lengths_sorted)
    csum = 0
    for length in lengths_sorted:
        csum += length
        if csum >= total / 2:
            return length
    return lengths_sorted[-1]


def parse_assembly_info(path):
    """Parse Flye's assembly_info.txt (tab-separated, header starting with
    '#') into a list of contig lengths (column 2)."""
    if not _exists_nonempty(path):
        return []
    lengths = []
    with open(path, newline="") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                lengths.append(int(parts[1].strip()))
            except ValueError:
                continue
    return lengths


def stats_from_lengths(lengths):
    if not lengths:
        return {}
    return {
        "contig_count": len(lengths),
        "total_length": sum(lengths),
        "n50": _n50(lengths),
        "largest_contig": max(lengths),
    }


# ---------------------------------------------------------------------------
# MetaQUAST/QUAST report.tsv parser
#
# Format: two columns, tab-separated, no reliable header ("Assembly" \t
# <assembly-name>), one metric per row, e.g.:
#   # contigs (>= 0 bp)   123
#   Largest contig        45000
#   Total length (>= 0 bp) 5000000
#   N50                   12345
# ---------------------------------------------------------------------------

_METAQUAST_KEYS = {
    "# contigs (>= 0 bp)": "contig_count",
    "Largest contig": "largest_contig",
    "Total length (>= 0 bp)": "total_length",
    "N50": "n50",
}


def parse_metaquast_report(path):
    if not _exists_nonempty(path):
        return {}
    stats = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            key = row[0].strip()
            field = _METAQUAST_KEYS.get(key)
            if not field:
                continue
            try:
                stats[field] = int(float(row[1].strip()))
            except ValueError:
                continue
    return stats


def main():
    ap = argparse.ArgumentParser(
        description="Extract assembly quality metrics (N50, etc.) for one sample."
    )
    ap.add_argument("--sample", required=True, help="Sample name")
    ap.add_argument("--assembly-info", default="",
                    help="Path to Flye assembly_info.txt")
    ap.add_argument("--metaquast-report", default="",
                    help="Path to MetaQUAST/QUAST report.tsv (optional)")
    ap.add_argument("--output", required=True, help="Output JSON path")
    args = ap.parse_args()

    stats = parse_metaquast_report(args.metaquast_report)
    source = "metaquast" if stats.get("n50") else "none"

    if not stats.get("n50"):
        lengths = parse_assembly_info(args.assembly_info)
        computed = stats_from_lengths(lengths)
        if computed:
            stats = computed
            source = "assembly_info"

    result = {"sample": args.sample, "source": source}
    result.update(stats)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"Assembly stats written to: {args.output}")


if __name__ == "__main__":
    main()
