#!/usr/bin/env python3
"""
parse_host_removal.py – Parse samtools flagstat output from host removal and
produce a JSON summary with host/non-host read counts and percentages.
"""
import argparse
import json
import os
import re


def parse_flagstat(flagstat_path):
    """Return (total_reads, mapped_reads) from a samtools flagstat file."""
    total = 0
    mapped = 0
    with open(flagstat_path) as fh:
        for line in fh:
            m_total = re.match(r"(\d+)\s+\+\s+\d+\s+in total", line)
            if m_total:
                total = int(m_total.group(1))
            m_mapped = re.match(r"(\d+)\s+\+\s+\d+\s+mapped", line)
            if m_mapped:
                mapped = int(m_mapped.group(1))
    return total, mapped


def main():
    ap = argparse.ArgumentParser(
        description="Parse samtools flagstat and generate host removal stats JSON."
    )
    ap.add_argument("--sample", required=True, help="Sample name")
    ap.add_argument("--flagstat", required=True, help="Path to samtools flagstat output")
    ap.add_argument("--output", required=True, help="Output JSON path")
    args = ap.parse_args()

    total, host_reads = parse_flagstat(args.flagstat)

    non_host_reads = max(0, total - host_reads)
    host_pct = round(100.0 * host_reads / total, 1) if total > 0 else 0.0
    non_host_pct = round(100.0 * non_host_reads / total, 1) if total > 0 else 0.0

    stats = {
        "sample": args.sample,
        "total_reads": total,
        "host_reads": host_reads,
        "host_percentage": host_pct,
        "non_host_reads": non_host_reads,
        "non_host_percentage": non_host_pct,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(stats, fh, indent=2)

    print(f"Host removal stats written to: {args.output}")


if __name__ == "__main__":
    main()
