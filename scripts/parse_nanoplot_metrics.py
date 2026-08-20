#!/usr/bin/env python3
"""
parse_nanoplot_metrics.py – extract summary metrics from NanoPlot NanoStats.txt
files for both raw and filtered reads.

Outputs a JSON file with the following structure:
{
  "raw": {"mean_read_length": ..., "median_read_length": ..., ...},
  "filtered": {"mean_read_length": ..., ...}
}
Values are numeric where possible; absent sections are represented as {}.
"""
import argparse
import json
import os
import re


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _exists_nonempty(path):
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


# ---------------------------------------------------------------------------
# NanoStats.txt parser
#
# NanoPlot writes a plain-text file with lines such as:
#   Mean read length:                 12,345.6
#   Median read length:               8,901.2
#   Read length N50:                  15,432
#   Mean read quality:                10.5
#   Median read quality:              11.0
#   Active channels:                  512
#   Number of reads:                  12345
#   Total bases:                      123456789
#   Number, percentage and megabases of reads above quality cutoffs ...
# ---------------------------------------------------------------------------

_FIELD_MAP = {
    "mean read length": "mean_read_length",
    "median read length": "median_read_length",
    "read length n50": "read_length_n50",
    "mean read quality": "mean_read_quality",
    "median read quality": "median_read_quality",
    "active channels": "active_channels",
    "number of reads": "number_of_reads",
    "total bases": "total_bases",
    "total bases (mb)": "total_bases_mb",
}


def parse_nanostats(path):
    """Return a dict of metrics parsed from a NanoStats.txt file."""
    if not _exists_nonempty(path):
        return {}
    metrics = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if ":" not in line:
                continue
            key_raw, _, value_raw = line.partition(":")
            key = key_raw.strip().lower()
            if key not in _FIELD_MAP:
                continue
            # Remove commas used as thousands separators, then convert
            value_str = value_raw.strip().replace(",", "")
            try:
                value = float(value_str)
                # Use int when value is a whole number for cleaner JSON
                if value == int(value):
                    value = int(value)
            except ValueError:
                value = value_str
            metrics[_FIELD_MAP[key]] = value
    return metrics


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Parse NanoPlot NanoStats.txt files into a summary JSON."
    )
    ap.add_argument("--raw-nanostats", default="",
                    help="Path to raw NanoStats.txt (e.g. QC/Raw/<sample>/<sample>_NanoStats.txt)")
    ap.add_argument("--filtered-nanostats", default="",
                    help="Path to filtered NanoStats.txt (e.g. QC/Filtered/<sample>/<sample>_NanoStats.txt)")
    ap.add_argument("--output", required=True, help="Output JSON file path")
    args = ap.parse_args()

    result = {
        "raw": parse_nanostats(args.raw_nanostats),
        "filtered": parse_nanostats(args.filtered_nanostats),
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2)


if __name__ == "__main__":
    main()
