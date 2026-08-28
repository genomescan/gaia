#!/usr/bin/env python3
"""Gate check run right after assembly, before binning starts.

Counts how many contigs in a MetaFlye ``assembly_info.txt`` file are at
least ``--min-length`` bp long, and compares that count against
``--min-contigs``. The result ("PASS" or "SKIP") is written to
``--output`` so that downstream binning rules can decide whether to run.

This script never fails the workflow: if the threshold is not met it
simply reports "SKIP" so that binning-related steps can be skipped
cleanly instead of erroring out on a too-fragmented/too-small assembly.
"""
import argparse
import sys


def count_long_contigs(assembly_info_path, min_length):
    """Count contigs >= min_length bp in a Flye assembly_info.txt file.

    The file is tab-separated with a header line starting with
    '#seq_name' and the contig length in the second column.
    """
    count = 0
    with open(assembly_info_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            try:
                length = int(fields[1])
            except ValueError:
                continue
            if length >= min_length:
                count += 1
    return count


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Count contigs >= a minimum length in a Flye assembly_info.txt "
            "file and decide whether binning should proceed."
        )
    )
    parser.add_argument(
        "--assembly-info", required=True, help="Path to MetaFlye assembly_info.txt"
    )
    parser.add_argument(
        "--min-length", type=int, default=1500, help="Minimum contig length in bp"
    )
    parser.add_argument(
        "--min-contigs",
        type=int,
        default=10,
        help="Minimum number of qualifying contigs required to proceed with binning",
    )
    parser.add_argument(
        "--output", required=True, help="Path to write the gate status (PASS/SKIP)"
    )
    args = parser.parse_args()

    count = count_long_contigs(args.assembly_info, args.min_length)
    status = "PASS" if count >= args.min_contigs else "SKIP"

    with open(args.output, "w") as out:
        out.write(f"{status}\n")
        out.write(f"contigs_at_or_above_{args.min_length}bp={count}\n")
        out.write(f"min_contigs_required={args.min_contigs}\n")

    if status == "SKIP":
        sys.stderr.write(
            f"WARNING: assembly has only {count} contig(s) >= {args.min_length} bp "
            f"(minimum required: {args.min_contigs}). "
            "Skipping downstream binning steps for this sample.\n"
        )
    else:
        sys.stderr.write(
            f"Assembly has {count} contig(s) >= {args.min_length} bp "
            f"(minimum required: {args.min_contigs}). Proceeding with binning.\n"
        )


if __name__ == "__main__":
    main()
