#!/usr/bin/env python3
"""
collect_versions.py – Collect tool versions from containers and write a
versions.json file to the pipeline run directory.

Usage (called by Snakemake rule):
    python collect_versions.py --tools nanoplot nanoqc minimap2 samtools chopper filtlong \
        --output /path/to/versions.json
"""
import argparse
import json
import os
import re
import subprocess


# Map tool name → version extraction command and regex pattern
TOOL_VERSION_CMDS = {
    "nanoplot": ("NanoPlot --version", r"NanoPlot\s+([\d.]+)"),
    "nanoqc": ("nanoQC --version", r"([\d.]+)"),
    "minimap2": ("minimap2 --version", r"([\d.]+(?:-r\d+)?)"),
    "samtools": ("samtools --version", r"samtools\s+([\d.]+)"),
    "chopper": ("chopper --version", r"([\d.]+)"),
    "filtlong": ("filtlong --version", r"Filtlong\s+v?([\d.]+)"),
    "kraken2": ("kraken2 --version", r"Kraken\s+version\s+([\d.]+)"),
    "centrifuger": ("centrifuger --version", r"([\d.]+)"),
    "metaflye": ("flye --version", r"([\d.]+)"),
    "metaquast": ("metaquast.py --version", r"([\d.]+)"),
    "jgi_summarize_bam_contig_depths": (
        "jgi_summarize_bam_contig_depths 2>&1 | head -1",
        r"([\d.]+)",
    ),
    "metabat2": ("metabat2 --help 2>&1 | head -1", r"([\d.]+)"),
    "semibin2": ("SemiBin2 --version", r"SemiBin2\s+([\d.]+)"),
    "comebin": ("run_comebin.sh --version 2>&1 | head -1", r"([\d.]+)"),
    "dastool": ("DAS_Tool --version 2>&1", r"DAS_Tool\s+([\d.]+)"),
    "checkm2": ("checkm2 --version", r"([\d.]+)"),
    "gtdbtk": ("gtdbtk --version", r"gtdbtk:\s*([\d.]+)"),
    "eukcc": ("eukcc --version", r"eukcc\s+([\d.]+)"),
    "bat": ("CAT_pack --version 2>&1 | head -1", r"([\d.]+)"),
    "drep": ("dRep --version", r"dRep\s+([\d.]+)"),
}


def get_version(tool_name):
    """Try to get the version string for a tool. Returns 'unknown' on failure."""
    if tool_name not in TOOL_VERSION_CMDS:
        return "unknown"
    cmd, pattern = TOOL_VERSION_CMDS[tool_name]
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=15,
        )
        output = result.stdout + result.stderr
        m = re.search(pattern, output)
        if m:
            return m.group(1)
    except Exception:
        pass
    return "unknown"


def main():
    ap = argparse.ArgumentParser(
        description="Collect tool versions and write versions.json."
    )
    ap.add_argument(
        "--tools",
        nargs="+",
        required=True,
        help="Tool names to collect versions for",
    )
    ap.add_argument("--output", required=True, help="Output versions.json path")
    args = ap.parse_args()

    versions = {}
    for tool in args.tools:
        ver = get_version(tool)
        versions[tool] = f"{tool} {ver}"
        print(f"  {tool}: {ver}")

    # Merge with existing versions.json if present (for incremental updates)
    if os.path.exists(args.output):
        try:
            with open(args.output) as fh:
                existing = json.load(fh)
            existing.update(versions)
            versions = existing
        except Exception:
            pass

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(versions, fh, indent=2)

    print(f"versions.json written to: {args.output}")


if __name__ == "__main__":
    main()
