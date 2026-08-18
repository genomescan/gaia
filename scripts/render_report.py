#!/usr/bin/env python3
"""
render_report.py – generate a standalone HTML report for a gaia pipeline run.

Combines:
  - Pipeline metadata (mode, sample, date)
  - Taxonomy top-10 species per tool (bar chart via Chart.js CDN)
  - Assembly / genome-bin summary tables
  - Genome inventory table

Inputs are all optional; sections are shown only when the relevant file
exists and is non-empty.

The HTML is rendered via a Jinja2 template located at
``scripts/templates/report.html.j2`` (resolved relative to this script).
"""
import argparse
import csv
import json
import os
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_TEMPLATES_DIR = Path(__file__).parent.parent / "templates" / "report"


def _exists_nonempty(path):
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def _read_table(path):
    """Return list of dicts from a tab/comma-delimited file."""
    if not _exists_nonempty(path):
        return []
    sample_bytes = open(path, "rb").read(4096)
    delim = "\t"
    try:
        delim = csv.Sniffer().sniff(
            sample_bytes.decode("utf-8", errors="replace"), delimiters="\t,;"
        ).delimiter
    except Exception:
        pass
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        return list(reader) if reader.fieldnames else []


def _read_json(path):
    if not _exists_nonempty(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

def render(sample, run_mode, taxonomy_json, genome_summary,
           genome_inventory, pipeline_summary, output):
    tax = _read_json(taxonomy_json)

    context = {
        "sample": sample,
        "run_mode": run_mode,
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "kraken2_species": tax.get("kraken2", []),
        "centrifuger_species": tax.get("centrifuger", []),
        "genome_summary": _read_table(genome_summary),
        "genome_inventory": _read_table(genome_inventory),
        "pipeline_summary": _read_table(pipeline_summary),
    }

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATES_DIR)),
        autoescape=select_autoescape(["html"]),
    )
    template = env.get_template("report.html.j2")
    html = template.render(**context)

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, "w", encoding="utf-8") as fh:
        fh.write(html)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Render a standalone HTML report for a gaia pipeline run."
    )
    ap.add_argument("--sample", required=True)
    ap.add_argument("--run-mode", required=True,
                    choices=["profiling_only", "assembly_binning_only", "both"])
    ap.add_argument("--taxonomy-json", default="",
                    help="JSON produced by parse_taxonomy.py")
    ap.add_argument("--genome-summary", default="")
    ap.add_argument("--genome-inventory", default="")
    ap.add_argument("--pipeline-summary", default="")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    render(
        sample=args.sample,
        run_mode=args.run_mode,
        taxonomy_json=args.taxonomy_json,
        genome_summary=args.genome_summary,
        genome_inventory=args.genome_inventory,
        pipeline_summary=args.pipeline_summary,
        output=args.output,
    )
    print(f"Report written to: {args.output}")


if __name__ == "__main__":
    main()
