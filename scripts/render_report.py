#!/usr/bin/env python3
"""
render_report.py – generate a standalone HTML report for a gaia pipeline run.

Combines:
  - Pipeline metadata (mode, sample, date)
  - Preprocessing summary (filtering method, host removal stats, versions)
  - Taxonomy top-10 species per tool (bar chart via Chart.js CDN)
  - Assembly / genome-bin summary tables
  - Genome inventory table
  - Embedded workflow diagram (_workflow_.png as base64)

Inputs are all optional; sections are shown only when the relevant file
exists and is non-empty.

The HTML is rendered via a Jinja2 template located at
``templates/report/report.html.j2`` (resolved relative to this script).
"""
import argparse
import base64
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


def _encode_image_base64(path):
    """Return base64-encoded data URI for an image, or empty string if not found."""
    if not path or not os.path.exists(path):
        return ""
    ext = os.path.splitext(path)[1].lstrip(".").lower()
    mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg",
            "svg": "image/svg+xml"}.get(ext, "image/png")
    with open(path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("ascii")
    return f"data:{mime};base64,{data}"


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

def render(
    sample,
    run_mode,
    taxonomy_json,
    genome_summary,
    genome_inventory,
    pipeline_summary,
    host_stats,
    versions_json,
    filtering_method,
    preprocessing_enabled,
    host_removal_enabled,
    host_ref,
    chopper_min_length,
    chopper_quality,
    filtlong_min_length,
    filtlong_keep_percent,
    workflow_png,
    output,
):
    tax = _read_json(taxonomy_json)
    host_stats_data = _read_json(host_stats)
    versions = _read_json(versions_json)

    # Build host removal Plotly chart data
    host_plot_data = None
    if host_stats_data:
        host_plot_data = {
            "labels": ["Host reads", "Non-host reads"],
            "values": [
                host_stats_data.get("host_reads", 0),
                host_stats_data.get("non_host_reads", 0),
            ],
            "percentages": [
                host_stats_data.get("host_percentage", 0),
                host_stats_data.get("non_host_percentage", 0),
            ],
        }

    context = {
        "sample": sample,
        "run_mode": run_mode,
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "kraken2_species": tax.get("kraken2", []),
        "centrifuger_species": tax.get("centrifuger", []),
        "genome_summary": _read_table(genome_summary),
        "genome_inventory": _read_table(genome_inventory),
        "pipeline_summary": _read_table(pipeline_summary),
        # Preprocessing info
        "preprocessing_enabled": preprocessing_enabled,
        "host_removal_enabled": host_removal_enabled,
        "host_ref": host_ref,
        "filtering_method": filtering_method,
        "chopper_min_length": chopper_min_length,
        "chopper_quality": chopper_quality,
        "filtlong_min_length": filtlong_min_length,
        "filtlong_keep_percent": filtlong_keep_percent,
        # Host removal stats
        "host_stats": host_stats_data,
        "host_plot_data": json.dumps(host_plot_data) if host_plot_data else None,
        # Versions
        "versions": versions,
        # Workflow diagram as base64
        "workflow_image": _encode_image_base64(workflow_png),
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
    ap.add_argument("--host-stats", default="",
                    help="JSON produced by parse_host_removal.py")
    ap.add_argument("--versions-json", default="",
                    help="versions.json produced by collect_versions.py")
    ap.add_argument("--filtering-method", default="chopper",
                    choices=["chopper", "filtlong"])
    ap.add_argument("--preprocessing-enabled", default="True")
    ap.add_argument("--host-removal-enabled", default="False")
    ap.add_argument("--host-ref", default="")
    ap.add_argument("--chopper-min-length", default="500")
    ap.add_argument("--chopper-quality", default="10")
    ap.add_argument("--filtlong-min-length", default="500")
    ap.add_argument("--filtlong-keep-percent", default="90")
    ap.add_argument("--workflow-png", default="",
                    help="Path to _workflow_.png for embedding as base64")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    render(
        sample=args.sample,
        run_mode=args.run_mode,
        taxonomy_json=args.taxonomy_json,
        genome_summary=args.genome_summary,
        genome_inventory=args.genome_inventory,
        pipeline_summary=args.pipeline_summary,
        host_stats=args.host_stats,
        versions_json=args.versions_json,
        filtering_method=args.filtering_method,
        preprocessing_enabled=args.preprocessing_enabled.lower() not in ("false", "0", ""),
        host_removal_enabled=args.host_removal_enabled.lower() not in ("false", "0", ""),
        host_ref=args.host_ref,
        chopper_min_length=args.chopper_min_length,
        chopper_quality=args.chopper_quality,
        filtlong_min_length=args.filtlong_min_length,
        filtlong_keep_percent=args.filtlong_keep_percent,
        workflow_png=args.workflow_png,
        output=args.output,
    )
    print(f"Report written to: {args.output}")


if __name__ == "__main__":
    main()
