#!/usr/bin/env python3
"""
render_report.py – generate a single run-level standalone HTML report for a
gaia pipeline run across one or more samples.

Combines:
  - Pipeline metadata (mode, samples, date)
  - Per-sample taxonomy top-10 species (bar chart + HTML table via Plotly)
  - Per-sample host removal and read QC (NanoPlot) summaries
  - A run-level assembly overview (N50 and other metrics, collapsed across
    samples rather than shown per sample)
  - A run-level binning overview (genome recovery counts, collapsed across
    samples; samples without binning results are called out explicitly)
  - Tool versions (from Reports/versions.json)
  - Embedded workflow diagram (_workflow_.png as base64)
  - Embedded Plotly runtime JS (injected from file path via argparse)

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
import statistics
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape
from markupsafe import Markup

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


def _read_text(path):
    if not path:
        raise ValueError("--plotly-js is required")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Plotly JS file not found: {path}")
    with open(path, encoding="utf-8") as fh:
        return fh.read()


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


def _split_paths(arg):
    """Split a space-separated (or empty) argument into a list of non-empty paths."""
    if not arg or not arg.strip():
        return []
    return [p for p in arg.split() if p.strip()]


def _round1(value):
    """Round a numeric value to 1 decimal place; pass through non-numeric values."""
    try:
        return round(float(value), 1)
    except (TypeError, ValueError):
        return value


def _round_species(rows):
    out = []
    for row in rows:
        row = dict(row)
        if "percent" in row:
            row["percent"] = _round1(row["percent"])
        out.append(row)
    return out


def _index_tables_by_sample(paths):
    """Read a list of TSV files and index their rows by the "sample" column
    found in the row itself, rather than relying on list position. This is
    needed because some per-sample TSVs (e.g. genome/pipeline summaries) are
    only produced for a subset of samples (those that passed the binning
    gate), so their paths list does not line up positionally with the full
    sample list.
    """
    by_sample = {}
    for path in paths:
        rows = _read_table(path)
        for row in rows:
            sample = row.get("sample")
            if sample:
                by_sample.setdefault(sample, []).append(row)
    return by_sample


def _index_json_by_sample(paths):
    by_sample = {}
    for path in paths:
        data = _read_json(path)
        sample = data.get("sample")
        if sample:
            by_sample[sample] = data
    return by_sample


def _build_assembly_overview(samples, assembly_stats_jsons, run_assembly):
    """Aggregate per-sample assembly QC (N50, contig count, ...) into a
    single run-level overview rather than per-sample sections."""
    stats_by_sample = _index_json_by_sample(assembly_stats_jsons)

    rows = []
    n50_values = []
    total_contigs = 0
    total_length = 0
    for sample in samples:
        st = stats_by_sample.get(sample, {})
        n50 = st.get("n50")
        contig_count = st.get("contig_count")
        length = st.get("total_length")
        rows.append({
            "sample": sample,
            "contig_count": contig_count if contig_count is not None else "—",
            "total_length": length if length is not None else "—",
            "n50": n50 if n50 is not None else "—",
            "largest_contig": st.get("largest_contig", "—"),
            "source": st.get("source", ""),
        })
        if isinstance(n50, (int, float)):
            n50_values.append(n50)
        if isinstance(contig_count, (int, float)):
            total_contigs += contig_count
        if isinstance(length, (int, float)):
            total_length += length

    return {
        "enabled": run_assembly,
        "rows": rows,
        "samples_total": len(samples),
        "samples_assembled": sum(1 for r in rows if r["n50"] != "—"),
        "total_contigs": total_contigs,
        "total_length": total_length,
        "median_n50": (
            round(statistics.median(n50_values)) if n50_values else None
        ),
        "min_n50": min(n50_values) if n50_values else None,
        "max_n50": max(n50_values) if n50_values else None,
    }


def _build_binning_overview(samples, genome_summaries, binning_enabled, run_assembly):
    """Aggregate per-sample binning/genome-recovery results into a single
    run-level overview, and flag samples for which binning was not done."""
    summary_by_sample = _index_tables_by_sample(genome_summaries)

    rows = []
    binned_count = 0
    not_binned_samples = []
    totals = {
        "total_final_genomes": 0,
        "total_bins_prok": 0,
        "total_bins_euk": 0,
        "high_quality_prok": 0,
        "medium_quality_prok": 0,
        "euk_bins_with_taxonomy": 0,
    }
    for sample in samples:
        sample_rows = summary_by_sample.get(sample)
        summary = sample_rows[0] if sample_rows else None
        if summary:
            binned_count += 1
            row = {"sample": sample, "binned": True}
            for key in totals:
                try:
                    val = int(summary.get(key, 0) or 0)
                except (TypeError, ValueError):
                    val = 0
                row[key] = val
                totals[key] += val
        else:
            not_binned_samples.append(sample)
            row = {"sample": sample, "binned": False}
            for key in totals:
                row[key] = "—"
        rows.append(row)

    return {
        "enabled": binning_enabled,
        "applicable": run_assembly,
        "rows": rows,
        "samples_total": len(samples),
        "samples_binned": binned_count,
        "not_binned_samples": not_binned_samples,
        **totals,
    }


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

def render(
    samples,
    run_mode,
    binning_enabled,
    taxonomy_jsons,
    nanoplot_jsons,
    assembly_stats_jsons,
    genome_summaries,
    # genome_inventories / pipeline_summaries are accepted for CLI/Snakemake
    # rule compatibility, but are intentionally not rendered per-sample
    # anymore: genome_summaries already feeds the collapsed binning
    # overview, and per-sample inventory/synthesis tables were removed to
    # keep multi-sample reports concise.
    genome_inventories,
    pipeline_summaries,
    host_stats_jsons,
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
    plotly_js,
    output,
):
    versions = _read_json(versions_json)
    run_assembly = run_mode in ("assembly_binning_only", "both")

    # Build per-sample data list (taxonomy, host removal, read QC only —
    # assembly and binning results are collapsed into run-level overviews
    # below instead of per-sample sections, to keep the report readable for
    # runs with many samples).
    samples_data = []
    for i, sample in enumerate(samples):
        tax = _read_json(taxonomy_jsons[i] if i < len(taxonomy_jsons) else "")
        host_stats_data = _read_json(host_stats_jsons[i] if i < len(host_stats_jsons) else "")
        nanoplot_data = _read_json(nanoplot_jsons[i] if i < len(nanoplot_jsons) else "")

        if host_stats_data:
            host_stats_data["host_percentage"] = _round1(host_stats_data.get("host_percentage", 0))
            host_stats_data["non_host_percentage"] = _round1(host_stats_data.get("non_host_percentage", 0))

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

        samples_data.append({
            "sample": sample,
            "kraken2_species": _round_species(tax.get("kraken2", [])),
            "centrifuger_species": _round_species(tax.get("centrifuger", [])),
            "host_stats": host_stats_data,
            "host_plot_data": host_plot_data,
            "nanoplot_raw": nanoplot_data.get("raw", {}),
            "nanoplot_filtered": nanoplot_data.get("filtered", {}),
        })

    context = {
        "samples": samples,
        "samples_data": samples_data,
        "run_mode": run_mode,
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        # Preprocessing info (run-level, same for all samples)
        "preprocessing_enabled": preprocessing_enabled,
        "host_removal_enabled": host_removal_enabled,
        "host_ref": host_ref,
        "filtering_method": filtering_method,
        "chopper_min_length": chopper_min_length,
        "chopper_quality": chopper_quality,
        "filtlong_min_length": filtlong_min_length,
        "filtlong_keep_percent": filtlong_keep_percent,
        # Versions
        "versions": versions,
        # Run-level overviews (collapsed across samples)
        "assembly_overview": _build_assembly_overview(samples, assembly_stats_jsons, run_assembly),
        "binning_overview": _build_binning_overview(samples, genome_summaries, binning_enabled, run_assembly),
        # Workflow diagram as base64
        "workflow_image": _encode_image_base64(workflow_png),
        # Inline Plotly runtime JS
        "plotly_js": Markup(f"<script>{plotly_js}</script>"),
    }

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATES_DIR)),
        autoescape=select_autoescape(["html"]),
    )
    env.filters["tojson"] = lambda value, **kwargs: Markup(json.dumps(value))
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
    ap.add_argument("--samples", required=True,
                    help="Comma-separated sample names")
    ap.add_argument("--run-mode", required=True,
                    choices=["profiling_only", "assembly_binning_only", "both"])
    ap.add_argument("--binning-enabled", default="False",
                    help="Whether binning was enabled for this run")
    ap.add_argument("--taxonomy-jsons", default="",
                    help="Space-separated JSON files produced by parse_taxonomy.py (one per sample, same order as --samples)")
    ap.add_argument("--nanoplot-jsons", default="",
                    help="Space-separated JSON files produced by parse_nanoplot_metrics.py (one per sample, same order as --samples)")
    ap.add_argument("--assembly-stats-jsons", default="",
                    help="Space-separated JSON files produced by parse_assembly_stats.py (one per sample)")
    ap.add_argument("--genome-summaries", default="")
    ap.add_argument("--genome-inventories", default="")
    ap.add_argument("--pipeline-summaries", default="")
    ap.add_argument("--host-stats-jsons", default="",
                    help="Space-separated host-removal stats JSON files (one per sample)")
    ap.add_argument("--versions-json", default=os.path.join("Reports", "versions.json"),
                    help="versions.json written by the run wrapper")
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
                    help="Path to _workflow_.png for embedding as base64 (should be pipeline repo root)")
    ap.add_argument("--plotly-js", required=True,
                    help="Path to Plotly JS source file to inject inline (e.g. templates/report/plotly-v1.58.5.js)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    samples = [s.strip() for s in args.samples.split(",") if s.strip()]

    render(
        samples=samples,
        run_mode=args.run_mode,
        binning_enabled=args.binning_enabled.lower() not in ("false", "0", ""),
        taxonomy_jsons=_split_paths(args.taxonomy_jsons),
        nanoplot_jsons=_split_paths(args.nanoplot_jsons),
        assembly_stats_jsons=_split_paths(args.assembly_stats_jsons),
        genome_summaries=_split_paths(args.genome_summaries),
        genome_inventories=_split_paths(args.genome_inventories),
        pipeline_summaries=_split_paths(args.pipeline_summaries),
        host_stats_jsons=_split_paths(args.host_stats_jsons),
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
        plotly_js=_read_text(args.plotly_js),
        output=args.output,
    )
    print(f"Report written to: {args.output}")


if __name__ == "__main__":
    main()
