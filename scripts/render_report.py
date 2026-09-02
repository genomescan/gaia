#!/usr/bin/env python3
"""
render_report.py – generate a single run-level standalone HTML report for a
gaia pipeline run across one or more samples.

Combines:
  - Pipeline metadata (mode, samples, date)
  - Run-level taxonomy relative-abundance charts
  - Run-level host removal and read QC (NanoPlot) summaries
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

# ---------------------------------------------------------------------------
# Interactive workflow diagram (replaces the static _workflow_.png in the
# report). This is a small, hand-maintained mapping of pipeline steps to the
# information shown in the on-hover tooltip in the report's "Pipeline
# overview" panel.
#
# To add/update a node:
#   1. Add/edit an entry in WORKFLOW_NODES below. "x"/"y" are percentages
#      (0-100) of the diagram viewBox and control where the dot is placed;
#      "category" must be a key in WORKFLOW_CATEGORIES (controls dot/line
#      colour and legend grouping); "description" is the tooltip text
#      (tool name + what it does in the pipeline).
#   2. If the new step should be visually connected to another step, add an
#      ("from_id", "to_id") tuple to WORKFLOW_EDGES.
# No other code changes are needed — render_report.py serialises these
# structures for report.html.j2 to draw as inline SVG.
# ---------------------------------------------------------------------------
WORKFLOW_CATEGORIES = {
    "preprocessing": {"label": "QC & Preprocessing", "color": "#6c8ebf"},
    "profiling": {"label": "Rapid Classification & Profiling", "color": "#9673a6"},
    "assembly": {"label": "Assembly & Binning", "color": "#82b366"},
    "prok": {"label": "Prok Bin Refinement & Classification", "color": "#b85450"},
    "euk": {"label": "Euk Bin Refinement & Classification", "color": "#d6b656"},
}

WORKFLOW_NODES = [
    {"id": "nanoplot_raw", "label": "NanoPlot + nanoQC (Raw)", "category": "preprocessing",
     "x": 8, "y": 20,
     "description": "Quality-checks the raw long reads (length, quality, yield) before any "
                    "preprocessing, giving a baseline to compare against after filtering."},
    {"id": "host_removal", "label": "minimap2 + samtools", "category": "preprocessing",
     "x": 24, "y": 20,
     "description": "Maps reads to the configured host reference genome and removes "
                    "host-derived alignments, keeping only non-host (microbial) reads."},
    {"id": "chopper", "label": "Chopper", "category": "preprocessing",
     "x": 40, "y": 20,
     "description": "Trims and filters reads by minimum length and quality score to remove "
                    "low-quality or overly short reads."},
    {"id": "filtlong", "label": "Filtlong", "category": "preprocessing",
     "x": 56, "y": 20,
     "description": "Performs additional read filtering, keeping the best-scoring reads by "
                    "length and quality (used as an alternative to Chopper)."},
    {"id": "nanoplot_filtered", "label": "NanoPlot + nanoQC (Filtered)", "category": "preprocessing",
     "x": 72, "y": 20,
     "description": "Re-runs quality control on the filtered reads to confirm the effect of "
                    "preprocessing before downstream analysis."},

    {"id": "kraken2", "label": "Kraken2", "category": "profiling",
     "x": 40, "y": 40,
     "description": "Fast k-mer based taxonomic classification run directly on reads, for "
                    "rapid compositional profiling."},
    {"id": "centrifuger", "label": "Centrifuger", "category": "profiling",
     "x": 56, "y": 40,
     "description": "Classifies reads against a compressed reference index to estimate "
                    "microbial community composition and abundance."},

    {"id": "metaflye", "label": "metaFlye", "category": "assembly",
     "x": 24, "y": 58,
     "description": "Assembles the filtered long reads into metagenomic contigs."},
    {"id": "metaquast", "label": "MetaQUAST", "category": "assembly",
     "x": 40, "y": 58,
     "description": "Evaluates assembly quality, reporting metrics such as N50, contig count "
                    "and total assembly length."},
    {"id": "coverage", "label": "Coverage / depth calculation", "category": "assembly",
     "x": 56, "y": 58,
     "description": "Maps reads back to the assembly (minimap2/samtools) and summarises "
                    "per-contig depth (jgi_summarize_bam_contig_depths) as input for binning."},
    {"id": "binners", "label": "SemiBin2 / MetaBAT2 / COMEBin", "category": "assembly",
     "x": 72, "y": 58,
     "description": "Groups assembled contigs into candidate genome bins using sequence "
                    "composition and coverage/depth signals."},

    {"id": "dastool", "label": "DAS Tool", "category": "prok",
     "x": 24, "y": 76,
     "description": "Integrates candidate bins from the multiple binners into a single "
                    "non-redundant, optimised set of prokaryotic bins."},
    {"id": "checkm2", "label": "CheckM2", "category": "prok",
     "x": 40, "y": 76,
     "description": "Estimates completeness and contamination of the refined prokaryotic "
                    "bins to assess genome quality."},
    {"id": "gtdbtk", "label": "GTDB-Tk", "category": "prok",
     "x": 56, "y": 76,
     "description": "Assigns taxonomy to prokaryotic bins using the GTDB reference database "
                    "and marker-gene phylogenetics."},

    {"id": "acr", "label": "ACR", "category": "euk",
     "x": 24, "y": 92,
     "description": "Identifies and refines candidate eukaryotic bins from the assembly for "
                    "downstream eukaryotic-specific processing."},
    {"id": "drep_eukcc", "label": "dRep / EukCC", "category": "euk",
     "x": 40, "y": 92,
     "description": "Dereplicates redundant eukaryotic bins (dRep) and assesses their "
                    "completeness/contamination (EukCC)."},
    {"id": "bat", "label": "BAT", "category": "euk",
     "x": 56, "y": 92,
     "description": "Assigns taxonomy to the refined eukaryotic bins/contigs using the Bin "
                    "Annotation Tool."},
]

# Visual connections drawn between node centres (in the order given).
WORKFLOW_EDGES = [
    ("nanoplot_raw", "host_removal"),
    ("host_removal", "chopper"),
    ("chopper", "filtlong"),
    ("filtlong", "nanoplot_filtered"),
    ("chopper", "kraken2"),
    ("kraken2", "centrifuger"),
    ("host_removal", "metaflye"),
    ("metaflye", "metaquast"),
    ("metaflye", "coverage"),
    ("coverage", "binners"),
    ("binners", "dastool"),
    ("dastool", "checkm2"),
    ("checkm2", "gtdbtk"),
    ("binners", "acr"),
    ("acr", "drep_eukcc"),
    ("drep_eukcc", "bat"),
]


def _build_workflow_diagram():
    """Assemble the interactive workflow diagram context: nodes (with resolved
    colours), edges (as endpoint coordinate pairs) and the category legend."""
    nodes_by_id = {n["id"]: n for n in WORKFLOW_NODES}
    nodes = [
        {**n, "color": WORKFLOW_CATEGORIES[n["category"]]["color"]}
        for n in WORKFLOW_NODES
    ]
    edges = []
    for src, dst in WORKFLOW_EDGES:
        a, b = nodes_by_id[src], nodes_by_id[dst]
        edges.append({
            "x1": a["x"], "y1": a["y"], "x2": b["x"], "y2": b["y"],
            "color": WORKFLOW_CATEGORIES[a["category"]]["color"],
        })
    legend = [
        {"label": cfg["label"], "color": cfg["color"]}
        for cfg in WORKFLOW_CATEGORIES.values()
    ]
    return {"nodes": nodes, "edges": edges, "legend": legend}


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


def _to_float(value):
    """Best-effort conversion to float; returns None for missing/invalid values
    (empty strings, "NA"/"N/A", None, NaN) instead of raising."""
    if value is None:
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    if f != f:  # NaN check without importing math for a single use
        return None
    return f


def _agg_stats(values):
    """Return {count, min, max, median} for a list of numeric values, filtering
    out anything non-numeric/NaN first. Returns None if nothing usable remains."""
    vals = [v for v in (_to_float(v) for v in values) if v is not None]
    if not vals:
        return None
    return {
        "count": len(vals),
        "min": min(vals),
        "max": max(vals),
        "median": statistics.median(vals),
    }


def _fmt_number(value, kind="float1"):
    """Human-friendly formatting for aggregate values. `kind` controls units:
    "int" (thousands-separated integer), "pct" (value already 0-100, shown with
    a % suffix), "pct_frac" (value is a 0-1 fraction, multiplied by 100), or
    "float1" (plain 1-decimal number, the default)."""
    if value is None:
        return "—"
    try:
        if kind == "int":
            return f"{int(round(value)):,}"
        if kind == "pct":
            return f"{value:.1f}%"
        if kind == "pct_frac":
            return f"{value * 100:.1f}%"
        return f"{value:.1f}"
    except (TypeError, ValueError):
        return "—"


def _fmt_range(stat, kind="float1"):
    """Render an aggregate stat dict (see _agg_stats) as a compact human
    readable range, e.g. "64.0 – 95.9% (median 71.3%, n=8)". Returns an em-dash
    when no data is available."""
    if not stat:
        return "—"
    lo = _fmt_number(stat["min"], kind)
    hi = _fmt_number(stat["max"], kind)
    med = _fmt_number(stat["median"], kind)
    if lo == hi:
        return f"{lo} (n={stat['count']})"
    return f"{lo} – {hi} (median {med}, n={stat['count']})"


def _round_species(rows):
    out = []
    for row in rows:
        row = dict(row)
        if "percent" in row:
            row["percent"] = _round1(row["percent"])
        out.append(row)
    return out


# Nanopore-oriented quality scale for the NanoPlot quality bars: modern ONT
# basecalling tops out around Q20-Q30, so bars are drawn on a fixed Q0-Q30
# axis with a Q20 ("high accuracy") reference marker instead of being scaled
# to the observed maximum.
_NANOPLOT_QUALITY_SCALE_MAX = 30.0
_NANOPLOT_QUALITY_TARGET = 20.0

# (metric key, label, unit suffix, _fmt_number kind, fixed scale maximum)
_NANOPLOT_BAR_METRICS = [
    ("median_read_length", "Median read length", " bp", "int", None),
    ("read_length_n50", "Read length N50", " bp", "int", None),
    ("median_read_quality", "Median read quality", "", "float1",
     _NANOPLOT_QUALITY_SCALE_MAX),
    ("number_of_reads", "Reads", "", "int", None),
    ("total_bases", "Total bases", " bp", "int", None),
]


def _nanoplot_bar(value, scale_max, unit, kind):
    """Render one thin-bar cell: percentage width plus a formatted label."""
    num = _to_float(value)
    if num is None or not scale_max:
        return {"value": None, "label": "—", "width": 0}
    width = max(0.0, min(100.0, num / scale_max * 100.0))
    prefix = "Q" if kind == "float1" and scale_max == _NANOPLOT_QUALITY_SCALE_MAX else ""
    return {
        "value": num,
        "label": f"{prefix}{_fmt_number(num, kind)}{unit}",
        "width": round(width, 1),
    }


def _build_nanoplot_bars(samples_data):
    """Build per-metric thin horizontal bar data (raw vs filtered per sample)."""
    groups = []
    for key, label, unit, kind, fixed_max in _NANOPLOT_BAR_METRICS:
        values = []
        for sd in samples_data:
            values.append(_to_float(sd["nanoplot_raw"].get(key)))
            values.append(_to_float(sd["nanoplot_filtered"].get(key)))
        observed = [v for v in values if v is not None]
        scale_max = fixed_max or (max(observed) if observed else None)
        rows = []
        for sd in samples_data:
            rows.append({
                "sample": sd["sample"],
                "raw": _nanoplot_bar(sd["nanoplot_raw"].get(key), scale_max, unit, kind),
                "filtered": _nanoplot_bar(sd["nanoplot_filtered"].get(key), scale_max, unit, kind),
            })
        groups.append({
            "key": key,
            "label": label,
            "has_data": bool(observed),
            "scale_max_label": (
                (f"Q{int(scale_max)}" if fixed_max else _fmt_number(scale_max, kind) + unit)
                if scale_max else "—"
            ),
            "scale_min_label": "Q0" if fixed_max else "0",
            "marker": (
                {"percent": round(_NANOPLOT_QUALITY_TARGET / scale_max * 100.0, 1),
                 "label": f"Q{int(_NANOPLOT_QUALITY_TARGET)}"}
                if fixed_max and scale_max else None
            ),
            "rows": rows,
        })
    return groups


def _build_nanoplot_overview(samples_data):
    """Build cross-sample nanoplot data for the overview table."""
    samples = [sd["sample"] for sd in samples_data]
    keys = [
        "median_read_length",
        "read_length_n50",
        "median_read_quality",
        "number_of_reads",
        "total_bases",
    ]
    raw = {key: [] for key in keys}
    filtered = {key: [] for key in keys}
    for sd in samples_data:
        for key in keys:
            raw[key].append(sd["nanoplot_raw"].get(key))
            filtered[key].append(sd["nanoplot_filtered"].get(key))
    return {
        "samples": samples,
        "raw": raw,
        "filtered": filtered,
        "bars": _build_nanoplot_bars(samples_data),
    }


def _build_taxonomy_overview(samples_data):
    """Build cross-sample taxonomy data for 100% stacked bar charts."""
    def _process(tool_key):
        species_reads = {}
        for sd in samples_data:
            for row in sd.get(tool_key + "_species", []):
                species_reads[row["name"]] = species_reads.get(row["name"], 0) + row["reads"]

        top_species = sorted(species_reads, key=lambda s: species_reads[s], reverse=True)[:5]
        samples = [sd["sample"] for sd in samples_data]
        has_data = any(sd.get(tool_key + "_species") for sd in samples_data)

        data = {}
        for sd in samples_data:
            smp = sd["sample"]
            rows = sd.get(tool_key + "_species", [])
            rows_by_name = {r["name"]: r for r in rows}
            total_reads = sum(r["reads"] for r in rows)
            data[smp] = {}
            assigned_pct = 0.0
            assigned_reads = 0
            for sp in top_species:
                if sp in rows_by_name:
                    r = rows_by_name[sp]
                    data[smp][sp] = {"reads": r["reads"], "percent": r["percent"]}
                    assigned_pct += r["percent"]
                    assigned_reads += r["reads"]
                else:
                    data[smp][sp] = {"reads": 0, "percent": 0.0}

            other_reads = max(0, total_reads - assigned_reads)
            other_pct = max(0.0, round(100.0 - assigned_pct, 1))
            if other_reads > 0 or other_pct > 0:
                data[smp]["Other"] = {"reads": other_reads, "percent": other_pct}

        has_other = any(
            "Other" in data[smp] and data[smp]["Other"]["reads"] > 0
            for smp in samples
        )
        species_list = top_species + (["Other"] if has_other else [])

        return {"has_data": has_data, "samples": samples, "species": species_list, "data": data}

    return {
        "kraken2": _process("kraken2"),
        "centrifuger": _process("centrifuger"),
    }


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


def _build_gtdbtk(bac120_paths, ar53_paths):
    """Merge bac120 + ar53 GTDB-Tk summaries across all samples into one list."""
    rows = []
    seen = set()
    for path in list(bac120_paths) + list(ar53_paths):
        for row in _read_table(path):
            gid = row.get("user_genome", "").strip()
            if not gid or gid in seen:
                continue
            seen.add(gid)
            classification = row.get("classification", "")
            parts = [part.strip() for part in classification.split(";") if part.strip()]
            species = parts[-1] if parts else classification
            if "__" in species:
                species = species.split("__", 1)[1]

            def _fmt(value, decimals=1):
                try:
                    return round(float(value), decimals)
                except (ValueError, TypeError):
                    return "N/A"

            ani = row.get("closest_genome_ani") or row.get("closest_placement_ani") or row.get("fastani_ani", "")
            af = row.get("closest_genome_af") or row.get("closest_placement_af") or row.get("fastani_af", "")
            red = row.get("red_value", "") or "N/A"
            warning = row.get("warnings", "") or "N/A"
            note = row.get("note", "") or ""
            rows.append({
                "user_genome": gid,
                "classification": species,
                "ani": _fmt(ani, 1),
                "af": _fmt(af, 1),
                "red": "N/A" if red in ("", "N/A", "none") else red,
                "warnings": "N/A" if warning in ("", "none") else warning,
                "notes": note,
            })
    return rows


# CheckM2 quality_report.tsv columns we aggregate for the run-level CheckM2
# summary, and how each should be formatted (see _fmt_number).
_CHECKM2_FIELDS = [
    ("completeness", "Completeness", "pct"),
    ("contamination", "Contamination", "pct"),
    ("coding_density", "Coding_Density", "pct_frac"),
    ("average_gene_length", "Average_Gene_Length", "float1"),
    ("genome_size", "Genome_Size", "int"),
    ("gc_content", "GC_Content", "pct_frac"),
    ("total_coding_sequences", "Total_Coding_Sequences", "int"),
]


def _checkm2_sample_from_path(path):
    """CheckM2 reports live at .../CheckM2/<sample>/quality_report.tsv, so the
    parent directory name identifies the sample."""
    parent = os.path.basename(os.path.dirname(os.path.abspath(path or "")))
    return parent or "—"


def _build_checkm2_overview(checkm2_paths):
    """Aggregate CheckM2 quality_report.tsv files (one per binned sample) into
    a single run-level summary plus a per-sample/per-genome table. Returns
    {"available": False} gracefully when no CheckM2 output could be read, so
    the template can hide the section."""
    rows = []
    genome_rows = []
    for path in checkm2_paths:
        sample = _checkm2_sample_from_path(path)
        for row in _read_table(path):
            rows.append(row)
            genome_row = {
                "sample": sample,
                "user_genome": (row.get("Name") or "").strip() or "—",
            }
            for key, column, kind in _CHECKM2_FIELDS:
                genome_row[key] = _fmt_number(_to_float(row.get(column)), kind)
            genome_rows.append(genome_row)

    if not rows:
        return {"available": False, "genome_count": 0, "metrics": {}, "rows": []}

    metrics = {}
    for key, column, kind in _CHECKM2_FIELDS:
        stat = _agg_stats(row.get(column) for row in rows)
        metrics[key] = {"stat": stat, "range": _fmt_range(stat, kind)}

    genome_rows.sort(key=lambda r: (r["sample"], r["user_genome"]))

    return {
        "available": True,
        "genome_count": len(rows),
        "sample_count": len({r["sample"] for r in genome_rows}),
        "metrics": metrics,
        "rows": genome_rows,
    }


def _build_top_taxa(samples_data, top_n=5):
    """Aggregate top taxa across all samples (reads summed per species) using
    whichever taxonomy tool has data (kraken2 preferred, else centrifuger)."""
    for tool_key in ("kraken2", "centrifuger"):
        species_reads = {}
        total_reads = 0
        for sd in samples_data:
            for row in sd.get(tool_key + "_species", []):
                reads = row.get("reads") or 0
                species_reads[row["name"]] = species_reads.get(row["name"], 0) + reads
                total_reads += reads
        if species_reads and total_reads:
            top = sorted(species_reads.items(), key=lambda kv: kv[1], reverse=True)[:top_n]
            return {
                "available": True,
                "tool": tool_key,
                "taxa": [
                    {"name": name, "reads": reads, "percent": _round1(reads / total_reads * 100)}
                    for name, reads in top
                ],
            }
    return {"available": False, "tool": None, "taxa": []}


def _build_general_stats(samples_data, checkm2_overview, host_removal_enabled):
    """Build the compact run-level "General Statistics" context: per-sample
    host-removal rings, run-level N50/CheckM2 ranges and top taxa. Every
    sub-part degrades gracefully (has_ring=False, stat=None, etc.) when its
    source data is missing, so the template never has to deal with broken or
    partial values."""
    n50_values = []
    for sd in samples_data:
        n50 = _to_float(sd["nanoplot_filtered"].get("read_length_n50"))
        if n50 is None:
            n50 = _to_float(sd["nanoplot_raw"].get("read_length_n50"))
        if n50 is not None:
            n50_values.append(n50)
    n50_stat = _agg_stats(n50_values)

    rings = []
    any_host_data = False
    for sd in samples_data:
        host = sd.get("host_stats") or {}
        host_pct = host.get("host_percentage")
        non_host_pct = host.get("non_host_percentage")
        has_ring = host_pct is not None and non_host_pct is not None
        if has_ring:
            any_host_data = True
        rings.append({
            "sample": sd["sample"],
            "host_percentage": host_pct,
            "non_host_percentage": non_host_pct,
            "has_ring": has_ring,
        })

    host_removal_available = bool(host_removal_enabled and any_host_data)
    top_taxa = _build_top_taxa(samples_data)

    completeness = checkm2_overview.get("metrics", {}).get("completeness", {}).get("stat")
    contamination = checkm2_overview.get("metrics", {}).get("contamination", {}).get("stat")

    return {
        "has_any_data": bool(n50_stat or checkm2_overview.get("available") or
                              host_removal_available or top_taxa["available"]),
        "n50": n50_stat,
        "n50_range": _fmt_range(n50_stat, "int"),
        "rings": rings,
        "host_removal_available": host_removal_available,
        "checkm2_available": checkm2_overview.get("available", False),
        "completeness_range": _fmt_range(completeness, "pct"),
        "contamination_range": _fmt_range(contamination, "pct"),
        "top_taxa": top_taxa,
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
    gtdbtk_bac120s,
    gtdbtk_ar53s,
    # genome_inventories / pipeline_summaries are accepted for CLI/Snakemake
    # rule compatibility, but are intentionally not rendered per-sample
    # anymore: genome_summaries already feeds the collapsed binning
    # overview, and per-sample inventory/synthesis tables were removed to
    # keep multi-sample reports concise.
    genome_inventories,
    pipeline_summaries,
    host_stats_jsons,
    checkm2_reports,
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

    # Build per-sample data used to assemble the run-level overview sections.
    samples_data = []
    for i, sample in enumerate(samples):
        tax = _read_json(taxonomy_jsons[i] if i < len(taxonomy_jsons) else "")
        host_stats_data = _read_json(host_stats_jsons[i] if i < len(host_stats_jsons) else "")
        nanoplot_data = _read_json(nanoplot_jsons[i] if i < len(nanoplot_jsons) else "")

        if host_stats_data:
            host_stats_data["host_percentage"] = _round1(host_stats_data.get("host_percentage", 0))
            host_stats_data["non_host_percentage"] = _round1(host_stats_data.get("non_host_percentage", 0))

        samples_data.append({
            "sample": sample,
            "kraken2_species": _round_species(tax.get("kraken2", [])),
            "centrifuger_species": _round_species(tax.get("centrifuger", [])),
            "host_stats": host_stats_data,
            "nanoplot_raw": nanoplot_data.get("raw", {}),
            "nanoplot_filtered": nanoplot_data.get("filtered", {}),
        })

    checkm2_overview = _build_checkm2_overview(checkm2_reports)

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
        "gtdbtk_rows": _build_gtdbtk(gtdbtk_bac120s, gtdbtk_ar53s),
        "nanoplot_overview": _build_nanoplot_overview(samples_data),
        "taxonomy_overview": _build_taxonomy_overview(samples_data),
        "checkm2_overview": checkm2_overview,
        "general_stats": _build_general_stats(samples_data, checkm2_overview, host_removal_enabled),
        # Interactive workflow diagram (SVG nodes/edges/legend). The static
        # _workflow_.png is still embedded as a base64 fallback for contexts
        # where the interactive SVG cannot be shown (e.g. printing).
        "workflow_diagram": _build_workflow_diagram(),
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
    ap.add_argument("--gtdbtk-bac120s", default="",
                    help="Space-separated GTDB-Tk bac120 TSV files (one per sample)")
    ap.add_argument("--gtdbtk-ar53s", default="",
                    help="Space-separated GTDB-Tk ar53 TSV files (one per sample)")
    ap.add_argument("--genome-inventories", default="")
    ap.add_argument("--pipeline-summaries", default="")
    ap.add_argument("--host-stats-jsons", default="",
                    help="Space-separated host-removal stats JSON files (one per sample)")
    ap.add_argument("--checkm2-reports", default="",
                    help="Space-separated CheckM2 quality_report.tsv files (one per binned sample)")
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
        gtdbtk_bac120s=_split_paths(args.gtdbtk_bac120s),
        gtdbtk_ar53s=_split_paths(args.gtdbtk_ar53s),
        genome_inventories=_split_paths(args.genome_inventories),
        pipeline_summaries=_split_paths(args.pipeline_summaries),
        host_stats_jsons=_split_paths(args.host_stats_jsons),
        checkm2_reports=_split_paths(args.checkm2_reports),
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
