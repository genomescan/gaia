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
"""
import argparse
import csv
import json
import os
from datetime import datetime

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _exists_nonempty(path):
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def _read_table(path):
    """Return list of dicts from a tab/comma-delimited file."""
    if not _exists_nonempty(path):
        return []
    sample_bytes = open(path, "rb").read(4096)
    delim = "\t"
    try:
        import csv as _csv
        delim = _csv.Sniffer().sniff(sample_bytes.decode("utf-8", errors="replace"),
                                     delimiters="\t,;").delimiter
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
# Chart.js horizontal bar chart snippet
# ---------------------------------------------------------------------------

_CHART_TEMPLATE = """
<canvas id="{chart_id}" style="max-height:320px;"></canvas>
<script>
(function(){{
  var ctx = document.getElementById('{chart_id}').getContext('2d');
  new Chart(ctx, {{
    type: 'bar',
    data: {{
      labels: {labels_json},
      datasets: [{{
        label: 'Reads',
        data: {values_json},
        backgroundColor: '{color}',
        borderColor: '{color}',
        borderWidth: 1
      }}]
    }},
    options: {{
      indexAxis: 'y',
      responsive: true,
      plugins: {{
        legend: {{ display: false }},
        tooltip: {{
          callbacks: {{
            label: function(context) {{
              var pct = {percents_json}[context.dataIndex];
              return 'Reads: ' + context.parsed.x.toLocaleString() + ' (' + pct + '%)';
            }}
          }}
        }}
      }},
      scales: {{
        x: {{ beginAtZero: true }}
      }}
    }}
  }});
}})();
</script>
"""


def _taxonomy_chart(species_list, chart_id, color):
    if not species_list:
        return "<p class='text-muted'>No species-level reads detected.</p>"
    labels = json.dumps([s["name"] for s in species_list])
    values = json.dumps([s["reads"] for s in species_list])
    percents = json.dumps([s["percent"] for s in species_list])
    return _CHART_TEMPLATE.format(
        chart_id=chart_id,
        labels_json=labels,
        values_json=values,
        percents_json=percents,
        color=color,
    )


# ---------------------------------------------------------------------------
# HTML template helpers
# ---------------------------------------------------------------------------

def _badge(text, cls="secondary"):
    return f'<span class="badge bg-{cls}">{text}</span>'


def _table_html(rows, columns=None):
    """Render a list-of-dicts as a Bootstrap table."""
    if not rows:
        return "<p class='text-muted'>No data available.</p>"
    cols = columns or list(rows[0].keys())
    header = "".join(f"<th>{c}</th>" for c in cols)
    body_rows = []
    for row in rows:
        cells = "".join(f"<td>{row.get(c, '')}</td>" for c in cols)
        body_rows.append(f"<tr>{cells}</tr>")
    body = "\n".join(body_rows)
    return f"""
<div class="table-responsive">
  <table class="table table-sm table-striped table-hover align-middle">
    <thead class="table-dark"><tr>{header}</tr></thead>
    <tbody>{body}</tbody>
  </table>
</div>"""


def _section(title, content, icon="📋"):
    return f"""
<div class="card mb-4 shadow-sm">
  <div class="card-header fw-semibold">{icon} {title}</div>
  <div class="card-body">{content}</div>
</div>
"""


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------

def _section_pipeline_info(sample, run_mode, generated_at):
    mode_badge = {
        "profiling_only":        ("Profiling only",        "info"),
        "assembly_binning_only": ("Assembly &amp; Binning", "warning"),
        "both":                  ("Both",                  "success"),
    }.get(run_mode, (run_mode, "secondary"))

    content = f"""
<dl class="row mb-0">
  <dt class="col-sm-3">Sample</dt>
  <dd class="col-sm-9"><code>{sample}</code></dd>
  <dt class="col-sm-3">Pipeline mode</dt>
  <dd class="col-sm-9">{_badge(mode_badge[0], mode_badge[1])}</dd>
  <dt class="col-sm-3">Generated</dt>
  <dd class="col-sm-9">{generated_at}</dd>
</dl>"""
    return _section("Pipeline information", content, "ℹ️")


def _section_taxonomy(taxonomy_json_path):
    data = _read_json(taxonomy_json_path)
    if not data:
        return ""

    parts = []

    kraken2 = data.get("kraken2", [])
    centrifuger = data.get("centrifuger", [])

    if not kraken2 and not centrifuger:
        return _section("Taxonomy – top species", "<p class='text-muted'>No taxonomy data available.</p>", "🔬")

    if kraken2:
        chart = _taxonomy_chart(kraken2, "chart_kraken2", "rgba(54, 162, 235, 0.75)")
        top_pct = sum(s["percent"] for s in kraken2)
        parts.append(f"""
<h6 class="mt-2">Kraken2 – top {len(kraken2)} species
  <small class="text-muted ms-2">(top-{len(kraken2)} account for {top_pct:.1f}% of all reads)</small>
</h6>
{chart}
{_table_html(kraken2, ["name", "reads", "percent"])}
""")

    if centrifuger:
        chart = _taxonomy_chart(centrifuger, "chart_centrifuger", "rgba(255, 159, 64, 0.75)")
        top_pct = sum(s["percent"] for s in centrifuger)
        parts.append(f"""
<h6 class="mt-4">Centrifuger – top {len(centrifuger)} species
  <small class="text-muted ms-2">(top-{len(centrifuger)} account for {top_pct:.4f} abundance)</small>
</h6>
{chart}
{_table_html(centrifuger, ["name", "reads", "percent"])}
""")

    return _section("Taxonomy – top species per classifier", "".join(parts), "🔬")


def _section_genome_summary(genome_summary_path):
    rows = _read_table(genome_summary_path)
    if not rows:
        return ""
    content = _table_html(rows)
    return _section("Genome summary", content, "🧬")


def _section_genome_inventory(genome_inventory_path):
    rows = _read_table(genome_inventory_path)
    if not rows:
        return ""
    content = _table_html(rows)
    return _section("Genome inventory", content, "📦")


def _section_pipeline_summary(pipeline_summary_path):
    rows = _read_table(pipeline_summary_path)
    if not rows:
        return ""
    content = _table_html(rows)
    return _section("Combined pipeline summary", content, "📊")


# ---------------------------------------------------------------------------
# Full HTML page
# ---------------------------------------------------------------------------

_HTML_SKELETON = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>GAIA report – {sample}</title>
  <link rel="stylesheet"
        href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css"
        integrity="sha384-QWTKZyjpPEjISv5WaRU9OFeRpok6YctnYmDr5pNlyT2bRjXh0JMhjY6hW+ALEwIH"
        crossorigin="anonymous"/>
  <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.3/dist/chart.umd.min.js"
          integrity="sha384-3hr+MFI1GBqiRkwfIiBlMH2wgFpKsVcBIi5FiLLnGBH25q/b1yx7sxkBNbJlBAC"
          crossorigin="anonymous"></script>
  <style>
    body {{ background: #f8f9fa; }}
    .navbar-brand {{ font-weight: 700; letter-spacing: .05em; }}
    canvas {{ background: #fff; border-radius: .25rem; padding: .5rem; }}
  </style>
</head>
<body>
<nav class="navbar navbar-dark bg-dark mb-4">
  <div class="container-fluid">
    <span class="navbar-brand">🧫 GAIA — Metagenomics Pipeline Report</span>
  </div>
</nav>
<div class="container pb-5">
{sections}
</div>
</body>
</html>
"""


def render(sample, run_mode, taxonomy_json, genome_summary,
           genome_inventory, pipeline_summary, output):
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    sections = []
    sections.append(_section_pipeline_info(sample, run_mode, generated_at))

    tax_section = _section_taxonomy(taxonomy_json)
    if tax_section:
        sections.append(tax_section)

    gs = _section_genome_summary(genome_summary)
    if gs:
        sections.append(gs)

    gi = _section_genome_inventory(genome_inventory)
    if gi:
        sections.append(gi)

    ps = _section_pipeline_summary(pipeline_summary)
    if ps:
        sections.append(ps)

    html = _HTML_SKELETON.format(
        sample=sample,
        sections="\n".join(sections),
    )
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
