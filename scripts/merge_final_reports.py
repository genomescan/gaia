#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
from pathlib import Path
from collections import defaultdict

FASTA_EXTS = (".fa", ".fna", ".fasta", ".fas", ".fa.gz", ".fna.gz", ".fasta.gz")


def detect_delimiter(path):
    with open(path, "r", newline="") as fh:
        sample = fh.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except Exception:
        return "\t"


def read_table(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    delim = detect_delimiter(path)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        return list(reader) if reader.fieldnames else []


def safe_float(x, default=""):
    try:
        return float(str(x).strip().replace("%", ""))
    except Exception:
        return default


def choose(d, keys, default=""):
    for k in keys:
        if k in d and d[k] not in ("", None):
            return d[k]
    return default


def fasta_files(folder):
    p = Path(folder)
    if not p.exists():
        return []
    out = []
    for f in sorted(p.iterdir()):
        if f.is_file() and any(f.name.lower().endswith(ext) for ext in FASTA_EXTS):
            out.append(f)
    return out


def strip_fasta_ext(name):
    lower = name.lower()
    for ext in FASTA_EXTS:
        if lower.endswith(ext):
            return name[: -len(ext)]
    return name


def fasta_stats(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    lengths = []
    gc = 0
    total = 0

    with opener(path, "rt") as fh:
        seq = []
        for line in fh:
            if line.startswith(">"):
                if seq:
                    s = "".join(seq).upper()
                    l = len(s)
                    lengths.append(l)
                    total += l
                    gc += s.count("G") + s.count("C")
                    seq = []
            else:
                seq.append(line.strip())

        if seq:
            s = "".join(seq).upper()
            l = len(s)
            lengths.append(l)
            total += l
            gc += s.count("G") + s.count("C")

    lengths_sorted = sorted(lengths, reverse=True)
    csum = 0
    n50 = 0
    for l in lengths_sorted:
        csum += l
        if csum >= total / 2:
            n50 = l
            break

    return {
        "contig_count": len(lengths),
        "total_length": total,
        "N50": n50,
        "longest_contig": max(lengths) if lengths else 0,
        "gc_content": round(gc / total, 5) if total else "",
    }


def load_checkm2(path):
    rows = read_table(path)
    out = {}
    for r in rows:
        gid = choose(r, ["Name", "Bin Id", "bin", "genome", "Genome", "Bin"], "")
        if gid:
            gid = os.path.basename(gid)
            gid_stripped = strip_fasta_ext(gid)
            payload = {
                "completeness": choose(r, ["Completeness", "completeness", "Completeness_General"], ""),
                "contamination": choose(r, ["Contamination", "contamination"], ""),
            }
            out[gid] = payload
            out[gid_stripped] = payload
    return out


def load_gtdb(*paths):
    out = {}
    for path in paths:
        rows = read_table(path)
        for r in rows:
            gid = choose(r, ["user_genome", "genome", "Genome"], "")
            if gid:
                gid = os.path.basename(gid)
                gid_stripped = strip_fasta_ext(gid)
                classification = choose(r, ["classification", "Classification"], "")
                domain = ""
                if classification:
                    parts = classification.split(";")
                    if parts:
                        domain = parts[0].replace("d__", "")
                payload = {"taxonomy": classification, "domain": domain}
                out[gid] = payload
                out[gid_stripped] = payload
    return out


def load_eukcc(path):
    rows = read_table(path)
    out = {}
    for r in rows:
        gid = choose(r, ["bin", "Bin", "genome", "Genome", "Name", "name"], "")
        if gid:
            gid = os.path.basename(gid)
            gid_stripped = strip_fasta_ext(gid)
            payload = {
                "completeness": choose(r, ["completeness", "Completeness"], ""),
                "contamination": choose(r, ["contamination", "Contamination"], ""),
            }
            out[gid] = payload
            out[gid_stripped] = payload
    return out


def load_bat(path):
    rows = read_table(path)
    out = {}
    for r in rows:
        gid = choose(r, ["Bin", "bin", "genome", "Genome"], "")
        if gid:
            gid = os.path.basename(gid)
            gid_stripped = strip_fasta_ext(gid)
            payload = {
                "taxonomy": choose(r, ["classification", "lineage", "taxon", "taxonomy"], ""),
                "domain": choose(r, ["superkingdom", "domain"], ""),
            }
            out[gid] = payload
            out[gid_stripped] = payload
    return out


def load_drep_clusters(path):
    rows = read_table(path)
    out = {}
    for r in rows:
        gid = choose(r, ["genome", "Genome", "representative", "cluster_member"], "")
        if gid:
            gid = os.path.basename(gid)
            gid_stripped = strip_fasta_ext(gid)
            cluster = choose(r, ["secondary_cluster", "cluster", "Secondary_cluster"], "")
            out[gid] = cluster
            out[gid_stripped] = cluster
    return out


def load_selected_manifest(path):
    rows = read_table(path)
    out = {}
    for r in rows:
        gid = choose(r, ["selected_bin", "bin", "genome"], "")
        if gid:
            gid = os.path.basename(gid)
            gid_stripped = strip_fasta_ext(gid)
            payload = {
                "cluster_id": choose(r, ["cluster_id"], ""),
                "selection_reason": choose(r, ["selection_reason"], ""),
            }
            out[gid] = payload
            out[gid_stripped] = payload
    return out


def load_dastool_map(path):
    counts = defaultdict(int)
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return counts
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                counts[parts[1]] += 1
    return counts


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--prok-bins", required=True)
    ap.add_argument("--checkm2", required=True)
    ap.add_argument("--gtdb-bac", required=True)
    ap.add_argument("--gtdb-ar", required=True)
    ap.add_argument("--dastool-map", required=True)
    ap.add_argument("--euk-bins", required=True)
    ap.add_argument("--eukcc", required=True)
    ap.add_argument("--bat", required=True)
    ap.add_argument("--drep-clusters", required=True)
    ap.add_argument("--kept-manifest", required=True)
    ap.add_argument("--selected-manifest", required=True)
    ap.add_argument("--inventory", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--trace", required=True)
    args = ap.parse_args()

    checkm2 = load_checkm2(args.checkm2)
    gtdb = load_gtdb(args.gtdb_bac, args.gtdb_ar)
    eukcc = load_eukcc(args.eukcc)
    bat = load_bat(args.bat)
    drep = load_drep_clusters(args.drep_clusters)
    selected = load_selected_manifest(args.selected_manifest)
    dastool_counts = load_dastool_map(args.dastool_map)

    inventory = []
    trace = []

    # ------------------------------------------------------------------
    # Prokaryotic final bins
    # ------------------------------------------------------------------
    for f in fasta_files(args.prok_bins):
        bid = f.name
        for ext in FASTA_EXTS:
            if bid.lower().endswith(ext):
                bid = bid[: -len(ext)]
                break

        stats = fasta_stats(f)
        q = checkm2.get(f.name, {}) or checkm2.get(bid, {})
        tax = gtdb.get(f.name, {}) or gtdb.get(bid, {})

        row = {
            "sample": args.sample,
            "bin_id": bid,
            "branch": "prok",
            "final_fasta_path": str(f),
            "source_binners": "metabat2,semibin2,comebin",
            "refinement_tool": "dastool",
            "dereplicated": "no",
            **stats,
            "completeness": q.get("completeness", ""),
            "contamination": q.get("contamination", ""),
            "taxonomy": tax.get("taxonomy", ""),
            "domain": tax.get("domain", ""),
        }
        inventory.append(row)

        trace.append({
            "sample": args.sample,
            "final_bin_id": bid,
            "branch": "prok",
            "origin_tool": "metabat2,semibin2,comebin",
            "origin_bin": "",
            "refined_bin": bid,
            "dereplicated_cluster": "",
            "notes": f"contigs_in_refined_bin={dastool_counts.get(bid, '')}",
        })

    # ------------------------------------------------------------------
    # Eukaryotic final selected bins
    # ------------------------------------------------------------------
    for f in fasta_files(args.euk_bins):
        bid = f.name
        for ext in FASTA_EXTS:
            if bid.lower().endswith(ext):
                bid = bid[: -len(ext)]
                break

        stats = fasta_stats(f)
        q = eukcc.get(f.name, {}) or eukcc.get(bid, {})
        tax = bat.get(f.name, {}) or bat.get(bid, {})
        sel = selected.get(f.name, {}) or selected.get(bid, {})
        cluster_id = sel.get("cluster_id", "") or drep.get(f.name, "") or drep.get(bid, "")

        row = {
            "sample": args.sample,
            "bin_id": bid,
            "branch": "euk",
            "final_fasta_path": str(f),
            "source_binners": "metabat2,semibin2,comebin",
            "refinement_tool": "drep_compare+eukcc_select",
            "dereplicated": "cluster-selected",
            **stats,
            "completeness": q.get("completeness", ""),
            "contamination": q.get("contamination", ""),
            "taxonomy": tax.get("taxonomy", ""),
            "domain": tax.get("domain", ""),
        }
        inventory.append(row)

        trace.append({
            "sample": args.sample,
            "final_bin_id": bid,
            "branch": "euk",
            "origin_tool": "metabat2,semibin2,comebin",
            "origin_bin": "",
            "refined_bin": bid,
            "dereplicated_cluster": cluster_id,
            "notes": sel.get("selection_reason", ""),
        })

    inv_fields = [
        "sample", "bin_id", "branch", "final_fasta_path",
        "source_binners", "refinement_tool", "dereplicated",
        "contig_count", "total_length", "N50", "longest_contig",
        "gc_content", "completeness", "contamination", "taxonomy", "domain"
    ]
    write_tsv(args.inventory, inventory, inv_fields)

    prok_rows = [r for r in inventory if r["branch"] == "prok"]
    euk_rows = [r for r in inventory if r["branch"] == "euk"]

    summary = [{
        "sample": args.sample,
        "total_final_genomes": len(inventory),
        "total_bins_prok": len(prok_rows),
        "total_bins_euk": len(euk_rows),
        "high_quality_prok": sum(
            1 for r in prok_rows
            if safe_float(r.get("completeness", ""), 0) >= 90 and safe_float(r.get("contamination", ""), 999) <= 5
        ),
        "medium_quality_prok": sum(
            1 for r in prok_rows
            if safe_float(r.get("completeness", ""), 0) >= 50 and safe_float(r.get("contamination", ""), 999) <= 10
        ),
        "euk_bins_with_taxonomy": sum(1 for r in euk_rows if r.get("taxonomy", "")),
    }]

    write_tsv(
        args.summary,
        summary,
        [
            "sample",
            "total_final_genomes",
            "total_bins_prok",
            "total_bins_euk",
            "high_quality_prok",
            "medium_quality_prok",
            "euk_bins_with_taxonomy",
        ],
    )

    write_tsv(
        args.trace,
        trace,
        [
            "sample",
            "final_bin_id",
            "branch",
            "origin_tool",
            "origin_bin",
            "refined_bin",
            "dereplicated_cluster",
            "notes",
        ],
    )


if __name__ == "__main__":
    main()
