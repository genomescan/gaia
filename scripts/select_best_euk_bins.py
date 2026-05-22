#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import shutil
from pathlib import Path

FASTA_EXTS = (".fa", ".fna", ".fasta", ".fas", ".fa.gz", ".fna.gz", ".fasta.gz")


def parse_args():
    p = argparse.ArgumentParser(
        description="Select one final eukaryotic representative per dRep cluster using EukCC metrics."
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument("--clusters", required=True)
    p.add_argument("--eukcc", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--manifest", required=True)
    return p.parse_args()


def detect_delimiter(path: Path) -> str:
    with open(path, "r", newline="") as fh:
        sample = fh.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except Exception:
        return "\t"


def read_table(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    delim = detect_delimiter(path)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        return list(reader) if reader.fieldnames else []


def choose(rec, keys, default=""):
    for key in keys:
        if key in rec and rec[key] not in ("", None):
            return rec[key]
    return default


def safe_float(x, default=None):
    try:
        return float(str(x).strip().replace("%", ""))
    except Exception:
        return default


def strip_fasta_ext(name: str) -> str:
    lower = name.lower()
    for ext in FASTA_EXTS:
        if lower.endswith(ext):
            return name[: -len(ext)]
    return name


def fasta_files(folder: Path):
    out = []
    if not folder.exists():
        return out
    for p in sorted(folder.iterdir()):
        if p.is_file() and any(p.name.lower().endswith(ext) for ext in FASTA_EXTS):
            out.append(p)
    return out


def fasta_stats(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    lengths = []
    total = 0

    with opener(path, "rt") as fh:
        seq = []
        for line in fh:
            if line.startswith(">"):
                if seq:
                    s = "".join(seq)
                    l = len(s)
                    lengths.append(l)
                    total += l
                    seq = []
            else:
                seq.append(line.strip())
        if seq:
            s = "".join(seq)
            l = len(s)
            lengths.append(l)
            total += l

    lengths_sorted = sorted(lengths, reverse=True)
    csum = 0
    n50 = 0
    for l in lengths_sorted:
        csum += l
        if csum >= total / 2:
            n50 = l
            break

    return {
        "total_length": total,
        "N50": n50,
    }


def load_eukcc(path: Path):
    out = {}
    for rec in read_table(path):
        genome = choose(rec, ["bin", "Bin", "genome", "Genome", "name", "Name"], "")
        if not genome:
            continue
        genome = os.path.basename(genome)
        payload = {
            "completeness": safe_float(choose(rec, ["completeness", "Completeness"], ""), None),
            "contamination": safe_float(choose(rec, ["contamination", "Contamination"], ""), None),
        }
        out[genome] = payload
        out[strip_fasta_ext(genome)] = payload
    return out


def load_clusters(path: Path):
    out = {}
    for rec in read_table(path):
        genome = choose(rec, ["genome", "Genome"], "")
        cluster = choose(rec, ["secondary_cluster", "Secondary_cluster", "cluster", "Cluster"], "")
        if not genome:
            continue
        genome = os.path.basename(genome)
        cluster = cluster or strip_fasta_ext(genome)
        out[genome] = cluster
        out[strip_fasta_ext(genome)] = cluster
    return out


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    manifest_path = Path(args.manifest)

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    files = fasta_files(input_dir)
    if not files:
        with manifest_path.open("w") as out:
            out.write("cluster_id\tselected_bin\tselected_path\tcompleteness\tcontamination\ttotal_length\tN50\tselection_reason\n")
        return

    cluster_map = load_clusters(Path(args.clusters))
    eukcc = load_eukcc(Path(args.eukcc))
    grouped = {}

    for fasta in files:
        name = fasta.name
        stem = strip_fasta_ext(name)
        cluster_id = cluster_map.get(name) or cluster_map.get(stem) or stem
        stats = fasta_stats(fasta)
        q = eukcc.get(name) or eukcc.get(stem) or {}

        grouped.setdefault(cluster_id, []).append({
            "path": fasta,
            "name": name,
            "completeness": q.get("completeness"),
            "contamination": q.get("contamination"),
            "total_length": stats["total_length"],
            "N50": stats["N50"],
        })

    rows = []
    for cluster_id, recs in sorted(grouped.items()):
        def sort_key(rec):
            has_q = 0 if rec["completeness"] is not None else 1
            completeness = -(rec["completeness"] if rec["completeness"] is not None else -1)
            contamination = rec["contamination"] if rec["contamination"] is not None else 1e9
            return (has_q, completeness, contamination, -rec["total_length"], -rec["N50"], rec["name"])

        best = sorted(recs, key=sort_key)[0]
        dest = output_dir / best["name"]

        try:
            if dest.exists() or dest.is_symlink():
                dest.unlink()
            dest.symlink_to(best["path"].resolve())
        except Exception:
            shutil.copy2(best["path"], dest)

        rows.append({
            "cluster_id": cluster_id,
            "selected_bin": best["name"],
            "selected_path": str(dest),
            "completeness": "" if best["completeness"] is None else best["completeness"],
            "contamination": "" if best["contamination"] is None else best["contamination"],
            "total_length": best["total_length"],
            "N50": best["N50"],
            "selection_reason": "best_eukcc_then_size",
        })

    with manifest_path.open("w", newline="") as out:
        fields = [
            "cluster_id",
            "selected_bin",
            "selected_path",
            "completeness",
            "contamination",
            "total_length",
            "N50",
            "selection_reason",
        ]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


if __name__ == "__main__":
    main()
