#!/usr/bin/env python3
import argparse, csv, os

def detect_delimiter(path):
    with open(path, "r", newline="") as fh:
        sample = fh.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except Exception:
        return "\t"

def read_table(path):
    if not path or not os.path.exists(path) or os.path.getsize(path)==0:
        return []
    delim=detect_delimiter(path)
    with open(path, newline="") as fh:
        reader=csv.DictReader(fh, delimiter=delim)
        return list(reader) if reader.fieldnames else []

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--kraken2-report", default="")
    ap.add_argument("--centrifuger-report", default="")
    ap.add_argument("--genome-summary", required=True)
    ap.add_argument("--inventory", required=True)
    ap.add_argument("--output", required=True)
    args=ap.parse_args()

    gen=read_table(args.genome_summary)
    inv=read_table(args.inventory)

    tools=[]
    if args.kraken2_report:
        tools.append("kraken2")
    if args.centrifuger_report:
        tools.append("centrifuger")

    gen0=gen[0] if gen else {}
    row={
        "sample": args.sample,
        "taxonomy_tools_run": ",".join(tools),
        "assembly_completed": "yes",
        "prok_final_bins": gen0.get("total_bins_prok","0"),
        "euk_final_bins": gen0.get("total_bins_euk","0"),
        "total_final_genomes": gen0.get("total_final_genomes","0"),
    }
    with open(args.output, "w", newline="") as fh:
        w=csv.DictWriter(fh, fieldnames=list(row.keys()), delimiter="\t")
        w.writeheader()
        w.writerow(row)

if __name__ == "__main__":
    main()
