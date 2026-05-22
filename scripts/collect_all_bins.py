#!/usr/bin/env python3
import argparse
import os
import shutil
import sys
from pathlib import Path

VALID_EXTS = (".fa", ".fasta", ".fna")

def detect_bins(bin_dir: Path):
    if not bin_dir.exists():
        return []
    return sorted(
        [p for p in bin_dir.iterdir() if p.is_file() and p.suffix.lower() in VALID_EXTS]
    )

def main():
    parser = argparse.ArgumentParser(description="Collect raw bins from multiple binners into one normalized directory.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--metabat2", default="")
    parser.add_argument("--semibin2", default="")
    parser.add_argument("--comebin", default="")
    args = parser.parse_args()

    sources = [
        ("MetaBAT2", Path(args.metabat2) if args.metabat2 else None),
        ("SemiBin2", Path(args.semibin2) if args.semibin2 else None),
        ("COMEBin", Path(args.comebin) if args.comebin else None),
    ]

    outdir = Path(args.outdir)
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest_path = Path(args.manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    collected = 0

    with manifest_path.open("w") as mf:
        mf.write("sample\tsource_tool\tsource_file\tcollected_name\tcollected_path\n")
        for source_name, source_dir in sources:
            if source_dir is None or not source_dir.exists():
                continue
            for i, src in enumerate(detect_bins(source_dir), start=1):
                dest_name = f"{source_name}.{src.stem}.fa"
                if dest_name in seen:
                    dest_name = f"{source_name}.{src.stem}.{i}.fa"
                seen.add(dest_name)
                dest = outdir / dest_name
                try:
                    os.symlink(src.resolve(), dest)
                except OSError:
                    shutil.copy2(src, dest)
                mf.write(f"{args.sample}\t{source_name}\t{src}\t{dest_name}\t{dest}\n")
                collected += 1

    if collected == 0:
        sys.stderr.write("ERROR: No bin FASTA files were found to collect.\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
