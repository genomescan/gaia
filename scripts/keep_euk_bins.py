#!/usr/bin/env python3
import argparse
import shutil
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser(description="Keep only ACR-refined bins labeled as eukaryotic.")
    p.add_argument("--input-dir", required=True, help="Directory containing ACR refined FASTA bins.")
    p.add_argument("--output-dir", required=True, help="Directory where retained Euk bins will be written.")
    p.add_argument("--manifest", required=True, help="Manifest TSV listing retained bins.")
    return p.parse_args()

def is_fasta(path: Path) -> bool:
    suffixes = "".join(path.suffixes).lower()
    return path.is_file() and (
        path.suffix.lower() in {".fa", ".fna", ".fasta", ".fas"} or
        suffixes.endswith(".fa.gz") or suffixes.endswith(".fna.gz") or suffixes.endswith(".fasta.gz")
    )

def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    manifest = Path(args.manifest)

    if not input_dir.exists():
        raise SystemExit(f"ERROR: input directory does not exist: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)

    retained = []
    for fasta in sorted(input_dir.iterdir()):
        if not is_fasta(fasta):
            continue
        name = fasta.name
        if ".Euk." not in name and not name.endswith(".Euk.fa") and ".Euk.fa." not in name:
            continue

        dest = output_dir / name
        try:
            if dest.exists() or dest.is_symlink():
                dest.unlink()
            dest.symlink_to(fasta.resolve())
        except Exception:
            shutil.copy2(fasta, dest)

        retained.append((name, str(fasta), str(dest)))

    with manifest.open("w") as out:
        out.write("bin_name\tsource_path\tretained_path\n")
        for row in retained:
            out.write("\t".join(row) + "\n")

    if len(retained) == 0:
        print("WARNING: no Euk bins found in ACR output — continuing with empty set")

if __name__ == "__main__":
    main()
