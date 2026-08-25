#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<USAGE
Usage: $(basename "$0") -i READS -r REFERENCE -q INSERT_FASTA -o OUTDIR [options]

Identify candidate insertion loci in PacBio HiFi reads using only the insert
sequence as bait. The workflow:
  1. converts supported inputs to FASTA when needed
  2. finds insert-positive reads by aligning reads to the insert sequence
  3. extracts left/right read flanks around the insert hit
  4. aligns flanks to the reference genome
  5. summarizes per-read and per-site insertion candidates

Required arguments:
  -i  Input reads: FASTQ, FASTQ.GZ, FQ, FQ.GZ, FASTA, FASTA.GZ, FA, FA.GZ or BAM
  -r  Reference genome FASTA
  -q  Insert FASTA (single insert sequence used as bait)
  -o  Output directory

Optional arguments:
  -t  Threads for minimap2/samtools (default: 4)
  -f  Flank size to extract on each side of the insert hit (default: 1000)
  -m  Minimum flank length to keep (default: 100)
  -b  Minimum insert-aligned bases required to keep a read (default: 100)
  -Q  Minimum MAPQ to call a flank unique enough for site inference (default: 20)
  -g  Maximum gap when merging insert hit segments on one read (default: 200)
  -a  Aligner for mapping steps: minimap2 or auto (default: auto)
      auto currently resolves to minimap2, which is the aligner already used by this repository.
  -k  Keep the normalized reads FASTA instead of deleting it at the end
  -h  Show this help message

Outputs written to OUTDIR:
  insert_sequence.fasta        normalized single-record insert FASTA
  reads.normalized.fasta       normalized reads FASTA (removed unless -k is set)
  insert_hits.paf              raw read-vs-insert alignments
  insert_hits.tsv              per-read insert hit summary
  flanks.fasta                 extracted left/right flanks
  flank_alignments.paf         raw flank-vs-genome alignments
  flank_hits.tsv               best flank mapping summary
  candidate_insertions.tsv     per-read candidate insertion calls
  candidate_sites.tsv          grouped site-level support summary
USAGE
}

THREADS=4
FLANK_SIZE=1000
MIN_FLANK_LEN=100
MIN_INSERT_BASES=100
MIN_MAPQ=20
MERGE_GAP=200
ALIGNER=auto
KEEP_INTERMEDIATE=0
INPUT=""
REFERENCE=""
INSERT_FASTA=""
OUTDIR=""

while getopts ":i:r:q:o:t:f:m:b:Q:g:a:kh" opt; do
    case "$opt" in
        i) INPUT="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        q) INSERT_FASTA="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        f) FLANK_SIZE="$OPTARG" ;;
        m) MIN_FLANK_LEN="$OPTARG" ;;
        b) MIN_INSERT_BASES="$OPTARG" ;;
        Q) MIN_MAPQ="$OPTARG" ;;
        g) MERGE_GAP="$OPTARG" ;;
        a) ALIGNER="$OPTARG" ;;
        k) KEEP_INTERMEDIATE=1 ;;
        h) usage; exit 0 ;;
        \?)
            echo "Error: invalid option -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Error: option -$OPTARG requires an argument" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$INPUT" || -z "$REFERENCE" || -z "$INSERT_FASTA" || -z "$OUTDIR" ]]; then
    echo "Error: -i, -r, -q and -o are required." >&2
    usage
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "Error: input reads not found: $INPUT" >&2
    exit 1
fi
if [[ ! -f "$REFERENCE" ]]; then
    echo "Error: reference FASTA not found: $REFERENCE" >&2
    exit 1
fi
if [[ ! -f "$INSERT_FASTA" ]]; then
    echo "Error: insert FASTA not found: $INSERT_FASTA" >&2
    exit 1
fi

case "$ALIGNER" in
    auto|minimap2) ;;
    *)
        echo "Error: unsupported aligner '$ALIGNER'. Use 'auto' or 'minimap2'." >&2
        exit 1
        ;;
esac

command -v minimap2 >/dev/null 2>&1 || {
    echo "Error: minimap2 is required but was not found on PATH." >&2
    exit 1
}
command -v python3 >/dev/null 2>&1 || {
    echo "Error: python3 is required but was not found on PATH." >&2
    exit 1
}

if [[ "$INPUT" == *.bam ]]; then
    command -v samtools >/dev/null 2>&1 || {
        echo "Error: samtools is required for BAM input but was not found on PATH." >&2
        exit 1
    }
fi

mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)

READS_FASTA="$OUTDIR/reads.normalized.fasta"
INSERT_NORMALIZED="$OUTDIR/insert_sequence.fasta"
INSERT_PAF="$OUTDIR/insert_hits.paf"
INSERT_TSV="$OUTDIR/insert_hits.tsv"
FLANKS_FASTA="$OUTDIR/flanks.fasta"
FLANKS_PAF="$OUTDIR/flank_alignments.paf"
FLANK_HITS_TSV="$OUTDIR/flank_hits.tsv"
CANDIDATE_TSV="$OUTDIR/candidate_insertions.tsv"
SITE_TSV="$OUTDIR/candidate_sites.tsv"

log() {
    echo "[INFO] $*"
}

log "Input reads: $INPUT"
log "Reference genome: $REFERENCE"
log "Insert FASTA: $INSERT_FASTA"
log "Output directory: $OUTDIR"
log "Threads: $THREADS"
log "Flank size: $FLANK_SIZE"
log "Min flank length: $MIN_FLANK_LEN"
log "Min insert-aligned bases: $MIN_INSERT_BASES"
log "Min flank MAPQ: $MIN_MAPQ"
log "Merge gap: $MERGE_GAP"

export INPUT REFERENCE INSERT_FASTA OUTDIR READS_FASTA INSERT_NORMALIZED INSERT_PAF INSERT_TSV FLANKS_FASTA \
    MIN_INSERT_BASES FLANK_SIZE MIN_FLANK_LEN MERGE_GAP

log "Normalizing insert FASTA and input reads"
python3 - <<'PY'
import gzip
import os
import sys
from pathlib import Path

input_path = Path(os.environ["INPUT"])
insert_path = Path(os.environ["INSERT_FASTA"])
reads_out = Path(os.environ["READS_FASTA"])
insert_out = Path(os.environ["INSERT_NORMALIZED"])


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def detect_text_format(path: Path):
    with open_text(path) as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"): 
                return "fasta"
            if line.startswith("@"):
                return "fastq"
            break
    raise SystemExit(f"Could not determine sequence format for {path}")


def iter_fasta(path: Path):
    name = None
    seq = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"): 
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def iter_fastq(path: Path):
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise SystemExit(f"Malformed FASTQ record in {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise SystemExit(f"Malformed FASTQ record in {path}")
            yield header[1:].split()[0], seq.strip()


def write_wrapped(handle, name, sequence, width=80):
    handle.write(f">{name}\n")
    for i in range(0, len(sequence), width):
        handle.write(sequence[i:i+width] + "\n")


insert_records = [(name, seq.upper()) for name, seq in iter_fasta(insert_path) if seq]
if not insert_records:
    raise SystemExit("Insert FASTA did not contain any sequence")

insert_name, insert_seq = insert_records[0]
if len(insert_records) > 1:
    sys.stderr.write(
        f"[WARN] Insert FASTA contains {len(insert_records)} records; using the first record '{insert_name}'.\n"
    )
if not set(insert_seq).issubset(set("ACGTRYSWKMBDHVN")):
    sys.stderr.write("[WARN] Insert sequence contains non-IUPAC characters; keeping sequence as provided.\n")

with open(insert_out, "w", encoding="utf-8") as handle:
    write_wrapped(handle, insert_name, insert_seq)

if input_path.suffix.lower() != ".bam":
    fmt = detect_text_format(input_path)
    iterator = iter_fastq(input_path) if fmt == "fastq" else iter_fasta(input_path)
    count = 0
    with open(reads_out, "w", encoding="utf-8") as handle:
        for name, seq in iterator:
            if not seq:
                continue
            write_wrapped(handle, name, seq.upper())
            count += 1

    if count == 0:
        raise SystemExit(f"No reads were found in {input_path}")
PY

if [[ "$INPUT" == *.bam ]]; then
    log "Converting BAM input to FASTA with samtools"
    samtools fasta -n "$INPUT" > "$READS_FASTA"
fi

if [[ ! -s "$READS_FASTA" ]]; then
    echo "Error: normalized reads FASTA was not created: $READS_FASTA" >&2
    exit 1
fi

log "Finding insert-positive reads by aligning reads to the insert sequence"
minimap2 -x map-hifi -c --secondary=yes -N 20 -t "$THREADS" "$INSERT_NORMALIZED" "$READS_FASTA" > "$INSERT_PAF"

log "Summarizing insert hits and extracting flanks"
python3 - <<'PY'
import csv
import os
from collections import defaultdict
from pathlib import Path

reads_fasta = Path(os.environ["READS_FASTA"])
insert_paf = Path(os.environ["INSERT_PAF"])
insert_tsv = Path(os.environ["INSERT_TSV"])
flanks_fasta = Path(os.environ["FLANKS_FASTA"])
min_insert_bases = int(os.environ["MIN_INSERT_BASES"])
flank_size = int(os.environ["FLANK_SIZE"])
min_flank_len = int(os.environ["MIN_FLANK_LEN"])
merge_gap = int(os.environ["MERGE_GAP"])


def read_fasta(path):
    name = None
    seq = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"): 
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def parse_tags(fields):
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def read_paf(path):
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            qname = fields[0]
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            strand = fields[4]
            tname = fields[5]
            tlen = int(fields[6])
            tstart = int(fields[7])
            tend = int(fields[8])
            matches = int(fields[9])
            block_len = int(fields[10])
            mapq = int(fields[11])
            tags = parse_tags(fields[12:])
            yield {
                "qname": qname,
                "qlen": qlen,
                "qstart": qstart,
                "qend": qend,
                "strand": strand,
                "tname": tname,
                "tlen": tlen,
                "tstart": tstart,
                "tend": tend,
                "matches": matches,
                "block_len": block_len,
                "mapq": mapq,
                "tags": tags,
            }


def cluster_hits(hits):
    by_strand = defaultdict(list)
    for hit in hits:
        aligned = hit["qend"] - hit["qstart"]
        if aligned <= 0:
            continue
        by_strand[hit["strand"]].append(hit)

    best = None
    for strand, strand_hits in by_strand.items():
        strand_hits = sorted(strand_hits, key=lambda h: (h["qstart"], h["qend"]))
        clusters = []
        current = None
        for hit in strand_hits:
            if current is None:
                current = {
                    "strand": strand,
                    "start": hit["qstart"],
                    "end": hit["qend"],
                    "hits": [hit],
                    "total_aligned": hit["qend"] - hit["qstart"],
                    "best_mapq": hit["mapq"],
                    "best_matches": hit["matches"],
                }
                continue
            if hit["qstart"] <= current["end"] + merge_gap:
                current["end"] = max(current["end"], hit["qend"])
                current["hits"].append(hit)
                current["total_aligned"] += hit["qend"] - hit["qstart"]
                current["best_mapq"] = max(current["best_mapq"], hit["mapq"])
                current["best_matches"] = max(current["best_matches"], hit["matches"])
            else:
                clusters.append(current)
                current = {
                    "strand": strand,
                    "start": hit["qstart"],
                    "end": hit["qend"],
                    "hits": [hit],
                    "total_aligned": hit["qend"] - hit["qstart"],
                    "best_mapq": hit["mapq"],
                    "best_matches": hit["matches"],
                }
        if current is not None:
            clusters.append(current)

        clusters.sort(key=lambda c: (c["total_aligned"], c["end"] - c["start"], c["best_matches"], c["best_mapq"]), reverse=True)
        chosen = clusters[0]
        chosen["cluster_count"] = len(clusters)
        chosen["strand_total_aligned"] = sum(c["total_aligned"] for c in clusters)
        if best is None or (
            chosen["total_aligned"],
            chosen["end"] - chosen["start"],
            chosen["best_matches"],
            chosen["best_mapq"],
        ) > (
            best["total_aligned"],
            best["end"] - best["start"],
            best["best_matches"],
            best["best_mapq"],
        ):
            best = chosen

    return best


reads = {name: seq for name, seq in read_fasta(reads_fasta)}
hits_by_read = defaultdict(list)
for hit in read_paf(insert_paf):
    hits_by_read[hit["qname"]].append(hit)

with open(insert_tsv, "w", encoding="utf-8", newline="") as insert_handle, \
     open(flanks_fasta, "w", encoding="utf-8") as flank_handle:
    writer = csv.writer(insert_handle, delimiter="\t")
    writer.writerow([
        "read_id",
        "read_length",
        "insert_strand",
        "insert_read_start",
        "insert_read_end",
        "insert_span_bases",
        "segment_count",
        "cluster_count",
        "total_aligned_bases",
        "best_segment_matches",
        "best_segment_mapq",
        "left_flank_length",
        "right_flank_length",
        "status",
    ])

    kept_reads = 0
    for read_id in sorted(reads):
        seq = reads[read_id]
        cluster = cluster_hits(hits_by_read.get(read_id, []))
        if cluster is None or cluster["total_aligned"] < min_insert_bases:
            continue

        start = max(0, cluster["start"])
        end = min(len(seq), cluster["end"])
        left = seq[max(0, start - flank_size):start]
        right = seq[end:min(len(seq), end + flank_size)]

        left_len = len(left)
        right_len = len(right)
        status = []
        if cluster["cluster_count"] > 1:
            status.append("multiple_insert_clusters")
        if left_len < min_flank_len:
            status.append("short_left_flank")
        if right_len < min_flank_len:
            status.append("short_right_flank")
        if not status:
            status.append("ok")

        writer.writerow([
            read_id,
            len(seq),
            cluster["strand"],
            start,
            end,
            end - start,
            len(cluster["hits"]),
            cluster["cluster_count"],
            cluster["total_aligned"],
            cluster["best_matches"],
            cluster["best_mapq"],
            left_len,
            right_len,
            ",".join(status),
        ])

        if left_len >= min_flank_len:
            flank_handle.write(f">{read_id}__left\n{left}\n")
        if right_len >= min_flank_len:
            flank_handle.write(f">{read_id}__right\n{right}\n")
        kept_reads += 1

if kept_reads == 0:
    raise SystemExit(
        f"No insert-positive reads met the minimum aligned-bases threshold ({min_insert_bases})."
    )
PY

if [[ ! -s "$FLANKS_FASTA" ]]; then
    echo "Error: no flanks were extracted; try lowering -m or -b." >&2
    exit 1
fi

log "Aligning extracted flanks to the reference genome with minimap2"
minimap2 -x map-hifi -c --secondary=yes -N 20 -t "$THREADS" "$REFERENCE" "$FLANKS_FASTA" > "$FLANKS_PAF"

log "Summarizing flank mappings and candidate insertion loci"
export INSERT_PAF FLANKS_PAF INSERT_TSV FLANK_HITS_TSV CANDIDATE_TSV SITE_TSV MIN_MAPQ
python3 - <<'PY'
import csv
import os
from collections import defaultdict
from pathlib import Path

insert_tsv = Path(os.environ["INSERT_TSV"])
flanks_paf = Path(os.environ["FLANKS_PAF"])
flank_hits_tsv = Path(os.environ["FLANK_HITS_TSV"])
candidate_tsv = Path(os.environ["CANDIDATE_TSV"])
site_tsv = Path(os.environ["SITE_TSV"])
min_mapq = int(os.environ["MIN_MAPQ"])


def parse_tags(fields):
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def read_paf(path):
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            yield {
                "qname": fields[0],
                "qlen": int(fields[1]),
                "qstart": int(fields[2]),
                "qend": int(fields[3]),
                "strand": fields[4],
                "tname": fields[5],
                "tlen": int(fields[6]),
                "tstart": int(fields[7]),
                "tend": int(fields[8]),
                "matches": int(fields[9]),
                "block_len": int(fields[10]),
                "mapq": int(fields[11]),
                "tags": parse_tags(fields[12:]),
            }


def choose_best_alignment(records, side):
    if not records:
        return {
            "side": side,
            "status": "unmapped",
            "mapq": 0,
            "is_unique": False,
            "is_anchored": False,
            "chrom": ".",
            "strand": ".",
            "junction_0based": "",
            "junction_1based": "",
            "query_start": "",
            "query_end": "",
            "query_length": "",
            "target_start": "",
            "target_end": "",
            "matches": "",
            "aligned_block_bases": "",
            "secondary_count": 0,
            "note": "no_alignment",
        }

    def rank_key(rec):
        tp = rec["tags"].get("tp", "")
        primary = 1 if tp == "P" else 0
        return (primary, rec["mapq"], rec["matches"], rec["block_len"])

    ordered = sorted(records, key=rank_key, reverse=True)
    best = ordered[0]
    secondary_count = max(0, len(ordered) - 1)
    tied_best = [
        rec for rec in ordered[1:]
        if rec["mapq"] == best["mapq"] and rec["matches"] == best["matches"] and rec["block_len"] == best["block_len"]
    ]

    if side == "left":
        anchored = best["qend"] == best["qlen"]
        junction0 = best["tend"] if best["strand"] == "+" else best["tstart"]
        overhang = best["qlen"] - best["qend"]
    else:
        anchored = best["qstart"] == 0
        junction0 = best["tstart"] if best["strand"] == "+" else best["tend"]
        overhang = best["qstart"]

    notes = []
    if not anchored:
        notes.append(f"junction_unanchored_overhang={overhang}")
    if tied_best:
        notes.append("best_alignment_tied")
    if best["mapq"] < min_mapq:
        notes.append("low_mapq")

    is_unique = anchored and best["mapq"] >= min_mapq and not tied_best
    status = "unique" if is_unique else "ambiguous"
    if not anchored:
        status = "unanchored"
    elif best["mapq"] < min_mapq:
        status = "low_mapq"

    return {
        "side": side,
        "status": status,
        "mapq": best["mapq"],
        "is_unique": is_unique,
        "is_anchored": anchored,
        "chrom": best["tname"],
        "strand": best["strand"],
        "junction_0based": junction0,
        "junction_1based": junction0 + 1,
        "query_start": best["qstart"],
        "query_end": best["qend"],
        "query_length": best["qlen"],
        "target_start": best["tstart"],
        "target_end": best["tend"],
        "matches": best["matches"],
        "aligned_block_bases": best["block_len"],
        "secondary_count": secondary_count,
        "note": ",".join(notes) if notes else "ok",
    }


def flip_strand(strand):
    return "+" if strand == "-" else "-"


def infer_insert_orientation(insert_strand, left, right):
    read_vs_genome = None
    if left["is_unique"]:
        read_vs_genome = left["strand"]
    elif right["is_unique"]:
        read_vs_genome = right["strand"]

    if read_vs_genome not in {"+", "-"} or insert_strand not in {"+", "-"}:
        return "."
    return insert_strand if read_vs_genome == "+" else flip_strand(insert_strand)


insert_rows = {}
with open(insert_tsv, "r", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        insert_rows[row["read_id"]] = row

by_flank = defaultdict(list)
for rec in read_paf(flanks_paf):
    by_flank[rec["qname"]].append(rec)

flank_summaries = {}
with open(flank_hits_tsv, "w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "read_id",
        "flank_side",
        "mapping_status",
        "chromosome",
        "junction_0based",
        "junction_1based",
        "genome_strand",
        "mapq",
        "query_start",
        "query_end",
        "query_length",
        "target_start",
        "target_end",
        "matches",
        "aligned_block_bases",
        "secondary_alignment_count",
        "note",
    ])

    for flank_name in sorted(by_flank):
        if flank_name.endswith("__left"):
            read_id = flank_name[:-6]
            side = "left"
        elif flank_name.endswith("__right"):
            read_id = flank_name[:-7]
            side = "right"
        else:
            read_id = flank_name
            side = "unknown"

        summary = choose_best_alignment(by_flank[flank_name], side)
        flank_summaries[(read_id, side)] = summary
        writer.writerow([
            read_id,
            side,
            summary["status"],
            summary["chrom"],
            summary["junction_0based"],
            summary["junction_1based"],
            summary["strand"],
            summary["mapq"],
            summary["query_start"],
            summary["query_end"],
            summary["query_length"],
            summary["target_start"],
            summary["target_end"],
            summary["matches"],
            summary["aligned_block_bases"],
            summary["secondary_count"],
            summary["note"],
        ])

# Ensure present even for missing sides.
for read_id in insert_rows:
    flank_summaries.setdefault((read_id, "left"), choose_best_alignment([], "left"))
    flank_summaries.setdefault((read_id, "right"), choose_best_alignment([], "right"))

site_counts = {}
with open(candidate_tsv, "w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "read_id",
        "insert_strand_in_read",
        "insert_orientation_vs_genome",
        "candidate_status",
        "candidate_chromosome",
        "left_junction_1based",
        "right_junction_1based",
        "left_mapq",
        "right_mapq",
        "left_flank_status",
        "right_flank_status",
        "flank_support",
        "junction_span_bp",
        "insert_read_start",
        "insert_read_end",
        "insert_span_bases",
        "insert_hit_status",
        "notes",
    ])

    for read_id in sorted(insert_rows):
        row = insert_rows[read_id]
        left = flank_summaries[(read_id, "left")]
        right = flank_summaries[(read_id, "right")]
        insert_strand = row["insert_strand"]
        orientation = infer_insert_orientation(insert_strand, left, right)

        status = "unresolved"
        chrom = "."
        left_j = left["junction_1based"] if left["is_unique"] else ""
        right_j = right["junction_1based"] if right["is_unique"] else ""
        span = ""
        notes = []

        flank_support = int(left["is_unique"]) + int(right["is_unique"])
        if left["is_unique"] and right["is_unique"]:
            if left["chrom"] == right["chrom"]:
                chrom = left["chrom"]
                status = "paired_same_chromosome"
                span = right["junction_0based"] - left["junction_0based"]
            else:
                status = "paired_discordant_chromosomes"
                chrom = f"{left['chrom']}|{right['chrom']}"
        elif left["is_unique"]:
            status = "left_only"
            chrom = left["chrom"]
        elif right["is_unique"]:
            status = "right_only"
            chrom = right["chrom"]

        if left["status"] not in {"unique", "unmapped"}:
            notes.append(f"left:{left['note']}")
        if right["status"] not in {"unique", "unmapped"}:
            notes.append(f"right:{right['note']}")
        if row["status"] != "ok":
            notes.append(f"insert:{row['status']}")

        writer.writerow([
            read_id,
            insert_strand,
            orientation,
            status,
            chrom,
            left_j,
            right_j,
            left["mapq"],
            right["mapq"],
            left["status"],
            right["status"],
            flank_support,
            span,
            row["insert_read_start"],
            row["insert_read_end"],
            row["insert_span_bases"],
            row["status"],
            ";".join(notes) if notes else "ok",
        ])

        if status == "paired_same_chromosome":
            key = ("paired", chrom, left_j, right_j, orientation)
            label = f"{chrom}:{left_j}-{right_j}"
        elif status == "left_only":
            key = ("single", chrom, "left", left_j, orientation)
            label = f"{chrom}:{left_j}"
        elif status == "right_only":
            key = ("single", chrom, "right", right_j, orientation)
            label = f"{chrom}:{right_j}"
        else:
            continue

        if key not in site_counts:
            site_counts[key] = {
                "locus": label,
                "status": status,
                "chrom": chrom,
                "left_junction_1based": left_j,
                "right_junction_1based": right_j,
                "orientation": orientation,
                "supporting_reads": 0,
                "left_support": 0,
                "right_support": 0,
                "best_left_mapq": 0,
                "best_right_mapq": 0,
            }
        site_counts[key]["supporting_reads"] += 1
        site_counts[key]["left_support"] += int(left["is_unique"])
        site_counts[key]["right_support"] += int(right["is_unique"])
        site_counts[key]["best_left_mapq"] = max(site_counts[key]["best_left_mapq"], left["mapq"])
        site_counts[key]["best_right_mapq"] = max(site_counts[key]["best_right_mapq"], right["mapq"])

with open(site_tsv, "w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "candidate_locus",
        "candidate_status",
        "chromosome",
        "left_junction_1based",
        "right_junction_1based",
        "insert_orientation_vs_genome",
        "supporting_reads",
        "left_flank_support",
        "right_flank_support",
        "best_left_mapq",
        "best_right_mapq",
    ])
    for record in sorted(site_counts.values(), key=lambda r: (-r["supporting_reads"], r["chrom"], str(r["left_junction_1based"]), str(r["right_junction_1based"]))):
        writer.writerow([
            record["locus"],
            record["status"],
            record["chrom"],
            record["left_junction_1based"],
            record["right_junction_1based"],
            record["orientation"],
            record["supporting_reads"],
            record["left_support"],
            record["right_support"],
            record["best_left_mapq"],
            record["best_right_mapq"],
        ])
PY

if [[ "$KEEP_INTERMEDIATE" -eq 0 ]]; then
    rm -f "$READS_FASTA"
fi

log "Done. Key outputs:"
log "  $INSERT_TSV"
log "  $FLANK_HITS_TSV"
log "  $CANDIDATE_TSV"
log "  $SITE_TSV"
