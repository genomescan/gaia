#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") -r HOST_REF -i READS -o OUT_PREFIX [-t THREADS]

Required arguments:
  -r  Host reference FASTA (e.g. ref/host.fa)
  -i  Input reads (FASTQ or FASTQ.GZ, ONT long reads)
  -o  Output prefix (e.g. clean/sample1)

Optional:
  -t  Threads (default: 4)
  -h  Show this help message

Outputs:
  OUT_PREFIX.host.bam           BAM with reads aligned to host
  OUT_PREFIX.host.flagstat.txt  samtools flagstat summary on host BAM
  OUT_PREFIX.nohost.fastq.gz    FASTQ of reads not mapping to host
EOF
}

THREADS=4
HOST_REF=""
READS=""
OUT_PREFIX=""

while getopts ":r:i:o:t:h" opt; do
    case $opt in
        r) HOST_REF="$OPTARG" ;;
        i) READS="$OPTARG" ;;
        o) OUT_PREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage; exit 0 ;;
        \?)
            echo "Error: Invalid option -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Error: Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$HOST_REF" || -z "$READS" || -z "$OUT_PREFIX" ]]; then
    echo "Error: -r, -i and -o are required." >&2
    usage
    exit 1
fi

if [[ ! -f "$HOST_REF" ]]; then
    echo "Error: Host reference FASTA not found: $HOST_REF" >&2
    exit 1
fi

if [[ ! -f "$READS" ]]; then
    echo "Error: Reads file not found: $READS" >&2
    exit 1
fi

HOST_INDEX="${HOST_REF}.mmi"

echo "[INFO] Using host reference: $HOST_REF"
echo "[INFO] Host index: $HOST_INDEX"
echo "[INFO] Input reads: $READS"
echo "[INFO] Output prefix: $OUT_PREFIX"
echo "[INFO] Threads: $THREADS"

if [[ ! -f "$HOST_INDEX" ]]; then
    echo "[INFO] Host index not found. Building index with minimap2..."
    minimap2 -d "$HOST_INDEX" "$HOST_REF"
    echo "[INFO] Index built."
else
    echo "[INFO] Host index already exists. Reusing."
fi

BAM="${OUT_PREFIX}.host.bam"
FLAGSTAT="${OUT_PREFIX}.host.flagstat.txt"
NOHOST_FASTQ_GZ="${OUT_PREFIX}.nohost.fastq.gz"

# Align reads to host and create BAM
echo "[INFO] Aligning reads to host with minimap2..."
minimap2 -t "$THREADS" -ax map-ont "$HOST_INDEX" "$READS" \
  | samtools view -b -o "$BAM" -
echo "[INFO] Host-aligned BAM written to: $BAM"

# Flagstat QC
echo "[INFO] Running samtools flagstat on host BAM..."
samtools flagstat "$BAM" > "$FLAGSTAT"
echo "[INFO] Flagstat summary written to: $FLAGSTAT"

# Extract unmapped reads to FASTQ.GZ
echo "[INFO] Extracting unmapped (non-host) reads to FASTQ..."
samtools fastq -f 4 "$BAM" \
  | pigz -p "$THREADS" -c > "$NOHOST_FASTQ_GZ"
echo "[INFO] Non-host FASTQ written to: $NOHOST_FASTQ_GZ"

echo "[INFO] Host removal completed successfully."
