#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tla_insert_softclip_msa.py - build a pileup / MSA of soft-clipped read ends
against a transgene insert (cassette) sequence for a list of TLA insertion sites.

Rationale
---------
Reads are mapped against the *host* genome only, so the insert-derived part of a
read that spans an insertion junction is soft-clipped. Collecting those clips per
insertion site and aligning them back to the insert sequence gives a per-site
pileup of the insert. Comparing the pileups across all sites shows whether a
candidate SNP in the insert/cassette is present at a single site only (i.e. one
integrated copy carries the variant) or at all sites.

Outputs (per --outdir)
----------------------
  pileup/site<N>_<chrom>_<pos>.pileup.txt  human readable MSA / pileup (wrapped)
  msa_txt/site<N>_..._msa.tsv              one aligned row per clip (machine readable)
  plots/site<N>_...png                     MSA heatmap + depth / non-ref fraction
  tla_insert_softclip_report.pdf           all plots + cross-site overview
  site_summary.tsv                         per-site QC summary
  insert_variants.tsv                      every candidate variant position per site
  variant_matrix.tsv                       sites x candidate positions, non-ref fraction
  pileup_counts.tsv                        (optional) full per-position base counts

Main improvements over the first version
----------------------------------------
* clips are only kept when the clip *boundary* is within the window of the site
  (instead of any read overlapping the window), which removes off-target clips;
* reads are properly filtered (secondary / duplicate / QC-fail, mapping quality)
  and clips are sub-sampled by reservoir sampling instead of stopping at the
  first N reads (which biased the pileup towards the left of the window);
* alignment is seed-anchored: k-mer diagonal voting gives an ungapped placement
  in O(len(clip)); the O(n*m) pure-Python Smith-Waterman is gone. If `edlib` is
  installed it is used for gapped (indel-aware) alignment;
* base qualities are used: low quality bases are shown in lower case and are not
  counted for variant calling, which is essential to trust a single SNP;
* per-position statistics keep strand information, so strand-biased artefacts can
  be filtered out;
* variant calling with depth / allele fraction / allele count / strand filters,
  plus a cross-site variant matrix and a cross-site heatmap that immediately show
  whether a variant is site specific;
* sites can be supplied in a TSV file instead of being hard-coded;
* parallel execution shares the insert sequence and k-mer index through a worker
  initializer instead of pickling them for every site.
"""

import argparse
import csv
import json
import os
import random
import sys
import time
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import pysam
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

try:  # optional, only used for gapped alignment
    import edlib

    HAVE_EDLIB = True
except ImportError:  # pragma: no cover - depends on the environment
    edlib = None
    HAVE_EDLIB = False


# --------------------------------------------------------------------------- #
# sites
# --------------------------------------------------------------------------- #

DEFAULT_SITES = [
    {"site": "1", "chrom": "NW_023276806.1", "pos": 25049554, "gene": "-"},
    {"site": "2", "chrom": "NW_023276806.1", "pos": 86062356, "gene": "Shox2"},
    {"site": "3", "chrom": "NW_023276806.1", "pos": 140526348, "gene": "-"},
    {"site": "4", "chrom": "NW_023276806.1", "pos": 144178550, "gene": "-"},
    {"site": "5", "chrom": "NW_023276806.1", "pos": 144179076, "gene": "-"},
    {"site": "7", "chrom": "NW_023276807.1", "pos": 188752223, "gene": "Nek4"},
    {"site": "8", "chrom": "NC_048595.1", "pos": 325802138, "gene": "Zswim6"},
    {"site": "9", "chrom": "NC_048595.1", "pos": 354895962, "gene": "-"},
    {"site": "10", "chrom": "NC_048596.1", "pos": 56258446, "gene": "-"},
    {"site": "11", "chrom": "NC_048596.1", "pos": 141603535, "gene": "Herc2"},
    {"site": "12", "chrom": "NC_048597.1", "pos": 64605473, "gene": "-"},
    {"site": "13", "chrom": "NC_048597.1", "pos": 79238191, "gene": "Wipi2"},
    {"site": "14", "chrom": "NC_048597.1", "pos": 98503120, "gene": "-"},
    {"site": "15", "chrom": "NC_048597.1", "pos": 148118184, "gene": "-"},
    {"site": "16", "chrom": "NC_048597.1", "pos": 161065561, "gene": "-"},
    {"site": "17", "chrom": "NC_048597.1", "pos": 164899642, "gene": "Ddx10"},
    {"site": "18", "chrom": "NC_048598.1", "pos": 123523511, "gene": "-"},
    {"site": "19", "chrom": "NC_048599.1", "pos": 1299538, "gene": "Tcfl5"},
    {"site": "20", "chrom": "NC_048599.1", "pos": 28422608, "gene": "-"},
    {"site": "21", "chrom": "NC_048599.1", "pos": 97786256, "gene": "Sestd1"},
    {"site": "22", "chrom": "NC_048604.1", "pos": 40407124, "gene": "-"},
    {"site": "23", "chrom": "NC_048604.1", "pos": 67788949, "gene": "Hdac8"},
]


def load_sites(path):
    """Read sites from a TSV/CSV with columns site, chrom, pos and optional gene."""
    if not path:
        return [dict(s) for s in DEFAULT_SITES]
    delim = "," if path.endswith(".csv") else "\t"
    sites = []
    with open(path, newline="") as fh:
        reader = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter=delim)
        missing = {"site", "chrom", "pos"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing column(s) {sorted(missing)}")
        for row in reader:
            sites.append(
                {
                    "site": str(row["site"]).strip(),
                    "chrom": str(row["chrom"]).strip(),
                    "pos": int(str(row["pos"]).strip()),
                    "gene": (row.get("gene") or "-").strip() or "-",
                }
            )
    if not sites:
        raise ValueError(f"{path}: no sites found")
    return sites


def site_sort_key(site_id):
    try:
        return (0, int(site_id), "")
    except ValueError:
        return (1, 0, str(site_id))


# --------------------------------------------------------------------------- #
# sequence helpers
# --------------------------------------------------------------------------- #

_COMP = str.maketrans("ACGTURYKMBDHVNacgturykmbdhvn", "TGCAAYRMKVHDBNtgcaayrmkvhdbn")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def read_first_fasta_seq(path):
    """Return the first record of a FASTA file as an upper case string."""
    seq, name, started = [], None, False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    break
                name = line[1:].strip().split()[0] if len(line) > 1 else "insert"
                started = True
                continue
            seq.append(line.strip())
    return name or "insert", "".join(seq).upper()


def parse_positions(s):
    if not s or not s.strip():
        return []
    return sorted({int(x.strip()) for x in s.split(",") if x.strip()})


def build_kmer_index(ref, k, max_hits):
    """k-mer -> list of start offsets; over-represented k-mers are dropped."""
    idx = defaultdict(list)
    for i in range(len(ref) - k + 1):
        idx[ref[i : i + k]].append(i)
    return {km: pos for km, pos in idx.items() if len(pos) <= max_hits}


# --------------------------------------------------------------------------- #
# soft clip extraction
# --------------------------------------------------------------------------- #

BAM_CSOFT_CLIP = 4


class Clip:
    """A single soft-clipped read end, always stored in read orientation."""

    __slots__ = ("read_name", "side", "seq", "qual", "mapq", "is_reverse", "ref_pos")

    def __init__(self, read_name, side, seq, qual, mapq, is_reverse, ref_pos):
        self.read_name = read_name
        self.side = side
        self.seq = seq
        self.qual = qual
        self.mapq = mapq
        self.is_reverse = is_reverse
        self.ref_pos = ref_pos


def extract_clips(
    bam_path,
    chrom,
    pos1,
    window=20,
    mapq=20,
    min_clip=15,
    max_clips=300,
    keep_duplicates=False,
    seed=0,
):
    """Collect soft clips whose clip boundary sits within `window` bp of pos1.

    Returns (clips, n_seen) where n_seen is the number of clips before
    sub-sampling to `max_clips` (reservoir sampling, so the sample is unbiased).
    """
    pos0 = pos1 - 1
    start = max(0, pos0 - window)
    end = pos0 + window + 1
    rng = random.Random(seed)
    kept, n_seen = [], 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        try:
            iterator = bam.fetch(chrom, start, end)
        except (ValueError, KeyError):
            sys.stderr.write(f"[WARN] contig {chrom} not in {bam_path}\n")
            return [], 0

        for read in iterator:
            if read.is_unmapped or read.is_secondary or read.is_qcfail:
                continue
            if read.is_duplicate and not keep_duplicates:
                continue
            if read.mapping_quality < mapq or read.cigartuples is None:
                continue
            query = read.query_sequence
            if not query:
                continue
            quals = read.query_qualities
            quals = list(quals) if quals is not None else [40] * len(query)

            cig = read.cigartuples
            candidates = []
            if cig[0][0] == BAM_CSOFT_CLIP and cig[0][1] >= min_clip:
                n = cig[0][1]
                candidates.append(("left", query[:n], quals[:n], read.reference_start))
            if cig[-1][0] == BAM_CSOFT_CLIP and cig[-1][1] >= min_clip:
                n = cig[-1][1]
                candidates.append(("right", query[-n:], quals[-n:], read.reference_end))

            for side, seq, qual, ref_pos in candidates:
                if ref_pos is None or abs(ref_pos - pos0) > window:
                    continue  # clip boundary is not the junction we are looking at
                n_seen += 1
                clip = Clip(
                    read.query_name,
                    side,
                    seq.upper(),
                    np.asarray(qual, dtype=np.int16),
                    read.mapping_quality,
                    bool(read.is_reverse),
                    ref_pos,
                )
                if len(kept) < max_clips:
                    kept.append(clip)
                else:  # reservoir sampling keeps the sample representative
                    j = rng.randrange(n_seen)
                    if j < max_clips:
                        kept[j] = clip
    return kept, n_seen


# --------------------------------------------------------------------------- #
# alignment of a clip against the insert
# --------------------------------------------------------------------------- #


class Alignment:
    __slots__ = ("q2r", "matches", "mismatches", "identity", "ref_start", "ref_end", "revcomp")

    def __init__(self, q2r, matches, mismatches, identity, ref_start, ref_end, is_rc):
        self.q2r = q2r  # query index -> insert index (0-based)
        self.matches = matches
        self.mismatches = mismatches
        self.identity = identity
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.revcomp = is_rc


def _ungapped_offsets(seq, insert_seq, kmer_index, k, max_offsets=3):
    """Vote for the best diagonal(s) (insert_index - query_index) using k-mers."""
    if len(seq) < k:
        return []
    votes = Counter()
    for i in range(len(seq) - k + 1):
        for hit in kmer_index.get(seq[i : i + k], ()):
            votes[hit - i] += 1
    if not votes:
        return []
    return [off for off, _ in votes.most_common(max_offsets)]


def _score_ungapped(seq, insert_seq, offset):
    """Score an ungapped placement of seq on insert_seq at the given offset."""
    L, n = len(insert_seq), len(seq)
    qs = max(0, -offset)
    qe = min(n, L - offset)
    if qe - qs <= 0:
        return None
    matches = mismatches = 0
    q2r = {}
    for qi in range(qs, qe):
        ri = qi + offset
        q2r[qi] = ri
        if seq[qi] == insert_seq[ri]:
            matches += 1
        else:
            mismatches += 1
    total = matches + mismatches
    return Alignment(q2r, matches, mismatches, matches / total, qs + offset, qe + offset, False)


def _align_edlib(seq, insert_seq):
    """Gapped (indel aware) infix alignment; returns query->insert mapping."""
    res = edlib.align(seq, insert_seq, mode="HW", task="path")
    if res is None or res.get("editDistance", -1) < 0 or not res.get("locations"):
        return None
    ref_start = res["locations"][0][0]
    qi, ri = 0, ref_start
    matches = mismatches = 0
    q2r = {}
    for length, op in _iter_cigar(res["cigar"]):
        if op in ("=", "X", "M"):
            for _ in range(length):
                q2r[qi] = ri
                if seq[qi] == insert_seq[ri]:
                    matches += 1
                else:
                    mismatches += 1
                qi += 1
                ri += 1
        elif op == "I":  # insertion in the clip relative to the insert
            qi += length
        elif op == "D":
            ri += length
    total = matches + mismatches
    if total == 0:
        return None
    return Alignment(q2r, matches, mismatches, matches / total, ref_start, ri, False)


def _iter_cigar(cigar):
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            yield int(num or 1), ch
            num = ""


def align_clip(seq, insert_seq, kmer_index, k, use_edlib=False, max_offsets=3):
    """Align a clip (both orientations) to the insert; returns the best Alignment."""
    best = None
    for is_rc, oriented in ((False, seq), (True, revcomp(seq))):
        cands = []
        if use_edlib and HAVE_EDLIB:
            aln = _align_edlib(oriented, insert_seq)
            if aln is not None:
                cands.append(aln)
        for off in _ungapped_offsets(oriented, insert_seq, kmer_index, k, max_offsets):
            aln = _score_ungapped(oriented, insert_seq, off)
            if aln is not None:
                cands.append(aln)
        for aln in cands:
            aln.revcomp = is_rc
            if best is None or (aln.matches, aln.identity) > (best.matches, best.identity):
                best = aln
    return best


# --------------------------------------------------------------------------- #
# pileup
# --------------------------------------------------------------------------- #

BASES = ("A", "C", "G", "T")
BASE_IDX = {b: i for i, b in enumerate(BASES)}

# cell codes used by the MSA heatmap
CELL_EMPTY, CELL_MATCH, CELL_MISMATCH, CELL_LOWQUAL = 0, 1, 2, 3


class Pileup:
    """Per-insert-position base counts, split by strand and base quality."""

    def __init__(self, length):
        self.length = length
        self.counts = np.zeros((length, 4), dtype=np.int32)  # high quality only
        self.counts_fwd = np.zeros((length, 4), dtype=np.int32)
        self.counts_rev = np.zeros((length, 4), dtype=np.int32)
        self.lowqual = np.zeros(length, dtype=np.int32)

    def add(self, ref_idx, base, is_reverse, high_quality):
        bi = BASE_IDX.get(base)
        if bi is None:
            return
        if not high_quality:
            self.lowqual[ref_idx] += 1
            return
        self.counts[ref_idx, bi] += 1
        if is_reverse:
            self.counts_rev[ref_idx, bi] += 1
        else:
            self.counts_fwd[ref_idx, bi] += 1

    @property
    def depth(self):
        return self.counts.sum(axis=1)

    def stats(self, insert_seq):
        """Return depth, non-ref count, non-ref fraction, top alt base per position."""
        depth = self.depth
        ref_idx = np.array([BASE_IDX.get(b, -1) for b in insert_seq])
        ref_cnt = np.zeros(self.length, dtype=np.int32)
        valid = ref_idx >= 0
        ref_cnt[valid] = self.counts[np.arange(self.length)[valid], ref_idx[valid]]
        nonref_cnt = depth - ref_cnt
        with np.errstate(invalid="ignore", divide="ignore"):
            nonref_frac = np.where(depth > 0, nonref_cnt / np.maximum(depth, 1), 0.0)
        alt = np.full(self.length, "-", dtype=object)
        alt_cnt = np.zeros(self.length, dtype=np.int32)
        for i in range(self.length):
            order = np.argsort(-self.counts[i])
            for bi in order:
                if bi != ref_idx[i] and self.counts[i, bi] > 0:
                    alt[i] = BASES[bi]
                    alt_cnt[i] = self.counts[i, bi]
                    break
        return depth, nonref_cnt, nonref_frac, alt, alt_cnt


def build_rows(aligned, insert_seq, min_bq):
    """Turn alignments into MSA rows (characters) and heat-map codes."""
    L = len(insert_seq)
    rows, codes = [], []
    pileup = Pileup(L)
    for rec in aligned:
        seq, qual, aln = rec["seq"], rec["qual"], rec["aln"]
        row = ["-"] * L
        code = np.zeros(L, dtype=np.int8)
        for qi, ri in aln.q2r.items():
            if not (0 <= qi < len(seq) and 0 <= ri < L):
                continue
            base = seq[qi]
            hq = qual[qi] >= min_bq
            pileup.add(ri, base, rec["is_reverse"], hq)
            row[ri] = base if hq else base.lower()
            if base == insert_seq[ri]:
                code[ri] = CELL_MATCH if hq else CELL_LOWQUAL
            else:
                code[ri] = CELL_MISMATCH if hq else CELL_LOWQUAL
        rows.append("".join(row))
        codes.append(code)
    return rows, codes, pileup


def call_variants(pileup, insert_seq, min_depth, min_alt_count, min_alt_frac, min_strand):
    """Return a list of candidate variant dicts for one site."""
    depth, nonref_cnt, nonref_frac, alt, alt_cnt = pileup.stats(insert_seq)
    out = []
    for i in range(pileup.length):
        if depth[i] < min_depth or alt_cnt[i] == 0:
            continue
        frac = alt_cnt[i] / depth[i]
        if alt_cnt[i] < min_alt_count or frac < min_alt_frac:
            continue
        bi = BASE_IDX[alt[i]]
        fwd = int(pileup.counts_fwd[i, bi])
        rev = int(pileup.counts_rev[i, bi])
        filters = []
        if min(fwd, rev) < min_strand:
            filters.append("strand_bias")
        out.append(
            {
                "insert_pos": i + 1,
                "insert_ref": insert_seq[i],
                "alt": alt[i],
                "depth": int(depth[i]),
                "alt_count": int(alt_cnt[i]),
                "alt_frac": round(float(frac), 4),
                "alt_fwd": fwd,
                "alt_rev": rev,
                "nonref_count": int(nonref_cnt[i]),
                "nonref_frac": round(float(nonref_frac[i]), 4),
                "lowqual_bases": int(pileup.lowqual[i]),
                "zygosity": "het" if 0.2 <= frac <= 0.8 else ("hom" if frac > 0.8 else "low"),
                "filter": ";".join(filters) if filters else "PASS",
            }
        )
    out.sort(key=lambda r: (-r["alt_frac"], -r["alt_count"]))
    return out


# --------------------------------------------------------------------------- #
# per-site worker
# --------------------------------------------------------------------------- #

_G = {}


def _worker_init(insert_seq, kmer_index, opts):
    _G["insert_seq"] = insert_seq
    _G["kmer_index"] = kmer_index
    _G["opts"] = opts


def process_site(site_obj):
    insert_seq = _G["insert_seq"]
    kmer_index = _G["kmer_index"]
    o = _G["opts"]
    t0 = time.time()

    clips, n_seen = extract_clips(
        o["bam"],
        site_obj["chrom"],
        site_obj["pos"],
        window=o["window"],
        mapq=o["mapq"],
        min_clip=o["min_clip"],
        max_clips=o["max_clips"],
        keep_duplicates=o["keep_duplicates"],
        seed=o["seed"],
    )

    aligned, n_rejected = [], 0
    for clip in clips:
        aln = align_clip(
            clip.seq,
            insert_seq,
            kmer_index,
            o["kmer_k"],
            use_edlib=o["use_edlib"],
            max_offsets=o["max_offsets"],
        )
        if aln is None or aln.matches < o["min_match"] or aln.identity < o["min_identity"]:
            n_rejected += 1
            continue
        seq = revcomp(clip.seq) if aln.revcomp else clip.seq
        qual = clip.qual[::-1] if aln.revcomp else clip.qual
        aligned.append(
            {
                "read_name": clip.read_name,
                "side": clip.side,
                "mapq": clip.mapq,
                "is_reverse": clip.is_reverse,
                "revcomp_used": aln.revcomp,
                "seq": seq,
                "qual": qual,
                "aln": aln,
            }
        )

    # sort like a pileup: by first aligned insert position, then length
    aligned.sort(key=lambda r: (r["aln"].ref_start, -(r["aln"].matches)))
    rows, codes, pileup = build_rows(aligned, insert_seq, o["min_bq"])
    depth, nonref_cnt, nonref_frac, alt, alt_cnt = pileup.stats(insert_seq)
    variants = call_variants(
        pileup,
        insert_seq,
        o["min_depth"],
        o["min_alt_count"],
        o["min_alt_frac"],
        o["min_strand_count"],
    )

    meta = [
        {
            "read_name": r["read_name"],
            "side": r["side"],
            "mapq": r["mapq"],
            "strand": "-" if r["is_reverse"] else "+",
            "revcomp_used": r["revcomp_used"],
            "insert_start": r["aln"].ref_start + 1,
            "insert_end": r["aln"].ref_end,
            "matches": r["aln"].matches,
            "mismatches": r["aln"].mismatches,
            "identity": round(r["aln"].identity, 4),
        }
        for r in aligned
    ]

    return {
        "site_obj": site_obj,
        "n_clips_seen": n_seen,
        "n_clips_used": len(clips),
        "n_aligned": len(aligned),
        "n_rejected": n_rejected,
        "rows": rows,
        "codes": np.array(codes, dtype=np.int8) if codes else np.zeros((0, len(insert_seq)), np.int8),
        "meta": meta,
        "depth": depth,
        "nonref_cnt": nonref_cnt,
        "nonref_frac": nonref_frac,
        "alt": alt,
        "alt_cnt": alt_cnt,
        "counts": pileup.counts,
        "counts_fwd": pileup.counts_fwd,
        "counts_rev": pileup.counts_rev,
        "lowqual": pileup.lowqual,
        "variants": variants,
        "elapsed": time.time() - t0,
    }


# --------------------------------------------------------------------------- #
# text pileup
# --------------------------------------------------------------------------- #


def write_pileup_text(path, title, insert_seq, result, meta, wrap, min_depth, highlight_frac,
                      annotate, max_rows):
    """Write a wrapped, human readable pileup: '.' = matches the insert."""
    rows = result["rows"][:max_rows]
    labels = [
        f"{m['read_name'][:28]:<28} {m['strand']}{m['side'][0]}"
        for m in meta[:max_rows]
    ]
    lab_w = max([len(l) for l in labels], default=10)
    lab_w = max(lab_w, 12)
    depth, nonref_frac = result["depth"], result["nonref_frac"]
    L = len(insert_seq)

    consensus = []
    for i in range(L):
        if depth[i] == 0:
            consensus.append("n")
            continue
        bi = int(np.argmax(result["counts"][i]))
        cons = BASES[bi]
        consensus.append(cons if cons != insert_seq[i] else ".")

    marks = []
    for i in range(L):
        if depth[i] >= min_depth and nonref_frac[i] >= highlight_frac:
            marks.append("*")
        elif (i + 1) in annotate:
            marks.append("|")
        else:
            marks.append(" ")

    with open(path, "w") as fh:
        fh.write(f"# {title}\n")
        fh.write(
            f"# clips seen={result['n_clips_seen']} used={result['n_clips_used']} "
            f"aligned={result['n_aligned']} rejected={result['n_rejected']}\n"
        )
        fh.write(f"# rows shown: {len(rows)} of {len(result['rows'])}\n")
        fh.write("# '.' in a read row = base identical to the insert; lower case = low base quality\n")
        fh.write(f"# '*' marks positions with depth >= {min_depth} and non-ref fraction >= {highlight_frac}\n")
        if result["variants"]:
            fh.write("# candidate variants: ")
            fh.write(
                ", ".join(
                    f"{v['insert_pos']}{v['insert_ref']}>{v['alt']}"
                    f"({v['alt_count']}/{v['depth']},{v['alt_frac']:.2f},{v['filter']})"
                    for v in result["variants"]
                )
            )
            fh.write("\n")
        else:
            fh.write("# candidate variants: none\n")
        fh.write("\n")

        for start in range(0, L, wrap):
            end = min(L, start + wrap)
            if depth[start:end].max(initial=0) == 0 and not any(start < p <= end for p in annotate):
                continue  # nothing aligned in this block
            ruler = _ruler(start, end)
            fh.write(f"{'':<{lab_w}}  {ruler}\n")
            fh.write(f"{'INSERT':<{lab_w}}  {insert_seq[start:end]}\n")
            fh.write(f"{'CONSENSUS':<{lab_w}}  {''.join(consensus[start:end])}\n")
            fh.write(f"{'VARIANT':<{lab_w}}  {''.join(marks[start:end])}\n")
            for label, row in zip(labels, rows):
                chunk = row[start:end]
                if chunk.strip("-") == "":
                    continue  # this read does not cover the block
                pretty = "".join(
                    "." if c.upper() == insert_seq[start + j] and c != "-" else c
                    for j, c in enumerate(chunk)
                )
                fh.write(f"{label:<{lab_w}}  {pretty}\n")
            fh.write("\n")


def _ruler(start, end):
    out = []
    i = start
    while i < end:
        pos = i + 1
        if pos % 10 == 0:
            label = str(pos)
            if i + len(label) <= end:
                out.append(label)
                i += len(label)
                continue
        out.append(".")
        i += 1
    return "".join(out)[: end - start]


# --------------------------------------------------------------------------- #
# plots
# --------------------------------------------------------------------------- #

MSA_CMAP = ListedColormap(["#f5f5f5", "#dddddd", "#e41a1c", "#fdd0a2"])


def plot_site(result, insert_seq, title, annotate, max_plot_rows, min_depth, highlight_frac):
    codes = result["codes"][:max_plot_rows]
    depth = result["depth"]
    nonref_frac = result["nonref_frac"]
    L = len(insert_seq)

    fig = plt.figure(figsize=(15, max(6.5, 3.2 + len(codes) * 0.07)))
    gs = fig.add_gridspec(2, 1, height_ratios=[max(2.4, len(codes) * 0.06), 2.0], hspace=0.28)

    ax1 = fig.add_subplot(gs[0, 0])
    if len(codes):
        ax1.imshow(
            codes,
            aspect="auto",
            interpolation="nearest",
            cmap=MSA_CMAP,
            vmin=-0.5,
            vmax=3.5,
            extent=[0.5, L + 0.5, len(codes) + 0.5, 0.5],
        )
        ax1.set_ylabel(f"Aligned clips (n={len(result['rows'])})")
    else:
        ax1.text(0.02, 0.5, "No aligned soft-clips", transform=ax1.transAxes)
        ax1.set_yticks([])
    ax1.set_xlim(0.5, L + 0.5)
    ax1.set_title(title)
    handles = [
        plt.Line2D([], [], marker="s", ls="", color=MSA_CMAP(1), label="matches insert"),
        plt.Line2D([], [], marker="s", ls="", color=MSA_CMAP(2), label="mismatch"),
        plt.Line2D([], [], marker="s", ls="", color=MSA_CMAP(3), label="low base quality"),
    ]
    ax1.legend(handles=handles, loc="upper right", frameon=False, fontsize=8, ncol=3)

    ax2 = fig.add_subplot(gs[1, 0])
    x = np.arange(1, L + 1)
    ax2.fill_between(x, depth, color="0.75", step="mid", label="Depth")
    ax2.set_xlim(0.5, L + 0.5)
    ax2.set_ylabel("Depth")
    ax2.set_xlabel("Insert position (1-based)")
    ax2.axhline(min_depth, color="0.4", lw=0.8, ls=":")

    ax3 = ax2.twinx()
    ax3.plot(x, nonref_frac, color="crimson", lw=1.0, label="Non-ref fraction")
    ax3.set_ylim(0, 1.02)
    ax3.set_ylabel("Non-ref fraction")
    ax3.axhline(highlight_frac, color="crimson", lw=0.8, ls=":")

    for v in result["variants"]:
        for ax in (ax1, ax2):
            ax.axvline(v["insert_pos"], color="black", ls="--", lw=0.9)
        ax3.annotate(
            f"{v['insert_pos']}{v['insert_ref']}>{v['alt']} ({v['alt_frac']:.2f})",
            xy=(v["insert_pos"], min(1.0, v["alt_frac"])),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            fontsize=8,
        )
    for p in annotate:
        if 1 <= p <= L:
            for ax in (ax1, ax2):
                ax.axvline(p, color="tab:blue", ls=":", lw=0.9)

    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax3.get_legend_handles_labels()
    ax2.legend(h1 + h2, l1 + l2, loc="upper right", frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.07, right=0.93, top=0.94, bottom=0.07)
    return fig


def plot_cross_site(results, insert_seq, min_depth, highlight_frac):
    """Heatmap sites x insert positions of the non-ref fraction (masked on depth)."""
    L = len(insert_seq)
    mat = np.full((len(results), L), np.nan)
    labels = []
    for i, r in enumerate(results):
        frac = r["nonref_frac"].astype(float).copy()
        frac[r["depth"] < min_depth] = np.nan
        mat[i] = frac
        s = r["site_obj"]
        labels.append(f"site {s['site']} ({s['gene']})")

    fig, ax = plt.subplots(figsize=(15, max(4.0, 0.32 * len(results) + 2.2)))
    cmap = plt.get_cmap("Reds").copy()
    cmap.set_bad("#eeeeee")
    im = ax.imshow(
        np.ma.masked_invalid(mat),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=1,
        extent=[0.5, L + 0.5, len(results) - 0.5, -0.5],
    )
    ax.set_yticks(range(len(results)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Insert position (1-based)")
    ax.set_title(
        f"Non-ref fraction per insert position and site (grey = depth < {min_depth}); "
        f"a column that is red for one site only is a site-specific variant"
    )
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="Non-ref fraction")
    fig.subplots_adjust(left=0.16, right=0.98, top=0.90, bottom=0.10)
    return fig


def plot_variant_matrix(matrix_df, highlight_frac):
    fig, ax = plt.subplots(figsize=(min(18, 3 + 0.6 * matrix_df.shape[1]),
                                    max(3.5, 0.32 * matrix_df.shape[0] + 2)))
    cmap = plt.get_cmap("Reds").copy()
    cmap.set_bad("#eeeeee")
    data = np.ma.masked_invalid(matrix_df.to_numpy(dtype=float))
    im = ax.imshow(data, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(matrix_df.shape[1]))
    ax.set_xticklabels(matrix_df.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(matrix_df.shape[0]))
    ax.set_yticklabels(matrix_df.index, fontsize=8)
    for i in range(matrix_df.shape[0]):
        for j in range(matrix_df.shape[1]):
            val = matrix_df.iat[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if val >= 0.6 else "black")
    ax.set_title("Alt fraction per candidate variant position (grey = insufficient depth)")
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="Alt fraction")
    fig.subplots_adjust(left=0.22, right=0.98, top=0.90, bottom=0.18)
    return fig


def text_page(text, figsize=(11.69, 8.27), fontsize=11):
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")
    ax.text(0.02, 0.98, text, va="top", ha="left", fontsize=fontsize, family="monospace")
    return fig


def table_page(df, title, max_rows=40):
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.axis("off")
    show = df.head(max_rows)
    tbl = ax.table(cellText=show.astype(str).values, colLabels=list(show.columns), loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7)
    tbl.scale(1, 1.2)
    ax.set_title(title if len(df) <= max_rows else f"{title} (first {max_rows} of {len(df)} rows)")
    return fig


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Pileup / MSA of soft-clipped read ends against a transgene insert, per TLA site.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--bam", required=True, help="Indexed BAM aligned to the host genome")
    ap.add_argument("--insert-fasta", default="insert_sequence.fasta", help="Insert/cassette FASTA")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sites", default=None, help="TSV/CSV with columns site,chrom,pos[,gene]")
    ap.add_argument("--focus-site", default=None, help="Only process this site id")
    ap.add_argument("--annotate-positions", default="", help="Comma separated insert positions to mark")

    g = ap.add_argument_group("clip selection")
    g.add_argument("--window", type=int, default=20, help="Max distance clip boundary - site (bp)")
    g.add_argument("--mapq", type=int, default=20)
    g.add_argument("--min-clip", type=int, default=15, help="Minimum soft-clip length")
    g.add_argument("--max-clips", type=int, default=300, help="Max clips per site (reservoir sampled)")
    g.add_argument("--keep-duplicates", action="store_true")
    g.add_argument("--seed", type=int, default=0, help="Seed for clip sub-sampling")

    g = ap.add_argument_group("alignment")
    g.add_argument("--kmer-k", type=int, default=12, help="Seed k-mer size")
    g.add_argument("--kmer-max-hits", type=int, default=50, help="Drop k-mers occurring more often")
    g.add_argument("--max-offsets", type=int, default=3, help="Candidate diagonals evaluated per clip")
    g.add_argument("--min-identity", type=float, default=0.85, help="Min identity of clip vs insert")
    g.add_argument("--min-match", type=int, default=20, help="Min matching bases of clip vs insert")
    g.add_argument("--no-edlib", action="store_true", help="Do not use edlib even if installed")

    g = ap.add_argument_group("variant calling")
    g.add_argument("--min-bq", type=int, default=20, help="Min base quality to count a base")
    g.add_argument("--min-depth", type=int, default=10)
    g.add_argument("--min-alt-count", type=int, default=5)
    g.add_argument("--min-alt-frac", type=float, default=0.2)
    g.add_argument("--min-strand-count", type=int, default=1, help="Min alt reads per strand")

    g = ap.add_argument_group("output")
    g.add_argument("--max-plot-rows", type=int, default=150)
    g.add_argument("--max-pileup-rows", type=int, default=500)
    g.add_argument("--wrap", type=int, default=100, help="Columns per block in the text pileup")
    g.add_argument("--emit-counts", action="store_true", help="Write full per-position base counts")
    g.add_argument("--workers", type=int, default=1)
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    t_all = time.time()

    annotate = parse_positions(args.annotate_positions)
    os.makedirs(args.outdir, exist_ok=True)
    plots_dir = os.path.join(args.outdir, "plots")
    msa_dir = os.path.join(args.outdir, "msa_txt")
    pileup_dir = os.path.join(args.outdir, "pileup")
    for d in (plots_dir, msa_dir, pileup_dir):
        os.makedirs(d, exist_ok=True)

    insert_name, insert_seq = read_first_fasta_seq(args.insert_fasta)
    if not insert_seq:
        raise SystemExit(f"[ERROR] no sequence found in {args.insert_fasta}")
    if args.kmer_k > len(insert_seq):
        raise SystemExit("[ERROR] --kmer-k is larger than the insert sequence")
    kmer_index = build_kmer_index(insert_seq, args.kmer_k, args.kmer_max_hits)

    sites = load_sites(args.sites)
    if args.focus_site is not None:
        sites = [s for s in sites if s["site"] == str(args.focus_site)]
    if not sites:
        raise SystemExit("[ERROR] no sites selected")
    sites.sort(key=lambda s: site_sort_key(s["site"]))

    use_edlib = HAVE_EDLIB and not args.no_edlib
    opts = {
        "bam": args.bam,
        "window": args.window,
        "mapq": args.mapq,
        "min_clip": args.min_clip,
        "max_clips": args.max_clips,
        "keep_duplicates": args.keep_duplicates,
        "seed": args.seed,
        "kmer_k": args.kmer_k,
        "max_offsets": args.max_offsets,
        "min_identity": args.min_identity,
        "min_match": args.min_match,
        "use_edlib": use_edlib,
        "min_bq": args.min_bq,
        "min_depth": args.min_depth,
        "min_alt_count": args.min_alt_count,
        "min_alt_frac": args.min_alt_frac,
        "min_strand_count": args.min_strand_count,
    }

    print(f"[INFO] insert '{insert_name}' length={len(insert_seq)} from {args.insert_fasta}")
    print(f"[INFO] sites={len(sites)} workers={args.workers} gapped_alignment={'edlib' if use_edlib else 'ungapped'}")

    results = []
    if args.workers > 1:
        from concurrent.futures import ProcessPoolExecutor, as_completed

        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=_worker_init,
            initargs=(insert_seq, kmer_index, opts),
        ) as ex:
            futs = {ex.submit(process_site, s): s for s in sites}
            for i, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                _log_site(i, len(sites), r)
                results.append(r)
    else:
        _worker_init(insert_seq, kmer_index, opts)
        for i, s in enumerate(sites, 1):
            r = process_site(s)
            _log_site(i, len(sites), r)
            results.append(r)

    results.sort(key=lambda r: site_sort_key(r["site_obj"]["site"]))

    summary_rows, variant_rows, count_rows = [], [], []
    for r in results:
        s = r["site_obj"]
        passing = [v for v in r["variants"] if v["filter"] == "PASS"]
        summary_rows.append(
            {
                "site": s["site"],
                "gene": s["gene"],
                "chrom": s["chrom"],
                "pos": s["pos"],
                "clips_seen": r["n_clips_seen"],
                "clips_used": r["n_clips_used"],
                "clips_aligned": r["n_aligned"],
                "clips_rejected": r["n_rejected"],
                "insert_bp_covered": int((r["depth"] > 0).sum()),
                "mean_depth": round(float(r["depth"].mean()), 2),
                "max_depth": int(r["depth"].max()) if len(r["depth"]) else 0,
                "n_variants_pass": len(passing),
                "variants": ";".join(
                    f"{v['insert_pos']}{v['insert_ref']}>{v['alt']}"
                    f"({v['alt_count']}/{v['depth']},{v['alt_frac']:.2f})"
                    for v in passing
                ),
                "elapsed_sec": round(r["elapsed"], 2),
            }
        )
        for v in r["variants"]:
            variant_rows.append(
                {"site": s["site"], "gene": s["gene"], "chrom": s["chrom"], "pos": s["pos"], **v}
            )
        if args.emit_counts:
            for i in range(len(insert_seq)):
                if r["depth"][i] == 0:
                    continue
                count_rows.append(
                    {
                        "site": s["site"],
                        "insert_pos": i + 1,
                        "insert_ref": insert_seq[i],
                        "depth": int(r["depth"][i]),
                        "A": int(r["counts"][i, 0]),
                        "C": int(r["counts"][i, 1]),
                        "G": int(r["counts"][i, 2]),
                        "T": int(r["counts"][i, 3]),
                        "lowqual": int(r["lowqual"][i]),
                        "nonref_frac": round(float(r["nonref_frac"][i]), 4),
                    }
                )

    summary_df = pd.DataFrame(summary_rows)
    variant_df = pd.DataFrame(variant_rows)
    matrix_df = build_variant_matrix(results, args.min_depth)

    # ------------------------------------------------------------------ files
    for r in results:
        s = r["site_obj"]
        tag = f"site{s['site']}_{s['chrom']}_{s['pos']}"
        title = (
            f"Site {s['site']} | {s['gene']} | {s['chrom']}:{s['pos']} | "
            f"aligned {r['n_aligned']}/{r['n_clips_used']} clips"
        )
        write_pileup_text(
            os.path.join(pileup_dir, f"{tag}.pileup.txt"),
            title,
            insert_seq,
            r,
            r["meta"],
            args.wrap,
            args.min_depth,
            args.min_alt_frac,
            annotate,
            args.max_pileup_rows,
        )
        with open(os.path.join(msa_dir, f"{tag}.msa.tsv"), "w") as fh:
            fh.write(f"# {title}\n")
            cols = ["read_name", "side", "strand", "mapq", "revcomp_used", "insert_start",
                    "insert_end", "matches", "mismatches", "identity"]
            fh.write("\t".join(cols + ["aligned_row"]) + "\n")
            for m, row in zip(r["meta"], r["rows"]):
                fh.write("\t".join(str(m[c]) for c in cols) + f"\t{row}\n")

    summary_tsv = os.path.join(args.outdir, "site_summary.tsv")
    summary_df.to_csv(summary_tsv, sep="\t", index=False)

    variant_tsv = os.path.join(args.outdir, "insert_variants.tsv")
    if variant_df.empty:
        variant_df = pd.DataFrame(
            columns=["site", "gene", "chrom", "pos", "insert_pos", "insert_ref", "alt", "depth",
                     "alt_count", "alt_frac", "alt_fwd", "alt_rev", "nonref_count", "nonref_frac",
                     "lowqual_bases", "zygosity", "filter"]
        )
    variant_df.to_csv(variant_tsv, sep="\t", index=False)

    matrix_tsv = os.path.join(args.outdir, "variant_matrix.tsv")
    matrix_df.to_csv(matrix_tsv, sep="\t")

    counts_tsv = None
    if args.emit_counts:
        counts_tsv = os.path.join(args.outdir, "pileup_counts.tsv")
        pd.DataFrame(count_rows).to_csv(counts_tsv, sep="\t", index=False)

    site_specific = site_specific_variants(matrix_df, variant_df, args.min_alt_frac)

    # -------------------------------------------------------------------- pdf
    pdf_path = os.path.join(args.outdir, "tla_insert_softclip_report.pdf")
    with PdfPages(pdf_path) as pdf:
        cover = (
            "TLA insert soft-clip pileup report\n\n"
            f"BAM            : {args.bam}\n"
            f"Insert         : {args.insert_fasta} ({insert_name}, {len(insert_seq)} bp)\n"
            f"Sites          : {len(sites)}\n"
            f"Clip filters   : window={args.window} mapq>={args.mapq} clip>={args.min_clip} "
            f"max_clips={args.max_clips}\n"
            f"Alignment      : k={args.kmer_k} min_identity={args.min_identity} "
            f"min_match={args.min_match} gapped={'edlib' if use_edlib else 'no'}\n"
            f"Variant filters: bq>={args.min_bq} depth>={args.min_depth} "
            f"alt>={args.min_alt_count} frac>={args.min_alt_frac} "
            f"strand>={args.min_strand_count}\n\n"
            + _conclusion_text(site_specific, matrix_df)
        )
        pdf.savefig(text_page(cover))
        plt.close("all")

        pdf.savefig(plot_cross_site(results, insert_seq, args.min_depth, args.min_alt_frac))
        plt.close("all")

        if not matrix_df.empty and matrix_df.shape[1] > 0:
            pdf.savefig(plot_variant_matrix(matrix_df, args.min_alt_frac))
            plt.close("all")

        if not summary_df.empty:
            show = summary_df.drop(columns=["chrom"]).copy()
            pdf.savefig(table_page(show, "Summary per site"))
            plt.close("all")

        for r in results:
            s = r["site_obj"]
            tag = f"site{s['site']}_{s['chrom']}_{s['pos']}"
            title = (
                f"Site {s['site']} | {s['gene']} | {s['chrom']}:{s['pos']} | "
                f"aligned {r['n_aligned']}/{r['n_clips_used']} clips"
            )
            fig = plot_site(r, insert_seq, title, annotate, args.max_plot_rows,
                            args.min_depth, args.min_alt_frac)
            fig.savefig(os.path.join(plots_dir, f"{tag}.png"), dpi=150)
            pdf.savefig(fig)
            plt.close(fig)

    with open(os.path.join(args.outdir, "run_parameters.json"), "w") as fh:
        json.dump(vars(args), fh, indent=2, sort_keys=True)

    print(f"[OK] PDF            : {pdf_path}")
    print(f"[OK] Pileups        : {pileup_dir}")
    print(f"[OK] Summary        : {summary_tsv}")
    print(f"[OK] Variants       : {variant_tsv}")
    print(f"[OK] Variant matrix : {matrix_tsv}")
    if counts_tsv:
        print(f"[OK] Base counts    : {counts_tsv}")
    print(_conclusion_text(site_specific, matrix_df))
    print(f"[DONE] total runtime {time.time() - t_all:.2f}s")
    return 0


def _log_site(i, n, r):
    s = r["site_obj"]
    v = ",".join(
        f"{x['insert_pos']}{x['insert_ref']}>{x['alt']}({x['alt_frac']:.2f})"
        for x in r["variants"]
        if x["filter"] == "PASS"
    )
    print(
        f"[{i}/{n}] site {s['site']:>3} {s['chrom']}:{s['pos']} "
        f"clips={r['n_clips_used']} aligned={r['n_aligned']} "
        f"max_depth={int(r['depth'].max()) if len(r['depth']) else 0} "
        f"variants={v or '-'} ({r['elapsed']:.2f}s)"
    )


def build_variant_matrix(results, min_depth):
    """sites x candidate positions matrix of the alt fraction (NaN if low depth)."""
    positions = {}
    for r in results:
        for v in r["variants"]:
            if v["filter"] == "PASS":
                positions.setdefault(v["insert_pos"], v["insert_ref"])
    if not positions:
        return pd.DataFrame()
    cols, index, data = sorted(positions), [], []
    for r in results:
        s = r["site_obj"]
        index.append(f"site {s['site']} ({s['gene']})")
        row = []
        for p in cols:
            i = p - 1
            if r["depth"][i] < min_depth:
                row.append(np.nan)
                continue
            alt_base = positions[p]
            # fraction of reads carrying any non-reference base at this position
            row.append(float(r["nonref_frac"][i]))
        data.append(row)
    labels = [f"{p}{positions[p]}" for p in cols]
    return pd.DataFrame(data, index=index, columns=labels)


def site_specific_variants(matrix_df, variant_df, min_alt_frac):
    """Positions where a variant is seen in exactly one site (with enough depth)."""
    out = []
    if matrix_df.empty:
        return out
    for col in matrix_df.columns:
        vals = matrix_df[col]
        covered = vals.dropna()
        positive = covered[covered >= min_alt_frac]
        if len(positive) == 1 and len(covered) > 1:
            out.append(
                {
                    "position": col,
                    "site": positive.index[0],
                    "alt_frac": float(positive.iloc[0]),
                    "n_sites_with_depth": int(len(covered)),
                }
            )
    return out


def _conclusion_text(site_specific, matrix_df):
    if matrix_df.empty:
        return "No candidate variants passed the filters in any site."
    lines = ["Site-specific variant positions (present in exactly one site):"]
    if site_specific:
        for v in site_specific:
            lines.append(
                f"  insert position {v['position']} -> {v['site']} "
                f"(alt fraction {v['alt_frac']:.2f}; {v['n_sites_with_depth']} sites had enough depth)"
            )
    else:
        lines.append("  none - every candidate variant is shared by several sites")
    shared = [c for c in matrix_df.columns if c not in {v["position"] for v in site_specific}]
    if shared:
        lines.append("Shared / recurrent positions: " + ", ".join(shared))
    return "\n".join(lines)


if __name__ == "__main__":
    sys.exit(main())
