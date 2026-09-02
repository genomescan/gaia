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
  pileup_index.html                        interactive overview, links to the pileups
  pileup_html/site<N>_...html              scrollable, searchable pileup viewer
  pileup/site<N>_<chrom>_<pos>.pileup.txt  human readable MSA / pileup (wrapped)
  msa_txt/site<N>_..._msa.tsv              one aligned row per clip (machine readable)
  plots/site<N>_...png                     MSA heatmap + depth / non-ref fraction
  tla_insert_softclip_report.pdf           all plots + cross-site overview
  site_summary.tsv                         per-site QC summary
  insert_variants.tsv                      every candidate variant position per site
  variant_classification.tsv               site-specific vs shared vs reference mismatch
  variant_matrix.tsv                       sites x candidate variants, alt fraction
  pileup_consensus.fasta                   consensus of all clips over all sites
  insert_vs_consensus.tsv                  positions where the insert FASTA and the reads differ
  pileup_counts.tsv                        (optional) full per-position base counts

Soft clips contain cassette sequence
------------------------------------
A clip of a junction-spanning read starts (or ends) in the cassette that flanks
the insert, so part of it does not belong to the insert at all. Two mechanisms
keep those bases out of the insert pileup:

* badly matching alignment ends are trimmed (--trim-window / --trim-identity);
* with --cassette-fasta the whole construct is used as alignment reference, the
  insert is located inside it and all positions are reported both in construct
  and in insert coordinates (variants are only called inside the insert unless
  --call-in-cassette is given).

Interpreting the variant table
------------------------------
Positions that are non-reference in *every* site, usually at ~100%, are not
integration-specific SNPs: they are differences between the supplied insert
FASTA and the construct that was actually integrated. They are labelled
`insert_reference_mismatch` in variant_classification.tsv and summarised in
insert_vs_consensus.tsv; a genuine copy-specific SNP shows up as
`site_specific`, i.e. present in one site while the other sites are covered but
reference at that position.

Main improvements over the first version
----------------------------------------
* clips are only kept when the clip *boundary* is within the window of the site
  (instead of any read overlapping the window), which removes off-target clips;
* reads are properly filtered (secondary / duplicate / QC-fail, mapping quality)
  and clips are sub-sampled by reservoir sampling instead of stopping at the
  first N reads (which biased the pileup towards the left of the window);
* alignment is seed-anchored: k-mer diagonal voting gives an ungapped placement
  in O(len(clip)); the O(n*m) pure-Python Smith-Waterman is gone. If `edlib` is
  installed it is used for gapped (indel-aware) alignment; alignment ends that do
  not match (cassette flanks, chimeric clips) are trimmed away;
* base qualities are used: low quality bases are shown in lower case and are not
  counted for variant calling, which is essential to trust a single SNP;
* per-position statistics keep strand information, so strand-biased artefacts can
  be filtered out;
* variant calling with depth / allele fraction / allele count / strand filters,
  and a cross-site comparison of the *same* alt base that separates site-specific
  SNPs from reference discrepancies;
* an interactive HTML pileup so long inserts can be scrolled instead of being
  squeezed into one figure; matplotlib labels are staggered and capped;
* sites can be supplied in a TSV file instead of being hard-coded;
* parallel execution shares the reference and k-mer index through a worker
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


def read_fasta_records(path):
    """Return [(name, sequence), ...] for every record in a FASTA file."""
    records, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq).upper()))
                name = line[1:].strip().split()[0] if len(line) > 1 else f"seq{len(records)}"
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq).upper()))
    return records


def locate_insert(insert_seq, construct):
    """Find the insert inside the construct; returns (offset, identity) or None."""
    idx = construct.find(insert_seq)
    if idx >= 0:
        return idx, 1.0
    # tolerate a few differences: vote for the best diagonal with k-mers
    k = min(25, max(12, len(insert_seq) // 20))
    index = build_kmer_index(construct, k, max_hits=20)
    votes = Counter()
    for i in range(0, len(insert_seq) - k + 1):
        for hit in index.get(insert_seq[i : i + k], ()):
            votes[hit - i] += 1
    if not votes:
        return None
    offset, _ = votes.most_common(1)[0]
    if offset < 0 or offset + len(insert_seq) > len(construct):
        return None
    window = construct[offset : offset + len(insert_seq)]
    identity = sum(a == b for a, b in zip(insert_seq, window)) / len(insert_seq)
    return (offset, identity) if identity >= 0.9 else None


def build_reference(insert_seq, cassette_fasta):
    """Alignment reference: the insert, optionally embedded in the cassette/construct.

    Returns (name, reference_sequence, insert_offset).
    """
    if not cassette_fasta:
        return "insert", insert_seq, 0
    best = None
    for name, seq in read_fasta_records(cassette_fasta):
        for oriented, label in ((seq, name), (revcomp(seq), name + "(revcomp)")):
            hit = locate_insert(insert_seq, oriented)
            if hit and (best is None or hit[1] > best[3]):
                best = (label, oriented, hit[0], hit[1])
    if best is None:
        raise SystemExit(
            f"[ERROR] the insert sequence was not found inside any record of {cassette_fasta}. "
            "The cassette FASTA must contain the construct *including* the insert."
        )
    name, seq, offset, identity = best
    print(
        f"[INFO] insert found in cassette record '{name}' at offset {offset} "
        f"(identity {identity:.4f}); reference length {len(seq)} bp"
    )
    return name, seq, offset


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


def trim_alignment(seq, ref_seq, aln, window=15, min_identity=0.7):
    """Trim badly matching ends of an alignment.

    Soft clips of reads spanning an insertion junction contain vector/cassette
    sequence that flanks the insert. Such bases must not be forced into insert
    coordinates, otherwise they show up as a wall of 100% "variants". Ends are
    trimmed while the identity in a sliding window stays below `min_identity`.
    """
    if not aln.q2r:
        return aln
    qs = sorted(aln.q2r)
    match_flags = [1 if seq[q] == ref_seq[aln.q2r[q]] else 0 for q in qs]
    n = len(qs)
    if n <= window:
        return aln

    start = 0
    while start + window <= n and sum(match_flags[start : start + window]) / window < min_identity:
        start += 1
    if start + window > n:  # nothing matches well enough
        return Alignment({}, 0, 0, 0.0, 0, 0, aln.revcomp)
    end = n
    while end - window >= start and sum(match_flags[end - window : end]) / window < min_identity:
        end -= 1
    # walk inwards over the last mismatching bases of each end
    while start < end and not match_flags[start]:
        start += 1
    while end > start and not match_flags[end - 1]:
        end -= 1
    if start == 0 and end == n:
        return aln
    if end <= start:
        return Alignment({}, 0, 0, 0.0, 0, 0, aln.revcomp)

    kept = qs[start:end]
    q2r = {q: aln.q2r[q] for q in kept}
    matches = sum(match_flags[start:end])
    total = end - start
    ref_positions = [q2r[q] for q in kept]
    return Alignment(
        q2r,
        matches,
        total - matches,
        matches / total,
        min(ref_positions),
        max(ref_positions) + 1,
        aln.revcomp,
    )


def align_clip(seq, ref_seq, kmer_index, k, use_edlib=False, max_offsets=3,
               trim_window=15, trim_identity=0.7):
    """Align a clip (both orientations) to the reference; returns the best Alignment."""
    best = None
    for is_rc, oriented in ((False, seq), (True, revcomp(seq))):
        cands = []
        if use_edlib and HAVE_EDLIB:
            aln = _align_edlib(oriented, ref_seq)
            if aln is not None:
                cands.append(aln)
        for off in _ungapped_offsets(oriented, ref_seq, kmer_index, k, max_offsets):
            aln = _score_ungapped(oriented, ref_seq, off)
            if aln is not None:
                cands.append(aln)
        for aln in cands:
            aln = trim_alignment(oriented, ref_seq, aln, trim_window, trim_identity)
            if not aln.q2r:
                continue
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


def call_variants(pileup, ref_seq, min_depth, min_alt_count, min_alt_frac, min_strand,
                  insert_offset=0, insert_len=None, include_flanks=False):
    """Return a list of candidate variant dicts for one site.

    Positions are reported both in reference (= construct) coordinates and in
    insert coordinates; positions outside the insert are labelled as cassette.
    """
    insert_len = pileup.length if insert_len is None else insert_len
    depth, nonref_cnt, nonref_frac, alt, alt_cnt = pileup.stats(ref_seq)
    out = []
    for i in range(pileup.length):
        region = region_label(i, insert_offset, insert_len)
        if region != "insert" and not include_flanks:
            continue
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
                "ref_pos": i + 1,
                "region": region,
                "insert_pos": i + 1 - insert_offset,
                "insert_ref": ref_seq[i],
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


def region_label(ref_idx, insert_offset, insert_len):
    if ref_idx < insert_offset:
        return "cassette_5p"
    if ref_idx >= insert_offset + insert_len:
        return "cassette_3p"
    return "insert"


# --------------------------------------------------------------------------- #
# per-site worker
# --------------------------------------------------------------------------- #

_G = {}


def _worker_init(ref_seq, kmer_index, opts):
    _G["ref_seq"] = ref_seq
    _G["kmer_index"] = kmer_index
    _G["opts"] = opts


def process_site(site_obj):
    ref_seq = _G["ref_seq"]
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
            ref_seq,
            kmer_index,
            o["kmer_k"],
            use_edlib=o["use_edlib"],
            max_offsets=o["max_offsets"],
            trim_window=o["trim_window"],
            trim_identity=o["trim_identity"],
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

    # sort like a pileup: by first aligned reference position, then length
    aligned.sort(key=lambda r: (r["aln"].ref_start, -(r["aln"].matches)))
    rows, codes, pileup = build_rows(aligned, ref_seq, o["min_bq"])
    depth, nonref_cnt, nonref_frac, alt, alt_cnt = pileup.stats(ref_seq)
    variants = call_variants(
        pileup,
        ref_seq,
        o["min_depth"],
        o["min_alt_count"],
        o["min_alt_frac"],
        o["min_strand_count"],
        insert_offset=o["insert_offset"],
        insert_len=o["insert_len"],
        include_flanks=o["call_in_cassette"],
    )

    off = o["insert_offset"]
    meta = []
    for r in aligned:
        aln = r["aln"]
        q_positions = aln.q2r.keys()
        q_first, q_last = min(q_positions), max(q_positions)
        meta.append(
            {
                "read_name": r["read_name"],
                "side": r["side"],
                "mapq": r["mapq"],
                "strand": "-" if r["is_reverse"] else "+",
                "revcomp_used": aln.revcomp,
                "ref_start": aln.ref_start + 1,
                "ref_end": aln.ref_end,
                "insert_start": aln.ref_start + 1 - off,
                "insert_end": aln.ref_end - off,
                "clip_len": len(r["seq"]),
                "aligned_len": len(aln.q2r),
                "unaligned_5p": q_first,
                "unaligned_3p": len(r["seq"]) - 1 - q_last,
                "matches": aln.matches,
                "mismatches": aln.mismatches,
                "identity": round(aln.identity, 4),
            }
        )

    return {
        "site_obj": site_obj,
        "n_clips_seen": n_seen,
        "n_clips_used": len(clips),
        "n_aligned": len(aligned),
        "n_rejected": n_rejected,
        "rows": rows,
        "codes": np.array(codes, dtype=np.int8) if codes else np.zeros((0, len(ref_seq)), np.int8),
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
# interactive HTML pileup
# --------------------------------------------------------------------------- #

HTML_CSS = """
:root { --cw: 8.4px; }
body { font-family: system-ui, sans-serif; margin: 0; padding: 0 0 1rem 0; color: #222; }
header { padding: .6rem 1rem; border-bottom: 1px solid #ccc; position: sticky; top: 0;
         background: #fff; z-index: 30; }
h1 { font-size: 1.05rem; margin: 0 0 .3rem 0; }
.meta { font-size: .8rem; color: #555; }
.chips { margin: .4rem 0 0 0; }
.chip { display: inline-block; margin: 0 .25rem .25rem 0; padding: .1rem .45rem; font-size: .78rem;
        border: 1px solid #bbb; border-radius: 10px; cursor: pointer; background: #fff; }
.chip.site_specific { background: #ffe082; border-color: #f0a500; font-weight: 600; }
.chip.insert_reference_mismatch { background: #e0e0e0; color: #666; }
.controls { margin-top: .4rem; font-size: .82rem; }
.controls input[type=number] { width: 6rem; }
#view { overflow: auto; max-height: calc(100vh - 210px); border-top: 1px solid #eee;
        font-family: ui-monospace, "DejaVu Sans Mono", Menlo, monospace; font-size: 12px;
        line-height: 1.15; white-space: pre; position: relative; }
.line { display: flex; }
.lbl { position: sticky; left: 0; z-index: 10; background: #fff; border-right: 1px solid #ddd;
       padding-right: .4rem; flex: 0 0 auto; width: 27ch; overflow: hidden; }
.hdr { position: sticky; z-index: 20; background: #fff; }
.hdr .lbl { z-index: 25; background: #fff; }
.ruler { color: #888; }
.ref { font-weight: 600; }
.cons { color: #00695c; }
.marks { color: #c62828; font-weight: 700; }
.seq .A { color: #2e7d32; background: #e8f5e9; }
.seq .C { color: #1565c0; background: #e3f2fd; }
.seq .G { color: #ef6c00; background: #fff3e0; }
.seq .T { color: #c62828; background: #ffebee; }
.seq .lq { color: #999; background: #f5f5f5; }
.row:hover .lbl { background: #fffde7; }
.row.hit .lbl { background: #ffe082; }
table.summary { border-collapse: collapse; font-size: .8rem; margin: 1rem; }
table.summary th, table.summary td { border: 1px solid #ddd; padding: .2rem .45rem; text-align: left; }
table.summary th { background: #f5f5f5; }
td.na { background: #eee; color: #999; }
"""

HTML_JS = """
function charWidth() {
  const probe = document.getElementById('probe');
  return probe.getBoundingClientRect().width / 100;
}
function goto(pos) {
  const view = document.getElementById('view');
  const lblW = document.querySelector('.lbl').getBoundingClientRect().width;
  view.scrollLeft = Math.max(0, (pos - 1) * charWidth() - view.clientWidth / 2 + lblW);
  highlight(pos);
}
function highlight(pos) {
  document.getElementById('posbox').value = pos;
  const view = document.getElementById('view');
  const bar = document.getElementById('cursor');
  bar.style.left = 'calc(27ch + .4rem + ' + (pos - 1) + ' * ' + charWidth() + 'px)';
  bar.style.height = view.scrollHeight + 'px';
  bar.style.display = 'block';
  const only = document.getElementById('onlyalt').checked;
  document.querySelectorAll('.row').forEach(function (row) {
    const alts = (row.dataset.alts || '').split(',');
    const hit = alts.indexOf(String(pos)) >= 0;
    row.classList.toggle('hit', hit);
    row.style.display = (only && !hit) ? 'none' : '';
  });
}
function applyPos() { goto(parseInt(document.getElementById('posbox').value || '1', 10)); }
function setFont(v) {
  document.getElementById('view').style.fontSize = v + 'px';
  applyPos();
}
window.addEventListener('DOMContentLoaded', function () {
  const p = new URLSearchParams(location.search).get('pos');
  if (p) { goto(parseInt(p, 10)); }
});
"""


def _html_seq_line(row, ref_seq):
    """Render one aligned row: '.' for matches, coloured span for mismatches."""
    out = []
    for i, ch in enumerate(row):
        if ch == "-":
            out.append(" ")
        elif ch.islower():  # low base quality
            out.append(f'<span class="lq">{ch}</span>')
        elif ch == ref_seq[i]:
            out.append(".")
        else:
            out.append(f'<span class="{ch}">{ch}</span>')
    return "".join(out)


def _html_ruler(length, step=10):
    out = []
    i = 0
    while i < length:
        pos = i + 1
        if pos % step == 0:
            label = str(pos)
            if i + len(label) <= length:
                out.append(label)
                i += len(label)
                continue
        out.append("." if pos % 5 == 0 else " ")
        i += 1
    return "".join(out)[:length]


def write_pileup_html(path, title, ref_seq, result, insert_offset, insert_len, min_depth,
                      min_alt_frac, max_rows, classes_by_pos, back_link="index.html"):
    """Self-contained, scrollable HTML pileup for one site."""
    import html as _html

    L = len(ref_seq)
    rows = result["rows"][:max_rows]
    meta = result["meta"][:max_rows]
    depth = result["depth"]
    counts = result["counts"]

    consensus = []
    marks = [" "] * L
    for i in range(L):
        if depth[i] == 0:
            consensus.append(" ")
        else:
            cons = BASES[int(np.argmax(counts[i]))]
            consensus.append("." if cons == ref_seq[i] else cons)
    for v in result["variants"]:
        marks[v["ref_pos"] - 1] = "*" if v["filter"] == "PASS" else "?"

    def line(cls, label, content, extra=""):
        return (
            f'<div class="line {cls}"{extra}><span class="lbl">{label}</span>'
            f'<span class="seq">{content}</span></div>'
        )

    parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>{_html.escape(title)}</title>",
        f"<style>{HTML_CSS}</style><script>{HTML_JS}</script></head><body>",
        "<header>",
        f"<h1>{_html.escape(title)}</h1>",
        '<div class="meta">'
        f"clips seen {result['n_clips_seen']} &middot; used {result['n_clips_used']} &middot; "
        f"aligned {result['n_aligned']} &middot; rejected {result['n_rejected']} &middot; "
        f"rows shown {len(rows)}/{len(result['rows'])} &middot; "
        f"reference {L} bp (insert {insert_len} bp at offset {insert_offset}) &middot; "
        f'<a href="{back_link}">all sites</a></div>',
    ]

    chips = []
    for v in result["variants"]:
        cls = classes_by_pos.get((v["ref_pos"], v["alt"]), "")
        label = f"{v['insert_pos']}{v['insert_ref']}&gt;{v['alt']} {v['alt_frac']:.2f}"
        chips.append(
            f'<span class="chip {cls}" title="{cls or v["filter"]}; '
            f'{v["alt_count"]}/{v["depth"]} reads" onclick="goto({v["ref_pos"]})">{label}</span>'
        )
    parts.append('<div class="chips">' + ("".join(chips) if chips else "no candidate variants") + "</div>")
    parts.append(
        '<div class="controls">'
        'position <input id="posbox" type="number" min="1" value="1" onchange="applyPos()"> '
        '<button onclick="applyPos()">go</button> '
        '<label><input id="onlyalt" type="checkbox" onchange="applyPos()"> only reads with the '
        "alt base at this position</label> "
        'font <input type="range" min="7" max="18" value="12" oninput="setFont(this.value)">'
        "</div>"
    )
    parts.append("</header>")

    parts.append('<div id="view">')
    parts.append('<span id="probe" style="position:absolute;visibility:hidden;">'
                 + "M" * 100 + "</span>")
    parts.append('<div id="cursor" style="position:absolute;top:0;width:1px;'
                 'background:#f00;display:none;z-index:5;"></div>')
    parts.append(line("hdr ruler", "position", _html_ruler(L), extra=' style="top:0"'))
    parts.append(line("hdr ref", "INSERT/REF", _html.escape(ref_seq), extra=' style="top:1.15em"'))
    parts.append(line("hdr cons", "CONSENSUS", "".join(consensus), extra=' style="top:2.3em"'))
    parts.append(line("hdr marks", "VARIANT", "".join(marks), extra=' style="top:3.45em"'))

    variant_positions = [v["ref_pos"] for v in result["variants"]]
    for m, row in zip(meta, rows):
        alts = [p for p in variant_positions if row[p - 1].upper() not in ("-", ref_seq[p - 1])]
        label = _html.escape(f"{m['read_name'][:24]} {m['strand']}{m['side'][0]}")
        parts.append(
            line(
                "row",
                label,
                _html_seq_line(row, ref_seq),
                extra=f' data-alts="{",".join(str(a) for a in alts)}"'
                f' title="{_html.escape(m["read_name"])} mapq={m["mapq"]} '
                f'identity={m["identity"]} clip={m["clip_len"]}bp aligned={m["aligned_len"]}bp"',
            )
        )
    parts.append("</div></body></html>")

    with open(path, "w") as fh:
        fh.write("\n".join(parts))


def variant_summary(variants):
    if not variants:
        return "-"
    return ", ".join(
        f"{v['insert_pos']}{v['insert_ref']}>{v['alt']} ({v['alt_frac']:.2f})" for v in variants
    )


def write_index_html(path, args, insert_name, ref_seq, insert_len, results,
                     classes, matrix_df, conclusion, site_files):
    import html as _html

    parts = [
        "<!doctype html><html><head><meta charset='utf-8'><title>TLA insert soft-clip pileups</title>",
        f"<style>{HTML_CSS}</style></head><body>",
        "<header><h1>TLA insert soft-clip pileups</h1>",
        f'<div class="meta">BAM {_html.escape(args.bam)} &middot; insert '
        f"{_html.escape(insert_name)} ({insert_len} bp) &middot; reference {len(ref_seq)} bp "
        f"&middot; {len(results)} sites</div></header>",
        "<pre style='margin:1rem;font-size:.82rem'>" + _html.escape(conclusion) + "</pre>",
        "<table class='summary'><tr><th>site</th><th>gene</th><th>locus</th><th>clips aligned</th>"
        "<th>max depth</th><th>variants (PASS)</th><th>pileup</th></tr>",
    ]
    for r, fname in zip(results, site_files):
        s = r["site_obj"]
        passing = [v for v in r["variants"] if v["filter"] == "PASS"]
        parts.append(
            "<tr>"
            f"<td>{_html.escape(s['site'])}</td><td>{_html.escape(s['gene'])}</td>"
            f"<td>{_html.escape(s['chrom'])}:{s['pos']}</td>"
            f"<td>{r['n_aligned']}/{r['n_clips_used']}</td>"
            f"<td>{int(r['depth'].max()) if len(r['depth']) else 0}</td>"
            f"<td>{_html.escape(variant_summary(passing))}</td>"
            f"<td><a href='pileup_html/{_html.escape(fname)}'>open</a></td></tr>"
        )
    parts.append("</table>")

    if not matrix_df.empty:
        parts.append("<h2 style='margin:1rem;font-size:.95rem'>Alt fraction per candidate variant</h2>")
        parts.append("<table class='summary'><tr><th>site</th>")
        for col in matrix_df.columns:
            parts.append(f"<th>{_html.escape(str(col))}</th>")
        parts.append("</tr>")
        for idx, row in matrix_df.iterrows():
            parts.append(f"<tr><td>{_html.escape(str(idx))}</td>")
            for val in row:
                if np.isnan(val):
                    parts.append("<td class='na'>n/a</td>")
                else:
                    shade = int(255 - 155 * min(1.0, val))
                    parts.append(
                        f"<td style='background:rgb(255,{shade},{shade})'>{val:.2f}</td>"
                    )
            parts.append("</tr>")
        parts.append("</table>")

    if classes:
        parts.append("<h2 style='margin:1rem;font-size:.95rem'>Variant classification</h2>")
        parts.append("<table class='summary'><tr><th>variant</th><th>class</th>"
                     "<th>sites with alt</th><th>sites covered</th><th>max alt fraction</th></tr>")
        for c in classes:
            parts.append(
                f"<tr><td>{_html.escape(c['label'])}</td><td>{c['class']}</td>"
                f"<td>{_html.escape(c['sites_with_alt'] or '-')}</td>"
                f"<td>{c['n_sites_covered']}</td><td>{c['max_alt_frac']:.2f}</td></tr>"
            )
        parts.append("</table>")

    parts.append("</body></html>")
    with open(path, "w") as fh:
        fh.write("\n".join(parts))


# --------------------------------------------------------------------------- #
# plots
# --------------------------------------------------------------------------- #

MSA_CMAP = ListedColormap(["#f5f5f5", "#dddddd", "#e41a1c", "#fdd0a2"])


def plot_site(result, ref_seq, title, annotate, max_plot_rows, min_depth, highlight_frac,
              max_labels=12):
    codes = result["codes"][:max_plot_rows]
    depth = result["depth"]
    nonref_frac = result["nonref_frac"]
    L = len(ref_seq)

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
    ax2.set_xlabel("Reference position (1-based; insert + cassette flanks if provided)")
    ax2.axhline(min_depth, color="0.4", lw=0.8, ls=":")

    ax3 = ax2.twinx()
    ax3.plot(x, nonref_frac, color="crimson", lw=1.0, label="Non-ref fraction")
    ax3.set_ylim(0, 1.02)
    ax3.set_ylabel("Non-ref fraction")
    ax3.axhline(highlight_frac, color="crimson", lw=0.8, ls=":")

    # every candidate is marked, but only the strongest few are labelled and the
    # labels are staggered vertically so that they do not overlap
    for v in result["variants"]:
        for ax in (ax1, ax2):
            ax.axvline(v["ref_pos"], color="black", ls="--", lw=0.6, alpha=0.6)
    labelled = sorted(result["variants"], key=lambda v: -v["alt_frac"])[:max_labels]
    for n, v in enumerate(sorted(labelled, key=lambda v: v["ref_pos"])):
        ax3.annotate(
            f"{v['insert_pos']}{v['insert_ref']}>{v['alt']} ({v['alt_frac']:.2f})",
            xy=(v["ref_pos"], min(1.0, v["alt_frac"])),
            xytext=(0, 8 + 11 * (n % 3)),
            textcoords="offset points",
            ha="center",
            fontsize=7,
            rotation=0,
            arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"),
        )
    if len(result["variants"]) > max_labels:
        ax3.text(
            0.005,
            0.96,
            f"{len(result['variants'])} candidate positions, {max_labels} labelled",
            transform=ax3.transAxes,
            fontsize=7,
            color="0.35",
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


def plot_cross_site(results, ref_seq, min_depth, highlight_frac):
    """Heatmap sites x reference positions of the non-ref fraction (masked on depth)."""
    L = len(ref_seq)
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
    show_values = matrix_df.shape[1] <= 25
    for i in range(matrix_df.shape[0]):
        for j in range(matrix_df.shape[1]) if show_values else ():
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
    ap.add_argument("--insert-fasta", default="insert_sequence.fasta", help="FASTA with the insert")
    ap.add_argument(
        "--cassette-fasta",
        default=None,
        help="FASTA of the full construct (cassette) that contains the insert. Soft clips of "
             "junction spanning reads also contain cassette sequence; giving it here lets those "
             "bases align instead of being forced into insert coordinates.",
    )
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
    g.add_argument("--trim-window", type=int, default=15,
                   help="Window used to trim badly matching alignment ends (0 disables)")
    g.add_argument("--trim-identity", type=float, default=0.7,
                   help="Minimum identity in the trim window to keep an alignment end")

    g = ap.add_argument_group("variant calling")
    g.add_argument("--min-bq", type=int, default=20, help="Min base quality to count a base")
    g.add_argument("--min-depth", type=int, default=10)
    g.add_argument("--min-alt-count", type=int, default=5)
    g.add_argument("--min-alt-frac", type=float, default=0.2)
    g.add_argument("--min-strand-count", type=int, default=1, help="Min alt reads per strand")
    g.add_argument("--call-in-cassette", action="store_true",
                   help="Also call variants in the cassette flanks (needs --cassette-fasta)")

    g = ap.add_argument_group("output")
    g.add_argument("--max-plot-rows", type=int, default=150)
    g.add_argument("--max-pileup-rows", type=int, default=500)
    g.add_argument("--wrap", type=int, default=100, help="Columns per block in the text pileup")
    g.add_argument("--emit-counts", action="store_true", help="Write full per-position base counts")
    g.add_argument("--max-html-rows", type=int, default=500,
                   help="Maximum read rows per interactive HTML pileup")
    g.add_argument("--no-html", action="store_true", help="Do not write the HTML pileup viewer")
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
    html_dir = os.path.join(args.outdir, "pileup_html")
    dirs = [plots_dir, msa_dir, pileup_dir] + ([] if args.no_html else [html_dir])
    for d in dirs:
        os.makedirs(d, exist_ok=True)

    insert_name, insert_seq = read_first_fasta_seq(args.insert_fasta)
    if not insert_seq:
        raise SystemExit(f"[ERROR] no sequence found in {args.insert_fasta}")
    if args.kmer_k > len(insert_seq):
        raise SystemExit("[ERROR] --kmer-k is larger than the insert sequence")
    ref_name, ref_seq, insert_offset = build_reference(insert_seq, args.cassette_fasta)
    insert_len = len(insert_seq)
    if args.call_in_cassette and not args.cassette_fasta:
        print("[WARN] --call-in-cassette has no effect without --cassette-fasta")
    kmer_index = build_kmer_index(ref_seq, args.kmer_k, args.kmer_max_hits)

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
        "trim_window": args.trim_window,
        "trim_identity": args.trim_identity,
        "insert_offset": insert_offset,
        "insert_len": insert_len,
        "call_in_cassette": bool(args.call_in_cassette and args.cassette_fasta),
    }

    print(f"[INFO] insert '{insert_name}' length={insert_len} from {args.insert_fasta}")
    print(f"[INFO] sites={len(sites)} workers={args.workers} gapped_alignment={'edlib' if use_edlib else 'ungapped'}")

    results = []
    if args.workers > 1:
        from concurrent.futures import ProcessPoolExecutor, as_completed

        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=_worker_init,
            initargs=(ref_seq, kmer_index, opts),
        ) as ex:
            futs = {ex.submit(process_site, s): s for s in sites}
            for i, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                _log_site(i, len(sites), r)
                results.append(r)
    else:
        _worker_init(ref_seq, kmer_index, opts)
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
                "ref_bp_covered": int((r["depth"] > 0).sum()),
                "insert_bp_covered": int(
                    (r["depth"][insert_offset : insert_offset + insert_len] > 0).sum()
                ),
                "mean_clip_len": (
                    round(float(np.mean([m["clip_len"] for m in r["meta"]])), 1) if r["meta"] else 0.0
                ),
                "mean_unaligned_flank": (
                    round(
                        float(np.mean([m["unaligned_5p"] + m["unaligned_3p"] for m in r["meta"]])), 1
                    )
                    if r["meta"]
                    else 0.0
                ),
                "mean_depth": round(float(r["depth"].mean()), 2),
                "max_depth": int(r["depth"].max()) if len(r["depth"]) else 0,
                "n_variants_pass": len(passing),
                "variants": ";".join(
                    f"{v['insert_pos']}{v['insert_ref']}>{v['alt']}"
                    f"({v['alt_count']}/{v['depth']},{v['alt_frac']:.2f})"
                    for v in passing[:20]
                ) + (f";+{len(passing) - 20} more" if len(passing) > 20 else ""),
                "elapsed_sec": round(r["elapsed"], 2),
            }
        )
        for v in r["variants"]:
            variant_rows.append(
                {"site": s["site"], "gene": s["gene"], "chrom": s["chrom"], "pos": s["pos"], **v}
            )
        if args.emit_counts:
            for i in range(len(ref_seq)):
                if r["depth"][i] == 0:
                    continue
                count_rows.append(
                    {
                        "site": s["site"],
                        "ref_pos": i + 1,
                        "region": region_label(i, insert_offset, insert_len),
                        "insert_pos": i + 1 - insert_offset,
                        "insert_ref": ref_seq[i],
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
    keys = variant_keys(results)
    matrix_df = build_variant_matrix(results, keys, args.min_depth)
    classes = classify_variants(
        results, keys, args.min_depth, args.min_alt_count, args.min_alt_frac
    )
    class_by_pos = {(c["ref_pos"], c["alt"]): c["class"] for c in classes}
    if not variant_df.empty:
        variant_df["class"] = [
            class_by_pos.get((row.ref_pos, row.alt), "not_confirmed")
            for row in variant_df.itertuples()
        ]
    consensus_seq, diffs = consensus_and_diff(
        results, ref_seq, insert_offset, insert_len, args.min_depth
    )

    # ------------------------------------------------------------------ files
    site_files = []
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
            ref_seq,
            r,
            r["meta"],
            args.wrap,
            args.min_depth,
            args.min_alt_frac,
            annotate,
            args.max_pileup_rows,
        )
        site_files.append(f"{tag}.html")
        if not args.no_html:
            write_pileup_html(
                os.path.join(html_dir, f"{tag}.html"),
                title,
                ref_seq,
                r,
                insert_offset,
                insert_len,
                args.min_depth,
                args.min_alt_frac,
                args.max_html_rows,
                class_by_pos,
                back_link="../pileup_index.html",
            )
        with open(os.path.join(msa_dir, f"{tag}.msa.tsv"), "w") as fh:
            fh.write(f"# {title}\n")
            cols = ["read_name", "side", "strand", "mapq", "revcomp_used", "ref_start", "ref_end",
                    "insert_start", "insert_end", "clip_len", "aligned_len", "unaligned_5p",
                    "unaligned_3p", "matches", "mismatches", "identity"]
            fh.write("\t".join(cols + ["aligned_row"]) + "\n")
            for m, row in zip(r["meta"], r["rows"]):
                fh.write("\t".join(str(m[c]) for c in cols) + f"\t{row}\n")

    summary_tsv = os.path.join(args.outdir, "site_summary.tsv")
    summary_df.to_csv(summary_tsv, sep="\t", index=False)

    variant_tsv = os.path.join(args.outdir, "insert_variants.tsv")
    if variant_df.empty:
        variant_df = pd.DataFrame(
            columns=["site", "gene", "chrom", "pos", "ref_pos", "region", "insert_pos",
                     "insert_ref", "alt", "depth", "alt_count", "alt_frac", "alt_fwd", "alt_rev",
                     "nonref_count", "nonref_frac", "lowqual_bases", "zygosity", "filter", "class"]
        )
    variant_df.to_csv(variant_tsv, sep="\t", index=False)

    classes_tsv = os.path.join(args.outdir, "variant_classification.tsv")
    pd.DataFrame(
        classes,
        columns=["label", "ref_pos", "region", "insert_pos", "insert_ref", "alt", "class",
                 "sites_with_alt", "n_sites_alt", "n_sites_covered", "max_alt_frac"],
    ).to_csv(classes_tsv, sep="\t", index=False)

    consensus_fasta = os.path.join(args.outdir, "pileup_consensus.fasta")
    with open(consensus_fasta, "w") as fh:
        fh.write(f">{ref_name}_pileup_consensus\n")
        for i in range(0, len(consensus_seq), 60):
            fh.write(consensus_seq[i : i + 60] + "\n")
    diff_tsv = os.path.join(args.outdir, "insert_vs_consensus.tsv")
    pd.DataFrame(
        diffs,
        columns=["ref_pos", "region", "insert_pos", "insert_ref", "consensus", "depth",
                 "consensus_count", "consensus_frac"],
    ).to_csv(diff_tsv, sep="\t", index=False)

    matrix_tsv = os.path.join(args.outdir, "variant_matrix.tsv")
    matrix_df.to_csv(matrix_tsv, sep="\t")

    counts_tsv = None
    if args.emit_counts:
        counts_tsv = os.path.join(args.outdir, "pileup_counts.tsv")
        pd.DataFrame(count_rows).to_csv(counts_tsv, sep="\t", index=False)

    conclusion = _conclusion_text(classes, len(diffs), insert_len)

    # -------------------------------------------------------------------- pdf
    pdf_path = os.path.join(args.outdir, "tla_insert_softclip_report.pdf")
    with PdfPages(pdf_path) as pdf:
        cover = (
            "TLA insert soft-clip pileup report\n\n"
            f"BAM            : {args.bam}\n"
            f"Insert         : {args.insert_fasta} ({insert_name}, {insert_len} bp)\n"
            f"Reference      : {ref_name} ({len(ref_seq)} bp, insert at offset {insert_offset})\n"
            f"Sites          : {len(sites)}\n"
            f"Clip filters   : window={args.window} mapq>={args.mapq} clip>={args.min_clip} "
            f"max_clips={args.max_clips}\n"
            f"Alignment      : k={args.kmer_k} min_identity={args.min_identity} "
            f"min_match={args.min_match} gapped={'edlib' if use_edlib else 'no'}\n"
            f"Variant filters: bq>={args.min_bq} depth>={args.min_depth} "
            f"alt>={args.min_alt_count} frac>={args.min_alt_frac} "
            f"strand>={args.min_strand_count}\n\n"
            + conclusion
        )
        pdf.savefig(text_page(cover))
        plt.close("all")

        pdf.savefig(plot_cross_site(results, ref_seq, args.min_depth, args.min_alt_frac))
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
            fig = plot_site(r, ref_seq, title, annotate, args.max_plot_rows,
                            args.min_depth, args.min_alt_frac)
            fig.savefig(os.path.join(plots_dir, f"{tag}.png"), dpi=150)
            pdf.savefig(fig)
            plt.close(fig)

    index_html = os.path.join(args.outdir, "pileup_index.html")
    if not args.no_html:
        write_index_html(index_html, args, insert_name, ref_seq, insert_len, results,
                         classes, matrix_df, conclusion, site_files)

    with open(os.path.join(args.outdir, "run_parameters.json"), "w") as fh:
        json.dump(vars(args), fh, indent=2, sort_keys=True)

    print(f"[OK] PDF            : {pdf_path}")
    print(f"[OK] Pileups        : {pileup_dir}")
    print(f"[OK] Summary        : {summary_tsv}")
    print(f"[OK] Variants       : {variant_tsv}")
    print(f"[OK] Variant matrix : {matrix_tsv}")
    print(f"[OK] Classification : {classes_tsv}")
    print(f"[OK] Consensus      : {consensus_fasta} ({len(diffs)} differences -> {diff_tsv})")
    if not args.no_html:
        print(f"[OK] HTML viewer    : {index_html}")
    if counts_tsv:
        print(f"[OK] Base counts    : {counts_tsv}")
    print(conclusion)
    print(f"[DONE] total runtime {time.time() - t_all:.2f}s")
    return 0


def _log_site(i, n, r, max_shown=6):
    s = r["site_obj"]
    passing = [x for x in r["variants"] if x["filter"] == "PASS"]
    v = ",".join(
        f"{x['insert_pos']}{x['insert_ref']}>{x['alt']}({x['alt_frac']:.2f})"
        for x in passing[:max_shown]
    )
    if len(passing) > max_shown:
        v += f",+{len(passing) - max_shown} more"
    print(
        f"[{i}/{n}] site {s['site']:>3} {s['chrom']}:{s['pos']} "
        f"clips={r['n_clips_used']} aligned={r['n_aligned']} "
        f"max_depth={int(r['depth'].max()) if len(r['depth']) else 0} "
        f"variants={v or '-'} ({r['elapsed']:.2f}s)"
    )


def variant_keys(results):
    """Unique (ref_pos, ref_base, alt_base) keys of all PASSing calls."""
    keys = {}
    for r in results:
        for v in r["variants"]:
            if v["filter"] == "PASS":
                keys[(v["ref_pos"], v["alt"])] = v
    return [keys[k] for k in sorted(keys)]


def build_variant_matrix(results, keys, min_depth):
    """sites x candidate variants matrix holding the fraction of the *alt* base.

    The fraction is NaN when a site has less than `min_depth` coverage at that
    position, so "not observed" and "not covered" cannot be confused.
    """
    if not keys:
        return pd.DataFrame()
    index, data = [], []
    for r in results:
        s = r["site_obj"]
        index.append(f"site {s['site']} ({s['gene']})")
        row = []
        for v in keys:
            i = v["ref_pos"] - 1
            depth = int(r["depth"][i])
            if depth < min_depth:
                row.append(np.nan)
            else:
                row.append(float(r["counts"][i, BASE_IDX[v["alt"]]]) / depth)
        data.append(row)
    labels = [variant_label(v) for v in keys]
    return pd.DataFrame(data, index=index, columns=labels)


def variant_label(v):
    pos = v["insert_pos"] if v["region"] == "insert" else f"{v['region']}:{v['insert_pos']}"
    return f"{pos}{v['insert_ref']}>{v['alt']}"


def classify_variants(results, keys, min_depth, min_alt_count, min_alt_frac):
    """Decide per candidate variant whether it is specific for a single site.

    A variant that is present in *every* site that has coverage - especially at
    ~100% - is not an integration-specific SNP but a difference between the
    supplied insert FASTA and the construct that was actually integrated.
    """
    out = []
    for v in keys:
        i = v["ref_pos"] - 1
        bi = BASE_IDX[v["alt"]]
        covered, positive, fracs = [], [], {}
        for r in results:
            depth = int(r["depth"][i])
            if depth < min_depth:
                continue
            site = r["site_obj"]["site"]
            count = int(r["counts"][i, bi])
            frac = count / depth
            covered.append(site)
            fracs[site] = frac
            if count >= min_alt_count and frac >= min_alt_frac:
                positive.append(site)
        if not covered:
            continue
        if len(positive) == 1 and len(covered) > 1:
            cls = "site_specific"
        elif positive and len(positive) == len(covered):
            median_frac = float(np.median([fracs[s] for s in positive]))
            cls = "insert_reference_mismatch" if median_frac >= 0.9 else "shared_all_sites"
        elif len(positive) > 1:
            cls = "shared_subset"
        else:
            cls = "not_confirmed"
        out.append(
            {
                "label": variant_label(v),
                "ref_pos": v["ref_pos"],
                "region": v["region"],
                "insert_pos": v["insert_pos"],
                "insert_ref": v["insert_ref"],
                "alt": v["alt"],
                "class": cls,
                "sites_with_alt": ",".join(positive),
                "n_sites_alt": len(positive),
                "n_sites_covered": len(covered),
                "max_alt_frac": round(max(fracs.values()), 4) if fracs else 0.0,
            }
        )
    return out


def consensus_and_diff(results, ref_seq, insert_offset, insert_len, min_depth):
    """Pooled consensus over all sites plus the positions where it differs.

    A long list of differences means the insert FASTA is not the sequence that
    was integrated (wrong construct version, vector backbone, orientation, ...),
    which is the usual reason for hundreds of "100% variants".
    """
    total = np.zeros((len(ref_seq), 4), dtype=np.int64)
    for r in results:
        total += r["counts"]
    depth = total.sum(axis=1)
    consensus, diffs = [], []
    for i, refb in enumerate(ref_seq):
        if depth[i] < min_depth:
            consensus.append("N" if depth[i] == 0 else refb.lower())
            continue
        bi = int(np.argmax(total[i]))
        cons = BASES[bi]
        consensus.append(cons)
        if cons != refb:
            diffs.append(
                {
                    "ref_pos": i + 1,
                    "region": region_label(i, insert_offset, insert_len),
                    "insert_pos": i + 1 - insert_offset,
                    "insert_ref": refb,
                    "consensus": cons,
                    "depth": int(depth[i]),
                    "consensus_count": int(total[i, bi]),
                    "consensus_frac": round(float(total[i, bi]) / int(depth[i]), 4),
                }
            )
    return "".join(consensus), diffs


def _conclusion_text(classes, n_diffs=0, insert_len=0):
    lines = []
    if n_diffs:
        pct = 100.0 * n_diffs / insert_len if insert_len else 0.0
        lines.append(
            f"NOTE: the pooled consensus of all sites differs from the supplied insert FASTA at "
            f"{n_diffs} position(s) ({pct:.1f}% of the insert). If this number is large the FASTA "
            f"is not the construct that was integrated; see insert_vs_consensus.tsv and "
            f"pileup_consensus.fasta."
        )
    if not classes:
        lines.append("No candidate variants passed the filters in any site.")
        return "\n".join(lines)

    site_specific = [c for c in classes if c["class"] == "site_specific"]
    ref_mismatch = [c for c in classes if c["class"] == "insert_reference_mismatch"]
    shared = [c for c in classes if c["class"] in ("shared_all_sites", "shared_subset")]

    lines.append("Site-specific variants (present in exactly one site):")
    if site_specific:
        for v in site_specific:
            lines.append(
                f"  {v['label']} -> site {v['sites_with_alt']} "
                f"(alt fraction {v['max_alt_frac']:.2f}; {v['n_sites_covered']} sites covered)"
            )
    else:
        lines.append("  none")
    if ref_mismatch:
        lines.append(
            f"  {len(ref_mismatch)} position(s) are non-reference in *all* covered sites at >=90%: "
            f"these are insert-FASTA/reference discrepancies, not integration-specific SNPs."
        )
    if shared:
        lines.append(
            "  " + str(len(shared)) + " position(s) are shared by several sites: "
            + ", ".join(c["label"] for c in shared[:20])
            + (" ..." if len(shared) > 20 else "")
        )
    return "\n".join(lines)


if __name__ == "__main__":
    sys.exit(main())
