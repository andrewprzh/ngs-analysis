#!/usr/bin/env python3
############################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

"""Phase 1 of single-cell / spatial variant calling on top of IsoQuant: raw SNVs and short indels counted per UMI.

Input is IsoQuant single-cell / spatial output: read_info (UMI-deduplicated reads with gene, strand, barcode and UMI),
umi_groups (all reads of each deduplication group, --large_output umi_groups) and the alignments IsoQuant was run on.
Only reads assigned to a single gene are used (--gene_assignment_types).

One molecule (UMI) is one piece of evidence. At each site, among the reads of a molecule that cover the site,
the molecule is alt if strictly more than --consensus of them carry the alt allele, ref if strictly more than
--consensus of them carry the reference, and is not counted otherwise (no call).
A variant is reported if at least --min_alt_umis molecules from at least --min_alt_barcodes barcodes are alt,
its VAF (alt / (alt + ref) molecules) is at least --min_vaf, and the alt molecules exceed sequencing errors:
the binomial p-value of alt molecules out of depth molecules at the error rate of the variant class is at most
--max_error_pvalue. Error classes are SNVs by transcript-strand change and indels by type, length (1, 2, 3+) and
homopolymer length (1-5, 6+). Error rates are estimated first, on random loci holding --error_sample_fraction of
molecules: per-read errors at all covered positions, from molecules with a single read, leaving out alleles seen in
at least 10% (and at least 3) of the reads at a position. A class with fewer than 100 errors uses the rate of its
closest pooled parent (e.g. DEL:3+:hp6+ -> DEL:3+ -> DEL). Rates can also be given with --error_rates.
Molecules with several reads need a read majority, so their error rate is lower and the test is conservative.

All reported sites are kept; the column `unreliable` lists reasons to doubt a site (or "."):
  read_end     in at least --end_fraction of the alt reads the site is within --end_distance aligned bases of the
               alignment end (adapter bases aligned past transcript ends, e.g. the TSO GGG or the polyA);
  editing      A>G or C>T on the transcript strand (A>G / C>T in + strand molecules, T>C / G>A in - strand ones) in
               at least half of the stranded alt molecules; the strand of a molecule is that of its assigned isoform
               (read_info);
  homopolymer  an indel overlapping or next to (1 bp) a reference run of at least --homopolymer_flag_len bases.

Each locus (overlapping genes merged) is read up to three times: the first pass counts molecules with at least one
alt read per allele (an upper bound for the consensus count, used to select candidates), the second pass computes the
molecule consensus at the candidate sites (totals and alt molecules per barcode) and applies the cut-offs, and the
third pass, for loci with reported sites, counts depth per barcode, alt read positions and alt molecule strands
at those sites only. Reads of a UMI group can extend past the locus of their molecule, so two loci can report the
same allele; such records count different molecules and are merged into one.

Outputs:
  PREFIX.sites.tsv              one line per variant allele;
  PREFIX.barcode_counts.tsv.gz  alt and depth UMIs per variant and barcode (every barcode with depth > 0);
  PREFIX.error_rates.tsv        error rate per class.
"""

import argparse
import gzip
import heapq
import logging
import os
import shutil
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy
import pysam
from scipy.stats import binom

logger = logging.getLogger("call_variants")

SINGLE_GENE_TYPES = "unique,unique_minor_difference,inconsistent,inconsistent_non_intronic,inconsistent_genic"
NO_BARCODE = {"", ".", "*", "None"}

# CIGAR operations
OP_M, OP_I, OP_D, OP_N, OP_S, OP_H, OP_P, OP_EQ, OP_X = range(9)
ALIGNED_OPS = (OP_M, OP_EQ, OP_X)
EDGE_OPS = (OP_N, OP_S, OP_H)

BASES = "ACGT"
BASE_CODE = numpy.full(256, 255, dtype=numpy.uint8)
for _i, _b in enumerate(BASES):
    BASE_CODE[ord(_b)] = _i
    BASE_CODE[ord(_b.lower())] = _i
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

# variant keys are int64: SNV pos * 64 + base (0-3), deletion pos * 64 + 4 + (len - 1);
# insertions are keyed by (pos, inserted sequence) tuples
KEY_SHIFT = 64
DEL_BASE = 4

# read status at a site
ST_NONE, ST_ALT, ST_REF, ST_OTHER = 0, 1, 2, 3

# compress accumulated (candidate, barcode) observations when the buffer grows above this size
ACCUMULATOR_BUFFER = 20000000

# error rate estimation: alleles in at least ESTIMATION_MAX_VAF of the reads (and at least ESTIMATION_MIN_EXCLUDED
# reads) at a position are not errors; classes with fewer than MIN_CLASS_ALT errors use a pooled parent rate
ESTIMATION_MAX_VAF = 0.1
ESTIMATION_MIN_EXCLUDED = 3
MIN_CLASS_ALT = 100

# reference runs at most this far from an indel count for the homopolymer flag
HOMOPOLYMER_FLAG_DISTANCE = 1

# genome change -> isoform strand on which it is A>G or C>T on the transcript (editing-like)
EDITING_CHANGES = {"A>G": "+", "C>T": "+", "T>C": "-", "G>A": "-"}


def open_text(fname):
    return gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname)


# ---------------------------------------------------------------- input splitting

def parse_exon_span(exons):
    """Span (0-based start, end) of an IsoQuant exon string 's1-e1,s2-e2' (1-based inclusive)."""
    start, end = None, None
    for block in exons.split(","):
        if "-" not in block:
            continue
        s, e = block.split("-")
        s, e = int(s), int(e)
        start = s if start is None else min(start, s)
        end = e if end is None else max(end, e)
    if start is None:
        return None
    return start - 1, end


def split_read_info(read_info, gene_types, tmp_dir):
    """Write kept molecules per chromosome; return chromosome -> (file, gene spans)."""
    columns = {"read_id": 0, "chr": 1, "strand": 2, "gene_id": 3, "gene_assignment_type": 4, "exons": 9,
               "barcode": 13, "umi": 14}
    chr_files, chr_handles = {}, {}
    gene_spans = defaultdict(dict)
    seen = set()
    stats = defaultdict(int)
    for line in open_text(read_info):
        if line.startswith("#"):
            continue
        if line.startswith("read_id\t"):
            header = line.rstrip("\n").split("\t")
            columns = {name: header.index(name) for name in columns}
            continue
        v = line.rstrip("\n").split("\t")
        read_id = v[columns["read_id"]]
        if read_id in seen:
            # reads assigned to several isoforms of the same gene occur several times
            continue
        seen.add(read_id)
        stats["reads"] += 1
        if v[columns["gene_assignment_type"]] not in gene_types:
            stats["not single-gene"] += 1
            continue
        barcode = v[columns["barcode"]]
        if barcode in NO_BARCODE:
            stats["no barcode"] += 1
            continue
        span = parse_exon_span(v[columns["exons"]])
        chr_id, gene_id = v[columns["chr"]], v[columns["gene_id"]]
        if span is None or chr_id == "." or gene_id == ".":
            stats["no exons/gene"] += 1
            continue
        stats["molecules kept"] += 1
        if chr_id not in chr_handles:
            chr_files[chr_id] = os.path.join(tmp_dir, "chr%d.molecules.tsv" % len(chr_files))
            chr_handles[chr_id] = open(chr_files[chr_id], "w")
        chr_handles[chr_id].write("%s\t%s\t%s\t%s\n" % (read_id, gene_id, v[columns["strand"]], barcode))
        spans = gene_spans[chr_id]
        if gene_id in spans:
            s, e, n = spans[gene_id]
            spans[gene_id] = (min(s, span[0]), max(e, span[1]), n + 1)
        else:
            spans[gene_id] = (span[0], span[1], 1)
    for f in chr_handles.values():
        f.close()
    logger.info("read_info: " + ", ".join("%s %d" % (k, v) for k, v in stats.items()))
    return chr_files, gene_spans


def split_umi_groups(umi_groups, chr_files, tmp_dir):
    """Write groups with more than one read per chromosome."""
    group_files = {chr_id: os.path.join(tmp_dir, os.path.basename(f).replace(".molecules", ".groups"))
                   for chr_id, f in chr_files.items()}
    handles = {}
    multi, total = 0, 0
    if umi_groups:
        for line in open_text(umi_groups):
            if line.startswith("#"):
                continue
            read_id, chr_id, group_size, read_ids = line.rstrip("\n").split("\t")
            total += 1
            if int(group_size) < 2 or chr_id not in group_files:
                continue
            multi += 1
            if chr_id not in handles:
                handles[chr_id] = open(group_files[chr_id], "w")
            handles[chr_id].write("%s\t%s\n" % (read_id, read_ids))
        logger.info("umi_groups: %d molecules, %d with more than one read" % (total, multi))
    for f in handles.values():
        f.close()
    return {chr_id: f for chr_id, f in group_files.items() if chr_id in handles}


def merge_loci(spans):
    """Merge overlapping gene spans into loci: list of (start, end, genes, molecule count)."""
    loci = []
    for gene_id, (s, e, n) in sorted(spans.items(), key=lambda x: x[1][0]):
        if loci and s < loci[-1][1]:
            ls, le, genes, ln = loci[-1]
            genes.append(gene_id)
            loci[-1] = (ls, max(le, e), genes, ln + n)
        else:
            loci.append((s, e, [gene_id], n))
    return loci


def make_jobs(chr_files, group_files, gene_spans, chr_order, tmp_dir, threads):
    """Split each chromosome into jobs of consecutive loci with similar molecule counts; write per-job inputs."""
    total = sum(n for spans in gene_spans.values() for _, _, n in spans.values())
    target = max(10000, total // (threads * 8))
    jobs = []
    for chr_id in sorted(chr_files, key=lambda c: chr_order.get(c, len(chr_order))):
        loci = merge_loci(gene_spans[chr_id])
        chr_jobs, current, count = [], [], 0
        for locus in loci:
            current.append(locus)
            count += locus[3]
            if count >= target:
                chr_jobs.append(current)
                current, count = [], 0
        if current:
            chr_jobs.append(current)

        gene_job = {}
        for j, job_loci in enumerate(chr_jobs):
            for locus_index, (_, _, genes, _) in enumerate(job_loci):
                for g in genes:
                    gene_job[g] = (j, locus_index)
        prefixes = [os.path.join(tmp_dir, "job%d" % (len(jobs) + j)) for j in range(len(chr_jobs))]
        mol_handles = [open(p + ".molecules.tsv", "w") for p in prefixes]
        rep_job = {}
        for line in open(chr_files[chr_id]):
            read_id, gene_id, strand, barcode = line.rstrip("\n").split("\t")
            j, locus_index = gene_job[gene_id]
            rep_job[read_id] = j
            mol_handles[j].write("%s\t%s\t%s\t%s\t%d\n" % (read_id, gene_id, strand, barcode, locus_index))
        for f in mol_handles:
            f.close()
        group_handles = [open(p + ".groups.tsv", "w") for p in prefixes]
        if chr_id in group_files:
            for line in open(group_files[chr_id]):
                read_id = line.split("\t", 1)[0]
                if read_id in rep_job:
                    group_handles[rep_job[read_id]].write(line)
        for f in group_handles:
            f.close()
        for p, job_loci in zip(prefixes, chr_jobs):
            jobs.append((p, chr_id, [(s, e) for s, e, _, _ in job_loci], [n for _, _, _, n in job_loci]))
    return jobs


# ---------------------------------------------------------------- read parsing

class ReadAlignment:
    """Aligned blocks, intron-free spans and variant events of one alignment."""
    __slots__ = ("start", "end", "m_starts", "m_ends", "span_starts", "span_ends",
                 "snv_pos", "snv_base", "indels", "indel_pos")


class Reference:
    def __init__(self, fasta, chr_id, start, end):
        self.start = max(0, start)
        self.seq = fasta.fetch(chr_id, self.start, end).upper()
        self.end = self.start + len(self.seq)
        self.codes = BASE_CODE[numpy.frombuffer(self.seq.encode(), dtype=numpy.uint8)]

    def base(self, pos):
        return self.seq[pos - self.start]


def left_align_deletion(ref, s0, length, min_start):
    while s0 - 2 >= min_start and s0 - 1 >= ref.start and ref.base(s0 - 1) == ref.base(s0 + length - 1):
        s0 -= 1
    return s0


def left_align_insertion(ref, s0, seq, min_start):
    while s0 - 2 >= min_start and s0 - 1 >= ref.start and ref.base(s0 - 1) == seq[-1]:
        seq = seq[-1] + seq[:-1]
        s0 -= 1
    return s0, seq


def parse_alignment(aln, ref, max_indel_len):
    """Extract mismatches and short indels; indels next to introns, clips or alignment ends are ignored."""
    r = ReadAlignment()
    ops = aln.cigartuples
    query = aln.query_sequence
    qcodes = BASE_CODE[numpy.frombuffer(query.encode(), dtype=numpy.uint8)]
    rpos, qpos = aln.reference_start, 0
    span_start = rpos
    m_starts, m_ends, span_starts, span_ends = [], [], [], []
    snv_pos, snv_base = [], []
    indels, indel_pos = set(), set()
    last = len(ops) - 1
    for i, (op, length) in enumerate(ops):
        if op in ALIGNED_OPS:
            q = qcodes[qpos:qpos + length]
            rc = ref.codes[rpos - ref.start:rpos - ref.start + length]
            mism = numpy.flatnonzero((q != rc) & (q < 4) & (rc < 4))
            if len(mism):
                snv_pos.append(mism + rpos)
                snv_base.append(q[mism])
            m_starts.append(rpos)
            m_ends.append(rpos + length)
            rpos += length
            qpos += length
        elif op == OP_I or op == OP_D:
            internal = 0 < i < last and ops[i - 1][0] not in EDGE_OPS and ops[i + 1][0] not in EDGE_OPS
            if op == OP_I:
                if internal and length <= max_indel_len:
                    s0, seq = left_align_insertion(ref, rpos, query[qpos:qpos + length].upper(), span_start)
                    if "N" not in seq:
                        indels.add((s0, seq))
                        indel_pos.add(s0)
                qpos += length
            else:
                if internal and length <= max_indel_len:
                    s0 = left_align_deletion(ref, rpos, length, span_start)
                    indels.add(s0 * KEY_SHIFT + DEL_BASE + length - 1)
                    indel_pos.add(s0)
                rpos += length
        elif op == OP_N:
            span_starts.append(span_start)
            span_ends.append(rpos)
            rpos += length
            span_start = rpos
        elif op == OP_S:
            qpos += length
    span_starts.append(span_start)
    span_ends.append(rpos)
    r.start, r.end = aln.reference_start, rpos
    r.m_starts, r.m_ends = m_starts, m_ends
    r.span_starts, r.span_ends = span_starts, span_ends
    if snv_pos:
        r.snv_pos = numpy.concatenate(snv_pos)
        r.snv_base = numpy.concatenate(snv_base)
    else:
        r.snv_pos = numpy.zeros(0, dtype=numpy.int64)
        r.snv_base = numpy.zeros(0, dtype=numpy.uint8)
    r.indels, r.indel_pos = indels, indel_pos
    return r


def read_event_keys(r):
    """Integer keys of SNVs and deletions, and insertion tuples, of one read."""
    keys = r.snv_pos.astype(numpy.int64) * KEY_SHIFT + r.snv_base
    dels = [k for k in r.indels if not isinstance(k, tuple)]
    if dels:
        keys = numpy.concatenate([keys, numpy.array(dels, dtype=numpy.int64)])
    return keys, [k for k in r.indels if isinstance(k, tuple)]


# ---------------------------------------------------------------- locus processing

class Molecules:
    def __init__(self, job_prefix):
        self.rep_ids, self.genes, self.strands, self.barcodes, self.locus = [], [], [], [], []
        index = {}
        for line in open(job_prefix + ".molecules.tsv"):
            read_id, gene_id, strand, barcode, locus_index = line.rstrip("\n").split("\t")
            index[read_id] = len(self.rep_ids)
            self.rep_ids.append(read_id)
            self.genes.append(gene_id)
            self.strands.append(strand)
            self.barcodes.append(barcode)
            self.locus.append(int(locus_index))
        self.read_to_mol = dict(index)
        self.group_size = [1] * len(self.rep_ids)
        for line in open(job_prefix + ".groups.tsv"):
            read_id, read_ids = line.rstrip("\n").split("\t")
            m = index[read_id]
            members = read_ids.split(",")
            self.group_size[m] = len(members)
            for member in members:
                self.read_to_mol[member] = m
        barcode_index = {}
        self.barcode_idx = numpy.array([barcode_index.setdefault(b, len(barcode_index)) for b in self.barcodes],
                                       dtype=numpy.int64)
        self.barcode_list = [None] * len(barcode_index)
        for b, i in barcode_index.items():
            self.barcode_list[i] = b


class BamSet:
    """Several indexed BAMs read as one (e.g. two runs of one library): fetch() merges them by position."""

    def __init__(self, fnames):
        self.bams = [pysam.AlignmentFile(f) for f in fnames]
        self.references = self.bams[0].references
        for f, bam in zip(fnames[1:], self.bams[1:]):
            if bam.references != self.references:
                raise ValueError("%s has other references than %s" % (f, fnames[0]))

    def fetch(self, chr_id, start, end):
        if len(self.bams) == 1:
            return self.bams[0].fetch(chr_id, start, end)
        return heapq.merge(*(b.fetch(chr_id, start, end) for b in self.bams), key=lambda aln: aln.reference_start)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        for bam in self.bams:
            bam.close()


def locus_alignments(bam, chr_id, start, end, molecules, locus_index):
    """Primary alignments of reads belonging to molecules of this locus, as (molecule, alignment)."""
    for aln in bam.fetch(chr_id, start, end):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.cigartuples is None:
            continue
        m = molecules.read_to_mol.get(aln.query_name)
        if m is None or molecules.locus[m] != locus_index:
            continue
        yield m, aln


def reference_window(bam, chr_id, start, end, molecules, locus_index, fasta, max_indel_len):
    lo, hi = start, end
    for _, aln in locus_alignments(bam, chr_id, start, end, molecules, locus_index):
        lo = min(lo, aln.reference_start)
        hi = max(hi, aln.reference_end)
    return Reference(fasta, chr_id, lo - 1000, hi + 1000)


class KeyCounter:
    """Counts of int64 keys accumulated in numpy batches."""
    def __init__(self):
        self.keys = numpy.zeros(0, dtype=numpy.int64)
        self.counts = numpy.zeros(0, dtype=numpy.int64)
        self.buffer = []
        self.buffered = 0

    def add(self, keys):
        if len(keys):
            self.buffer.append(keys)
            self.buffered += len(keys)
            if self.buffered > ACCUMULATOR_BUFFER:
                self.compress()

    def compress(self):
        if not self.buffer:
            return
        keys = numpy.concatenate([self.keys] + self.buffer)
        counts = numpy.concatenate([self.counts, numpy.ones(self.buffered, dtype=numpy.int64)])
        self.keys, inverse = numpy.unique(keys, return_inverse=True)
        self.counts = numpy.bincount(inverse, weights=counts, minlength=len(self.keys)).astype(numpy.int64)
        self.buffer, self.buffered = [], 0


def count_alt_molecules(bam, chr_id, start, end, molecules, locus_index, ref, args):
    """Pass 1: number of molecules with at least one read carrying each allele."""
    counter = KeyCounter()
    ins_counts = defaultdict(int)
    pending = {}  # multi-read molecule -> [reads seen, int keys, insertion keys]
    for m, aln in locus_alignments(bam, chr_id, start, end, molecules, locus_index):
        r = parse_alignment(aln, ref, args.max_indel_len)
        keys, ins = read_event_keys(r)
        if molecules.group_size[m] == 1:
            counter.add(numpy.unique(keys))
            for k in set(ins):
                ins_counts[k] += 1
            continue
        state = pending.setdefault(m, [0, [], set()])
        state[0] += 1
        state[1].append(keys)
        state[2].update(ins)
        if state[0] == molecules.group_size[m]:
            finish_pending_count(pending.pop(m), counter, ins_counts)
    for state in pending.values():
        # reads of the group that were not aligned within this locus
        finish_pending_count(state, counter, ins_counts)
    counter.compress()
    keep = counter.counts >= args.min_alt_umis
    candidates = [int(k) for k in counter.keys[keep]]
    candidates += [k for k, c in ins_counts.items() if c >= args.min_alt_umis]
    return candidates


def finish_pending_count(state, counter, ins_counts):
    _, key_arrays, ins = state
    counter.add(numpy.unique(numpy.concatenate(key_arrays)))
    for k in ins:
        ins_counts[k] += 1


class Candidates:
    """Candidate alleles of a locus, sorted by position."""
    def __init__(self, keys):
        def position(k):
            return k[0] if isinstance(k, tuple) else k // KEY_SHIFT
        self.keys = sorted(keys, key=lambda k: (position(k), str(k)))
        self.index = {k: i for i, k in enumerate(self.keys)}
        self.pos = numpy.array([position(k) for k in self.keys], dtype=numpy.int64)
        self.is_snv = numpy.array([not isinstance(k, tuple) and k % KEY_SHIFT < DEL_BASE for k in self.keys],
                                  dtype=bool)
        self.snv_base = numpy.array([k % KEY_SHIFT if s else 255 for k, s in zip(self.keys, self.is_snv)],
                                    dtype=numpy.int64)
        self.indel_idx = numpy.flatnonzero(~self.is_snv)
        # window of reference bases a read has to cover: SNV - the base; indels - the anchor and the next base
        self.cover_start = self.pos.copy()
        self.cover_end = self.pos + 1
        for i in self.indel_idx:
            k = self.keys[i]
            length = (k % KEY_SHIFT - DEL_BASE + 1) if not isinstance(k, tuple) else 0
            self.cover_start[i] = self.pos[i] - 1
            self.cover_end[i] = self.pos[i] + length + 1
        self.max_cover = int((self.cover_end - self.pos).max()) if len(self.keys) else 0


def in_blocks(starts, ends, p, p_end=None):
    """Whether positions p lie within one of the sorted blocks [starts, ends), or, with p_end, whether the intervals
    [p, p_end) do."""
    if len(starts) == 0:
        return numpy.zeros(len(p), dtype=bool)
    b = numpy.searchsorted(starts, p, "right") - 1
    block_end = ends[numpy.maximum(b, 0)]
    return (b >= 0) & ((p < block_end) if p_end is None else (p_end <= block_end))


def read_statuses(r, cand):
    """Status of a read at the candidates it overlaps: (candidate indices, statuses)."""
    lo = numpy.searchsorted(cand.pos, r.start - 1, "left")
    hi = numpy.searchsorted(cand.pos, r.end, "right")
    if lo >= hi:
        return None, None
    idx = numpy.arange(lo, hi)
    pos = cand.pos[lo:hi]
    status = numpy.zeros(hi - lo, dtype=numpy.int8)
    m_starts = numpy.asarray(r.m_starts, dtype=numpy.int64)
    m_ends = numpy.asarray(r.m_ends, dtype=numpy.int64)
    s_starts = numpy.asarray(r.span_starts, dtype=numpy.int64)
    s_ends = numpy.asarray(r.span_ends, dtype=numpy.int64)

    snv = cand.is_snv[lo:hi]
    if snv.any():
        p = pos[snv]
        in_m = in_blocks(m_starts, m_ends, p)
        in_span = in_blocks(s_starts, s_ends, p)
        if len(r.snv_pos):
            order = numpy.argsort(r.snv_pos, kind="stable")
            sp, sb = r.snv_pos[order], r.snv_base[order]
            j = numpy.minimum(numpy.searchsorted(sp, p), len(sp) - 1)
            is_mm = sp[j] == p
            is_alt = is_mm & (sb[j] == cand.snv_base[lo:hi][snv])
        else:
            is_mm = numpy.zeros(len(p), dtype=bool)
            is_alt = is_mm
        st = numpy.zeros(len(p), dtype=numpy.int8)
        st[in_span] = ST_OTHER  # inside a deletion
        st[in_m & ~is_mm] = ST_REF
        st[in_m & is_mm] = ST_OTHER
        st[in_m & is_alt] = ST_ALT
        status[snv] = st

    indel = ~snv
    if indel.any():
        # an indel candidate is covered if its window lies in one intron-free span; then it is alt if the read has
        # this indel, other if the read has another indel at the position, ref if both window ends are aligned
        start, end = cand.cover_start[lo:hi][indel], cand.cover_end[lo:hi][indel]
        st = numpy.where(in_blocks(m_starts, m_ends, start) & in_blocks(m_starts, m_ends, end - 1),
                         ST_REF, ST_OTHER).astype(numpy.int8)
        if r.indel_pos:
            read_pos = numpy.sort(numpy.fromiter(r.indel_pos, dtype=numpy.int64, count=len(r.indel_pos)))
            p = pos[indel]
            j = numpy.minimum(numpy.searchsorted(read_pos, p), len(read_pos) - 1)
            st[read_pos[j] == p] = ST_OTHER
        st[~in_blocks(s_starts, s_ends, start, end)] = ST_NONE
        status[indel] = st
        for key in r.indels:
            c = cand.index.get(key)
            if c is not None and lo <= c < hi and status[c - lo] != ST_NONE:
                status[c - lo] = ST_ALT

    keep = status != ST_NONE
    return idx[keep], status[keep]


class SiteAccumulator:
    """Per-candidate molecule totals, per (candidate, barcode) alt UMIs and, with report, per (candidate, barcode)
    depth UMIs, alt reads near alignment ends and alt molecules by strand. At deep loci nearly every
    (candidate, barcode) pair is covered, so these are only collected for the sites that are reported."""
    def __init__(self, n_candidates, n_barcodes, report=False):
        self.n_barcodes = n_barcodes
        self.alt = numpy.zeros(n_candidates, dtype=numpy.int64)
        self.ref = numpy.zeros(n_candidates, dtype=numpy.int64)
        self.no_call = numpy.zeros(n_candidates, dtype=numpy.int64)
        self.alt_counter = KeyCounter()
        self.depth_counter = KeyCounter() if report else None
        if report:
            self.alt_reads = numpy.zeros(n_candidates, dtype=numpy.int64)
            self.alt_reads_near_end = numpy.zeros(n_candidates, dtype=numpy.int64)
            self.alt_plus = numpy.zeros(n_candidates, dtype=numpy.int64)
            self.alt_minus = numpy.zeros(n_candidates, dtype=numpy.int64)

    def add_molecule(self, barcode, cand_idx, mol_status):
        """mol_status: ST_ALT / ST_REF / ST_OTHER (no call) per candidate."""
        alt = mol_status == ST_ALT
        ref = mol_status == ST_REF
        self.alt[cand_idx[alt]] += 1
        self.ref[cand_idx[ref]] += 1
        self.no_call[cand_idx[~(alt | ref)]] += 1
        self.alt_counter.add(cand_idx[alt] * self.n_barcodes + barcode)
        if self.depth_counter is not None:
            self.depth_counter.add(cand_idx[alt | ref] * self.n_barcodes + barcode)


def end_distance(r, p):
    """Aligned reference bases (introns skipped) from positions p to the nearer end of the alignment."""
    starts = numpy.asarray(r.span_starts, dtype=numpy.int64)
    lengths = numpy.asarray(r.span_ends, dtype=numpy.int64) - starts
    before = numpy.cumsum(lengths) - lengths
    k = numpy.maximum(numpy.searchsorted(starts, p, "right") - 1, 0)
    left = before[k] + p - starts[k]
    return numpy.minimum(left, lengths.sum() - 1 - left)


def molecule_consensus(read_results, consensus):
    """Combine read statuses of one molecule: (candidate indices, molecule statuses)."""
    idx = numpy.concatenate([i for i, _ in read_results])
    st = numpy.concatenate([s for _, s in read_results])
    cands, inverse = numpy.unique(idx, return_inverse=True)
    cover = numpy.bincount(inverse, minlength=len(cands))
    alt = numpy.bincount(inverse, weights=(st == ST_ALT), minlength=len(cands))
    ref = numpy.bincount(inverse, weights=(st == ST_REF), minlength=len(cands))
    mol_status = numpy.full(len(cands), ST_OTHER, dtype=numpy.int8)
    mol_status[alt > consensus * cover] = ST_ALT
    mol_status[ref > consensus * cover] = ST_REF
    return cands, mol_status


def genotype_candidates(bam, chr_id, start, end, molecules, locus_index, ref, cand, args, report=False):
    """Molecule consensus at the candidate alleles: pass 2 for all candidates, pass 3 (report=True) for the
    selected ones. The consensus of a candidate does not depend on the other candidates."""
    acc = SiteAccumulator(len(cand.keys), len(molecules.barcode_list), report)
    alt_genes = [defaultdict(int) for _ in cand.keys]
    pending = {}  # multi-read molecule -> [reads seen, list of read results]

    def add(m, cand_idx, mol_status):
        acc.add_molecule(molecules.barcode_idx[m], cand_idx, mol_status)
        alt_idx = cand_idx[mol_status == ST_ALT]
        for c in alt_idx:
            alt_genes[c][molecules.genes[m]] += 1
        if report:
            if molecules.strands[m] == "+":
                acc.alt_plus[alt_idx] += 1
            elif molecules.strands[m] == "-":
                acc.alt_minus[alt_idx] += 1

    for m, aln in locus_alignments(bam, chr_id, start, end, molecules, locus_index):
        r = parse_alignment(aln, ref, args.max_indel_len)
        cand_idx, status = read_statuses(r, cand)
        if report and cand_idx is not None:
            alt_idx = cand_idx[status == ST_ALT]
            if len(alt_idx):
                acc.alt_reads[alt_idx] += 1
                acc.alt_reads_near_end[alt_idx[end_distance(r, cand.pos[alt_idx]) < args.end_distance]] += 1
        if molecules.group_size[m] == 1:
            if cand_idx is not None and len(cand_idx):
                # a single read is its own consensus; a read with another allele is a no call
                add(m, cand_idx, status)
            continue
        state = pending.setdefault(m, [0, []])
        state[0] += 1
        if cand_idx is not None and len(cand_idx):
            state[1].append((cand_idx, status))
        if state[0] == molecules.group_size[m]:
            finish_pending_molecule(m, pending.pop(m), add, args)
    for m, state in pending.items():
        finish_pending_molecule(m, state, add, args)
    acc.alt_counter.compress()
    if report:
        acc.depth_counter.compress()
    return acc, alt_genes


def finish_pending_molecule(m, state, add, args):
    if state[1]:
        cand_idx, mol_status = molecule_consensus(state[1], args.consensus)
        add(m, cand_idx, mol_status)


def homopolymer_length(ref, start):
    """Length of the reference run of identical bases starting at start."""
    if start < ref.start or start >= ref.end:
        return 0
    b = ref.base(start)
    end = start
    while end < ref.end and ref.base(end) == b:
        end += 1
    return end - start


def homopolymer_length_left(ref, end):
    """Length of the reference run of identical bases ending at end (exclusive)."""
    if end <= ref.start or end > ref.end:
        return 0
    b = ref.base(end - 1)
    start = end - 1
    while start > ref.start and ref.base(start - 1) == b:
        start -= 1
    return end - start


def nearby_homopolymer(ref, lo, hi):
    """Length of the longest reference run of identical bases (not N) overlapping [lo, hi)."""
    best = 0
    p = max(lo, ref.start)
    while p < min(hi, ref.end):
        run = homopolymer_length(ref, p)
        if ref.base(p) in BASES:
            best = max(best, homopolymer_length_left(ref, p + 1) + run - 1)
        p += run
    return best


def describe_allele(key, ref):
    """VCF-style (pos, ref, alt, type) and homopolymer length for an allele key."""
    if isinstance(key, tuple):
        s0, seq = key
        anchor = ref.base(s0 - 1)
        return s0, anchor, anchor + seq, "INS", homopolymer_length(ref, s0)
    pos, code = divmod(key, KEY_SHIFT)
    if code < DEL_BASE:
        hp = max(homopolymer_length_left(ref, pos), homopolymer_length(ref, pos + 1))
        return pos + 1, ref.base(pos), BASES[code], "SNV", hp
    length = code - DEL_BASE + 1
    ref_bases = ref.seq[pos - 1 - ref.start:pos + length - ref.start]
    return pos, ref_bases, ref_bases[0], "DEL", homopolymer_length(ref, pos)


def indel_class(vtype, length, hp):
    return "%s:%s:hp%s" % (vtype, "3+" if length >= 3 else length, "6+" if hp >= 6 else hp)


def error_class(vtype, ref_bases, alt_bases, tx_change, hp):
    """SNV: transcript-strand change; indels: type, length (1, 2, 3+) and homopolymer length (1-5, 6+)."""
    if vtype == "SNV":
        return "SNV:" + (tx_change if tx_change != "." else "%s>%s" % (ref_bases, alt_bases))
    return indel_class(vtype, abs(len(alt_bases) - len(ref_bases)), hp)


def class_rate(rates, cls):
    """Rate of the class, or of its closest pooled parent (DEL:3+:hp6+ -> DEL:3+ -> DEL -> ALL)."""
    while cls not in rates and ":" in cls:
        cls = cls.rsplit(":", 1)[0]
    return rates.get(cls, rates["ALL"])


def describe_candidate(key, c, alt_genes, molecules, ref, locus_strand):
    """Allele description, the gene with most alt molecules, its strand, transcript-strand change and error class."""
    pos, ref_bases, alt_bases, vtype, hp = describe_allele(key, ref)
    gene_counts = alt_genes[c]
    gene_id = max(sorted(gene_counts), key=lambda g: gene_counts[g]) if gene_counts else "."
    strand = molecules_gene_strand(molecules, gene_id) if gene_counts else locus_strand
    tx_change = "."
    if vtype == "SNV" and strand in ("+", "-"):
        if strand == "+":
            tx_change = "%s>%s" % (ref_bases, alt_bases)
        else:
            tx_change = "%s>%s" % (COMPLEMENT[ref_bases], COMPLEMENT[alt_bases])
    cls = error_class(vtype, ref_bases, alt_bases, tx_change, hp)
    return pos, ref_bases, alt_bases, vtype, hp, gene_id, gene_counts, strand, tx_change, cls


def locus_strands(molecules):
    """Majority gene strand of the molecules of each locus."""
    counts = defaultdict(lambda: defaultdict(int))
    for loc, s in zip(molecules.locus, molecules.strands):
        counts[loc][s] += 1
    return {loc: max(sorted(c), key=lambda s: c[s]) for loc, c in counts.items()}


def homopolymer_arrays(ref):
    """Run length starting at each position (to the run end), and the full run length at run starts (0 elsewhere)."""
    n = len(ref.codes)
    is_start = numpy.ones(n, dtype=bool)
    is_start[1:] = ref.codes[1:] != ref.codes[:-1]
    starts = numpy.flatnonzero(is_start)
    ends = numpy.append(starts[1:], n)
    run_end = ends[numpy.cumsum(is_start) - 1]
    hp_from = run_end - numpy.arange(n)
    hp_start = numpy.zeros(n, dtype=numpy.int64)
    hp_start[starts] = ends - starts
    hp_start[ref.codes > 3] = 0
    return hp_from, hp_start


def count_read_errors(bam, chr_id, start, end, molecules, locus_index, ref, args, class_counts):
    """Per-read errors at all positions of a locus, from molecules with a single read (no consensus).

    class_counts[class] = [alt, depth, excluded sites]. Alleles seen in at least ESTIMATION_MAX_VAF of the reads at
    their position (and in at least ESTIMATION_MIN_EXCLUDED reads) are likely variants and are left out, with their
    depth. SNV depth is the number of aligned reads
    per transcript-strand reference base; indel depth is the number of reads spanning homopolymer run starts of each
    length (left-aligned indels sit at run starts).
    """
    n = ref.end - ref.start
    m_diff = numpy.zeros((2, n + 1), dtype=numpy.int64)  # aligned bases, by gene strand (+, -)
    span_diff = numpy.zeros(n + 1, dtype=numpy.int64)  # intron-free spans, including deletions
    snv_counter = KeyCounter()
    indel_counts = defaultdict(int)
    for m, aln in locus_alignments(bam, chr_id, start, end, molecules, locus_index):
        strand = molecules.strands[m]
        if molecules.group_size[m] != 1 or strand not in ("+", "-"):
            continue
        s = int(strand == "-")
        r = parse_alignment(aln, ref, args.max_indel_len)
        numpy.add.at(m_diff[s], numpy.asarray(r.m_starts, dtype=numpy.int64) - ref.start, 1)
        numpy.add.at(m_diff[s], numpy.asarray(r.m_ends, dtype=numpy.int64) - ref.start, -1)
        numpy.add.at(span_diff, numpy.asarray(r.span_starts, dtype=numpy.int64) - ref.start, 1)
        numpy.add.at(span_diff, numpy.asarray(r.span_ends, dtype=numpy.int64) - ref.start, -1)
        snv_counter.add(((r.snv_pos.astype(numpy.int64) - ref.start) * 4 + r.snv_base) * 2 + s)
        for k in r.indels:
            indel_counts[k] += 1
    m_cov = numpy.cumsum(m_diff, axis=1)[:, :n]
    span_cov = numpy.cumsum(span_diff)[:n]
    codes = ref.codes.astype(numpy.int64)
    valid = codes < 4

    # SNVs, in transcript-strand bases (the complement of code c is 3 - c)
    depth_by_base = numpy.zeros(4)
    for s in (0, 1):
        tx_codes = codes[valid] if s == 0 else 3 - codes[valid]
        depth_by_base += numpy.bincount(tx_codes, weights=m_cov[s][valid], minlength=4)
    snv_counter.compress()
    keys, counts = snv_counter.keys, snv_counter.counts
    s, base, pos = keys % 2, (keys // 2) % 4, keys // 8
    rc = codes[pos]
    cov = m_cov[s, pos]
    on_ref = rc < 4
    tx_class = numpy.where(s == 1, 3 - rc, rc) * 4 + numpy.where(s == 1, 3 - base, base)
    keep = on_ref & ((counts < ESTIMATION_MAX_VAF * cov) | (counts < ESTIMATION_MIN_EXCLUDED))
    drop = on_ref & ~keep
    alt_by_class = numpy.bincount(tx_class[keep], weights=counts[keep], minlength=16)
    dropped_depth = numpy.bincount(tx_class[drop], weights=cov[drop], minlength=16)
    dropped_sites = numpy.bincount(tx_class[drop], minlength=16)
    for r_code in range(4):
        for a_code in range(4):
            if r_code != a_code:
                c = class_counts["SNV:%s>%s" % (BASES[r_code], BASES[a_code])]
                i = r_code * 4 + a_code
                c[0] += int(alt_by_class[i])
                c[1] += int(depth_by_base[r_code] - dropped_depth[i])
                c[2] += int(dropped_sites[i])

    # indels
    hp_from, hp_start = homopolymer_arrays(ref)
    depth_by_hp = numpy.bincount(numpy.minimum(hp_start, 6), weights=span_cov, minlength=7)
    dropped_depth = defaultdict(int)
    for k, count in indel_counts.items():
        if isinstance(k, tuple):
            s0, length, vtype = k[0], len(k[1]), "INS"
        else:
            s0, code = divmod(k, KEY_SHIFT)
            length, vtype = code - DEL_BASE + 1, "DEL"
        p = s0 - ref.start
        cls = indel_class(vtype, length, int(hp_from[p]))
        if count >= ESTIMATION_MAX_VAF * span_cov[p] and count >= ESTIMATION_MIN_EXCLUDED:
            class_counts[cls][2] += 1
            dropped_depth[cls] += int(span_cov[p])
        else:
            class_counts[cls][0] += count
    for vtype in ("DEL", "INS"):
        for length in (1, 2, 3):
            for hp in range(1, 7):
                cls = indel_class(vtype, length, hp)
                class_counts[cls][1] += int(depth_by_hp[hp]) - dropped_depth.get(cls, 0)


def estimate_error_rates(class_counts):
    """Rate per class; a class with fewer than MIN_CLASS_ALT errors gets the rate of its closest pooled parent
    with enough errors. Returns table rows (class, alt, depth, excluded sites, rate, source)."""
    pooled = defaultdict(lambda: [0, 0, 0])
    for cls, counts in class_counts.items():
        parts = cls.split(":")
        for group in [":".join(parts[:i]) for i in range(1, len(parts))] + ["ALL"]:
            for i in range(3):
                pooled[group][i] += counts[i]
    rates = {group: alt / depth if depth > 0 else 0.0 for group, (alt, depth, _) in pooled.items()}
    rates.setdefault("ALL", 0.0)
    rows = []
    for cls in sorted(class_counts):
        alt, depth, excluded = class_counts[cls]
        if alt >= MIN_CLASS_ALT and depth > 0:
            rate, source = alt / depth, "class"
        else:
            group = cls
            while ":" in group:
                group = group.rsplit(":", 1)[0]
                if pooled[group][0] >= MIN_CLASS_ALT:
                    break
            else:
                group = "ALL"
            rate, source = rates[group], group
        rows.append((cls, alt, depth, excluded, rate, source))
    for group in sorted(pooled):
        alt, depth, excluded = pooled[group]
        rows.append((group, alt, depth, excluded, rates[group], "pooled"))
    return rows


def select_sites(cand, acc, alt_genes, molecules, ref, locus_strand, args, stats):
    """Candidates passing the cut-offs and the error test, in candidate order: (key, description, rate, p-value)."""
    alt_cand = acc.alt_counter.keys // acc.n_barcodes
    selected = []
    for c in numpy.flatnonzero(acc.alt >= args.min_alt_umis).tolist():
        alt_umis = int(acc.alt[c])
        a_lo, a_hi = numpy.searchsorted(alt_cand, [c, c + 1])
        if a_hi - a_lo < args.min_alt_barcodes:
            continue
        depth = alt_umis + int(acc.ref[c])
        if alt_umis / depth < args.min_vaf:
            continue
        stats["tested"] += 1
        description = describe_candidate(cand.keys[c], c, alt_genes, molecules, ref, locus_strand)
        rate = class_rate(args.rates, description[-1])
        pvalue = float(binom.sf(alt_umis - 1, depth, rate))
        if pvalue <= args.max_error_pvalue:
            selected.append((cand.keys[c], description, rate, pvalue))
    return selected


def near_homopolymer(ref, pos, ref_bases, vtype, args):
    # the window holds the anchor base, the deleted bases (none for insertions) and the next base
    return vtype != "SNV" and nearby_homopolymer(ref, pos - HOMOPOLYMER_FLAG_DISTANCE, pos + len(ref_bases) - 1 +
                                                 HOMOPOLYMER_FLAG_DISTANCE) >= args.homopolymer_flag_len


def unreliable_reasons(counts, homopolymer, ref_bases, alt_bases, vtype, tx_change, args):
    """Fraction of alt reads near alignment ends, editing-like (0 / 1) and the reasons to doubt the site.
    counts: alt reads, alt reads near alignment ends, alt molecules of + and of - strand isoforms."""
    alt_reads, near_end, alt_plus, alt_minus = counts
    reasons = []
    end_fraction = near_end / alt_reads if alt_reads else 0.0
    if alt_reads and end_fraction >= args.end_fraction:
        reasons.append("read_end")
    editing = 0
    if vtype == "SNV":
        stranded = alt_plus + alt_minus
        strand = EDITING_CHANGES.get(ref_bases + ">" + alt_bases)
        if stranded == 0:
            editing = int(tx_change in ("A>G", "C>T"))
        elif strand == "+":
            editing = int(alt_plus >= 0.5 * stranded)
        elif strand == "-":
            editing = int(alt_minus >= 0.5 * stranded)
        if editing:
            reasons.append("editing")
    elif homopolymer:
        reasons.append("homopolymer")
    return end_fraction, editing, ",".join(reasons) if reasons else "."


def report_sites(chr_id, sites, acc, molecules, ref, args, sites_out, barcodes_out, stats):
    """Write the selected sites; acc is from genotyping the selected candidates only (same order), with report."""
    n_bc = acc.n_barcodes
    alt_keys, alt_counts = acc.alt_counter.keys, acc.alt_counter.counts
    depth_keys, depth_counts = acc.depth_counter.keys, acc.depth_counter.counts
    alt_cand = alt_keys // n_bc
    depth_cand = depth_keys // n_bc
    for c, (_, description, rate, pvalue) in enumerate(sites):
        pos, ref_bases, alt_bases, vtype, hp, gene_id, gene_counts, strand, tx_change, cls = description
        alt_umis = int(acc.alt[c])
        depth = alt_umis + int(acc.ref[c])
        a_lo, a_hi = numpy.searchsorted(alt_cand, [c, c + 1])
        d_lo, d_hi = numpy.searchsorted(depth_cand, [c, c + 1])
        variant_id = "%s:%d:%s:%s" % (chr_id, pos, ref_bases, alt_bases)
        counts = [int(x[c]) for x in (acc.alt_reads, acc.alt_reads_near_end, acc.alt_plus, acc.alt_minus)]
        end_fraction, editing, unreliable = unreliable_reasons(
            counts, near_homopolymer(ref, pos, ref_bases, vtype, args), ref_bases, alt_bases, vtype, tx_change, args)
        # the raw counts (HIDDEN_COLUMNS) are kept for merging records of the same allele and dropped in main
        sites_out.write("\t".join(map(str, [
            chr_id, pos, variant_id, ref_bases, alt_bases, vtype, gene_id, ",".join(sorted(gene_counts)), strand,
            tx_change, editing, hp, alt_umis, depth, int(acc.no_call[c]), "%.4f" % (alt_umis / depth),
            a_hi - a_lo, d_hi - d_lo, cls, "%.3g" % rate, "%.3g" % pvalue, "%.3f" % end_fraction,
            unreliable] + counts)) + "\n")
        alt_by_bc = dict(zip((alt_keys[a_lo:a_hi] % n_bc).tolist(), alt_counts[a_lo:a_hi].tolist()))
        for bc, d in zip((depth_keys[d_lo:d_hi] % n_bc).tolist(), depth_counts[d_lo:d_hi].tolist()):
            barcodes_out.write("%s\t%s\t%d\t%d\n" % (variant_id, molecules.barcode_list[bc], alt_by_bc.get(bc, 0), d))
        stats["reported"] += 1


def molecules_gene_strand(molecules, gene_id):
    if not hasattr(molecules, "gene_strand"):
        molecules.gene_strand = {}
        for g, s in zip(molecules.genes, molecules.strands):
            molecules.gene_strand.setdefault(g, s)
    return molecules.gene_strand.get(gene_id, ".")


def process_job(job):
    """Call variants in the loci of a job; with estimate=True only count read errors in the selected loci."""
    job_prefix, chr_id, loci, _, selected, estimate, args = job
    molecules = Molecules(job_prefix)
    stats = defaultdict(int)
    class_counts = defaultdict(lambda: [0, 0, 0])  # error class -> alt, depth, excluded sites
    strands = locus_strands(molecules)
    sites_file, barcodes_file = (os.devnull, os.devnull) if estimate else \
        (job_prefix + ".sites.tsv", job_prefix + ".barcodes.tsv")
    with BamSet(args.bam) as bam, pysam.FastaFile(args.reference) as fasta, \
            open(sites_file, "w") as sites_out, open(barcodes_file, "w") as bc_out:
        for locus_index, (start, end) in enumerate(loci):
            if selected is not None and locus_index not in selected:
                continue
            ref = reference_window(bam, chr_id, start, end, molecules, locus_index, fasta, args.max_indel_len)
            stats["loci"] += 1
            if estimate:
                count_read_errors(bam, chr_id, start, end, molecules, locus_index, ref, args, class_counts)
                continue
            candidates = count_alt_molecules(bam, chr_id, start, end, molecules, locus_index, ref, args)
            stats["candidates"] += len(candidates)
            if not candidates:
                continue
            cand = Candidates(candidates)
            acc, alt_genes = genotype_candidates(bam, chr_id, start, end, molecules, locus_index, ref, cand, args)
            sites = select_sites(cand, acc, alt_genes, molecules, ref, strands.get(locus_index, "."), args, stats)
            if not sites:
                continue
            # pass 3: depth per barcode, alt read positions and alt molecule strands for the selected sites only
            site_cand = Candidates([key for key, _, _, _ in sites])
            assert site_cand.keys == [key for key, _, _, _ in sites]
            acc, _ = genotype_candidates(bam, chr_id, start, end, molecules, locus_index, ref, site_cand, args,
                                         report=True)
            report_sites(chr_id, sites, acc, molecules, ref, args, sites_out, bc_out, stats)
    stats["molecules"] = len(molecules.rep_ids)
    return job_prefix, chr_id, dict(stats), dict(class_counts)


# ---------------------------------------------------------------- main

SITES_HEADER = ["#chrom", "pos", "variant_id", "ref", "alt", "type", "gene_id", "alt_genes", "gene_strand",
                "tx_change", "editing_like", "homopolymer_len", "alt_umis", "depth_umis", "no_call_umis", "vaf",
                "alt_barcodes", "covered_barcodes", "error_class", "error_rate", "error_pvalue", "alt_end_fraction",
                "unreliable"]
# raw counts after the site columns of the job outputs, for merging records of the same allele
HIDDEN_COLUMNS = ["alt_reads", "alt_reads_near_end", "alt_plus", "alt_minus"]


def merge_site_records(records, barcode_counts, args):
    """One record from the records of an allele reported by several loci, which count different molecules."""
    col = {name.lstrip("#"): i for i, name in enumerate(SITES_HEADER + HIDDEN_COLUMNS)}
    best = max(records, key=lambda r: int(r[col["alt_umis"]]))
    merged = list(best)

    def total(name):
        return sum(int(r[col[name]]) for r in records)

    alt, depth = total("alt_umis"), total("depth_umis")
    merged[col["alt_genes"]] = ",".join(sorted({g for r in records for g in r[col["alt_genes"]].split(",") if g}))
    merged[col["alt_umis"]], merged[col["depth_umis"]] = alt, depth
    merged[col["no_call_umis"]] = total("no_call_umis")
    merged[col["vaf"]] = "%.4f" % (alt / depth)
    merged[col["alt_barcodes"]] = sum(1 for a, _ in barcode_counts.values() if a)
    merged[col["covered_barcodes"]] = len(barcode_counts)
    rate = class_rate(args.rates, best[col["error_class"]])
    merged[col["error_pvalue"]] = "%.3g" % binom.sf(alt - 1, depth, rate)
    # the homopolymer flag depends on the reference only, so it is the same in every record
    homopolymer = any("homopolymer" in r[col["unreliable"]] for r in records)
    end_fraction, editing, unreliable = unreliable_reasons(
        [total(name) for name in HIDDEN_COLUMNS], homopolymer, best[col["ref"]], best[col["alt"]],
        best[col["type"]], best[col["tx_change"]], args)
    merged[col["alt_end_fraction"]] = "%.3f" % end_fraction
    merged[col["editing_like"]], merged[col["unreliable"]] = editing, unreliable
    return merged[:len(SITES_HEADER)]


def write_outputs(jobs, args):
    """Concatenate the job outputs in job order, merging the records of alleles reported by several loci.
    Returns the number of merged alleles."""
    seen, duplicated = set(), set()
    for job_prefix, _, _, _ in jobs:
        for line in open(job_prefix + ".sites.tsv"):
            vid = line.split("\t", 3)[2]
            (duplicated if vid in seen else seen).add(vid)
    del seen
    records = defaultdict(list)  # duplicated variant id -> its site records (split lines)
    dup_jobs = set()
    for j, (job_prefix, _, _, _) in enumerate(jobs):
        for line in open(job_prefix + ".sites.tsv"):
            v = line.rstrip("\n").split("\t")
            if v[2] in duplicated:
                records[v[2]].append(v)
                dup_jobs.add(j)
    barcode_counts = defaultdict(lambda: defaultdict(lambda: [0, 0]))  # variant -> barcode -> [alt, depth]
    for j in sorted(dup_jobs):
        for line in open(jobs[j][0] + ".barcodes.tsv"):
            vid, barcode, alt, depth = line.rstrip("\n").split("\t")
            if vid in duplicated:
                counts = barcode_counts[vid][barcode]
                counts[0] += int(alt)
                counts[1] += int(depth)

    sites_written, barcodes_written = set(), set()
    n_columns = len(SITES_HEADER)
    with open(args.output + ".sites.tsv", "w") as sites_out, \
            gzip.open(args.output + ".barcode_counts.tsv.gz", "wt", compresslevel=1) as bc_out:
        sites_out.write("\t".join(SITES_HEADER) + "\n")
        bc_out.write("#variant_id\tbarcode\talt_umis\tdepth_umis\n")
        for j, (job_prefix, _, _, _) in enumerate(jobs):
            for line in open(job_prefix + ".sites.tsv"):
                v = line.rstrip("\n").split("\t")
                if v[2] not in duplicated:
                    sites_out.write("\t".join(v[:n_columns]) + "\n")
                elif v[2] not in sites_written:
                    sites_written.add(v[2])
                    sites_out.write("\t".join(map(str, merge_site_records(records[v[2]], barcode_counts[v[2]], args)))
                                    + "\n")
            with open(job_prefix + ".barcodes.tsv") as inf:
                if j not in dup_jobs:
                    shutil.copyfileobj(inf, bc_out)
                    continue
                for line in inf:
                    vid = line.split("\t", 1)[0]
                    if vid not in duplicated:
                        bc_out.write(line)
                    elif vid not in barcodes_written:
                        # all barcodes of a merged allele at its first line, so its lines stay consecutive
                        barcodes_written.add(vid)
                        for barcode, (alt, depth) in barcode_counts[vid].items():
                            bc_out.write("%s\t%s\t%d\t%d\n" % (vid, barcode, alt, depth))
    return len(duplicated)


def load_error_rates(fname):
    rates = {}
    for line in open(fname):
        if line.startswith("#"):
            continue
        v = line.rstrip("\n").split("\t")
        rates[v[0]] = float(v[4])
    rates.setdefault("ALL", 0.0)
    return rates


def run_error_estimation(pool, jobs, args):
    """Count errors at the candidates of random loci holding --error_sample_fraction of molecules; rate table rows."""
    loci = [(j, locus_index, n) for j, (_, _, _, counts) in enumerate(jobs) for locus_index, n in enumerate(counts)]
    total = sum(n for _, _, n in loci)
    rng = numpy.random.default_rng(args.seed)
    selected = defaultdict(set)
    sampled = 0
    for i in rng.permutation(len(loci)):
        if sampled >= args.error_sample_fraction * total:
            break
        j, locus_index, n = loci[i]
        selected[j].add(locus_index)
        sampled += n
    logger.info("Estimating error rates on %d loci with %d of %d molecules" %
                (sum(len(s) for s in selected.values()), sampled, total))
    class_counts = defaultdict(lambda: [0, 0, 0])
    estimation_jobs = [jobs[j] + (selected[j], True, args) for j in sorted(selected)]
    for _, _, _, counts in pool.imap_unordered(process_job, estimation_jobs):
        for cls, values in counts.items():
            for i, v in enumerate(values):
                class_counts[cls][i] += v
    rows = estimate_error_rates(class_counts)
    if not any(r[0] == "ALL" and r[2] > 0 for r in rows):
        logger.warning("No reads for error rate estimation, the error test is disabled")
    logger.info("Error rates: " + ", ".join("%s %.2g" % (r[0], r[4]) for r in rows if r[5] == "pooled"))
    return rows


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--read_info", required=True, help="IsoQuant read_info (.tsv or .tsv.gz)")
    parser.add_argument("--umi_groups", help="IsoQuant UMI groups (--large_output umi_groups); "
                                             "without it every molecule is its representative read")
    parser.add_argument("--bam", required=True, nargs="+",
                        help="indexed BAM(s) IsoQuant was run on (all reads, not deduplicated)")
    parser.add_argument("--reference", required=True, help="reference genome FASTA (indexed)")
    parser.add_argument("-o", "--output", required=True, help="output prefix")
    parser.add_argument("-t", "--threads", type=int, default=16)
    parser.add_argument("--min_alt_umis", type=int, default=10, help="minimal number of alt molecules [%(default)s]")
    parser.add_argument("--min_alt_barcodes", type=int, default=10,
                        help="minimal number of barcodes with an alt molecule [%(default)s]")
    parser.add_argument("--min_vaf", type=float, default=0.005,
                        help="minimal fraction of alt molecules among alt + ref molecules [%(default)s]")
    parser.add_argument("--max_error_pvalue", type=float, default=1e-4,
                        help="maximal binomial p-value of the alt molecules under the class error rate [%(default)s]")
    parser.add_argument("--error_rates", help="error rates per class (PREFIX.error_rates.tsv of an earlier run) "
                                              "instead of estimating them")
    parser.add_argument("--error_sample_fraction", type=float, default=0.05,
                        help="fraction of molecules (random loci) used to estimate error rates [%(default)s]")
    parser.add_argument("--seed", type=int, default=0, help="random seed for the error rate sample [%(default)s]")
    parser.add_argument("--consensus", type=float, default=0.5,
                        help="a molecule is alt (ref) if strictly more than this fraction of its reads covering "
                             "the site are alt (ref) [%(default)s]")
    parser.add_argument("--max_indel_len", type=int, default=10, help="maximal indel length [%(default)s]")
    parser.add_argument("--end_distance", type=int, default=20,
                        help="unreliable read_end: alt within this many aligned bases of the alignment end "
                             "[%(default)s] ...")
    parser.add_argument("--end_fraction", type=float, default=0.5,
                        help="... in at least this fraction of the alt reads [%(default)s]")
    parser.add_argument("--homopolymer_flag_len", type=int, default=5,
                        help="unreliable homopolymer: indel overlapping or next to a reference run of at least "
                             "this length [%(default)s]")
    parser.add_argument("--gene_assignment_types", default=SINGLE_GENE_TYPES,
                        help="read_info gene_assignment_type values to keep [%(default)s]")
    parser.add_argument("--keep_tmp", action="store_true", help="keep temporary files")
    return parser.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", stream=sys.stderr)
    tmp_dir = args.output + ".tmp"
    os.makedirs(tmp_dir, exist_ok=True)

    with BamSet(args.bam) as bam:
        chr_order = {c: i for i, c in enumerate(bam.references)}
    gene_types = set(args.gene_assignment_types.split(","))
    chr_files, gene_spans = split_read_info(args.read_info, gene_types, tmp_dir)
    group_files = split_umi_groups(args.umi_groups, chr_files, tmp_dir)
    jobs = make_jobs(chr_files, group_files, gene_spans, chr_order, tmp_dir, args.threads)
    logger.info("Processing %d loci in %d jobs" % (sum(len(j[2]) for j in jobs), len(jobs)))

    with Pool(args.threads) as pool:
        if args.error_rates:
            args.rates = load_error_rates(args.error_rates)
            logger.info("Loaded error rates for %d classes from %s" % (len(args.rates), args.error_rates))
        else:
            args.rates = None
            rows = run_error_estimation(pool, jobs, args)
            args.rates = {cls: rate for cls, _, _, _, rate, _ in rows}
        with open(args.output + ".error_rates.tsv", "w") as outf:
            outf.write("#error_class\talt\tdepth\texcluded_sites\trate\tsource\n")
            if args.error_rates:
                for cls in sorted(args.rates):
                    outf.write("%s\t.\t.\t.\t%.6g\tloaded\n" % (cls, args.rates[cls]))
            else:
                for row in rows:
                    outf.write("%s\t%d\t%d\t%d\t%.6g\t%s\n" % row)

        totals = defaultdict(int)
        done = 0
        for _, chr_id, stats, _ in pool.imap_unordered(process_job, [j + (None, False, args) for j in jobs]):
            done += 1
            for k, v in stats.items():
                totals[k] += v
            logger.info("Job %d / %d done (%s): %s" % (done, len(jobs), chr_id,
                                                        ", ".join("%s %d" % kv for kv in stats.items())))

    merged = write_outputs(jobs, args)
    logger.info("Total: " + ", ".join("%s %d" % kv for kv in totals.items()))
    logger.info("Alleles reported by two or more loci, merged into one record: %d" % merged)
    logger.info("Saved %s.sites.tsv and %s.barcode_counts.tsv.gz" % (args.output, args.output))
    if not args.keep_tmp:
        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    main()
