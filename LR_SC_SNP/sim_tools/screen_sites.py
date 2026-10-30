#!/usr/bin/env python3
"""Screen phase 1 calls on a simulated dataset: false positives (FP), true positives (TP) and missed true variants
(FN, at least MIN_FN_TEMPLATES alt templates). Reads at each site are traced to their simulated templates (read
names hold the template id and its transcript).

Per site: annotation (expressed transcripts at the site, distance to splice sites and transcript ends, overlapping
genes, nearest truth variant) and read features (allele per primary read, MAPQ, distance to alignment ends and to
introns of the read, whether the read's source transcript has the site in an exon), the molecule counts the caller
would see (read_info / umi_groups), and for true variants the loss of alt templates step by step.

Usage: screen_sites.py DATASET_DIR OUT_PREFIX [THREADS] [SITES_PER_GROUP]
"""
import sys
import os
import glob
import gzip
import random
import bisect
from collections import defaultdict
from multiprocessing import Pool

import numpy
import pysam
from scipy.stats import binom

sys.path.insert(0, "/home/andreyp/ngs-analysis/LR_SC_SNP")
import call_variants as cv

REF_DIR = "/abga/work/andreyp/data/reference/mouse"
GENOME = REF_DIR + "/GRCm39.primary_assembly.genome.fa"
GTF = REF_DIR + "/gencode.vM36.basic.annotation.gtf"
MAX_READS = 20000
MIN_FN_TEMPLATES = 5
SINGLE_GENE_TYPES = set(cv.SINGLE_GENE_TYPES.split(","))
BIN = 100000
FAR = 10 ** 9

dataset = sys.argv[1].rstrip("/")
out_prefix = sys.argv[2]
threads = int(sys.argv[3]) if len(sys.argv) > 3 else 12
per_group = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
P = os.path.basename(dataset)
LR = dataset + "/lr_sc_snp"
BAM_PATH = dataset + "/" + P + ".bam"


def log(msg):
    print(msg, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- annotation

def load_gtf(fname):
    tx = {}
    genes = {}
    for line in open(fname):
        if line.startswith("#"):
            continue
        v = line.rstrip("\n").split("\t")
        if v[2] not in ("transcript", "exon"):
            continue
        attrs = {}
        for a in v[8].split(";"):
            a = a.strip()
            if a and " " in a:
                k, val = a.split(" ", 1)
                attrs.setdefault(k, val.strip('"'))
        tid = attrs["transcript_id"]
        if v[2] == "transcript":
            tx[tid] = [v[0], v[6], [], attrs["gene_id"]]
            genes[attrs["gene_id"]] = (attrs.get("gene_name", "."), attrs.get("gene_type", "."), v[6])
        else:
            tx[tid][2].append((int(v[3]) - 1, int(v[4])))
    for t in tx.values():
        t[2].sort()
    return tx, genes


def load_expressed(fname):
    expr = set()
    for line in open(fname):
        if line.startswith("#"):
            continue
        tid, count = line.split("\t")[:2]
        if float(count) > 0:
            expr.add(tid)
    return expr


log("loading annotation")
TX, GENES = load_gtf(GTF)
EXPRESSED = load_expressed(dataset + "/" + P + ".counts.tsv")
TX_BINS = defaultdict(list)
for tid in EXPRESSED:
    if tid not in TX:
        continue
    chrom, _, exons, _ = TX[tid]
    for b in range(exons[0][0] // BIN, exons[-1][1] // BIN + 1):
        TX_BINS[(chrom, b)].append(tid)


def exon_index(exons, p):
    for i, (s, e) in enumerate(exons):
        if s <= p < e:
            return i
    return -1


def tx_offsets(exons, strand, p, i):
    """Transcript-coordinate distances from p (in exon i) to the 5' and 3' ends."""
    left = sum(e - s for s, e in exons[:i]) + (p - exons[i][0])
    right = sum(e - s for s, e in exons[i + 1:]) + (exons[i][1] - p - 1)
    return (left, right) if strand == "+" else (right, left)


def splice_distance(exons, p):
    """Genomic distance from p to the nearest splice site (internal exon boundary) of a transcript."""
    d = FAR
    for i, (s, e) in enumerate(exons):
        if i > 0:
            d = min(d, abs(p - s))
        if i < len(exons) - 1:
            d = min(d, abs(p - (e - 1)))
    return d


def annotate(chrom, p):
    exonic, spanning = [], []
    for tid in TX_BINS.get((chrom, p // BIN), []):
        _, strand, exons, gene = TX[tid]
        if exons[0][0] <= p < exons[-1][1]:
            spanning.append(tid)
            i = exon_index(exons, p)
            if i >= 0:
                exonic.append((tid, i))
    res = {"expr_tx_exonic": len(exonic), "expr_tx_spanning": len(spanning),
           "genes_spanning": len({TX[t][3] for t in spanning}),
           "strands_spanning": "".join(sorted({TX[t][1] for t in spanning}))}
    d_splice, d5, d3 = FAR, FAR, FAR
    for tid, i in exonic:
        _, strand, exons, _ = TX[tid]
        d_splice = min(d_splice, splice_distance(exons, p))
        a, b = tx_offsets(exons, strand, p, i)
        d5, d3 = min(d5, a), min(d3, b)
    res.update({"dist_splice": d_splice, "dist_tx5": d5, "dist_tx3": d3})
    return res


# ---------------------------------------------------------------- sites

def allele_key(pos, ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return (pos - 1) * cv.KEY_SHIFT + cv.BASES.index(alt), "SNV"
    if len(ref) > len(alt):
        return pos * cv.KEY_SHIFT + cv.DEL_BASE + len(ref) - 2, "DEL"
    return (pos, alt[1:]), "INS"


def read_table(fname):
    rows = []
    header = None
    for line in open(fname):
        v = line.rstrip("\n").split("\t")
        if header is None:
            header = [h.lstrip("#") for h in v]
            continue
        rows.append(dict(zip(header, v)))
    return rows


log("loading sites")
calls = {r["variant_id"]: r for r in read_table(LR + "/calls.sites.tsv")}
eval_sites = read_table(LR + "/calls.eval.sites.tsv")
truth = read_table(LR + "/calls.eval.truth.tsv")
truth_pos = defaultdict(list)
for t in truth:
    truth_pos[t["chrom"]].append(int(t["pos"]))
for c in truth_pos:
    truth_pos[c].sort()
rates = cv.load_error_rates(LR + "/calls.error_rates.tsv")

sites = []  # dicts
for e in eval_sites:
    c = calls[e["variant_id"]]
    chrom, pos = c["chrom"], int(c["pos"])
    key, vtype = allele_key(pos, c["ref"], c["alt"])
    sites.append({"site": e["variant_id"], "group": e["label"], "type": vtype, "chrom": chrom, "pos": pos,
                  "ref": c["ref"], "alt": c["alt"], "key": key, "truth_id": e["truth_id"], "gene_id": c["gene_id"],
                  "tx_change": c["tx_change"], "editing_like": c["editing_like"], "hp": c["homopolymer_len"],
                  "alt_umis": int(c["alt_umis"]), "depth_umis": int(c["depth_umis"]), "vaf": float(c["vaf"]),
                  "alt_barcodes": int(c["alt_barcodes"]), "error_class": c["error_class"],
                  "alt_templates": None, "covering_templates": None})
for t in truth:
    if t["called"] == "1" or int(t["alt_templates"]) < MIN_FN_TEMPLATES:
        continue
    pos = int(t["pos"])
    key, vtype = allele_key(pos, t["ref"], t["alt"])
    sites.append({"site": "%s:%d:%s:%s" % (t["chrom"], pos, t["ref"], t["alt"]), "group": "FN", "type": vtype,
                  "chrom": t["chrom"], "pos": pos, "ref": t["ref"], "alt": t["alt"], "key": key,
                  "truth_id": t["truth_id"], "gene_id": ".", "alt_templates": int(t["alt_templates"]),
                  "covering_templates": int(t["covering_templates"])})
truth_by_id = {t["truth_id"]: t for t in truth}
for s in sites:
    if s["group"] == "TP":
        t = truth_by_id[s["truth_id"]]
        s["alt_templates"], s["covering_templates"] = int(t["alt_templates"]), int(t["covering_templates"])

# annotation for all sites
for s in sites:
    p0 = s["pos"] - 1 if s["type"] == "SNV" else s["pos"]  # SNV base / first base after the indel anchor
    s.update(annotate(s["chrom"], p0))
    tp = truth_pos[s["chrom"]]
    i = bisect.bisect_left(tp, s["pos"])
    d = FAR
    for j in (i - 1, i, i + 1):
        if 0 <= j < len(tp) and not (s["group"] in ("TP", "FN") and tp[j] == s["pos"]):
            d = min(d, abs(tp[j] - s["pos"]))
    s["dist_truth"] = d
    g = GENES.get(s["gene_id"])
    s["gene_type"] = g[1] if g else "."

# read-level screening on a sample per group and type
random.seed(0)
by_group = defaultdict(list)
for i, s in enumerate(sites):
    by_group[(s["group"], s["type"])].append(i)
sample = []
for k, idx in sorted(by_group.items()):
    n = per_group * 3 if k[0] == "FN" else per_group
    sample += idx if len(idx) <= n else random.sample(idx, n)
log("sites: " + ", ".join("%s/%s %d" % (k[0], k[1], len(v)) for k, v in sorted(by_group.items())) +
    "; read-level sample %d" % len(sample))

# alt templates of the true variants in the sample
need_var = {sites[i]["truth_id"]: i for i in sample if sites[i]["group"] in ("TP", "FN")}
alt_tpl = defaultdict(set)
with gzip.open(dataset + "/" + P + ".template_variants.tsv.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        tpl, vid = line.split("\t", 2)[:2]
        if vid in need_var:
            alt_tpl[vid].add(tpl)
log("alt templates loaded for %d variants" % len(alt_tpl))


# ---------------------------------------------------------------- reads

def parse_name(name):
    left = name.split("_aligned_")[0] if "_aligned_" in name else name.split("_unaligned_")[0]
    tpl = left.rsplit("_", 1)[0]
    orient = name.split("_")[-4] if "_aligned_" in name else "."
    return tpl, tpl.split("_")[2], orient


def locate(aln, p):
    """Query distance from p to the nearer alignment end (aligned bases), soft clip on that side, and genomic
    distance from p to the nearest intron boundary of the read."""
    rpos, qpos = aln.reference_start, 0
    ops = aln.cigartuples
    clip_l = ops[0][1] if ops[0][0] == cv.OP_S else 0
    clip_r = ops[-1][1] if ops[-1][0] == cv.OP_S else 0
    qp = None
    junc = FAR
    for op, length in ops:
        if op in cv.ALIGNED_OPS:
            if rpos <= p < rpos + length:
                qp = qpos + p - rpos
            rpos += length
            qpos += length
        elif op == cv.OP_I or op == cv.OP_S:
            qpos += length
        elif op == cv.OP_D:
            if rpos <= p < rpos + length:
                qp = qpos
            rpos += length
        elif op == cv.OP_N:
            junc = min(junc, abs(p - rpos), abs(p - (rpos + length)))
            rpos += length
    if qp is None:
        return FAR, 0, junc
    qlen = aln.query_length
    dl, dr = qp - clip_l, qlen - clip_r - 1 - qp
    return (dl, clip_l, junc) if dl <= dr else (dr, clip_r, junc)


def worker_init():
    global BAM, FASTA
    BAM = pysam.AlignmentFile(BAM_PATH)
    FASTA = pysam.FastaFile(GENOME)


def screen(i):
    s = sites[i]
    cand = cv.Candidates([s["key"]])
    p = int(cand.pos[0])
    cs, ce = int(cand.cover_start[0]), int(cand.cover_end[0])
    tpls = alt_tpl.get(s["truth_id"], set())
    alns, n_sec, n_supp = [], 0, 0
    for aln in BAM.fetch(s["chrom"], cs, ce):
        if aln.is_unmapped or aln.cigartuples is None:
            continue
        if aln.is_secondary:
            n_sec += 1
            continue
        if aln.is_supplementary:
            n_supp += 1
            continue
        alns.append(aln)
    capped = len(alns) > MAX_READS
    if capped:
        keep_alt = [a for a in alns if parse_name(a.query_name)[0] in tpls]
        rest = [a for a in alns if parse_name(a.query_name)[0] not in tpls]
        alns = keep_alt + random.Random(i).sample(rest, min(len(rest), MAX_READS))
    if not alns:
        return i, [], n_sec, n_supp, capped
    lo = min(a.reference_start for a in alns)
    hi = max(a.reference_end for a in alns)
    ref = cv.Reference(FASTA, s["chrom"], lo - 1000, hi + 1000)
    reads = []
    for aln in alns:
        r = cv.parse_alignment(aln, ref, 10)
        idx, st = cv.read_statuses(r, cand)
        status = int(st[0]) if idx is not None and len(idx) else cv.ST_NONE
        tpl, tid, orient = parse_name(aln.query_name)
        qdist, clip, junc = locate(aln, p)
        t = TX.get(tid)
        if t is not None and t[0] == s["chrom"]:
            src_exonic = int(exon_index(t[2], p) >= 0)
            src_gene = t[3]
        else:
            src_exonic, src_gene = 0, "."
        nm = aln.get_tag("NM") if aln.has_tag("NM") else -1
        aligned = aln.query_alignment_length
        reads.append((aln.query_name, status, aln.mapping_quality, tpl, int(tpl in tpls), orient, qdist, clip, junc,
                      src_exonic, src_gene, nm / aligned if aligned and nm >= 0 else -1))
    return i, reads, n_sec, n_supp, capped


FASTA_MAIN = pysam.FastaFile(GENOME)
log("screening reads at %d sites" % len(sample))
results = {}
with Pool(threads, initializer=worker_init) as pool:
    for n, (i, reads, n_sec, n_supp, capped) in enumerate(pool.imap_unordered(screen, sample, chunksize=4)):
        results[i] = (reads, n_sec, n_supp, capped)
        if (n + 1) % 500 == 0:
            log("  %d / %d sites" % (n + 1, len(sample)))

# molecules: read -> representative (umi_groups), representative -> gene, assignment type, barcode (read_info)
needed = {r[0] for reads, _, _, _ in results.values() for r in reads}
log("reads at sampled sites: %d" % len(needed))
rep_of = {}
umi_groups = glob.glob(dataset + "/isoquant/*/*.UMI_filtered.ED*.umi_groups.tsv.gz")[0]
with gzip.open(umi_groups, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        rep, _, size, members = line.rstrip("\n").split("\t")
        if size == "1":
            if rep in needed:
                rep_of[rep] = rep
            continue
        for m in members.split(","):
            if m in needed:
                rep_of[m] = rep
need_rep = set(rep_of.values())
rep_info = {}
read_info = glob.glob(dataset + "/isoquant/*/*.read_info.tsv.gz")[0]
with gzip.open(read_info, "rt") as f:
    cols = None
    for line in f:
        if line.startswith("#"):
            continue
        v = line.rstrip("\n").split("\t")
        if cols is None:
            cols = {name: v.index(name) for name in ("read_id", "gene_id", "gene_assignment_type", "barcode")}
            continue
        rid = v[cols["read_id"]]
        if rid in need_rep and rid not in rep_info:
            rep_info[rid] = (v[cols["gene_id"]], v[cols["gene_assignment_type"]], v[cols["barcode"]])
log("representatives %d, in read_info %d" % (len(need_rep), len(rep_info)))


def kept_state(read_id):
    """'kept' (molecule used by the caller), 'no_molecule' (not in umi_groups), 'no_read_info', 'multi_gene',
    'no_barcode'."""
    rep = rep_of.get(read_id)
    if rep is None:
        return "no_molecule", None
    info = rep_info.get(rep)
    if info is None:
        return "no_read_info", rep
    gene, gtype, barcode = info
    if gtype not in SINGLE_GENE_TYPES:
        return "multi_gene", rep
    if barcode in cv.NO_BARCODE:
        return "no_barcode", rep
    return "kept", rep


def med(x):
    return float(numpy.median(x)) if len(x) else float("nan")


def frac(x):
    return float(numpy.mean(x)) if len(x) else float("nan")


out_reads = open(out_prefix + ".alt_reads.tsv", "w")
out_reads.write("site\tgroup\ttype\tread_id\tmapq\ttemplate_alt\torient\tqdist_end\tclip\tdist_intron\tsrc_exonic\t"
                "src_gene_match\tread_err\tkept\n")
for i, (reads, n_sec, n_supp, capped) in results.items():
    s = sites[i]
    s.update({"screened": 1, "reads": len(reads), "secondary": n_sec, "supplementary": n_supp, "capped": int(capped)})
    by_status = defaultdict(list)
    for r in reads:
        by_status[r[1]].append(r)
    alt, refr = by_status[cv.ST_ALT], by_status[cv.ST_REF]
    s["alt_reads"], s["ref_reads"], s["other_reads"] = len(alt), len(refr), len(by_status[cv.ST_OTHER])
    if s["gene_id"] == ".":
        # missed variant: the gene of most alt-template reads
        genes = [r[10] for r in reads if r[4] and r[10] != "."]
        s["gene_id"] = max(sorted(set(genes)), key=genes.count) if genes else "."
        g = GENES.get(s["gene_id"])
        s["gene_type"] = g[1] if g else "."
    gene = s["gene_id"]
    for name, rr in (("alt", alt), ("ref", refr)):
        s[name + "_mapq_lt60"] = frac([r[2] < 60 for r in rr])
        s[name + "_mapq_lt20"] = frac([r[2] < 20 for r in rr])
        s[name + "_qdist_med"] = med([r[6] for r in rr])
        s[name + "_near_end20"] = frac([r[6] < 20 for r in rr])
        s[name + "_near_intron10"] = frac([r[8] < 10 for r in rr])
        s[name + "_src_exonic"] = frac([r[9] for r in rr])
        s[name + "_src_gene"] = frac([r[10] == gene for r in rr]) if gene != "." else float("nan")
        s[name + "_fwd"] = frac([r[5] == "F" for r in rr])
        s[name + "_read_err"] = med([r[11] for r in rr if r[11] >= 0])
        s[name + "_from_alt_tpl"] = frac([r[4] for r in rr])
    for r in alt:
        st, _ = kept_state(r[0])
        out_reads.write("\t".join(map(str, [s["site"], s["group"], s["type"], r[0], r[2], r[4], r[5], r[6], r[7], r[8],
                                            r[9], int(r[10] == gene), "%.4f" % r[11], st])) + "\n")

    # molecules the caller would count (kept reads grouped by representative, >50% consensus)
    mol = defaultdict(lambda: [0, 0, 0, None])  # rep -> [covering reads, alt, ref, barcode]
    lost = defaultdict(int)
    for r in reads:
        if r[1] == cv.ST_NONE:
            continue
        st, rep = kept_state(r[0])
        if st != "kept":
            lost[st] += 1
            continue
        m = mol[rep]
        m[0] += 1
        m[1] += r[1] == cv.ST_ALT
        m[2] += r[1] == cv.ST_REF
        m[3] = rep_info[rep][2]
    alt_m = [rep for rep, m in mol.items() if m[1] > 0.5 * m[0]]
    ref_m = [rep for rep, m in mol.items() if m[2] > 0.5 * m[0]]
    s["est_alt_umis"], s["est_depth_umis"] = len(alt_m), len(alt_m) + len(ref_m)
    s["est_alt_barcodes"] = len({mol[rep][3] for rep in alt_m})
    for k in ("no_molecule", "no_read_info", "multi_gene", "no_barcode"):
        s["covering_reads_" + k] = lost[k]
    cls = s.get("error_class")
    if cls is None:
        # class of a missed variant: SNV by transcript-strand change of the source gene, indels by length and hp
        genes = [r[10] for r in alt if r[10] != "."]
        strand = GENES[max(set(genes), key=genes.count)][2] if genes else "."
        if s["type"] == "SNV":
            ch = "%s>%s" % (s["ref"], s["alt"])
            if strand == "-":
                ch = "%s>%s" % (cv.COMPLEMENT[s["ref"]], cv.COMPLEMENT[s["alt"]])
            cls = "SNV:" + ch
        else:
            ref_win = cv.Reference(FASTA_MAIN, s["chrom"], s["pos"] - 50, s["pos"] + 50)
            hp_len = cv.homopolymer_length(ref_win, s["pos"])
            cls = cv.indel_class(s["type"], abs(len(s["alt"]) - len(s["ref"])), hp_len)
        s["error_class"] = cls
    d = s["est_depth_umis"]
    a = s["est_alt_umis"]
    s["est_pvalue"] = float(binom.sf(a - 1, d, cv.class_rate(rates, cls))) if d else float("nan")

    if s["group"] in ("TP", "FN"):
        # alt templates: covered at the site, alt by template consensus, in a kept molecule that is alt
        tpl_reads = defaultdict(list)
        for r in reads:
            if r[4] and r[1] != cv.ST_NONE:
                tpl_reads[r[3]].append(r)
        s["tpl_covered"] = len(tpl_reads)
        s["tpl_alt"] = sum(sum(r[1] == cv.ST_ALT for r in rr) > 0.5 * len(rr) for rr in tpl_reads.values())
        s["tpl_kept"] = sum(any(kept_state(r[0])[0] == "kept" for r in rr) for rr in tpl_reads.values())
        alt_set = set(alt_m)
        s["tpl_in_alt_mol"] = sum(any(kept_state(r[0])[1] in alt_set for r in rr) for rr in tpl_reads.values())
out_reads.close()

columns = []
for s in sites:
    for k in s:
        if k not in columns and k != "key":
            columns.append(k)
with open(out_prefix + ".sites.tsv", "w") as f:
    f.write("\t".join(columns) + "\n")
    for s in sites:
        f.write("\t".join(str(s.get(k, "")) for k in columns) + "\n")
log("saved %s.sites.tsv and %s.alt_reads.tsv" % (out_prefix, out_prefix))
