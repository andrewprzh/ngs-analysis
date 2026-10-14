############################################################################
# Copyright (c) 2026 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

# Generates error-free barcoded templates (as simulate_barcoded.py) with somatic variants
# (SNVs, insertions, deletions) introduced into the cDNA part.
# Variants are placed in genomic coordinates within the exons of expressed transcripts, so all isoforms covering
# a variant carry it. Each variant is assigned to a random subset of barcodes (carriers); every template from a carrier
# barcode that covers the variant carries the alternative allele with probability --allele_fraction.
#
# Notes:
#  - indels are reported at the sampled position, run `bcftools norm -f genome` before comparing to a caller;
#  - molecule counts in the ground truth are template-level, read simulators may truncate or drop templates;
#  - sites are uniform over the expressed exonic territory, no enrichment for homopolymers / STRs.

import sys
import re
import argparse
import bisect
import random
from collections import defaultdict
from traceback import print_exc

import numpy
import pysam

from simulate_barcoded import (create_template_10x, create_template_visiumhd, load_counts, load_barcodes,
                               get_random_seq, reverese_complement, BARCODE_LEN_10X, UMI_LEN_10X, UMI_LEN_10X_V2,
                               VISIUM_HD_BC_LEN, VISIUM_HD_UMI_LEN)


MODES = {"10x": (BARCODE_LEN_10X, UMI_LEN_10X, create_template_10x),
         "10x_v2": (BARCODE_LEN_10X, UMI_LEN_10X_V2, create_template_10x),
         "visium_hd": (VISIUM_HD_BC_LEN, VISIUM_HD_UMI_LEN, create_template_visiumhd)}

NUCLS = "ACGT"
SNV, INS, DEL = "SNV", "INS", "DEL"
BIN_SIZE = 10000
MAX_ATTEMPTS_PER_VARIANT = 1000


class Transcript:
    def __init__(self, transcript_id, chrom, strand, gene_id, gene_name):
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.strand = strand
        self.gene_id = gene_id
        self.gene_name = gene_name
        # 0-based half-open, sorted by start after finalize()
        self.exons = []
        self.cumulative = []
        self.length = 0

    def finalize(self):
        self.exons.sort()
        self.cumulative = []
        self.length = 0
        for s, e in self.exons:
            self.cumulative.append(self.length)
            self.length += e - s

    def genomic_offset(self, exon_index, pos):
        # offset in the genomic-orientation transcript sequence
        return self.cumulative[exon_index] + pos - self.exons[exon_index][0]


class Variant:
    def __init__(self, chrom, pos, vtype, ref, alt):
        # pos is 0-based: the substituted base, the first deleted base, or the base before which the insertion goes
        self.chrom = chrom
        self.pos = pos
        self.vtype = vtype
        self.ref = ref
        self.alt = alt
        self.index = -1
        self.variant_id = ""
        self.carrier_fraction = 0.0
        self.carriers = None
        self.transcripts = []
        self.alt_templates = 0
        self.covering_templates = 0

    def length(self):
        if self.vtype == SNV:
            return 1
        elif self.vtype == DEL:
            return len(self.ref)
        return len(self.alt)

    def check_span(self):
        # region that must be fully inside a single exon of every transcript overlapping it;
        # for indels one flanking base is required on each side
        if self.vtype == SNV:
            return self.pos, self.pos + 1
        elif self.vtype == DEL:
            return self.pos - 1, self.pos + len(self.ref) + 1
        return self.pos - 1, self.pos + 1

    def vcf_record(self, genome):
        if self.vtype == SNV:
            return self.pos + 1, self.ref, self.alt
        anchor = genome.fetch(self.chrom, self.pos - 1, self.pos).upper()
        if self.vtype == DEL:
            return self.pos, anchor + self.ref, anchor
        return self.pos, anchor, anchor + self.alt

    def transcript_pos(self, transcript, genomic_offset):
        # 0-based position of the first affected reference base in transcript orientation
        # (for insertions: the transcript base before which the sequence is inserted)
        if transcript.strand == '+':
            return genomic_offset
        if self.vtype == SNV:
            return transcript.length - genomic_offset - 1
        elif self.vtype == DEL:
            return transcript.length - genomic_offset - len(self.ref)
        return transcript.length - genomic_offset


def get_attribute(attributes, name):
    m = re.search(name + ' "([^"]+)"', attributes)
    return m.group(1) if m else None


def load_transcripts(gtf, transcript_ids):
    transcripts = {}
    for l in open(gtf):
        if l.startswith("#"):
            continue
        v = l.rstrip("\n").split("\t")
        if len(v) < 9 or v[2] != "exon":
            continue
        t_id = get_attribute(v[8], "transcript_id")
        if t_id not in transcript_ids:
            continue
        if t_id not in transcripts:
            gene_id = get_attribute(v[8], "gene_id")
            gene_name = get_attribute(v[8], "gene_name")
            transcripts[t_id] = Transcript(t_id, v[0], v[6], gene_id, gene_name if gene_name else gene_id)
        transcripts[t_id].exons.append((int(v[3]) - 1, int(v[4])))
    for t in transcripts.values():
        t.finalize()
    return transcripts


def build_exon_index(transcripts):
    exon_bins = defaultdict(lambda: defaultdict(list))
    intervals = defaultdict(list)
    for t in transcripts.values():
        for i, (s, e) in enumerate(t.exons):
            intervals[t.chrom].append((s, e))
            for b in range(s // BIN_SIZE, (e - 1) // BIN_SIZE + 1):
                exon_bins[t.chrom][b].append((s, e, t.transcript_id, i))
    return exon_bins, intervals


def merge_intervals(intervals):
    territory = []
    for chrom in sorted(intervals.keys()):
        merged = []
        for s, e in sorted(intervals[chrom]):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        territory += [(chrom, s, e) for s, e in merged]
    return territory


def overlapping_exons(exon_bins, chrom, start, end):
    found = set()
    chrom_bins = exon_bins.get(chrom, {})
    for b in range(start // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
        for exon in chrom_bins.get(b, []):
            if exon[0] < end and start < exon[1]:
                found.add(exon)
    return found


class VariantSampler:
    def __init__(self, args, genome, template_counts, exon_bins, territory):
        self.args = args
        self.genome = genome
        self.template_counts = template_counts
        self.exon_bins = exon_bins
        self.territory = territory
        self.territory_cumulative = []
        total = 0
        for chrom, s, e in territory:
            self.territory_cumulative.append(total)
            total += e - s
        self.territory_length = total
        # chrom -> sorted list of occupied (start, end) spans
        self.occupied = defaultdict(list)

    def random_position(self):
        offset = random.randrange(self.territory_length)
        i = bisect.bisect_right(self.territory_cumulative, offset) - 1
        chrom, s, e = self.territory[i]
        return chrom, s + offset - self.territory_cumulative[i]

    def indel_length(self):
        while True:
            length = int(numpy.random.geometric(self.args.indel_geom_p))
            if length <= self.args.max_indel_len:
                return length

    def is_isolated(self, chrom, start, end):
        pad = self.args.max_indel_len + 1
        spans = self.occupied[chrom]
        i = bisect.bisect_left(spans, (start, end))
        if i > 0 and spans[i - 1][1] + pad > start:
            return False
        if i < len(spans) and spans[i][0] < end + pad:
            return False
        return True

    def covering_transcripts(self, variant):
        # returns list of (transcript_id, exon_index) or None if the span crosses an exon boundary of some transcript
        start, end = variant.check_span()
        covering = []
        for s, e, t_id, exon_index in overlapping_exons(self.exon_bins, variant.chrom, start, end):
            if not (s <= start and end <= e):
                return None
            covering.append((t_id, exon_index))
        return sorted(covering)

    def try_variant(self, vtype):
        chrom, pos = self.random_position()
        if vtype == SNV:
            ref = self.genome.fetch(chrom, pos, pos + 1).upper()
            if ref not in NUCLS:
                return None
            variant = Variant(chrom, pos, SNV, ref, random.choice([n for n in NUCLS if n != ref]))
        elif vtype == INS:
            variant = Variant(chrom, pos, INS, "", get_random_seq(self.indel_length()))
        else:
            length = self.indel_length()
            variant = Variant(chrom, pos, DEL, self.genome.fetch(chrom, pos, pos + length).upper(), "")

        start, end = variant.check_span()
        if start < 0 or not self.is_isolated(chrom, start, end):
            return None
        covering = self.covering_transcripts(variant)
        if not covering:
            return None
        if any(c not in NUCLS for c in self.genome.fetch(chrom, start, end).upper()):
            return None
        if sum(self.template_counts[t_id] for t_id, _ in covering) < self.args.min_site_templates:
            return None

        variant.transcripts = covering
        bisect.insort(self.occupied[chrom], (start, end))
        return variant

    def sample(self, vtype, count):
        variants = []
        attempts = 0
        while len(variants) < count and attempts < count * MAX_ATTEMPTS_PER_VARIANT:
            attempts += 1
            variant = self.try_variant(vtype)
            if variant is not None:
                variants.append(variant)
        if len(variants) < count:
            print("Warning: placed only %d of %d %s variants" % (len(variants), count, vtype))
        return variants


def assign_carriers(variants, barcode_count, args):
    for v in variants:
        v.carrier_fraction = random.uniform(args.min_carrier_fraction, args.max_carrier_fraction)
        carrier_count = max(1, round(v.carrier_fraction * barcode_count))
        v.carriers = numpy.zeros(barcode_count, dtype=bool)
        v.carriers[random.sample(range(barcode_count), carrier_count)] = True


def transcript_sequence(genome, transcript):
    seq = "".join(genome.fetch(transcript.chrom, s, e) for s, e in transcript.exons).upper()
    return re.sub("[^ACGTN]", "N", seq)


def apply_variants(seq, edits):
    # edits are (genomic offset, variant), applied right to left so that offsets stay valid
    for offset, v in sorted(edits, key=lambda x: x[0], reverse=True):
        if v.vtype == SNV:
            seq = seq[:offset] + v.alt + seq[offset + 1:]
        elif v.vtype == DEL:
            seq = seq[:offset] + seq[offset + len(v.ref):]
        else:
            seq = seq[:offset] + v.alt + seq[offset:]
    return seq


def write_vcf(fname, variants, transcripts, genome, args):
    contig_order = {c: i for i, c in enumerate(genome.references)}
    with open(fname, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=simulate_barcoded_mutations.py\n")
        vcf.write("##reference=%s\n" % args.genome)
        for c, l in zip(genome.references, genome.lengths):
            vcf.write("##contig=<ID=%s,length=%d>\n" % (c, l))
        vcf.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: SNV, INS or DEL">\n')
        vcf.write('##INFO=<ID=LEN,Number=1,Type=Integer,Description="Number of substituted, inserted or deleted bases">\n')
        vcf.write('##INFO=<ID=CF,Number=1,Type=Float,Description="Fraction of barcodes carrying the variant">\n')
        vcf.write('##INFO=<ID=NCARRIERS,Number=1,Type=Integer,Description="Number of barcodes carrying the variant">\n')
        vcf.write('##INFO=<ID=TRANSCRIPTS,Number=.,Type=String,Description="Expressed transcripts covering the variant">\n')
        vcf.write('##INFO=<ID=GENES,Number=.,Type=String,Description="Genes of the covering transcripts">\n')
        vcf.write('##INFO=<ID=ALT_MOL,Number=1,Type=Integer,Description="Templates carrying the variant">\n')
        vcf.write('##INFO=<ID=COV_MOL,Number=1,Type=Integer,Description="Templates covering the variant">\n')
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in sorted(variants, key=lambda x: (contig_order[x.chrom], x.pos)):
            pos, ref, alt = v.vcf_record(genome)
            t_ids = [t_id for t_id, _ in v.transcripts]
            genes = sorted(set(transcripts[t_id].gene_name for t_id in t_ids))
            info = "TYPE=%s;LEN=%d;CF=%.4f;NCARRIERS=%d;TRANSCRIPTS=%s;GENES=%s;ALT_MOL=%d;COV_MOL=%d" % \
                   (v.vtype, v.length(), v.carrier_fraction, int(v.carriers.sum()), ",".join(t_ids), ",".join(genes),
                    v.alt_templates, v.covering_templates)
            vcf.write("%s\t%d\t%s\t%s\t%s\t.\tPASS\t%s\n" % (v.chrom, pos, v.variant_id, ref, alt, info))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="Barcoded templates with somatic SNVs and indels in the cDNA part. "
                                                 "Indels are reported at the sampled position: run bcftools norm "
                                                 "before comparing to a caller. Molecule counts are template-level.")
    parser.add_argument("--genome", help="reference genome (indexed FASTA)", required=True, type=str)
    parser.add_argument("--gtf", help="gene annotation in GTF format", required=True, type=str)
    parser.add_argument("--counts", "-c", help="expression table", type=str, required=True)
    parser.add_argument("--barcodes", "-b", help="barcode whitelist; for Visium HD pass a subset", type=str, required=True)
    parser.add_argument("--output", "-o", help="output prefix", type=str, required=True)
    parser.add_argument("--mode", help="[%s]" % " | ".join(MODES.keys()), type=str, default="10x")
    parser.add_argument("--template_count", help="total number of templates to generate", type=int, default=100000)
    parser.add_argument("--sub_rate", help="per-base substitution rate over expressed exons", type=float, default=0.0)
    parser.add_argument("--ins_rate", help="per-base insertion rate over expressed exons", type=float, default=0.0)
    parser.add_argument("--del_rate", help="per-base deletion rate over expressed exons", type=float, default=0.0)
    parser.add_argument("--indel_geom_p", help="parameter of geometric indel length distribution", type=float, default=0.5)
    parser.add_argument("--max_indel_len", help="maximal indel length", type=int, default=10)
    parser.add_argument("--min_carrier_fraction", help="minimal fraction of barcodes carrying a variant", type=float, default=0.01)
    parser.add_argument("--max_carrier_fraction", help="maximal fraction of barcodes carrying a variant", type=float, default=0.5)
    parser.add_argument("--allele_fraction", help="probability that a template from a carrier barcode has the variant",
                        type=float, default=0.5)
    parser.add_argument("--min_site_templates", help="minimal number of templates covering a variant site",
                        type=int, default=0)
    parser.add_argument("--seed", help="random seed", type=int, default=0)

    args = parser.parse_args()
    if args.mode not in MODES:
        print("Unknown mode %s" % args.mode)
        exit(-1)
    if not 0 < args.indel_geom_p <= 1 or args.max_indel_len < 1:
        print("Invalid indel length parameters")
        exit(-1)
    if not 0 <= args.min_carrier_fraction <= args.max_carrier_fraction <= 1:
        print("Invalid carrier fractions")
        exit(-1)
    return args


def main():
    args = parse_args()
    random.seed(args.seed)
    numpy.random.seed(args.seed)
    bc_len, umi_len, template_func = MODES[args.mode]

    genome = pysam.FastaFile(args.genome)
    barcodes = load_barcodes(args.barcodes)
    print("Using %d reference barcodes from %s" % (len(barcodes), args.barcodes))
    count_dict = {t_id: c for t_id, c in load_counts(args.counts, args.template_count).items() if c > 0}
    print("Loaded %d counts from %s" % (len(count_dict), args.counts))
    transcripts = load_transcripts(args.gtf, set(count_dict.keys()))
    missing = [t_id for t_id in count_dict if t_id not in transcripts]
    if missing:
        print("Warning: %d transcripts from the counts table are missing in %s, e.g. %s" %
              (len(missing), args.gtf, missing[0]))
        for t_id in missing:
            del count_dict[t_id]
    print("Loaded %d expressed transcripts from %s" % (len(transcripts), args.gtf))

    exon_bins, intervals = build_exon_index(transcripts)
    territory = merge_intervals(intervals)
    sampler = VariantSampler(args, genome, count_dict, exon_bins, territory)
    print("Expressed exonic territory: %d bp" % sampler.territory_length)

    variants = []
    for vtype, rate in ((SNV, args.sub_rate), (INS, args.ins_rate), (DEL, args.del_rate)):
        count = round(rate * sampler.territory_length)
        sampled = sampler.sample(vtype, count)
        print("Placed %d %s variants" % (len(sampled), vtype))
        variants += sampled
    for i, v in enumerate(variants):
        v.index = i
        v.variant_id = "VAR_%d" % i
    assign_carriers(variants, len(barcodes), args)

    # transcript -> list of (genomic offset, variant)
    transcript_variants = defaultdict(list)
    for v in variants:
        for t_id, exon_index in v.transcripts:
            transcript_variants[t_id].append((transcripts[t_id].genomic_offset(exon_index, v.pos), v))

    scale_factor = 1000000.0 / sum(count_dict.values())
    counts_file = args.output + ".counts.tsv"
    with open(counts_file, "w") as expr_table:
        expr_table.write("#transcript_id\tcounts\ttpm\n")
        for t_id in count_dict:
            expr_table.write("%s\t%d\t%.4f\n" % (t_id, count_dict[t_id], scale_factor * count_dict[t_id]))
    print("Saved templates counts to %s" % counts_file)

    # (variant index, barcode index) -> [alt templates, covering templates]
    cell_counts = defaultdict(lambda: [0, 0])
    template_fname = args.output + ".templates.fasta"
    output_fasta = open(template_fname, "w")
    output_counts = open(args.output + ".template_counts.tsv", "w")
    output_template_variants = open(args.output + ".template_variants.tsv", "w")
    output_template_variants.write("#template_id\tvariant_id\ttranscript_id\ttranscript_pos\ttype\n")
    mutated_templates = 0
    for t_id in count_dict:
        transcript = transcripts[t_id]
        reference_seq = transcript_sequence(genome, transcript)
        assert len(reference_seq) == transcript.length
        for i in range(count_dict[t_id]):
            bc_index = random.randrange(len(barcodes))
            bc = barcodes[bc_index]
            umi = get_random_seq(umi_len, 0)
            seq_id = "READ_%d_%s_%s_%s_+" % (i, t_id, bc, umi)

            edits = []
            for offset, v in transcript_variants.get(t_id, []):
                counts = cell_counts[(v.index, bc_index)]
                counts[1] += 1
                v.covering_templates += 1
                if v.carriers[bc_index] and random.random() < args.allele_fraction:
                    counts[0] += 1
                    v.alt_templates += 1
                    edits.append((offset, v))
                    output_template_variants.write("%s\t%s\t%s\t%d\t%s\n" % (seq_id, v.variant_id, t_id,
                                                                            v.transcript_pos(transcript, offset), v.vtype))
            seq = apply_variants(reference_seq, edits)
            if transcript.strand == '-':
                seq = reverese_complement(seq)
            if edits:
                mutated_templates += 1

            output_fasta.write(">%s\n%s\n" % (seq_id, template_func(seq, bc, umi)))
            output_counts.write("%s\t%d\t%.4f\n" % (seq_id, 1, scale_factor))
    output_fasta.close()
    output_counts.close()
    output_template_variants.close()
    print("Saved %d templates (%d with variants) to %s" % (sum(count_dict.values()), mutated_templates, template_fname))

    vcf_fname = args.output + ".variants.vcf"
    write_vcf(vcf_fname, variants, transcripts, genome, args)
    print("Saved %d variants to %s" % (len(variants), vcf_fname))

    with open(args.output + ".variant_carriers.tsv", "w") as outf:
        outf.write("#variant_id\tbarcodes\n")
        for v in variants:
            outf.write("%s\t%s\n" % (v.variant_id, ",".join(barcodes[i] for i in numpy.flatnonzero(v.carriers))))

    with open(args.output + ".variant_cell_counts.tsv", "w") as outf:
        outf.write("#variant_id\tbarcode\tis_carrier\talt_templates\tcovering_templates\n")
        for (v_index, bc_index), (alt, cov) in sorted(cell_counts.items()):
            v = variants[v_index]
            outf.write("%s\t%s\t%d\t%d\t%d\n" % (v.variant_id, barcodes[bc_index], int(v.carriers[bc_index]), alt, cov))
    print("Saved ground truth to %s.{variant_carriers,variant_cell_counts,template_variants}.tsv" % args.output)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
