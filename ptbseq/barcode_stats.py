import sys
from collections import defaultdict
import numpy


def read_barcodes(inf_name, gene_prefix=""):
    umi_dict = defaultdict(set)
    gene_dict = defaultdict(lambda: defaultdict(set))
    for l in open(inf_name):
        v = l.strip().split('\t')
        umi = v[2]
        bc = v[3]
        gene_id = v[4]
        umi_dict[bc].add(umi)

        if gene_prefix and not gene_id.startswith(gene_prefix):
            continue
        gene_dict[bc][gene_id].add(umi)

    return umi_dict, gene_dict


def count_barcode_stats(umi_dict):
    umi_counts = [len(x) for x in umi_dict.values()]
    bins = [i * 1000 for i in range(200)]
    hist = numpy.histogram(umi_counts, bins=bins)
    print(hist[1])


def count_gene_stats(gene_dict):
    gene_counts = defaultdict(int)
    gene_counts_filt = defaultdict(int)
    for bc in gene_dict.keys():
        gene2umi = gene_dict[bc]
        distinct_genes = len(gene2umi.keys())
        gene_counts[distinct_genes] += 1
        if distinct_genes > 1:
            max_umis = 0
            supported_by_2 = 0
            for gene_id in gene2umi.keys():
                max_umis = max(max_umis, len(gene2umi[gene_id]))
                if len(gene2umi[gene_id]) > 1:
                    supported_by_2 += 1
            if max_umis > 1:
                gene_counts_filt[supported_by_2] += 1
            else:
                gene_counts_filt[distinct_genes] += 1
        else:
            gene_counts_filt[distinct_genes] += 1

    print(gene_counts)
    print(gene_counts_filt)

    print(sum(gene_counts.values()))
    print(sum(gene_counts_filt.values()))


umi_dict, gene_dict = read_barcodes(sys.argv[1], "CROP")
count_barcode_stats(umi_dict)

