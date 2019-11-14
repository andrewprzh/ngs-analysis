import sys
from traceback import print_exc
from lr_snp_caller import *
from filter_snps import *


MIN_TOTAL_FRAC = 0.1
MIN_GROUP_FRAC = 0.2
MIN_FREQ = 0.3

groups = [(0, 141), (141, 271), (272, 376)]

def freqs(snp):
    freq = []
    for i in range(len(snp.sample_coverage)):
        if snp.sample_coverage[i] == 0:
            freq.append(0.0)
        else:
            freq.append(float(snp.sample_counts[i]) / float(snp.sample_coverage[i]))
    return freq


def frac_within_range(freq_list, r1, r2):
    total_abundancy = 0
    for i in range(r1, r2):
        if freq_list[i] >= MIN_FREQ:
            total_abundancy += 1
    return float(total_abundancy) / float(r2 - r1)


def is_snp_ok(snp):
    freq = freqs(snp)

    res = []
    res.append(frac_within_range(freq, 0, len(freq)))
#    if total_frac >= MIN_TOTAL_FRAC:
#        return True

    for g in groups:
        res.append(frac_within_range(freq, g[0], g[1] + 1))
#        if group_frac >= MIN_GROUP_FRAC:
#            return True
    return res


def snp_to_tsv(snp, chr_id, pos):
    freq = freqs(snp)
    res = chr_id + "_" + str(pos) + "_" + snp.reference_nucl + "_" + snp.alternative_nucl
    for i in range(len(snp.sample_coverage)):
        if snp.sample_coverage[i] == 0:
            res += '\t' + 'NA'
        else:
            res += '\t' + '{:.2f}'.format(freq[i])
    return res

def snp_to_tsv_avg(snp, chr_id, pos):
    freq = freqs(snp)
    res = chr_id + "_" + str(pos) + "_" + snp.reference_nucl + "_" + snp.alternative_nucl
    mean = 0.0
    count = 0
    for i in range(len(snp.sample_coverage)):
        if snp.sample_coverage[i] != 0:
            mean += freq[i]
            count += 1
    mean /= count

    for i in range(len(snp.sample_coverage)):
        if snp.sample_coverage[i] == 0:
            res += '\t' + '{:.2f}'.format(mean)
        else:
            res += '\t' + '{:.2f}'.format(freq[i])
    return res


def main():
    f = sys.argv[1]
#    print("Reading from " + f)
    reader = TSVParser(f, [], None)
#    print("Detected " + str(len(reader.sample_ids)) + " samples")
    snp_storage = SNPStorage()
    reader.fill_map(snp_storage)

    for chr_id in snp_storage.snp_map:
        for pos in snp_storage.snp_map[chr_id]:
            for snp in snp_storage.snp_map[chr_id][pos]:
                print(snp_to_tsv(snp, chr_id, pos))
#                frac = is_snp_ok(snp)
#                print(chr_id + "\t" + str(pos) + "\t" + snp.reference_nucl + "\t" + snp.alternative_nucl + "\t" + "\t".join(map('{:.2f}'.format, frac)))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)


