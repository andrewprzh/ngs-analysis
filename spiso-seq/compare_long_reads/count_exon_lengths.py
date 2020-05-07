############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import gffutils
import numpy

def process_db(db):
    processed_exons = set()
    lengths = []
    for e in db.features_of_type('exon'):
        etuple = (e.seqid, e.start, e.end)
        if etuple not in processed_exons:
            processed_exons.add(etuple)
            lengths.append(e.end + 1 - e.start)
    return lengths


def main():
    genedb = sys.argv[1]
    if not os.path.isfile(genedb):
        raise Exception("Gene database " + genedb + " does not exist")
    db = gffutils.FeatureDB(genedb, keep_order=True)
    lengths = process_db(db)

    count, exon_len = numpy.histogram(lengths, bins=[10 * i for i in range(61)] + [10000])
    for i in range(len(count)):
        print(str(exon_len[i]) + '\t' + str(count[i]))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
