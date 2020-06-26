############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import gffutils
import random
from traceback import print_exc

def process_gene_db(db, main_gtf_fname, excluded_gtf_fname, probability = 0.05):
    main_gtf = open(main_gtf_fname, "w")
    excl_gtf = open(excluded_gtf_fname, "w")
    for g in db.features_of_type('gene', order_by=('seqid', 'start')):
        exon_count = {}
        for t in db.children(g, featuretype='transcript', order_by='start'):
            exon_count[t.id] = 0
            for e in db.children(t, featuretype='exon', order_by='start'):
                exon_count[t.id] += 1

        to_remove = set()
        for t_id in exon_count.keys():
            if exon_count[t_id] == 1:
                to_remove.add(t_id)
        for t_id in to_remove:
            del exon_count[t_id]

        ignored_transcripts = set()
        if len(exon_count) > 5:
            for t_id in exon_count.keys():
                r = random.random()
                if r < probability:
                    ignored_transcripts.add(t_id)

        gene_str = '%s\t.\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s";\n' % (g.seqid, g.start, g.end, g.strand, g.id)
        main_gtf.write(gene_str)
        if len(ignored_transcripts):
            excl_gtf.write(gene_str)

        for t in db.children(g, featuretype='transcript', order_by='start'):
            transcript_str = '%s\t.\ttranscript\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s";\n' % \
                             (t.seqid, t.start, t.end, t.strand, g.id, t.id)
            prefix_columns = "%s\t.\texon\t" % t.seqid
            suffix_columns = '.\t%s\t.\tgene_id "%s"; transcript_id "%s"; \n' % \
                             (t.strand, g.id, t.id)

            if t.id in ignored_transcripts:
                current_file = excl_gtf
            else:
                current_file = main_gtf

            current_file.write(transcript_str)
            for e in db.children(t, featuretype='exon', order_by='start'):
                current_file.write(prefix_columns + "%d\t%d\t" % (e.start, e.end) + suffix_columns)

    main_gtf.close()
    excl_gtf.close()


def main():
    random.seed(11)
    if len(sys.argv) != 4:
        print("Usage %s <Input gffutils DB> <Remaining GTF> <Excluded GTF>" % sys.argv[0])
        exit(-1)
    gffutils_db = gffutils.FeatureDB(sys.argv[1], keep_order=True)
    process_gene_db(gffutils_db, sys.argv[2], sys.argv[3])


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
