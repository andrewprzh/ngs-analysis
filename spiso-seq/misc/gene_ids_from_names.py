import os
import sys
import gffutils
from traceback import print_exc

def get_codon_pair(transcript, db):
    start_codon = None
    stop_codon = None
    for s in db.children(transcript, featuretype='start_codon', order_by='start'):
        start_codon = s.start
        if s.strand == "+":
            break
    for s in db.children(transcript, featuretype='stop_codon', order_by='start'):
        stop_codon = s.start
        if s.strand == "-":
            break

    if start_codon is None:
        for s in db.children(transcript, featuretype='CDS', order_by='start'):
            if s.strand == "+":
                start_codon = s.start
                break
            else:
                start_codon = s.end
    if stop_codon is None:
        for s in db.children(transcript, featuretype='CDS', order_by='start'):
            if s.strand == "+":
                stop_codon = s.end + 1
            else:
                stop_codon = s.start - 2
                break

    return start_codon, stop_codon


def count_codons(gene_db, db):
    start_codons = set()
    stop_codons = set()
    for t in db.children(gene_db, featuretype = 'transcript', order_by='start'):
        start_codon, stop_codon = get_codon_pair(t, db)
        if stop_codon is not None:
            stop_codons.add(stop_codon)
        if start_codon is not None:
            start_codons.add(start_codon)
    return len(start_codons), len(stop_codons)

def main():

    if len(sys.argv) < 3:
        sys.stderr.write("Usage: " + sys.argv[0] + " <GFF> <gene list> \n")
        exit(0)

    gene_names = set()
    for l in open(sys.argv[2]):
        gene_names.add(l.strip())

    if not os.path.isfile(sys.argv[1]):
        raise Exception("Gene database " + sys.argv[1] + " does not exist")
    db = gffutils.FeatureDB(sys.argv[1], keep_order=True)

    found_genes = set()
    for g in db.features_of_type('gene', order_by=('seqid', 'start')):
        gene_db = db[g.id]
        gene_name = db[g.id].attributes['gene_name'][0]
        if gene_name in gene_names:
            print(g.id + ": " + str(count_codons(gene_db, db)))
            found_genes.add(gene_name)


    print("Total number " + str(len(gene_names)))
    print("Not found:")
    for g in gene_names:
        if g not in found_genes:
            print(g)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
