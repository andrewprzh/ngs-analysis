import os
import sys
import gffutils



def main():
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: " + sys.argv[0] + " <GFF> <gene list> \n")
        exit(0)

    gene_names = []
    for l in open(sys.argv[2]):
        gene_names.append(l.strip())

    if not os.path.isfile(sys.argv[1]):
        raise Exception("Gene database " + sys.argv[1] + " does not exist")
    self.db = gffutils.FeatureDB(sys.argv[1], keep_order=True)

    for gene in gene_names:
        genedb = db[gene]
        print(genedb.id)

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)