import sys

# CM038463.1      IsoQuant        transcript      475899  496537  .       -       .       gene_id "novel_gene_CM038463.1_8"; transcript_id "transcript7.CM038463.1.nnic"; exons "6";


def get_gene_id(attrs: str):
    v = attrs.split(" ")
    try:
        g_pos = v.index("gene_id")
    except ValueError:
        return ""
    return v[g_pos + 1]


processed_genes = set()
include_current = True

for l in open(sys.argv[1]):
    if l.startswith("#"):
        sys.stdout.write(l)

    v = l.split('\t')
    if v[2] == "gene":
        sys.stdout.write(l)
        include_current = True
    elif v[2] == "transcript":
        gene_id = get_gene_id(v[8])
        if gene_id in processed_genes:
            include_current = False
            continue
        else:
            processed_genes.add(gene_id)
            include_current = True
            sys.stdout.write(l)
    elif include_current:
        sys.stdout.write(l)
