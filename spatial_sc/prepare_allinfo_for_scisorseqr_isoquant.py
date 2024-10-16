import sys

# return  "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s" % (self.read_id, self.gene_id, cell_type,
#                                                                         self.barcode, self.umi, introns_str, TSS, polyA,
#                                                                         exons_str, read_type, len(intron_blocks),
#                                                                         self.transcript_id, self.transcript_type)
# if(ncol(all_info) == 9){
# 	all_info <- all_info[,c(1:4,6)]
# 	colnames(all_info) <- c("Readname","Gene","Celltype","Barcode","Isoform")
# } else if(ncol(all_info) >= 11){
# 	all_info <- all_info[,c(1:4,6:8)]
# 	colnames(all_info) <- c("Readname","Gene","Celltype","Barcode","Isoform","TSS","PolyA")
#


for l in open(sys.argv[1]):
    if l.startswith("#"):
        sys.stdout.write(l)
    v = l.strip().split('\t')
    assert len(v) == 13
    polyA = "NA" if v[7].startswith("No") else v[7]
    tss = "NA" if v[6].startswith("No") else v[6]
    sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (v[0], v[1], v[2], v[3], v[4],
                                                                       v[11], tss, polyA, v[5], v[8], v[9]))