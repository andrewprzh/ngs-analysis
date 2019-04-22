library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

#Change the next two lines!!!
GTFfile = "~/ablab/analysis/RNA_nanopores/genes/Illumina_ONT_FL.all_pasa.gtf" #CHANGE ME!!!
FASTAfile = "~/ablab/analysis/RNA_nanopores/genes/CC_2C_reference_pseudomolecules.fasta" #CHANGE ME!!!

#Load the annotation and reduce it
GTF <- import.gff(GTFfile, format="gtf", genome="GRCm38.71", feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementLengths(grl))

#Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(output) <- c("Length", "GC")

write.table(output, file="input/GC_lengths.tsv", sep="\t")
