library(DropletUtils)
library(dplyr)
library(R.utils)
library(plyr)
library(parallel)

# args <- commandArgs(trailing=TRUE)
bcFile <- "barcodes.tsv"  # make sure no -1 is there!
h5File <- "./NovaPertrubseqShortMay/outs/molecule_info.h5"
outName <- "short_gene/out"

bc <- read.table(bcFile)

## read in 10x HDF5 file
mol.info <- read10xMolInfo(h5File, get.umi=TRUE, get.gene = TRUE, 
		get.cell = TRUE, barcode.length = 16)

molecule_df <- mol.info$data %>% as.data.frame() %>% 
	filter(cell %in% bc$V1)

molecule_df <- molecule_df %>% 
	mutate(Gene = mol.info$genes[gene])

## convert from 2-bit encoding to nucleotide sequences
twoBit_toSeq <- function(deci,uLen){
	bin_str <- intToBin(deci)
	len <- length(strsplit(bin_str,"")[[1]])
	leading <- 2*uLen - len
	if(leading > 0){
		padded <- paste(strrep(0,leading), bin_str, sep = "")
	} else {padded = bin_str}
	split_bin <- sapply(seq(1, 2*uLen, 2),
		function(ix) substr(padded,ix,ix+1))
	mapped_toSeq <- mapvalues(split_bin, from = c("00","01","10","11"), 
				to = c("A","C","G","T"),warn_missing = FALSE)
	mapped_str <- paste(mapped_toSeq, collapse ="")
	return(mapped_str)
}


par_umis <- unlist(mclapply(molecule_df$umi, function(umi) twoBit_toSeq(umi,12),mc.cores=12))

molecule_df$UMI <- par_umis

## preserving structure reqd for downstream processing
# molecule_df <- molecule_df[,c(4, 6,2, 1, 7)]

write.table(molecule_df[,c(4, 6,8, 1, 7)],  paste0(outName,".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
