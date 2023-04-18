library(gplots)
library(dplyr)

# Read normalized counts
  data = read.table("/home/andreyp/ablab/analysis/pertrubseq/dee/SCv3.UMI_filtered.exon_filtered.21.03/full_guides/oneVScontrol.dPSI.defaults.tsv", header=T, sep="\t", as.is=TRUE)

genes = data[,1]
mat = as.matrix(data[,2:ncol(data)])
row.names(mat) <- genes

mat[is.na(mat)] <- -0.01
my_palette <- colorRampPalette(c("red", "white", "blue")) (n=20)
breaks <- seq(min(mat, na.rm = T), max(mat, na.rm = T), length.out = 21)
heatmap.2(mat, trace="none", na.color = "black", scale="none", 
          col = my_palette)



subs_m <- as.matrix(data) ## Main matrix

DR <- dist(subs_m) ## Get the euclidean distance
DR <- as.data.frame(as.matrix(DR))

DR2 <- DR %>% mutate(numNA = ncol(DR) - rowMeans(is.na(.))*ncol(DR)) ## find number of NAs per row
tooManyNAs_R <- unname(which(DR2$numNA <= ncol(DR)-200 )) ## I set it to 200 for a 6000 row matrix but you will have to adjust the parameter
length(tooManyNAs_R) ## was like 1200 


subs_m <- as.matrix(subs)
subs_m <- subs_m[-c(tooManyNAs_R),] ## Remove rows with too many NAs
dim(subs_m) ## was like 4500


