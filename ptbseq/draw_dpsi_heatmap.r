library(gplots)

# Read normalized counts
data = read.table("/home/andreyp/ablab/analysis/pertrubseq/dee/heatmaps/allVSnegative.dPSI_heatmap.tsv", header=T, sep="\t", as.is=TRUE)

genes = data[,1]
mat = as.matrix(data[,2:ncol(data)])
row.names(mat) <- genes

mat[is.na(mat)] <- -0.01
my_palette <- colorRampPalette(c("red", "white", "blue")) (n=20)
breaks <- seq(min(mat, na.rm = T), max(mat, na.rm = T), length.out = 21)
heatmap.2(mat, trace="none", na.color = "black", scale="none", 
          col = my_palette)

