library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library (NMF)
library(gplots)
library( "RColorBrewer" )
library(pheatmap)

theme_set(theme_bw())

setwd("~/ablab/analysis/UK/")
countData = read.table("Novogene.HiFi.transcript_grouped_counts.tsv", header=TRUE, sep="\t", row.names=1 )
unexpressed_genes <- as.data.frame(rowMeans(countData))
colnames(unexpressed_genes) <- c("avg_count")
gene_id <- row.names(unexpressed_genes)
unexpressed_genes <- cbind(gene_id, unexpressed_genes)
unexpressed_genes <- subset(unexpressed_genes, unexpressed_genes[, "avg_count"] > 0)
countData <- countData[unexpressed_genes$gene_id,]

samples = names(countData)
samplesData = read.table("pheno_data.tsv", header=TRUE, sep="\t", row.names=1 )
head(samplesData)

# Build the dataframe from the conditions
condition = factor(samplesData[row.names(samplesData) %in% samples,][[4]])

colData = data.frame(samples=samples, condition=condition)
summary(colData)

# create DeSeq object
ddsMat = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)

#Set the reference to be compared
ddsMat$condition = relevel(ddsMat$condition,"51")

#ddsMat$batch <- factor(ddsMat$batch, levels=c('A', 'B', 'C'))
# Run DESeq2
#ddsMat <- estimateSizeFactors(ddsMat)
#ddsMat <- estimateDispersions(ddsMat)
ddsMat <- DESeq(ddsMat)
#ddsMat <- DESeq(ddsMat, full=design(ddsMat), test="LRT", reduced = ~ batch)

# normalize data
counts_normalized <- counts(ddsMat, normalized=TRUE)
counts_stabilized <-getVarianceStabilizedData(ddsMat)
norm_df <- mutate(as.data.frame(counts_normalized), gene_id = row.names(counts_normalized)) 
write_tsv(norm_df, "deseq_normalization.tsv")
logNormCounts <- log2(counts_normalized + 1)

#log count pearson corrleation
logDist <- as.dist (1 - cor( logNormCounts , method = "pearson" ) )
plot(hclust(logDist),labels=colnames(logNormCounts),main=" log transformed read counts distance : Pearson correlation ")

#rlog count pearson correlation
#ddsRlog <- rlog ( ddsMat , blind = TRUE )
#rlogNormCounts <- assay ( ddsRlog )
#rlogDist <- as.dist (1 - cor( rlogNormCounts , method = "pearson" ) )
#plot(hclust(rlogDist),labels=colnames(rlogNormCounts),main=" rlog transformed read counts distance : Pearson correlation ")


# plot PCA
pca <- prcomp(t(counts_stabilized))
pca_df <- tibble(sample = row.names(pca$x),
                 PC1 = pca$x[, "PC1"], 
                 PC2 = pca$x[, "PC2"], 
                 PC3 = pca$x[, "PC3"]) 

require(gridExtra)
plot1 <- ggplot(pca_df) +
  geom_point(aes(PC1, PC2, colour = condition), size = 3) +
  labs(x = paste0("PC1 (", summary(pca)$importance[2, 1]*100, "%)"),
       y = paste0("PC2 (", summary(pca)$importance[2, 2]*100, "%)")) +
  scale_colour_brewer(palette = "Set1") +
  geom_text_repel(data = pca_df, 
                  aes(PC1, PC2, label = substr(sample,1,5) ),
                  min.segment.length = unit(0.5, "lines")) +
  theme(legend.position = "top")
grid.arrange(plot1, ncol = 1)


# get results
res <- results(ddsMat)
gene_id <- row.names(res)
res <- as.data.frame(res)
res <- cbind(gene_id, res)
res <- res %>%
  filter(!is.na(stat)) %>%
  arrange(desc(stat)) %>%
  select(gene_id, everything())
write_tsv(res,  "deseq_results.tsv")

p_adj_cutoff <- 0.05
log2_cutoff <- 1

res %>%
  filter(!is.na(padj)) %>%
  ggplot() +
  geom_point(aes(log2FoldChange, -log10(padj), colour = (abs(log2FoldChange) > log2_cutoff & padj < p_adj_cutoff))) +
  theme_bw() +
  scale_color_manual(values = c("black", "red")) +
  geom_text_repel(data = filter(res, padj < p_adj_cutoff),
                  aes(log2FoldChange, -log10(padj), label = ""),
                  min.segment.length = unit(0.1, "lines"))


res_sorted <- res[order(res$padj), ]
DGEgenes <-  subset ( res_sorted , padj < p_adj_cutoff & abs(log2FoldChange) > log2_cutoff)$gene_id 
matDGEgenes <- logNormCounts[DGEgenes, ]
low_expressed_genes <- as.matrix(rowMeans(matDGEgenes))
matDGEgenes <- matDGEgenes - rowMeans(matDGEgenes)

aheatmap(matDGEgenes, Rowv = TRUE , Colv = TRUE , distfun = "euclidean" , hclustfun = "average" )
write_tsv(subset ( res_sorted , gene_id %in% DGEgenes) %>% 
            select(gene_id,log2FoldChange,padj), paste(exp_name, "_DE_transcritps.tsv")) 

# Get normalized counts and write this to a file
normalized_counts = counts(ddsMat,normalized=TRUE)

# Turn it into a dataframe to have proper column names.
normalized_counts_dt = data.frame("gene_id"=rownames(normalized_counts),normalized_counts)

# Save the normalize data matrix.
write.table(normalized_counts_dt, file=paste(exp_name, "_norm-matrix-deseq2.txt"), sep="\t",  row.name=TRUE, col.names=TRUE,quote=FALSE)

gene = normalized_counts_dt[DGEgenes,1]
vals = as.matrix(normalized_counts_dt[DGEgenes,2:ncol(normalized_counts_dt)])

# Adds a little noise to each element
# To avoid the clusteing function failing on zero
# variance datalines.
vals = jitter(vals, factor = 1, amount=0.00001)

# Calculate zscore
score = NULL
for (i in 1:nrow(vals)) {
  row=vals[i,]
  zscore=(row-mean(row))/sd(row)
  score =rbind(score,zscore)
}

row.names(score) = gene
zscore=score

# Generate new heatmap
mat = as.matrix(zscore)

colors = colorRampPalette(c("green","black","red"),space="rgb")(256)
heatmap.2(mat,col=colors,density.info="none",trace="none", margins=c(14,14),lhei=c(1,5))

write_tsv(data.frame(gene, mat),  paste(exp_name, "_heatmap.tsv"))


#other heat maps

rld <- rlog( ddsMat )
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL 
heatmap.2( sampleDistMatrix, trace="none")

heatmap.2( assay(rld)[ DGEgenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[
             colData(rld)$condition ] )

matrix <- assay(rld)[ DGEgenes, ]
matrix <- matrix - rowMeans(matrix)
annotation_data <- as.data.frame(colData(rld))
pheatmap(matrix, annotation_col=annotation_data)
