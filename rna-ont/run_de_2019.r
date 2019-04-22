library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library (NMF)
library(gplots)
library(RColorBrewer)
library(pheatmap)

theme_set(theme_bw())

#read data
countData = read.table("input/gene_count_matrix.csv", header=TRUE, sep=",", row.names=1 )
samples = names(countData)
samplesData = read.table("input/pheno_data.tsv", header=TRUE, sep="\t", row.names=1 )
head(samplesData)

# Build the dataframe from the conditions
condition = factor(samplesData[row.names(samplesData) %in% samples,][[4]])
colData = data.frame(samples=samples, condition=condition)
summary(colData)

# create DeSeq object
ddsMat = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)

#Set the reference to be compared
ddsMat$condition = relevel(ddsMat$condition,"14")

ddsMat <- DESeq(ddsMat)

# normalize data
counts_normalized <- counts(ddsMat, normalized=TRUE)
counts_stabilized <-getVarianceStabilizedData(ddsMat)
norm_df <- mutate(as.data.frame(counts_normalized), gene_id = row.names(counts_normalized)) 
#write_tsv(norm_df, "deseq_normalization.tsv")
logNormCounts <- log2(counts_normalized + 1)

#log count pearson corrleation
logDist <- as.dist (1 - cor( logNormCounts , method = "pearson" ) )
plot(hclust(logDist),labels=colnames(logNormCounts),main=" log transformed read counts distance : Pearson correlation ")

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
res2 <- results(ddsMat,name = "condition_82_vs_14")
gene_id <- row.names(res)
res <- as.data.frame(res)
res <- cbind(gene_id, res)
res <- res %>%
  filter(!is.na(stat)) %>%
  arrange(desc(stat)) %>%
  select(gene_id, everything())
write_tsv(res,  "deseq_results.tsv")

#cutoffs
p_adj_cutoff <- 0.01
log2_cutoff <- 3

#volcano plot
res %>%
  filter(!is.na(padj)) %>%
  ggplot() +
  geom_point(aes(log2FoldChange, -log10(padj), colour = (abs(log2FoldChange) > log2_cutoff & padj < p_adj_cutoff))) +
  theme_bw() +
  scale_color_manual(values = c("black", "red")) +
  geom_text_repel(data = filter(res, padj < p_adj_cutoff),
                  aes(log2FoldChange, -log10(padj), label = ""),
                  min.segment.length = unit(0.1, "lines"))

#select DE genes
res_sorted <- res[order(res$padj), ]
DGEgenes <-  subset ( res_sorted , padj < p_adj_cutoff & abs(log2FoldChange) > log2_cutoff)$gene_id 
matDGEgenes <- subset(logNormCounts, gene_id %in% DGEgenes)
write_tsv(subset ( res_sorted , gene_id %in% DGEgenes) %>% 
            select(gene_id,log2FoldChange,padj), "DE_genes.tsv") 

aheatmap(matDGEgenes, Rowv = TRUE , Colv = TRUE , distfun = "euclidean" , hclustfun = "average" )

# Get normalized counts and write this to a file
normalized_counts = counts(ddsMat,normalized=TRUE)

# Turn it into a dataframe to have proper column names.
normalized_counts_dt = data.frame("gene_id"=rownames(normalized_counts),normalized_counts)
gene = subset(normalized_counts_dt, gene_id %in% DGEgenes)[1]
vals = as.matrix(subset(normalized_counts_dt, gene_id %in% DGEgenes)[2:ncol(normalized_counts_dt)])
#TOO SLOW
#gene = normalized_counts_dt[1]
#vals = as.matrix(normalized_counts_dt[2:ncol(normalized_counts_dt)])

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
zscore=score
row.names(zscore) = row.names(gene)

# Generate new heatmap
mat = as.matrix(zscore)
aheatmap(mat, Colv = NA)

plot(hclust(logDist),labels=colnames(logNormCounts),main=" log transformed read counts distance : Pearson correlation ")
colors = colorRampPalette(c("yellow","black","red"),space="rgb")(256)
heatmap.2( mat, col=colors,density.info="none",trace="none", margins=c(4,4),lhei=c(1,7), Colv=NA)
write_tsv(data.frame(gene,mat),  "heatmap.tsv")


# == heatmap with gene clusters ==
#read description file
frequent_terms <- function(descr) {
  res <- NULL
  uni_descr <- sapply(unique(descr[,2]), function(x) which(descr[,2] == x))
  
  for (i in 1:length(uni_descr)) {
    index_list <- uni_descr[[i]]
    if (length(index_list) > 5) {
      for (j in 1:length(index_list)) {
        res <- rbind(res, descr[index_list[[j]],])
      }
    }
  }
  res
}


merge_matrix <- function(m) {
  short_col_names <- unique(sapply(strsplit(colnames(m), "_"), function(x) x[1]))
  newm = NULL
  cnames <- colnames(m)
  for (i in 1:nrow(m)) {
    r = c()
    for (j in 1:length(short_col_names)) {
      r <- c(r, mean(m[i, cnames[startsWith(cnames, short_col_names[j])]]))
    }
    newm <- rbind(newm, r)
  }
  colnames(newm) <- short_col_names
  row.names(newm) <- row.names(m)
  newm
}


color41 <- c("#A52A2A", "#D2691E", "#1122F5", "#FFFF00", "#FFE4C4", "#3CB371", "#A9A9A9", "#4169E1", "#008B8B", "#6B8E23", "#F00FF0", "#FF4500",
             "#11FF50", "#333333", "#7FFFD4", "#ADFF2F", "#FFA500", "#9F1AF0", "#000080", "#FFE4B5", "#DAA520", "#F4A460", "#BDB76B", "#DB7093",
             "#7FFF00", "#9932CC", "#B8860B", "red", "#006400", "#483D8B", "#FA8072", "#C0C0C0", "#D3D3D3", "#40E0D0", "#F08080", "#FAF0E6", "#DDA0DD",
             "#AFEEEE", "#8B0000", "#FFEFD5", "#0000CD")

scolors <- c("#A52A2A", "#C2291E", "#B14215", "#EEDE00", "#12F424", "#15D434", "#22C941", "#29BA55", "#3423E1")


gcolors <- function(group_colors, uni, colors) {
  for (i in 1:length(uni)) {
    index_list <- uni[[i]]
    for (j in 1:length(index_list)) {
      group_colors[index_list[[j]]] <- colors[i] 
    }
  }
  group_colors
}


# ================= FOR GO  =========================


description <-read.table("all_genes.go_terms_annot.txt", sep="\t", row.names = NULL, header = TRUE)
description[,1] <- sapply(strsplit(sapply(description[,1],toString), "-"), function(x) x[1])
description <- distinct(description)
#description <- frequent_terms(description)

go_terms <- read.table("enriched_go_terms_categories.txt", sep="\t", row.names = NULL, header = TRUE)
description <- subset(description, description$GO_term %in% go_terms$GO_term)
go_mat <- as.matrix(go_terms)
row.names(go_mat) <- go_terms$GO_term

new_desc = NULL
for (i in 1:length(description[,1])) {
  new_desc <- rbind(new_desc, c(description[i,]$Gene_id, go_mat[toString(description[i,]$GO_term),]))
}

description <- as.data.frame(new_desc[order(new_desc[,2], new_desc[,3]),])
colnames(description) <- c("Gene_id", "GO_term", "Domain", "Description")
short_desc <- as.data.frame(subset(description, description$Gene_id %in% rownames(mat)))


sub_mat <- subset(mat, rownames(mat) %in% short_desc$Gene_id)
dev.off()
clust_map <- heatmap.2(sub_mat)
gene_order <- rownames(sub_mat[clust_map$rowInd,])

short_desc <- short_desc[ order(match(short_desc[,1], gene_order)), ] 
short_desc <- short_desc[order(short_desc[,3], short_desc[,4]),]

new_mat = NULL
for (i in 1:length(short_desc[,1])) {
  new_mat <- rbind(new_mat, mat[toString(short_desc[i,]$Gene_id),])
}

new_mat <- merge_matrix(new_mat)
len <- length(rownames(new_mat))

colnames(new_mat) <- sapply(strsplit(colnames(new_mat), "_"), function(x) x[1])
uni_desc <- sapply(unique(short_desc[,4]), function(x) which(short_desc[,4] == x))

group_colors <- c(rep("red", len))
sample_colors <- c(rep("green", length(samples_desc)))

group_colors <- gcolors(group_colors, uni_desc, color41[4:41])
sample_colors <- gcolors(sample_colors, uni_samples, scolors)

#heatmap
legend_name <- unique(short_desc[,4])
colors = colorRampPalette(c("blue","#CCCCCC","red"),space="rgb")(256)

png("figures/GO2.png", width = 2000, height = 2500, res = 300)
heatmap.2(new_mat, dendrogram = "none", col=colors,  RowSideColors = group_colors,Rowv = NA, trace = "none",density.info="none", 
          margins=c(10,6),lhei=c(1,6), Colv=NA, labCol = colnames(new_mat), labRow = FALSE, sepwidth=c(0.0001,0.001))
par(lend = 1)
legend("top",  inset = 0, legend = legend_name, col = color41[4:41], lty= 1, lwd = 4, ncol=3, cex=0.5)
dev.off()




# ================= FOR KEGG  =========================

#description <-read.table("selected_kegg_genes.txt", sep="\t", row.names = NULL, header = TRUE)
#description <-read.table("kegg+oil_short.txt", sep="\t", row.names = NULL, header = TRUE)
description <-read.table("flowering.genes.txt", sep="\t", row.names = NULL, header = TRUE)
description <- distinct(description)
short_desc <- as.data.frame(subset(description, description$Gene_id %in% rownames(mat)))

sub_mat <- subset(mat, rownames(mat) %in% short_desc$Gene_id)
dev.off()
clust_map <- heatmap.2(sub_mat)
gene_order <- rownames(sub_mat[clust_map$rowInd,])

short_desc <- short_desc[ order(match(short_desc[,1], gene_order)), ] 
short_desc <- short_desc[order(short_desc[,2]),]

new_mat = NULL
genes = NULL
for (i in 1:length(short_desc[,1])) {
  new_mat <- rbind(new_mat, mat[toString(short_desc[i,]$Gene_id),])
  genes <- rbind(genes, toString(short_desc[i,]$Gene_id))
}

length(short_desc[,1]) == length(new_mat[,1])

new_mat <- merge_matrix(new_mat)
len <- length(new_mat[,1])

gene_property_column = 2
start_color = 1

colnames(new_mat) <- sapply(strsplit(colnames(new_mat), "_"), function(x) x[1])
uni_desc <- sapply(unique(short_desc[,gene_property_column]), function(x) which(short_desc[,gene_property_column] == x))

group_colors <- c(rep("red", len))
sample_colors <- c(rep("green", length(samples_desc)))

group_colors <- gcolors(group_colors, uni_desc, color41[start_color:41])
sample_colors <- gcolors(sample_colors, uni_samples, scolors)

#heatmap
legend_name <- unique(short_desc[,gene_property_column])
colors = colorRampPalette(c("blue","#CCCCCC","red"),space="rgb")(256)

png("figures/keggg.DE.3.png", width = 2500, height = 2500, res = 300)
heatmap.2(new_mat, dendrogram = "none", col=colors,  RowSideColors = group_colors,Rowv = NA, trace = "none",density.info="none", 
          margins=c(10,10),lhei=c(1,6), Colv=NA, labCol = colnames(new_mat), labRow = genes, sepwidth=c(0.0001,0.001),
          sepcolor="white",
          colsep=1:ncol(mat),
          rowsep=1:nrow(mat))
par(lend = 1)
legend("top",  inset = 0, legend = legend_name, col = color41[start_color:41], lty= 1, lwd = 4, ncol=3, cex=0.5)
dev.off()


# ======= ONLY FOR TF =========
description <-read.table("deg_tf_blast2go_annot.txt", sep="\t", row.names = NULL)
description[,1] <- sapply(strsplit(description[,1], "-"), function(x) x[1])
description <- unique(description)
description <- frequent_terms(description)

new_name <- description[,1]
len <- length(new_name)
new_mat <- mat[new_name,]
new_mat <- merge_matrix(mat[new_name,])
colnames(new_mat) <- sapply(strsplit(colnames(new_mat), "_"), function(x) x[1])
uni_desc <- sapply(unique(description[,2]), function(x) which(description[,2] == x))
#samples_desc <- colnames(new_mat)
#uni_samples <- sapply(unique(samples_desc), function(x) which(samples_desc == x))

#colors

group_colors = c(rep("red", len))
sample_colors = c(rep("green", length(samples_desc)))

group_colors <- gcolors(group_colors, uni_desc, color41)
sample_colors <- gcolors(sample_colors, uni_samples, scolors)

#heatmap

legend_name <- unique(description[,2])


colors = colorRampPalette(c("blue","#CCCCCC","red"),space="rgb")(256)
#plot(hclust(logDist),labels=colnames(logNormCounts),main=" log transformed read counts distance : Pearson correlation ")
png("figures/TF_heatmap_4.png", width = 2000, height = 2500, res = 300)
heatmap.2(new_mat, dendrogram = "none", col=colors,  RowSideColors = group_colors,Rowv = NA, trace = "none",density.info="none", 
          margins=c(6,6),lhei=c(1,6), Colv=NA, labCol = colnames(new_mat), labRow = FALSE, sepwidth=c(0.0001,0.001),
          sepcolor="white",
          colsep=1:ncol(mat),
          rowsep=1:nrow(mat))
par(lend = 1)
legend("left",  inset = 0.002, legend = legend_name, col = color41, lty= 1, lwd = 6, ncol=1, cex=0.5)
dev.off()
#legend("top",  inset = 0.002, legend = unique(samples_desc), col = scolors, lty= 1, lwd = 5, ncol=9, cex=0.5)


