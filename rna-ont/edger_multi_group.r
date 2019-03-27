library(edgeR)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(NMF)
library(gplots)
library(RColorBrewer)
library(pheatmap)

# === DE analysis ===

#read count data
countData = read.table("input/gene_count_matrix.csv", header=TRUE, sep=",", row.names=1 )
head(countData)
samples = names(countData)

#samples grouping
samplesData = read.table("input/pheno_data.tsv", header=TRUE, sep="\t", row.names=1 )
head(samplesData)
#create condition accroding 
condition = factor(samplesData[row.names(samplesData) %in% samples,][[4]])
designMat <- model.matrix(~condition)


# DGE list for EdgeR
dgList <- DGEList(counts=countData, genes=rownames(countData), group = condition)
#filter unexpressed genes
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]

#normalization
dgList <- calcNormFactors(dgList, method="TMM")

#dispersion
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)

#de analysis
fit <- glmQLFit(dgList, designMat)
qlf <- glmQLFTest(fit, coef=2:9)

#OR
#fit2 <- glmFit(dgList, designMat)
#lrt <- glmLRT(fit2, coef=2:9)

# === Check plots ===
plot_adj_pvalue_cutoff <- 0.01 #~ same as FDR
plot_lfc_cutoff <- 2

simpleTable <- NULL
log_fc_cols <- which(substr(colnames(qlf$table), 1, 5) == "logFC")
for (i in 1:nrow(qlf$table)) {
  lfcs <- qlf$table[i, log_fc_cols]
  lfc <- max(0, max(lfcs)) - min(0, min(lfcs))
  simpleTable <- rbind(simpleTable,c(lfc, qlf$table[i,]$PValue, qlf$table[i,]$logCPM))  
}
row.names(simpleTable) <- row.names(qlf$table)
colnames(simpleTable) <- c("LFC", "pvalue", "logCPM")
simpleTable <- as.matrix(simpleTable)
summary(simpleTable)

#plot volcano
plot(simpleTable[,"LFC"], -log10(simpleTable[,"pvalue"]),
     xlim=c(0, 15), ylim=c(0, 25),
     xlab="log2 fold change", ylab="-log10 p-value",
     type="n")
# then add the points
sel <- which(simpleTable[,"LFC"] <= plot_lfc_cutoff | simpleTable[,"pvalue"] >= plot_adj_pvalue_cutoff) 
points(simpleTable[sel,"LFC"], -log10(simpleTable[sel,"pvalue"]),col="black")
sel <- which(simpleTable[,"LFC"] > plot_lfc_cutoff & simpleTable[,"pvalue"] < plot_adj_pvalue_cutoff) 
points(simpleTable[sel,"LFC"], -log10(simpleTable[sel,"pvalue"]),col="red")

#plot smear
plot(simpleTable[,"logCPM"], simpleTable[,"LFC"],
     xlim=c(-10, 10), ylim=c(0, 15),
     xlab="logCPM", ylab="log2 fold change",
     type="n")
# then add the points
sel <- which(simpleTable[,"pvalue"] >= plot_adj_pvalue_cutoff) 
points(simpleTable[sel,"logCPM"], simpleTable[sel,"LFC"],col="black")
sel <- which(simpleTable[,"pvalue"] < plot_adj_pvalue_cutoff)
points(simpleTable[sel,"logCPM"], simpleTable[sel,"LFC"],col="red")

# === Choosing DE genes ===
adj_pvalue_cutoff <- 0.01 #~ same as FDR
lfc_cutoff <- 2

sigTable <- topTags(qlf, p.value = adj_pvalue_cutoff, n = length(row.names(qlf$table)))
sigTable <- sigTable$table

deTable <- NULL
log_fc_cols <- which(substr(colnames(sigTable), 1, 5) == "logFC")
for (i in 1:nrow(sigTable)) {
  lfcs <- sigTable[i, log_fc_cols]
  lfc <- max(0, max(lfcs)) - min(0, min(lfcs))
  if (lfc > lfc_cutoff) {
    deTable <- rbind(deTable, sigTable[i,])  
  }
}

summary(deTable$FDR)
summary(deTable$PValue)

deGenes <- row.names(deTable)
length(deGenes)

#output DE genes
write.table(deTable, file="output/DE_genes_data.tsv", quote=FALSE, sep='\t')
write.table(deGenes, file="output/DE_genes.list.txt", quote=FALSE, sep='\t')


# === HEAT MAP ===


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

#normalized counts
normalized_counts <- cpm(dgList, normalized.lib.sizes=FALSE)
#select DE genes
de_normalized_counts <- normalized_counts[deGenes,]
vals = as.matrix(de_normalized_counts)

#heatmap
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
row.names(zscore) = row.names(de_normalized_counts)

# Generate heatmap
mat = as.matrix(zscore)

#nice sample names


dev.off()
clust_map <- heatmap.2(mat)
sorted_mat <- mat[clust_map$rowInd,]
sorted_mat <- merge_matrix(sorted_mat)

#colnames(mat) <- sapply(strsplit(colnames(mat), "_"), function(x) x[1])

#color scheme
colors = colorRampPalette(c("blue","white","red"),space="rgb")(256)
png("output/general_heat_map.png", width = 2000, height = 2500, res = 300)
#heatmap.2(mat, col=colors,density.info="none",trace="none", margins=c(4,4),lhei=c(1,7), Colv=NA)
heatmap.2(sorted_mat, dendrogram = "none", col=colors,  Rowv = NA, trace = "none", density.info="none", 
          margins=c(4,4),lhei=c(1,7), Colv=NA, labCol = colnames(sorted_mat), labRow = FALSE)
dev.off()


#=== Supplementary functions for heatmaps ===
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



# ================= ANNOTATION HEATMAPS  =========================

description <-read.table("input/new_ann/new_degs_tf.gene_id.txt", sep="\t", row.names = NULL, header = FALSE)
colnames(description) <- c("gene_id", "description")
row.names(description) <- description$gene_id

sub_mat <- subset(mat, rownames(mat) %in% description$gene_id)
dev.off()
clust_map <- heatmap.2(sub_mat)
gene_order <- rownames(sub_mat[clust_map$rowInd,])

sorted_descr <- description[ order(match(description[,1], gene_order)), ] 
sorted_descr <- sorted_descr[order(sorted_descr[,2]),]

new_mat = NULL
genes = NULL
for (i in 1:length(sorted_descr[,1])) {
  new_mat <- rbind(new_mat, mat[toString(sorted_descr[i,]$gene_id),])
  genes <- rbind(genes, toString(sorted_descr[i,]$gene_id))
}

length(sorted_descr[,1]) == length(new_mat[,1])

new_mat <- merge_matrix(new_mat)
len <- length(new_mat[,1])

gene_property_column = 2
start_color = 1

colnames(new_mat) <- sapply(strsplit(colnames(new_mat), "_"), function(x) x[1])
uni_desc <- sapply(unique(sorted_descr[,gene_property_column]), function(x) which(sorted_descr[,gene_property_column] == x))

group_colors <- c(rep("red", len))
sample_colors <- c(rep("green", length(samples_desc)))

group_colors <- gcolors(group_colors, uni_desc, color41[start_color:41])
sample_colors <- gcolors(sample_colors, uni_samples, scolors)

#heatmap
legend_name <- unique(sorted_descr[,gene_property_column])
colors = colorRampPalette(c("blue","#CCCCCC","red"),space="rgb")(256)

png("output/new_TF_all.png", width = 2500, height = 2500, res = 300)
heatmap.2(new_mat, dendrogram = "none", col=colors,  RowSideColors = group_colors,Rowv = NA, trace = "none",density.info="none", 
          margins=c(10,10),lhei=c(1,6), Colv=NA, labCol = colnames(new_mat), labRow = NA, sepwidth=c(0.0001,0.001))
par(lend = 1)
legend("top",  inset = 0, legend = legend_name, col = color41[start_color:41], lty= 1, lwd = 4, ncol=3, cex=0.5)
dev.off()

