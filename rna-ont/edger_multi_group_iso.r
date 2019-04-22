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
countData = read.table("input/transcript_counts.tsv", header=TRUE, sep="\t" ,  row.names=1)

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
total_conditions <- length(unique(condition))
fit <- glmQLFit(dgList, designMat)
qlf <- glmQLFTest(fit, coef=2:total_conditions)

#OR
#fit2 <- glmFit(dgList, designMat)
#lrt <- glmLRT(fit2, coef=2:total_conditions)

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
     xlim=c(0, 25), ylim=c(0, 25),
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
#write.table(deTable, file="output/DE_cap_genes_data.tsv", quote=FALSE, sep='\t')
write.table(deGenes, file="output/DE_isoforms.list.txt", quote=FALSE, sep='\t', row.names = FALSE)


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

get_color_gradient <- function(m) {
  #color scheme
  def_stop <- 5.0
  mmax <- abs(max(m))
  mmin <- abs(min(m))
  absmax <- max(15.0, mmax)
  absmin <- max(15.0, mmin)
  color_resolution <- 0.1

  low_gradient <- colorpanel( (absmin - def_stop) / color_resolution , "#000080", "blue" )
  mid_low_gradient <- colorpanel( def_stop / color_resolution, "blue", "white" )
  #mid_low_gradient <- mid_low_gradient[2:length(mid_low_gradient)]
  mid_top_gradient <- colorpanel( def_stop / color_resolution, "white", "red" )
  #mid_top_gradient <- mid_top_gradient[2:length(mid_top_gradient)]
  top_gradient <- colorpanel( (absmax - def_stop) / color_resolution, "red", "#800000" )
  #top_gradient <- top_gradient[2:length(top_gradient)]
  
  full_gradient <- c(low_gradient, mid_low_gradient, mid_top_gradient, top_gradient)
  
  center_coord <- absmin / color_resolution
  lower_bound <- center_coord - (mmin / color_resolution)
  upper_bound <- center_coord + (mmax / color_resolution)
 
  full_gradient[lower_bound:upper_bound]
}

#normalized counts
#normalized_counts <- cpm(dgList, normalized.lib.sizes=TRUE, log = TRUE)

gene_lengths <- read.table("input/gene_lenths.tsv", header=TRUE, sep="\t" )
gene_lengths <- read.table("input/new_ann/isoforms/isoform_lengths.tsv" , header=TRUE, sep="\t")

gene_lengths$Length <- pmax(1, gene_lengths$Length - 200)
row.names(gene_lengths) <- gene_lengths$Geneid
gene_lengths <- gene_lengths[row.names(dgList),]

raw_normalized_counts <- rpkm(dgList, normalized.lib.sizes=TRUE, gene.length = gene_lengths$Length)
raw_normalized_counts <- na.omit(raw_normalized_counts, cols=colnames(raw_normalized_counts))
norm_fact <- colSums(raw_normalized_counts)
norm_fact <- log2(norm_fact / 10^6)

normalized_counts <- rpkm(dgList, normalized.lib.sizes=TRUE, gene.length = gene_lengths$Length, log = TRUE)
for (c in colnames(normalized_counts)) {
  normalized_counts[,c] <- normalized_counts[,c] - norm_fact[c]
}

#select DE genes -- select a way
#de_normalized_counts <- subset(normalized_counts, startsWith(row.names(normalized_counts), "novel"))
#mod_genes <-read.table("input/new_ann/updated_genes.txt", sep="\t", row.names = NULL, header = FALSE)
#de_normalized_counts <- subset(normalized_counts, row.names(normalized_counts) %in% mod_genes[,1])
de_normalized_counts <- subset(normalized_counts, row.names(normalized_counts) %in% deGenes)
de_normalized_counts <- normalized_counts

#use this only for taking CC and OC taking subsamples
oc_samples <- samples[startsWith(samples, "OC")]
gs_samples <- samples[startsWith(samples, "GS")]
yl_samples <- samples[startsWith(samples, "YL")]
cc_samples <- samples[startsWith(samples, "CC")]
sub_samples <- c(gs_samples, yl_samples, cc_samples, oc_samples)
de_normalized_counts <- de_normalized_counts[, sub_samples]
head(de_normalized_counts)

merged_counts <- merge_matrix(de_normalized_counts)
head(merged_counts)
vals = as.matrix(merged_counts)
#heatmap
# Adds a little noise to each element
# To avoid the clusteing function failing on zero
# variance datalines.
vals = jitter(vals, factor = 1, amount=0.00001)

# Calculate zscore
score = NULL
for (i in 1:nrow(vals)) {
  row=vals[i, ]
  zscore = (row - mean(row)) #/sd(row)
  score = rbind(score, zscore)
}
row.names(score) = row.names(merged_counts)

# Generate heatmap
mat = as.matrix(score)
mat <- na.omit(mat, cols=colnames(mat))

head(mat)
sum(mat[1,])

draw_tree <- TRUE #NA
sorted_mat <- mat
if (is.na(draw_tree)) {
  clust_map <- aheatmap(mat, Colv = NA, scale = "none")
  sorted_mat <- mat[clust_map$rowInd,]
}
sum(sorted_mat[3,])

#color scheme
h_colors <- get_color_gradient(sorted_mat)

#colors = colorRampPalette(c("blue","white","red"),space="rgb")(256)
aheatmap(sorted_mat, color = h_colors, labCol = colnames(sorted_mat), Colv=NA, labRow = NA, Rowv = draw_tree,
         legend=TRUE, fontsize = 7, scale="none", filename = "isoforms/isoforms_heatmap.pdf")


#=== Supplementary functions for heatmaps ===
frequent_terms <- function(descr, freq = 10) {
  res <- NULL
  uni_descr <- sapply(unique(descr[,2]), function(x) which(descr[,2] == x))
  
  for (i in 1:length(uni_descr)) {
    index_list <- uni_descr[[i]]
    if (length(index_list) > freq) {
      for (j in 1:length(index_list)) {
        res <- rbind(res, descr[index_list[[j]],])
      }
    }
  }
  res
}


side_colors <- c("#A52A2A", "#D2691E", "#1122F5", "#FFFF00", "#FFE4C4", "#3CB371", "#A9A9A9", "#4169E1", "#008B8B", "#6B8E23", "#F00FF0", "#FF4500",
             "#11FF50", "#333333", "#7FFFD4", "#ADFF2F", "#FFA500", "#9F1AF0", "#000080?", "#FFE4B5", "#DAA520", "#F4A460", "#BDB76B", "#DB7093",
             "#7FFF00", "#9932CC", "#B8860B", "red", "#006400", "#483D8B", "#FA8072", "#C0C0C0", "#D3D3D3", "#40E0D0", "#F08080", "#FAF0E6", "#DDA0DD",
             "#AFEEEE", "#8B0000", "#FFEFD5", "#0000CD", "#A52A2A", "#D2691E", "#1122F5", "#FFFF00", "#FFE4C4", "#3CB371", "#A9A9A9", "#4169E1", "#008B8B", "#6B8E23", "#F00FF0", "#FF4500",
             "#11FF50", "#333333", "#7FFFD4", "#ADFF2F", "#FFA500", "#9F1AF0", "#000080", "#FFE4B5", "#DAA520", "#F4A460", "#BDB76B", "#DB7093",
             "#7FFF00", "#9932CC", "#B8860B", "red", "#006400", "#483D8B", "#FA8072", "#C0C0C0", "#D3D3D3", "#40E0D0", "#F08080", "#FAF0E6", "#DDA0DD",
             "#AFEEEE", "#8B0000", "#FFEFD5", "#0000CD", "#A52A2A", "#D2691E", "#1122F5", "#FFFF00", "#FFE4C4", "#3CB371", "#A9A9A9", "#4169E1", "#008B8B", "#6B8E23", "#F00FF0", "#FF4500",
             "#11FF50", "#333333", "#7FFFD4", "#ADFF2F", "#FFA500", "#9F1AF0", "#000080", "#FFE4B5", "#DAA520", "#F4A460", "#BDB76B", "#DB7093",
             "#7FFF00", "#9932CC", "#B8860B", "red", "#006400", "#483D8B", "#FA8072", "#C0C0C0", "#D3D3D3", "#40E0D0", "#F08080", "#FAF0E6", "#DDA0DD",
             "#AFEEEE", "#8B0000", "#FFEFD5", "#0000CD", "#A52A2A", "#D2691E", "#1122F5", "#FFFF00", "#FFE4C4", "#3CB371", "#A9A9A9", "#4169E1", "#008B8B", "#6B8E23", "#F00FF0", "#FF4500",
             "#11FF50", "#333333", "#7FFFD4", "#ADFF2F", "#FFA500", "#9F1AF0", "#000080", "#FFE4B5", "#DAA520", "#F4A460", "#BDB76B", "#DB7093",
             "#7FFF00", "#9932CC", "#B8860B", "red", "#006400", "#483D8B", "#FA8072", "#C0C0C0", "#D3D3D3", "#40E0D0", "#F08080", "#FAF0E6", "#DDA0DD",
             "#AFEEEE", "#8B0000", "#FFEFD5", "#0000CD")

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

description <-read.table("input/new_ann/isoforms/terpenes_for_heatmap.iso.selected.txt", sep="\t", row.names = NULL, header = FALSE)
description <-read.table("input/new_ann/isoforms/terpenes_for_heatmap.iso.txt", sep="\t", row.names = NULL, header = FALSE)
description <-read.table("input/new_ann/isoforms/phenilpropanoids_for_heatmap.iso.selected.txt", sep="\t", row.names = NULL, header = FALSE)
description <-read.table("input/new_ann/isoforms/phenilpropanoids_for_heatmap.iso.txt", sep="\t", row.names = NULL, header = FALSE)

colnames(description) <- c("gene_id", "description")
row.names(description) <- description$gene_id
length(unique(description$description))
#description <- frequent_terms(description, 30)
length(unique(description$description))

sub_mat <- subset(mat, rownames(mat) %in% description$gene_id)
clust_map <- aheatmap(sub_mat, Colv = NA, scale = "none")
gene_order <- rownames(sub_mat[clust_map$rowInd,])

sorted_descr <- description[ order(match(description[,1], gene_order)), ] 
sorted_descr <- sorted_descr[order(sorted_descr[,2]),]

new_mat = NULL
genes = NULL
genes_term = NULL
for (i in 1:length(sorted_descr[,1])) {
  if (toString(sorted_descr[i,]$gene_id) %in% row.names(mat)) {
    new_mat <- rbind(new_mat, mat[toString(sorted_descr[i,]$gene_id),])
    genes <- rbind(genes, toString(sorted_descr[i,]$gene_id))
    genes_term <- rbind(genes_term, toString(sorted_descr[i,]$description))
  }
}
row.names(new_mat) <- genes
length(sorted_descr[,1])
length(new_mat[,1])
length(unique(genes_term))

annotation = data.frame(Description = factor(genes_term))
ann_colors = list(Description = side_colors)

h_colors <- get_color_gradient(new_mat)

aheatmap(new_mat, color = h_colors, annRow = annotation, annColors = ann_colors, labRow = row.names(new_mat), labCol = colnames(new_mat),
         Colv=NA, Rowv = NA, legend=TRUE, border_color = "white", cellwidth = 15, cellheight = 8, fontsize = 8, scale="none", 
         filename = "isoforms/terpenes_selected.pdf")

# ================= GO HEATMAPS WITHOUT DESCTIPTION  =========================

description <-read.table("input/new_ann/go_annot_1to1.txt", sep="\t", row.names = NULL, header = TRUE)
colnames(description) <- c("gene_id", "description")

go_genes <- unique(description$gene_id)
length(go_genes)

sub_mat <- subset(mat, rownames(mat) %in% go_genes)
clust_map <- aheatmap(sub_mat, Colv = NA, scale = "none")
genes <- rownames(sub_mat[clust_map$rowInd,])

new_mat = NULL
for (i in 1:length(genes)) {
  if (toString(genes[i]) %in% row.names(mat)) {
    new_mat <- rbind(new_mat, mat[genes[i],])
  }
}
row.names(new_mat) <- genes
length(new_mat[,1])

new_mat <- merge_matrix(new_mat)

top_colors <- 100
bottom_colors <- top_colors * abs(min(new_mat)/max(new_mat))
breaks = seq(min(new_mat),max(new_mat),length.out=top_colors+bottom_colors)
gradient_low = colorpanel( sum( breaks[-1]<=0 ), "blue", "white" )
gradient_top = colorpanel( sum( breaks[-1]>0 ), "white", "red" )
sub_colors = c(gradient_low,gradient_top)

aheatmap(new_mat, color = sub_colors, labRow = row.names(new_mat), labCol = colnames(new_mat),
         Colv=NA, Rowv = NA, legend=TRUE, border_color = "white", cellwidth = 18, cellheight = 3, fontsize = 8, scale="none", 
         filename = "output/go_1to1.pdf")
dev.off()


