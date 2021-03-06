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

#use this only for taking CC and OC taking subsamples
oc_samples <- samples[startsWith(samples, "OC")]
c_samples <- samples[startsWith(samples, "CC")]
sub_samples <- c(cc_samples, oc_samples)
countData <- countData[, sub_samples]

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

# === Choosing DE genes 1 TO 1===
res1vs2 <- glmLRT(fit, coef = 2)
res2vs3 <- glmLRT(fit, contrast=c(0,-1,1,0,0,0,0))
res3vs4 <- glmLRT(fit, contrast=c(0,0,-1,1,0,0,0))
res4vs5 <- glmLRT(fit, contrast=c(0,0,0,-1,1,0,0))
res5vs6 <- glmLRT(fit, contrast=c(0,0,0,0,-1,1,0))
res6vs7 <- glmLRT(fit, contrast=c(0,0,0,0,0,-1,1))
all_res <- c(res1vs2, res2vs3, res3vs4, res4vs5, res5vs6, res6vs7)

deGenes1to1 <- NULL

adj_pvalue_cutoff <- 0.01 #~ same as FDR
lfc_cutoff <- 2

for (res in list(res1vs2, res2vs3, res3vs4, res4vs5, res5vs6, res6vs7)) {
  sigTable <- topTags(res, p.value = adj_pvalue_cutoff,  n = length(row.names(res$table)))
  sigTable <- as.data.frame(sigTable$table)
  sigTable <- sigTable[abs(sigTable$logFC) > lfc_cutoff, ]

  summary(sigTable$FDR)
  summary(sigTable$PValue)
  summary(abs(sigTable$logFC))
  
  deGenes <-sigTable[, c("genes","FDR")]
  deGenes1to1 <- rbind(deGenes1to1, deGenes)
  print(length(deGenes[,1]))
}

length(deGenes1to1[,1])
unique_genes <- unique(deGenes1to1$genes)
for (g in unique_genes) {
  sub_genes <- subset(deGenes1to1, deGenes1to1$genes == toString(g))
  if (length(sub_genes[,2]) > 1) {
    fdr <- min(sub_genes[,2])
    for (i in row.names(sub_genes)) {
      deGenes1to1[i,2] <- fdr
    }
  }
}
deGenes1to1 <- deGenes1to1[unique_genes,]
length(deGenes1to1[,1])


#output DE genes
#write.table(deTable, file="output/DE_cap_genes_data.tsv", quote=FALSE, sep='\t')
write.table(deGenes1to1$genes, file="output_1to1/DE_genes_1to1.list.txt", quote=FALSE, sep='\t', row.names = FALSE)

all_DE_genes <- row.names(deGenes1to1)
top_1to1 <- deGenes1to1[order(deGenes1to1$FDR),]

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
gene_lengths$Length <- pmax(1, gene_lengths$Length - 200)
row.names(gene_lengths) <- gene_lengths$Geneid
gene_lengths <- gene_lengths[row.names(dgList),]

raw_normalized_counts <- rpkm(dgList, normalized.lib.sizes=TRUE, gene.length = gene_lengths$Length)
norm_fact <- colSums(raw_normalized_counts)
norm_fact <- log2(norm_fact / 10^6)

normalized_counts <- rpkm(dgList, normalized.lib.sizes=TRUE, gene.length = gene_lengths$Length, log = TRUE)
for (c in colnames(normalized_counts)) {
  normalized_counts[,c] <- normalized_counts[,c] - norm_fact[c]
}

#select DE genes
#de_normalized_counts <- subset(normalized_counts, startsWith(row.names(normalized_counts), "novel"))
de_normalized_counts <- subset(normalized_counts, row.names(normalized_counts) %in% all_DE_genes)


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
h_colors = get_color_gradient(sorted_mat)

#colors = colorRampPalette(c("blue","white","red"),space="rgb")(256)
aheatmap(sorted_mat, color = h_colors, labCol = colnames(sorted_mat), Colv=NA, labRow = NA, Rowv = draw_tree,
         legend=TRUE, fontsize = 6, scale="none", filename = "TPM/cap_heatmap.pdf")


#=== Supplementary functions for heatmaps ===
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

description <-read.table("input/new_ann/TF_blast_list.txt", sep="\t", row.names = NULL, header = TRUE)
#description <-read.table("input/new_ann/flowering_blast_list.txt", sep="\t", row.names = NULL, header = TRUE)
description <- description[,c("SeqName", "Description")]
colnames(description) <- c("gene_id", "description")

selected_genes <- subset(deGenes1to1, deGenes1to1$genes %in% description$gene_id)
top_selected <- selected_genes[order(selected_genes$FDR),]
#topN <- 100
topN <- length(top_selected$genes)
top_selected <- top_selected[c(1:topN),]
top_selected_genes <- top_selected$genes

row.names(description) <- description$gene_id
description <- description[top_selected_genes,]
head(description)
length(description[,1])
length(unique(description$description))
#description <- frequent_terms(description, 30)
#length(unique(description$description))

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

h_colors = get_color_gradient(new_mat)

aheatmap(new_mat, color = h_colors, annRow = annotation, annColors = ann_colors, 
         labRow = row.names(new_mat), labCol = colnames(new_mat),
         Colv=NA, Rowv = NA, legend=TRUE, fontsize = 3, scale="none",
         filename = "updated_description_out/TF.all.pdf")



# ================= UP-DOWN REGULATION  =========================
stats <- NULL
tf_families <- c("bHLH", "NAC", "WRKY", "MYB", "C2H2", "B3", "MADS")
description <- as.matrix(description)
cutoff <- 0
for (f in tf_families) {
  fam_genes <- row.names(subset(description, str_detect(description[,"description"], f)))
  print(paste(f,toString(length(fam_genes))))
  submat <- mat[fam_genes,]
  fam_stat_up <- NULL
  fam_stat_down <-NULL
  for (s in colnames(mat)) {
    fam_stat_up <- c(fam_stat_up, length(which(submat[,s] >= cutoff)))
    fam_stat_down <- c(fam_stat_down, length(which(submat[,s] <= -cutoff)))
  }
  sub_stat <- rbind(fam_stat_up, fam_stat_down)
  row.names(sub_stat) <-c(paste(f,"up",sep ="_"), paste(f, "down",sep ="_"))
  stats <- rbind(stats, sub_stat)
}
colnames(stats) <- colnames(mat)
write.table(stats, file="updated_description_out/TF_regulated.2.tsv", quote=FALSE, sep='\t', row.names = TRUE)


#TODO
stats1v1 <- NULL
tf_families <- c("bHLH", "NAC", "WRKY", "MYB", "C2H2", "B3", "MADS")
description <- as.matrix(description)
cutoff <- 1
for (f in tf_families) {
  fam_genes <- row.names(subset(description, str_detect(description[,"description"], f)))
  submat <- mat[fam_genes,]
  fam_stat_up <- NULL
  fam_stat_down <-NULL
  for (s in colnames(mat)) {
    fam_stat_up <- c(fam_stat_up, length(which(submat[,s] >= cutoff)))
    fam_stat_down <- c(fam_stat_down, length(which(submat[,s] <= -cutoff)))
  }
  sub_stat <- rbind(fam_stat_up, fam_stat_down)
  row.names(sub_stat) <-c(paste(f,"up",sep ="_"), paste(f, "down",sep ="_"))
  stats <- rbind(stats, sub_stat)
}
colnames(stats) <- colnames(mat)
write.table(stats, file="updated_description_out/TF_regulated.2.tsv", quote=FALSE, sep='\t', row.names = TRUE)
