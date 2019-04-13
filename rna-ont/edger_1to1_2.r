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


deGenes1to1 <- NULL

adj_pvalue_cutoff <- 0.001 #~ same as FDR
lfc_cutoff <- 2

for (n in 1:6) {
  countData = read.table("input/gene_count_matrix.csv", header=TRUE, sep=",", row.names=1 )
  countData <- countData[,c(1:12,17:32)]
  head(countData)
  samples = names(countData)
  
  sample_range <- c(((n-1)*4+1):((n+1)*4))
  print(sample_range)
  countData <- countData[, sample_range]
  
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
  fit <- glmFit(dgList, designMat)
  qlf <- glmLRT(fit, coef = 2)
  
  sigTable <- topTags(qlf, p.value = adj_pvalue_cutoff,  n = length(row.names(qlf$table)))
  sigTable <- as.data.frame(sigTable$table)
  sigTable <- sigTable[abs(sigTable$logFC) > lfc_cutoff, ]
  
  summary(sigTable$FDR)
  summary(sigTable$PValue)
  summary(abs(sigTable$logFC))
  
  deGenes <- row.names(sigTable)
  deGenes1to1 <- c(deGenes1to1, deGenes)
  print(length(deGenes))
}

length(deGenes1to1)
deGenes1to1 <- unique(deGenes1to1)
length(deGenes1to1)

#output DE genes
#write.table(deTable, file="output/DE_cap_genes_data.tsv", quote=FALSE, sep='\t')
write.table(deGenes1to1, file="output_1to1_2/DE_genes.top.txt", quote=FALSE, sep='\t', row.names = FALSE)


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

#normalized counts
normalized_counts <- cpm(dgList, normalized.lib.sizes=TRUE)
#select DE genes
#de_normalized_counts <- subset(normalized_counts, startsWith(row.names(normalized_counts), "novel"))
de_normalized_counts <- subset(normalized_counts, row.names(normalized_counts) %in% deGenes1to1)

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

clust_map <- aheatmap(mat, Colv = NA, scale = "none")
sorted_mat <- mat[clust_map$rowInd,]
sorted_mat <- merge_matrix(sorted_mat)

#colnames(mat) <- sapply(strsplit(colnames(mat), "_"), function(x) x[1])

#color scheme
top_colors <- 100
bottom_colors <- top_colors * abs(min(sorted_mat)/max(sorted_mat))
breaks = seq(min(sorted_mat),max(sorted_mat),length.out=top_colors+bottom_colors)
gradient_low = colorpanel( sum( breaks[-1]<=0 ), "blue", "white" )
gradient_top = colorpanel( sum( breaks[-1]>0 ), "white", "red" )
h_colors = c(gradient_low,gradient_top)

#colors = colorRampPalette(c("blue","white","red"),space="rgb")(256)
aheatmap(sorted_mat, color = h_colors, labCol = colnames(sorted_mat), Colv=NA, Rowv = NA, labRow = NA, 
         legend=TRUE, fontsize = 6, scale="none", filename = "output_1to1_2/cap_heatmap.pdf")
dev.off()


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

description <-read.table("input/new_ann/new_degs_tf.gene_id.txt", sep="\t", row.names = NULL, header = FALSE)
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

new_mat <- merge_matrix(new_mat)
length(unique(genes_term))

annotation = data.frame(Descritpion = factor(genes_term))
ann_colors = list(Descritpion = side_colors)

top_colors <- 100
bottom_colors <- top_colors * abs(min(new_mat)/max(new_mat))
breaks = seq(min(new_mat),max(new_mat),length.out=top_colors+bottom_colors)
gradient_low = colorpanel( sum( breaks[-1]<=0 ), "blue", "white" )
gradient_top = colorpanel( sum( breaks[-1]>0 ), "white", "red" )
sub_colors = c(gradient_low,gradient_top)

aheatmap(new_mat, color = sub_colors, annRow = annotation, annColors = ann_colors, labRow = row.names(new_mat), labCol = colnames(new_mat),
         Colv=NA, Rowv = NA, legend=TRUE, fontsize = 6, scale="none", 
         filename = "output_1to1/TF_1to1.top100.pdf")
dev.off()


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


