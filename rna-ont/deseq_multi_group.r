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
samplesData = read.table("input/pheno_data.tsv", header=TRUE, sep="\t", row.names=1 )
head(samplesData)
samplesData <- samplesData[order(samplesData$pc),]

countData = read.table("input/gene_count_matrix.csv", header=TRUE, sep=",", row.names=1 )
countData <- countData[,row.names(samplesData)]
samples = names(countData)
condition = factor(samplesData[row.names(samplesData) %in% samples,][[4]])

# Build the dataframe from the conditions


colData = data.frame(samples=samples, condition=condition)
summary(colData)

# create DeSeq object
ddsMat = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)

#Set the reference to be compared
ddsMat$condition = relevel(ddsMat$condition,"9")

ddsMat <- DESeq(ddsMat,  test = "LRT", reduced = ~ 1)

# get results
unique_conds <- unique(condition)
total_conditions <- length(unique_conds)
pairwise_res <- NULL
for (i in 2:length(unique_conds)) {
  pairwise_res <- rbind(pairwise_res, c(results = results(ddsMat, contrast=c("condition","9",toString(unique_conds[[i]])))))
}
row.names(pairwise_res) <- unique_conds[c(2:total_conditions)]

all_res = NULL
def_cond <- "14"
def_res <- as.data.frame(pairwise_res[def_cond,])
gene_list <- row.names(def_res)
cmp_conds <- sapply(unique_conds[c(2:total_conditions)], function(x) toString(x))
for (i in 1:length(gene_list)) {
  max_lfc <- 0
  min_lfc <- 0
  for (cond in cmp_conds) {
    lfc <- as.data.frame(pairwise_res[cond,])$log2FoldChange[i]
    max_lfc <- max(max_lfc, lfc)
    min_lfc <- min(min_lfc, lfc)
  }
  res_lfc <- max_lfc - min_lfc
  all_res <- rbind(all_res, c(def_res$baseMean[i],
                              res_lfc, 
                              def_res$pvalue[i],
                              def_res$padj[i]))
  if (i %% 100 == 0) {
    print(i)
  }
}
colnames(all_res) <- c("baseMean","log2FoldChange","pvalue","padj")
row.names(all_res) <- gene_list
all_res <- as.data.frame(all_res)

# === Choosing DE genes ===
adj_pvalue_cutoff <- 0.01 
lfc_cutoff <- 2

DEres <- subset(all_res, !is.na(padj) & padj < adj_pvalue_cutoff & log2FoldChange > lfc_cutoff)
summary(DEres$pvalue)
summary(DEres$padj)

deGenesDeseq <- row.names(DEres)
deGenesDeseq <- row.names(DEres)
length(deGenesDeseq)

#output DE genes
#write.table(deTable, file="output/DE_cap_genes_data.tsv", quote=FALSE, sep='\t')
write.table(deGenesDeseq, file="output/deseq/DE_genes_stringtie.txt", quote=FALSE, sep='\t', row.names = FALSE)
