#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")
library(edgeR)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
#if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

#BiocManager::install("Biobase")
library(NMF)
library(gplots)
library(RColorBrewer)
library(pheatmap)

# === DE analysis ===


#read count data
for (i in 1:3) {
  for (j in (i+1):4) {
    print(i)
    print(j)
    countData = read.table("~/ablab/analysis/italy/dipsacus/ST_new/Dipsacus.SUP.gene_grouped_counts.tsv", header=TRUE, sep="\t", row.names=1 )
    samples = names(countData)
    all_cons = c("LCD", "LAL","LAD","LCL")
    
    sub_samples = samples[c(i,j)]
    countData <- countData[, sub_samples]
    
    
    #use this only for taking subsamples
    #samples1 <- samples[startsWith(samples, "sample.1")]
    #samples2 <- samples[startsWith(samples, "sample.2")]
    #samples3 <- samples[startsWith(samples, "sample.3")]
    #samples4 <- samples[startsWith(samples, "sample.4")]
    #sub_samples <- c(samples3, samples4)
    #
    
    # head(countData)
    #samples = names(countData)
    
    #create condition according 
    condition = factor(all_cons[c(i,j)])
    
    designMat <- model.matrix(~condition)
    
    # DGE list for EdgeR
    dgList <- DGEList(counts=countData, genes=rownames(countData), group = condition)
    dgList <- calcNormFactors(dgList, method="TMM")
    
    #dispersion
    housekeeping = as.vector(read.table("~/ablab/analysis/italy/dipsacus/ST_new/housekeeping_genes.tsv", header=FALSE, sep="\t"))
    head(housekeeping)
    
    dgListD <- dgList
    dgListD$samples$group <- 1
    hk_index <- which(row.names(dgListD) %in% housekeeping[[1]])
    # head(dgListD[hk_index,])
    dgListD <- estimateDisp(dgListD[hk_index,], trend="none", tagwise=FALSE)
    dgList$common.dispersion <- dgListD$common.dispersion
    
    #de analysis
    total_conditions <- length(unique(condition))
    fit <- glmFit(dgList, designMat)
    lrt <- glmLRT(fit)
    
    
    # === Check plots ===
    plot_adj_pvalue_cutoff <- 0.05 #~ same as FDR. Originariamente 0.01
    plot_lfc_cutoff <- 2
    
    simpleTable <- lrt$table
    simpleTable <- as.matrix(simpleTable)
    summary(simpleTable)
    
    #plot volcano
    plot(simpleTable[,"logFC"], -log10(simpleTable[,"PValue"]),
         xlim=c(-5, 5), ylim=c(0, 5),
         xlab="log2 fold change", ylab="-log10 p-value",
         type="n")
    # then add the points
    sel <- which(abs(simpleTable[,"logFC"]) <= plot_lfc_cutoff | simpleTable[,"PValue"] >= plot_adj_pvalue_cutoff) 
    points(simpleTable[sel,"logFC"], -log10(simpleTable[sel,"PValue"]),col="black")
    sel <- which(abs(simpleTable[,"logFC"]) > plot_lfc_cutoff & simpleTable[,"PValue"] < plot_adj_pvalue_cutoff) 
    points(simpleTable[sel,"logFC"], -log10(simpleTable[sel,"PValue"]),col="red")
    
    #plot smear
    plot(simpleTable[,"logCPM"], simpleTable[,"logFC"],
         xlim=c(-10, 10), ylim=c(0, 15),
         xlab="logCPM", ylab="log2 fold change",
         type="n")
    # then add the points
    sel <- which(simpleTable[,"PValue"] >= plot_adj_pvalue_cutoff) 
    points(simpleTable[sel,"logCPM"], simpleTable[sel,"logFC"],col="black")
    sel <- which(simpleTable[,"PValue"] < plot_adj_pvalue_cutoff)
    points(simpleTable[sel,"logCPM"], simpleTable[sel,"logFC"],col="red")
    
    # === Choosing DE genes ===
    adj_pvalue_cutoff <- 0.05 #~ same as FDR. Originariamente 0.01
    lfc_cutoff <- 2
    
    sigTable <- topTags(lrt, p.value = adj_pvalue_cutoff, n = length(row.names(lrt$table)), adjust.method = "BH")
    sigTable <- sigTable$table
    deTable <- sigTable[which(abs(sigTable["logFC"]) > lfc_cutoff),]
    # head(deTable)
    summary(deTable$PValue)
    deGenes <- row.names(deTable)
    length(deGenes)
    #output DE table
    write.table(deTable, file=paste0("~/ablab/analysis/italy/dipsacus/ST_new/1vs1/DEG_table.list.", i, "_vs_", j, ".lfc2_p005.txt"), quote=FALSE, sep='\t', row.names = FALSE)
    #output DE genes
    write.table(deGenes, file=paste0("~/ablab/analysis/italy/dipsacus/ST_new/1vs1/DE_genes.list.", i, "_vs_", j, ".lfc2_p005.txt"), quote=FALSE, sep='\t', row.names = FALSE)
  }
}


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

de_normalized_counts <- subset(countData, row.names(countData) %in% deGenes)
vals = as.matrix(de_normalized_counts)
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
row.names(score) = row.names(de_normalized_counts)

# Generate heatmap
mat = as.matrix(score)
head(mat)
nrow(mat)
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
         legend=TRUE, fontsize = 7, scale="none",
         filename = "all_conditions_heatmap_p001.pdf")
dev.off()
