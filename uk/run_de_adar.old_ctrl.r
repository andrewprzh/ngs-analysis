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
library(sva)
library(gridExtra)

theme_set(theme_bw())


#  B1A_1_HiFi      B1A_2_HiFi      B1A_3_HiFi      B1A
#  id1b_1_HiFi     id1b_2_HiFi     id1b_3_HiFi     1B
#  id72_1_HiFi     id72_2_HiFi     id72_3_HiFi     OLD coltrol
#   id7_1_HiFi      id7_2_HiFi      id7_3_HiFi     ID7  
#   B1A_1_Revio     B1A_2_Revio     B1A_3_Revio    B1A for batch effect  
#  C_1_Revio       C_2_Revio          C_3_Revio New control


plot_pca <- function(countData, samplesData, out_prefix) { 
  # general PCA
  samples = colnames(countData)
  condition = factor(samplesData[][[1]])
  colData = data.frame(samples=samples, condition=condition)
  summary(colData)
  ddsMatAll = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)
  ddsMatAll <- estimateSizeFactors(ddsMatAll)
  ddsMatAll <- estimateDispersions(ddsMatAll)
  normalized_counts_all <- counts(ddsMatAll, normalized = TRUE)
  stabilized_counts_all <-getVarianceStabilizedData(ddsMatAll)
  
  norm_df <- mutate(as.data.frame(normalized_counts_all), gene_id = row.names(normalized_counts_all)) 
  logNormCountsAll <- log2(normalized_counts_all + 1)
  
  logDist <- as.dist (1 - cor( logNormCountsAll , method = "pearson" ) )
  pdf(file = paste0(out_prefix, "_hclust.pdf"), width = pdfWidth, height = pdfHeight) 
  plot(hclust(logDist),labels=colnames(logNormCountsAll),main=" log transformed read counts distance : Pearson correlation ")
  dev.off()
  
  # plot PCA
  pca <- prcomp(t(stabilized_counts_all))
  pca_df <- tibble(sample = row.names(pca$x),
                   PC1 = pca$x[, "PC1"], 
                   PC2 = pca$x[, "PC2"], 
                   PC3 = pca$x[, "PC3"]) 
  
  #  require(gridExtra)
  pdf(file = paste0(out_prefix, "_PCA.pdf"), width = pdfWidth, height = pdfHeight) 
  #plotPCA(vst(ddsMatAll, blind = FALSE), intgroup = "condition")
  
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
    
  dev.off()
}

pdfWidth=11
pdfHeight=8;
for (ann in c("Ensemble", "RefSeq")) {
  folder_prefix = tolower(ann)
  for (folder in c(folder_prefix, paste0(folder_prefix, "_novel"))) {
    wd <- paste0("~/ablab/analysis/UK/DIE/revio_old_control//", folder, "/")
    setwd(wd)
    samplesData = read.table("../new_samples.tsv", header=TRUE, sep="\t", row.names=1 )
    
    for (feature_name in c("gene", "transcript")) {
      rawData = read.table(paste0("ADAR.PB.", ann, ".", feature_name, "_grouped_counts.tsv.header.tsv"), header=TRUE, sep="\t", row.names=1 )
      plot_pca(rawData, samplesData, paste0(folder, "_", feature_name, "_raw_data"))
      combat_counts <- rawData
      
      #combat_counts <- ComBat_seq(
      #  as.matrix(rawData),
      #  batch = samplesData$batch,
      #  group = samplesData$group 
      #)
      #plot_pca(combat_counts, samplesData, paste0(folder, "_", feature_name, "_batch_corrected"))

      # === RUN DE analysis ===
      for (experiment in c(1, 2, 4, 5)) {
        if (experiment == 1) {
          exp_name = paste0("B1A_vs_WT.", folder)
          countData <- combat_counts[,c(1,2,3,7,8,9)]
        }
        if (experiment == 2) {
          exp_name = paste0("1B_vs_WT.", folder)
          countData <- combat_counts[,c(4,5,6,7,8,9)]
        }
        if (experiment == 3) {
          exp_name = paste0("WT_vs_Ctrl.", folder)
          countData <- combat_counts[,c(7,8,9,16,17,18)]
        }
        if (experiment == 4) {
          exp_name = paste0("ID7_vs_WT.", folder)
          countData <- combat_counts[,c(10,11,12,7,8,9)]
        }
        if (experiment == 5) {
          exp_name = paste0("B1A_new_vs_WT.", folder)
          countData <- combat_counts[,c(13,14,15,7,8,9)]
        }

        samples = colnames(countData)
        head(samplesData)
        
        # Build the data frame from the conditions
        condition = factor(samplesData[row.names(samplesData) %in% samples,][[2]])
        colData = data.frame(samples=samples, condition=condition)
        summary(colData)
        
        # create DeSeq object
        ddsMat = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)
        #Set the reference to be compared
        ddsMat$condition = relevel(ddsMat$condition,"CT")
        ddsMat <- DESeq(ddsMat)
        
        # normalize data
        counts_normalized <- counts(ddsMat, normalized=TRUE)
        counts_stabilized <-getVarianceStabilizedData(ddsMat)
        norm_df <- mutate(as.data.frame(counts_normalized), gene_id = row.names(counts_normalized)) 
        #write_tsv(norm_df, "deseq_normalization.tsv")
        logNormCounts <- log2(counts_normalized + 1)

        # get DE results
        res <- results(ddsMat)
        gene_id <- row.names(res)
        res <- as.data.frame(res)
        res <- cbind(gene_id, res)
        res <- res %>%
          filter(!is.na(stat)) %>%
          arrange(desc(stat)) %>%
          select(gene_id, everything())
        write_tsv(res,  paste0(exp_name, "_deseq_", feature_name, "s_results.tsv"))
        
        #cutoffs
        p_adj_cutoff <- 0.05
        log2_cutoff <- 1
        
        #volcano plot
        pdf(file = paste0(exp_name, "_", feature_name, "s_volcano.pdf"), width = pdfWidth, height = pdfHeight) 
        res %>%
          filter(!is.na(padj)) %>%
          ggplot() +
          geom_point(aes(log2FoldChange, -log10(padj), colour = (abs(log2FoldChange) > log2_cutoff & padj < p_adj_cutoff))) +
          theme_bw() +
          scale_color_manual(values = c("black", "red")) +
          geom_text_repel(data = filter(res, padj < p_adj_cutoff),
                          aes(log2FoldChange, -log10(padj), label = ""),
                          min.segment.length = unit(0.1, "lines"))
        dev.off()
        
        #select DE genes
        res_sorted <- res[order(res$padj), ]
        DGEgenes <-  subset ( res_sorted , padj < p_adj_cutoff & abs(log2FoldChange) > log2_cutoff)$gene_id 
        matDGEgenes <- subset(logNormCounts, gene_id %in% DGEgenes)
        write_tsv(subset ( res_sorted , gene_id %in% DGEgenes) %>% 
                    select(gene_id,log2FoldChange,padj), paste0(exp_name, "_DE_", feature_name, "s.tsv")) 
        
        pdf(file = paste0(exp_name, "_", feature_name, "s_raw_heatmap.pdf"), width = pdfWidth, height = pdfHeight) 
        aheatmap(matDGEgenes, Rowv = TRUE , Colv = TRUE , distfun = "euclidean" , hclustfun = "average" )
        dev.off()
        
        # Get normalized counts and write this to a file
        normalized_counts = counts(ddsMat,normalized=TRUE)
        
        # Turn it into a dataframe to have proper column names.
        normalized_counts_dt = data.frame("gene_id"=rownames(normalized_counts),normalized_counts)
        gene = subset(normalized_counts_dt, gene_id %in% DGEgenes)[1]
        vals = as.matrix(subset(normalized_counts_dt, gene_id %in% DGEgenes)[2:ncol(normalized_counts_dt)])
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
        pdf(file = paste0(exp_name, "_zscore_", feature_name, "s_heatmap.pdf"), width = pdfWidth, height = pdfHeight) 
        aheatmap(mat, Colv = NA)
        #plot(hclust(logDist),labels=colnames(logNormCounts),main=" log transformed read counts distance : Pearson correlation ")
        #colors = colorRampPalette(c("blue","black","red"),space="rgb")(256)
        #heatmap.2( mat, col=colors,density.info="none",trace="none", margins=c(10,19), lhei=c(1,7), Colv=NA)
        dev.off()
        write_tsv(data.frame(gene,mat),  paste0(exp_name, "_zscore_", feature_name, "s_heatmap.tsv"))
      }
    }
  }
}  


