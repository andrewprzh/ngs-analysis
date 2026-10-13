library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library (NMF)
library(RColorBrewer)
library(gridExtra)

theme_set(theme_bw())


#  B1A_1_HiFi      B1A_2_HiFi      B1A_3_HiFi      B1A
#  id1b_1_HiFi     id1b_2_HiFi     id1b_3_HiFi     1B
#  id72_1_HiFi     id72_2_HiFi     id72_3_HiFi     HiFi control
#   id7_1_HiFi      id7_2_HiFi      id7_3_HiFi     ID7
#   B1A_1_Revio     B1A_2_Revio     B1A_3_Revio    B1A Revio
#  C_1_Revio       C_2_Revio          C_3_Revio    Revio control
#
# No batch correction; every comparison is done within a single platform.


plot_pca <- function(countData, samplesData, out_prefix) {
  # general PCA
  samples = colnames(countData)
  condition = factor(samplesData[samples, "group"])
  batch = factor(samplesData[samples, "batch"])
  colData = data.frame(samples=samples, condition=condition, batch=batch)
  ddsMatAll = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)
  ddsMatAll <- estimateSizeFactors(ddsMatAll)
  ddsMatAll <- estimateDispersions(ddsMatAll)
  normalized_counts_all <- counts(ddsMatAll, normalized = TRUE)
  stabilized_counts_all <-getVarianceStabilizedData(ddsMatAll)

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
                   PC3 = pca$x[, "PC3"],
                   condition = condition,
                   batch = batch)

  pdf(file = paste0(out_prefix, "_PCA.pdf"), width = pdfWidth, height = pdfHeight)
  plot1 <- ggplot(pca_df) +
    geom_point(aes(PC1, PC2, colour = condition, shape = batch), size = 3) +
    labs(x = paste0("PC1 (", round(summary(pca)$importance[2, 1]*100, 2), "%)"),
         y = paste0("PC2 (", round(summary(pca)$importance[2, 2]*100, 2), "%)")) +
    scale_colour_brewer(palette = "Set1") +
    geom_text_repel(data = pca_df,
                    aes(PC1, PC2, label = sample),
                    min.segment.length = unit(0.5, "lines")) +
    theme(legend.position = "top")
  grid.arrange(plot1, ncol = 1)
  dev.off()
}

pdfWidth=11
pdfHeight=8;

base_dir <- "~/ablab/analysis/UK/DIE/oct26/"
samplesData = read.table(paste0(base_dir, "samples.tsv"), header=TRUE, sep="\t", row.names=1, check.names=FALSE)

# folder -> count table prefix (followed by "<feature>_grouped_file_name_counts.tsv")
tables <- list(
  "ensemble94"           = "ADAR.PB_Revio.Ensemble94.",
  "ensemble94_extended"  = "ADAR.PB_Revio.Ensemble94.discovered_",
  "ensemble113"          = "ADAR.PB_Revio.Ensemble113.",
  "ensemble113_extended" = "ADAR.PB_Revio.Ensemble113.discovered_"
)

# comparison name -> list(KO samples, control samples)
reps <- 1:3
experiments <- list(
  "B1A_Revio_vs_C_Revio"   = list(paste0("B1A_", reps, "_Revio"),  paste0("C_", reps, "_Revio")),
  "B1A_HiFi_vs_id72_HiFi"  = list(paste0("B1A_", reps, "_HiFi"),   paste0("id72_", reps, "_HiFi")),
  "id1b_HiFi_vs_id72_HiFi" = list(paste0("id1b_", reps, "_HiFi"),  paste0("id72_", reps, "_HiFi")),
  "id7_HiFi_vs_id72_HiFi"  = list(paste0("id7_", reps, "_HiFi"),   paste0("id72_", reps, "_HiFi"))
)

#cutoffs
p_adj_cutoff <- 0.05
log2_cutoff <- 1

for (folder in names(tables)) {
  wd <- paste0(base_dir, folder, "/")
  setwd(wd)

  for (feature_name in c("gene", "transcript")) {
    rawData = read.table(paste0(tables[[folder]], feature_name, "_grouped_file_name_counts.tsv"),
                         header=TRUE, sep="\t", row.names=1, check.names=FALSE)
    rawData <- round(as.matrix(rawData))
    storage.mode(rawData) <- "integer"
    stopifnot(all(colnames(rawData) %in% row.names(samplesData)))

    hifi_samples  <- colnames(rawData)[samplesData[colnames(rawData), "batch"] == "HiFi"]
    revio_samples <- colnames(rawData)[samplesData[colnames(rawData), "batch"] == "Revio"]
    plot_pca(rawData, samplesData, paste0(folder, "_", feature_name, "_all_raw_data"))
    plot_pca(rawData[, hifi_samples], samplesData, paste0(folder, "_", feature_name, "_HiFi_raw_data"))
    plot_pca(rawData[, revio_samples], samplesData, paste0(folder, "_", feature_name, "_Revio_raw_data"))

    # === RUN DE analysis ===
    for (experiment in names(experiments)) {
      exp_name = paste0(experiment, ".", folder)
      samples = c(experiments[[experiment]][[1]], experiments[[experiment]][[2]])
      countData <- rawData[, samples]

      # Build the data frame from the conditions
      condition = factor(samplesData[samples, "type"])
      colData = data.frame(samples=samples, condition=condition)
      message("=== ", exp_name, " / ", feature_name)
      print(colData)

      # create DeSeq object
      ddsMat = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)
      #Set the reference to be compared
      ddsMat$condition = relevel(ddsMat$condition,"CT")
      ddsMat <- DESeq(ddsMat)

      # normalize data
      counts_normalized <- counts(ddsMat, normalized=TRUE)
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

      #volcano plot
      pdf(file = paste0(exp_name, "_", feature_name, "s_volcano.pdf"), width = pdfWidth, height = pdfHeight)
      p <- res %>%
        filter(!is.na(padj)) %>%
        ggplot() +
        geom_point(aes(log2FoldChange, -log10(padj), colour = (abs(log2FoldChange) > log2_cutoff & padj < p_adj_cutoff))) +
        theme_bw() +
        scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), name = "DE")
      print(p)
      dev.off()

      #select DE genes
      res_sorted <- res[order(res$padj), ]
      DGEgenes <-  subset ( res_sorted , padj < p_adj_cutoff & abs(log2FoldChange) > log2_cutoff)$gene_id
      message(exp_name, " ", feature_name, "s: ", length(DGEgenes), " DE")
      write_tsv(subset ( res_sorted , gene_id %in% DGEgenes) %>%
                  select(gene_id,log2FoldChange,padj), paste0(exp_name, "_DE_", feature_name, "s.tsv"))

      matDGEgenes <- logNormCounts[rownames(logNormCounts) %in% DGEgenes, , drop = FALSE]
      if (nrow(matDGEgenes) >= 2) {
        pdf(file = paste0(exp_name, "_", feature_name, "s_raw_heatmap.pdf"), width = pdfWidth, height = pdfHeight)
        aheatmap(matDGEgenes, Rowv = TRUE , Colv = TRUE , distfun = "euclidean" , hclustfun = "average" )
        dev.off()
      }

      # z-scores of normalized counts for DE features
      vals = counts_normalized[rownames(counts_normalized) %in% DGEgenes, , drop = FALSE]

      # Adds a little noise to each element
      # To avoid the clusteing function failing on zero
      # variance datalines.
      vals = jitter(vals, factor = 1, amount=0.00001)
      vals = vals[apply(vals, 1, sd) > 0, , drop = FALSE]
      mat = if (nrow(vals) > 0) t(scale(t(vals))) else vals
      gene = data.frame(gene_id = rownames(mat))

      if (nrow(mat) >= 2) {
        pdf(file = paste0(exp_name, "_zscore_", feature_name, "s_heatmap.pdf"), width = pdfWidth, height = pdfHeight)
        aheatmap(mat, Colv = NA)
        dev.off()
      }
      write_tsv(data.frame(gene, mat, check.names = FALSE),  paste0(exp_name, "_zscore_", feature_name, "s_heatmap.tsv"))
    }
  }
}
