library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(NMF)
library(gplots)
library(RColorBrewer)
library(pheatmap)

theme_set(theme_bw())

# Ensure dplyr functions are used
rename <- dplyr::rename
select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter
arrange <- dplyr::arrange

# Function to plot PCA
plot_pca <- function(transformed_data, samplesData, out_prefix, title_suffix = "") {
  samples <- colnames(transformed_data)
  group <- factor(samplesData[samples, ]$group)
  
  # Hierarchical clustering
  logDist <- as.dist(1 - cor(transformed_data, method = "pearson"))
  
  pdf(file = paste0(out_prefix, "_hclust.pdf"), width = pdfWidth, height = pdfHeight)
  plot(hclust(logDist), labels = colnames(transformed_data),
       main = paste0("log transformed counts distance: Pearson correlation", title_suffix))
  dev.off()
  
  # PCA
  pca <- prcomp(t(transformed_data))
  pca_df <- tibble(
    sample = row.names(pca$x),
    PC1 = pca$x[, "PC1"],
    PC2 = pca$x[, "PC2"],
    PC3 = pca$x[, "PC3"],
    condition = group
  )
  
  require(gridExtra)
  pdf(file = paste0(out_prefix, "_PCA.pdf"), width = pdfWidth, height = pdfHeight)
  plot1 <- ggplot(pca_df) +
    geom_point(aes(PC1, PC2, colour = condition), size = 3) +
    labs(
      title = title_suffix,
      x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")
    ) +
    scale_colour_brewer(palette = "Set1") +
    geom_text_repel(
      data = pca_df,
      aes(PC1, PC2, label = sample),
      min.segment.length = unit(0.5, "lines")
    ) +
    theme(legend.position = "top")
  grid.arrange(plot1, ncol = 1)
  dev.off()
}

pdfWidth <- 11
pdfHeight <- 8

for (ann in c("Ensemble", "RefSeq")) {
  folder_prefix <- tolower(ann)
  
  for (folder in c(folder_prefix, paste0(folder_prefix, "_novel"))) {
    wd <- paste0("~/ablab/analysis/UK/DIE/revio/", folder, "/")
    setwd(wd)
    
    samplesData <- read.table("../new_samples.tsv", 
                              header = TRUE, sep = "\t", row.names = 1)
    
    for (feature_name in c("gene", "transcript")) {
      cat("\n=== Processing:", ann, folder, feature_name, "===\n")
      
      # Read raw data - ALL samples
      rawData <- read.table(
        paste0("ADAR.PB.", ann, ".", feature_name, "_grouped_counts.tsv.header.tsv"),
        header = TRUE, sep = "\t", row.names = 1
      )
      
      # Create colData with ALL samples
      samples <- colnames(rawData)
      group <- factor(samplesData[samples, ]$group)
      batch <- factor(samplesData[samples, ]$batch)
      
      colData <- data.frame(
        row.names = samples,
        group = group,
        batch = batch
      )
      
      # Create DESeq2 dataset with batch in design
      dds <- DESeqDataSetFromMatrix(
        countData = rawData,
        colData = colData,
        design = ~ batch + group
      )
      
      # Pre-filtering: remove low count features
      keep <- rowSums(counts(dds)) >= 10
      dds <- dds[keep, ]
      cat("Kept", sum(keep), "out of", length(keep), "features after filtering\n")
      
      # Run DESeq2
      dds <- DESeq(dds)
      
      # Get transformed data for visualization
      vsd <- vst(dds, blind = FALSE)
      
      # Plot PCA with batch in model
      plot_pca(assay(vsd), samplesData,
               paste0(folder, "_", feature_name, "_raw_data"),
               " (batch in model)")
      
      # Remove batch effect for visualization only
      mat <- assay(vsd)
      mat_corrected <- limma::removeBatchEffect(mat, batch = batch, 
                                                design = model.matrix(~ group))
      plot_pca(mat_corrected, samplesData,
               paste0(folder, "_", feature_name, "_batch_corrected"),
               " (batch corrected)")
      
      # === RUN DE analysis ===
      
      # Define comparisons
      comparisons_info <- list(
        list(name = "B1A_vs_Ctrl", test = "B1A", control = "CTRL"),
        list(name = "1B_vs_Ctrl", test = "1B", control = "CTRL"),
        list(name = "ID7_vs_Ctrl", test = "ID7", control = "CTRL")
      )
      
      # Process each comparison
      for (comp_info in comparisons_info) {
        comp_name <- comp_info$name
        test_group <- comp_info$test
        control_group <- comp_info$control
        exp_name <- paste0(comp_name, ".", folder)
        
        cat("\n--- Processing comparison:", comp_name, "---\n")
        
        # Check if both groups exist
        if (!(test_group %in% levels(group) && control_group %in% levels(group))) {
          cat("Warning: Groups", test_group, "or", control_group, "not found. Skipping.\n")
          next
        }
        
        # Get results for this contrast
        res <- results(dds, contrast = c("group", test_group, control_group))
        
        # Convert to data frame and clean up
        res_df <- as.data.frame(res) %>%
          mutate(gene_id = row.names(res)) %>%
          rename(
            log2FoldChange = log2FoldChange,
            pvalue = pvalue,
            padj = padj,
            stat = stat
          ) %>%
          select(gene_id, log2FoldChange, baseMean, stat, pvalue, padj, lfcSE) %>%
          filter(!is.na(stat)) %>%
          arrange(pvalue)
        
        write_tsv(res_df, paste0(exp_name, "_deseq2_", feature_name, "s_results.tsv"))
        
        # Cutoffs
        p_adj_cutoff <- 0.05
        log2_cutoff <- 1
        
        # Count DE features
        n_de <- sum(res_df$padj < p_adj_cutoff & abs(res_df$log2FoldChange) > log2_cutoff, 
                    na.rm = TRUE)
        cat("Found", n_de, "DE", feature_name, "s\n")
        
        # Volcano plot
        pdf(file = paste0(exp_name, "_", feature_name, "s_volcano.pdf"), 
            width = pdfWidth, height = pdfHeight)
        print(
          res_df %>%
            filter(!is.na(padj)) %>%
            ggplot() +
            geom_point(aes(log2FoldChange, -log10(padj), 
                           colour = (abs(log2FoldChange) > log2_cutoff & padj < p_adj_cutoff))) +
            theme_bw() +
            scale_color_manual(values = c("black", "red")) +
            labs(title = paste0(comp_name, " - ", feature_name, "s")) +
            geom_text_repel(
              data = filter(res_df, padj < p_adj_cutoff & abs(log2FoldChange) > log2_cutoff),
              aes(log2FoldChange, -log10(padj), label = gene_id),
              min.segment.length = unit(0.1, "lines"),
              max.overlaps = 20
            )
        )
        dev.off()
        
        # Select DE genes
        DGEgenes <- res_df %>%
          filter(padj < p_adj_cutoff & abs(log2FoldChange) > log2_cutoff) %>%
          pull(gene_id)
        
        if (length(DGEgenes) == 0) {
          cat("No significant DE", feature_name, "s found\n")
          next
        }
        
        # Get normalized counts for ALL samples
        normalized_counts <- counts(dds, normalized = TRUE)
        logNormCounts <- log2(normalized_counts + 1)
        matDGEgenes <- logNormCounts[DGEgenes, , drop = FALSE]
        
        # Write DE genes
        res_df %>%
          filter(gene_id %in% DGEgenes) %>%
          select(gene_id, log2FoldChange, padj) %>%
          write_tsv(paste0(exp_name, "_DE_", feature_name, "s.tsv"))
        
        # Raw heatmap (using all samples)
        if (length(DGEgenes) > 1) {
          pdf(file = paste0(exp_name, "_", feature_name, "s_raw_heatmap.pdf"), 
              width = pdfWidth, height = pdfHeight)
          aheatmap(matDGEgenes, Rowv = TRUE, Colv = TRUE, 
                   distfun = "euclidean", hclustfun = "average",
                   main = paste0(comp_name, " - All samples"))
          dev.off()
          
          # Z-score heatmap (using all samples)
          vals <- matDGEgenes
          vals <- jitter(vals, factor = 1, amount = 0.00001)
          
          zscore <- t(scale(t(vals)))
          
          pdf(file = paste0(exp_name, "_zscore_", feature_name, "s_heatmap.pdf"), 
              width = pdfWidth, height = pdfHeight)
          aheatmap(zscore, Colv = NA, 
                   main = paste0(comp_name, " - Z-score (all samples)"))
          dev.off()
          
          zscore_df <- data.frame(gene_id = rownames(zscore), zscore, check.names = FALSE)
          write_tsv(zscore_df, paste0(exp_name, "_zscore_", feature_name, "s_heatmap.tsv"))
        } else if (length(DGEgenes) == 1) {
          cat("Only 1 DE gene, skipping heatmap\n")
        }
      }
    }
  }
}

cat("\n=== Analysis complete ===\n")