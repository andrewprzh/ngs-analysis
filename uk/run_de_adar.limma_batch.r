library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(limma)
library(edgeR)
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
plot_pca <- function(logExpr, samplesData, out_prefix, title_suffix = "") {
  samples <- colnames(logExpr)
  group <- factor(samplesData[samples, ]$group)
  
  # Hierarchical clustering
  logDist <- as.dist(1 - cor(logExpr, method = "pearson"))
  
  pdf(file = paste0(out_prefix, "_hclust.pdf"), width = pdfWidth, height = pdfHeight)
  plot(hclust(logDist), labels = colnames(logExpr),
       main = paste0("log transformed expression distance: Pearson correlation", title_suffix))
  dev.off()
  
  # PCA on log-expression data
  pca <- prcomp(t(logExpr))
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
    wd <- paste0("~/ablab/analysis/UK/DIE/revio_limma//", folder, "/")
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
      
      # Create DGEList with ALL samples
      group <- factor(samplesData$group)
      batch <- factor(samplesData$batch)
      
      dge <- DGEList(counts = rawData, group = group)
      
      # Filter low expressed features
      keep <- filterByExpr(dge, group = group)
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      cat("Kept", sum(keep), "out of", length(keep), "features after filtering\n")
      
      # Calculate normalization factors
      dge <- calcNormFactors(dge)
      
      # Design matrix with batch effect - ALL samples
      design <- model.matrix(~0 + group + batch)
      colnames(design) <- gsub("group", "", colnames(design))
      
      # voom transformation on ALL samples
      v <- voom(dge, design, plot = FALSE)
      
      # Plot PCA for raw data (with batch in model but not removed)
      plot_pca(v$E, samplesData, 
               paste0(folder, "_", feature_name, "_raw_data"),
               " (raw, batch in model)")
      
      # Remove batch effect for visualization only
      v_corrected <- removeBatchEffect(v$E, batch = batch, design = design)
      plot_pca(v_corrected, samplesData, 
               paste0(folder, "_", feature_name, "_batch_corrected"),
               " (batch corrected)")
      
      # === RUN DE analysis on ALL samples ===
      
      # Fit linear model on ALL samples
      fit <- lmFit(v, design)
      
      # Define all contrasts
      # Note: R prefixes group names starting with numbers with 'X'
      group_levels <- levels(group)
      cat("Group levels:", paste(group_levels, collapse = ", "), "\n")
      
      # Build contrast names handling the X prefix for numeric group names
      contrast_list <- list()
      comparisons_info <- list(
        list(name = "B1A_vs_Ctrl", test = "B1A", control = "CTRL"),
        list(name = "1B_vs_Ctrl", test = "ID1B", control = "CTRL"),
        list(name = "ID7_vs_Ctrl", test = "ID7", control = "CTRL")
      )
      
      # Create contrast strings
      contrast_strings <- c()
      comparison_names <- c()
      
      for (comp_info in comparisons_info) {
        test_name <- comp_info$test
        control_name <- comp_info$control
        
        # Check if group names need X prefix (if they start with a number)
        if (grepl("^[0-9]", test_name)) {
          test_name <- paste0("X", test_name)
        }
        if (grepl("^[0-9]", control_name)) {
          control_name <- paste0("X", control_name)
        }
        
        # Check if these groups exist in the design
        if (test_name %in% colnames(design) && control_name %in% colnames(design)) {
          contrast_str <- paste0(test_name, " - ", control_name)
          contrast_strings <- c(contrast_strings, contrast_str)
          comparison_names <- c(comparison_names, comp_info$name)
        }
      }
      
      if (length(contrast_strings) == 0) {
        cat("Warning: No valid contrasts found. Available columns:", 
            paste(colnames(design), collapse = ", "), "\n")
        next
      }
      
      # Make contrasts
      contrast.matrix <- makeContrasts(
        contrasts = contrast_strings,
        levels = design
      )
      colnames(contrast.matrix) <- comparison_names
      
      # Fit contrasts
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      
      # Process each comparison
      for (i in 1:length(comparison_names)) {
        comp_name <- comparison_names[i]
        exp_name <- paste0(comp_name, ".", folder)
        
        cat("\n--- Processing comparison:", comp_name, "---\n")
        
        # Get results
        res <- topTable(fit2, coef = i, number = Inf, sort.by = "none")
        res <- res %>%
          mutate(gene_id = row.names(res)) %>%
          rename(
            log2FoldChange = logFC,
            pvalue = P.Value,
            padj = adj.P.Val,
            stat = t
          ) %>%
          select(gene_id, log2FoldChange, AveExpr, stat, pvalue, padj, B) %>%
          arrange(pvalue)
        
        write_tsv(res, paste0(exp_name, "_limma_", feature_name, "s_results.tsv"))
        
        # Cutoffs
        p_adj_cutoff <- 0.05
        log2_cutoff <- 1
        
        # Count DE features
        n_de <- sum(res$padj < p_adj_cutoff & abs(res$log2FoldChange) > log2_cutoff, na.rm = TRUE)
        cat("Found", n_de, "DE", feature_name, "s\n")
        
        # Volcano plot
        pdf(file = paste0(exp_name, "_", feature_name, "s_volcano.pdf"), 
            width = pdfWidth, height = pdfHeight)
        print(
          res %>%
            filter(!is.na(padj)) %>%
            ggplot() +
            geom_point(aes(log2FoldChange, -log10(padj), 
                           colour = (abs(log2FoldChange) > log2_cutoff & padj < p_adj_cutoff))) +
            theme_bw() +
            scale_color_manual(values = c("black", "red")) +
            labs(title = paste0(comp_name, " - ", feature_name, "s")) +
            geom_text_repel(
              data = filter(res, padj < p_adj_cutoff & abs(log2FoldChange) > log2_cutoff),
              aes(log2FoldChange, -log10(padj), label = gene_id),
              min.segment.length = unit(0.1, "lines"),
              max.overlaps = 20
            )
        )
        dev.off()
        
        # Select DE genes
        DGEgenes <- res %>%
          filter(padj < p_adj_cutoff & abs(log2FoldChange) > log2_cutoff) %>%
          pull(gene_id)
        
        if (length(DGEgenes) == 0) {
          cat("No significant DE", feature_name, "s found\n")
          next
        }
        
        # Get normalized counts (logCPM) for ALL samples
        logCPM <- cpm(dge, log = TRUE)
        matDGEgenes <- logCPM[DGEgenes, , drop = FALSE]
        
        # Write DE genes
        res %>%
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

cat("\n=== Analysis complete ===\n