
library("Matrix")

processTable = function(table, uncompressed_table, gene_id, sample_name_prefix) {
  
  if ( sum(table) < 50)  {
    return(0)
  }
  if ( nrow(table) < 2)  {
    return(0)
  }
  if ( ncol(table) < 2)  {
    return(0)
  }
  if ( nnzero(table) < 2)  {
    return(0)
  }
  if ( (min(colSums(table))*min(rowSums(table))/sum(table)) < 5)  {
    return(0)
  }
  write.table(table, file = paste(sample_name_prefix, "tables.tsv", sep =""), row.names = FALSE, col.names = FALSE, append = TRUE)
  write("=====", file = paste(sample_name_prefix, "tables.tsv", sep =""), append = TRUE)
  matrix <<- table
#  if ( all(rowSums(table == table[1,][col(table)]) == ncol(table)) )  {
#    return(0)
#  }

#  if ( length(unique(table[,2])) == 1 )  {
#    return(0)
#  }
#  if ( length(unique(table[,1])) == 1 )  {
#    return(0)
#  }
  matrices <<- unlist(list(matrices, list(table)), recursive=FALSE)
  uncompressed_matrices <<- unlist(list(uncompressed_matrices, list(uncompressed_table)), recursive=FALSE)
  gene_names <<- c(gene_names, gene_id)
  print(table)
  fisher <- chisq.test(table)
  print(fisher$p.value)
  p <- fisher$p.value
  if (p < 0.01) {
    write(sub("====", "", gene_name), file = paste(sample_name_prefix, "genes.tsv", sep =""), append = TRUE)
  }
  pvalues <<- c(pvalues, fisher$p.value)
}

preprocessTable = function(table) {
  if ( nrow(table) < 2)  {
    return(table)
  }
  if ( ncol(table) < 2)  {
    return(table)
  }
  print(table)
  while ( nrow(table) > 2) {
    smallest_index <- which.min(rowSums(table))
    smallest_row <- table[smallest_index,]
    table <- table[-smallest_index,]
    smallest_index2 <- which.min(rowSums(table))
    smallest_row2 <- table[smallest_index2,]
    combined_row <- smallest_row + smallest_row2
    table <- table[-smallest_index2,]
    table <- rbind(table, combined_row)    
  }
  while ( ncol(table) > 2) {
    smallest_index <- which.min(colSums(table))
    smallest_col <- table[,smallest_index]
    table <- table[,-smallest_index]
    smallest_index2 <- which.min(colSums(table))
    smallest_col2 <- table[,smallest_index2]
    combined_col <- smallest_col + smallest_col2
    table <- table[,-smallest_index2]
    table <- cbind(table, combined_col)    
  }
  
  return(table)
}

processFile = function(filepath, output_prefix) {
  con = file(filepath, "r")
  m = matrix(ncol = 2, nrow = 0)
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    if ( substring(line, 1, 1) == "=" || substring(line, 1, 1) == "E") {
      try(processTable(preprocessTable(m), m, gene_name, output_prefix))
      m = matrix(ncol = 2, nrow = 0)
      gene_name <<- line
      next
    }
    if (nrow(m) == 0) {
      m = matrix(ncol = length(unlist(strsplit(line, '\t'))), nrow = 0)
    }
    #print(m)
    m <- rbind(m, as.integer(unlist(strsplit(line, '\t'))))
  }
  
  close(con)
}

output_prefix = "/home/andrey/ablab/analysis/RNA_10x/codon_dependency/human/contigs_10x/15.11.19/human.10x_contigs.d2.ann_"
close( file( paste(output_prefix, "genes.tsv", sep =""), open="w" ) )
close( file( paste(output_prefix, "tables.tsv", sep =""), open="w" ) )
close( file( paste(output_prefix, "results.tsv", sep =""), open="w" ) )
pvalues <- c()
matrix <- matrix()
matrices <- list()
uncompressed_matrices <- list()
gene_names <- c()
gene_name <- ""

processFile("/home/andrey/ablab/analysis/RNA_10x/codon_dependency/human/contigs_10x/15.11.19/human.10x_contigs.d2.annotated_codons.codon_stats.raw.tsv", output_prefix)
sorted_pvalues_2 <- NULL
sorted_pvalues_2 <- sort(pvalues, index.return=TRUE, decreasing=FALSE)
#CTCF
a <- p.adjust(sorted_pvalues_2$x, method = "BH")
df <- data.frame(gene_names = gene_names[sorted_pvalues_2$ix])
df$pvalues <- a
write.table(df, file = paste(output_prefix, "results.tsv", sep =""), row.names=FALSE, col.names = FALSE)
