setwd("~/Desktop/ngs-analysis/rna_fischer/")
library("Matrix")
pvalues <- c()
matrix <- matrix()
matrices <- list()
uncompressed_matrices <- list()
gene_names <- c()

gene_name <- ""
processTable = function(table, uncompressed_table, gene_id) {
  
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
  write.table(table, file = "tables.csv", row.names = FALSE, col.names = FALSE, append = TRUE)
  write("=====", file = "tables.csv", append = TRUE)
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
    write(gene_name, file = "genes.csv", append = TRUE)
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

processFile = function(filepath) {
  con = file(filepath, "r")
  m = matrix(ncol = 2, nrow = 0)
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    if ( substring(line, 1, 1) == "=" ) {
      try(processTable(preprocessTable(m), m, gene_name))
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

processFile("reads.wd.codon_stats.tsv")
sorted_pvalues_2 <- sort(pvalues, index.return=TRUE, decreasing=FALSE)
#CTCF
a <- p.adjust(sorted_pvalues_2$x, method = "BH")
df <- data.frame(gene_names = gene_names[sorted_pvalues_2$ix])
df$pvalues <- a
write.table(df, file = "reads.csv",row.names=FALSE, col.names = FALSE)
