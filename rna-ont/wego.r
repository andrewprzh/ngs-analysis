
library(ggplot2)
library(forcats)
library(dplyr)

WEGO <- read.table("./wego.tsv", header = FALSE, sep = "\t")
head(WEGO)
colnames(WEGO) <- c("ID", "Domain", "Number of genes", "Percentage of genes", "Ml", "Description")
head(WEGO)

WEGO2 <- WEGO %>%
 ## Group the entries by "Domain"
 group_by(Domain) %>%
 ## Take the top 5 entries per "Domain" according to "Percentage of genes"
 top_n(5, `Percentage of genes`) %>% 
 ## Ungroup the entries
 ungroup() %>% 
 ## Arrange the entries by "Domain", then by "Percentage of genes"
 arrange(Domain, `Percentage of genes`) %>% 
 ## Take note of the arrangement by creating a "Position" column
 mutate(Position = n():1)  

head(WEGO2)

png("WEGO.png", width = 2000, height = 1500, res = 300)

normalizer <- max(c(WEGO2$`Number of genes`))/max(c(WEGO2$`Percentage of genes`))
 ## Plot "Description" in the x-axis following the order stated in the "Position" column
## vs "Percentage of genes" in the first y-axis
ggplot(data = WEGO2, aes(x = fct_reorder(Description, dplyr::desc(Position)), y = `Percentage of genes`, fill = Domain)) +
 ## Plot "Description" in the x-axis following the order stated in the "Position" column
 ## vs normalized "Number of genes" in the second y-axis
 geom_col(data = WEGO2, aes(x = fct_reorder(Description, dplyr::desc(Position)), y = c(`Number of genes`)/normalizer)) +
 ## Add a second y-axis based on the transformation of "Percentage of genes" to "Number of genes".
 ## Notice that the transformation undoes the normalization for the earlier geom_col.
 scale_y_continuous(sec.axis = sec_axis(trans = ~.*normalizer, name = "Number of genes")) +
 ## Modify the aesthetic of the theme
 theme(axis.text.x = element_text(angle = 70, hjust = 1), axis.title.y = element_text(size = 8),
         legend.text = element_text(size = 7), legend.title = element_text(size = 8),
         legend.key.size =  unit(0.2, "in"), plot.title = element_text(size = 11, hjust = 0.5)) +
 ## Add a title to the plot
 labs(x = NULL, title = "Gene Ontology (GO) Annotation")
dev.off()



