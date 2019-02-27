library(ggplot2)
library(scales)
library(reshape2)

len_dis1 <- read.table("/home/andrey/ablab/analysis/RNA_nanopores/assembly_august/ctg_lengths.illumina.txt", header=FALSE)
len_dis2 <- read.table("/home/andrey/ablab/analysis/RNA_nanopores/assembly_august/ctg_lengths.ont.txt", header=FALSE)
len_dis1 <- cbind(len_dis1, "SR-seq")
len_dis2 <- cbind(len_dis2, "Hybrid-seq")
colnames(len_dis1) <- c('Length', "Assembly")
colnames(len_dis2) <- c('Length', "Assembly")

dat <- rbind(len_dis1, len_dis2)

p <- ggplot()+
  geom_density(data = dat, aes(x=Length, fill=Assembly), alpha = 0.3) + 
  labs(title = "Assembled contigs distribution")+
  labs(y="Density")+
  labs(x="Assembled contigs length (bp)")  

png("contig_density_plot.png", width = 1900, height = 1200, res = 300)
plot(p) + scale_x_continuous(trans='log10') + annotation_logticks(sides = "b", base = 10, alpha = 0.75) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

len_dis1 <- read.table("/home/andrey/ablab/analysis/RNA_nanopores/assembly_august/tr_length.illumina.txt", header=FALSE)
len_dis2 <- read.table("/home/andrey/ablab/analysis/RNA_nanopores/assembly_august/tr_length.ont.txt", header=FALSE)
len_dis1 <- cbind(len_dis1, "SR-seq")
len_dis2 <- cbind(len_dis2, "Hybrid-seq")
colnames(len_dis1) <- c('Length', "Assembly")
colnames(len_dis2) <- c('Length', "Assembly")

dat <- rbind(len_dis1, len_dis2)

p <- ggplot()+
  geom_density(data = dat, aes(x=Length, fill=Assembly), alpha = 0.3) + 
  labs(title = "Annotated transcripts length distribution")+
  labs(y="Density")+
  labs(x="Annotated transcript length (bp)")  

png("annotation_density_plot.png", width = 1900, height = 1200, res = 300)
plot(p) + scale_x_continuous(trans='log10') + annotation_logticks(sides = "b", base = 10, alpha = 0.75) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()