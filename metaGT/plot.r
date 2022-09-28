library(ggplot2)
library(fields)
library(RColorBrewer)
buylrd = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
           "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026") 
Lab.palette <- colorRampPalette(c(buylrd))

exp_table <- read.csv("fixed/counts/countsVSsalmon/RAW.tpm.values.tsv", sep="\t")
real<-exp_table[,2]
isoq<-exp_table[,3]

fudgeit <- function(){
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  fields::image.plot(xm,ym,z, axis.args=list(cex.axis=1.5), col = colramp(256), legend.only = T, add =F)
}

cor.test(isoq, real, method = "spearman")

par(mar = c(4, 5, 2, 6) + .1) 
smoothScatter(isoq,real,colramp = Lab.palette, 
              xlim = c(0,400),ylim = c(0,400), xaxs="i", yaxs="i",ylab="Salmon",
              xlab = "MetaGT",cex.lab=1.5,cex.axis=1.5, postPlotHook = fudgeit) #+abline(coef = c(0,1),lty="dashed")
text(122,380, "Spearman’s Rho: 0.981", col="yellow", cex=1.4, font=1)
text(90,350, "P-value < 2.2×10", col="yellow", cex=1.4, font=1)
text(182, 360, "-16", col="yellow", cex=0.8, font=1)

