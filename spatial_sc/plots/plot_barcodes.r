data = read.table("/home/andreyp/ablab/analysis/curio/plots/precision_recall.tsv", header=TRUE, sep="\t", row.names=1 )
data  <- t(as.matrix(data))
pdfWidth=11
pdfHeight=8;
pdf(file = "/home/andreyp/ablab/analysis/curio/plots/precision_recall.pdf", width = pdfWidth, height = pdfHeight) 
barplot(data, beside = TRUE, xlab = "Minimal score", width = 2, ylim = c(0, 100), col = c("darkgray", "darkred"), cex.axis=2, cex.names=2, cex.lab = 2, lwd=2)
dev.off()

pbp = read.table("/home/andreyp/ablab/analysis/curio/plots/per_barcode_precision.tsv", header=TRUE, sep="\t", row.names=1 )
pbp  <- t(as.matrix(pbp))
zpbp  <- t(as.matrix(pbp))
zpbp <- zpbp[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),]
pdf(file = "/home/andreyp/ablab/analysis/curio/plots/precision_per_bc_full.pdf", width = pdfWidth, height = pdfHeight) 

par(mar=c(5, 6, 4, 2) + 0.1, oma=c(0,1,0,0))
barplot(pbp, beside = TRUE, xlab = "Precision", ylab = "# Barcodes", width = 2, ylim = c(0, 80000), col = "darkred",  cex.axis=2, cex.names=2, cex.lab = 2, lwd=2)
par(fig = c(0.2,0.8,0.3,0.9), new = T)  
barplot(zpbp, beside = TRUE, width = 2, ylim = c(0, 40), col = "darkred", cex.axis=2, cex.names=2, cex.lab = 2, lwd=2)

dev.off()



pdfWidth=11
pdfHeight=8;
#outGraphFileNamePdf=paste(paste(outGraphFileName,lastFigure,sep="_"),"pdf",sep=".");
#pdf("Layers_Sig_Density_ISO.pdf")
cexMain=1.8;
cexLab=1.6;
legendscale=1.8
marVector=c(4, 4, 4, 2) +0.1;
mypadj=-5.3;
mycex=1.6;
cexText=1.7
myadj=2;
omaVector=c(1,1,1,3);
omaVector=c(0,0,0,0);
mysubfigureLabelCex=2.5;
mgpVector=c(2.4,0.8,0)
mylwd=2.5;
exampleY=0.495;
box1X=c(0.4,1.1);
box1Y=c(0.8,10);
box2X=c(8.25,9.5);
box2Y=c(0.8,10);
xMtextLine=2.3
bpw <- .8
colorU = "darkred"
colorD = "darkred"
L1med <- median(cortex_avg[[7]], na.rm = T)
L4med <- median(WM_avg[[7]], na.rm = T)
bp <- boxplot(ps, cex.axis=1.6, cexLab=1.6, names=c("Cortex","WM"),ylab="P56-P28 Mouse Exons",cex.lab=cexLab,cex.main=1.8,outline=FALSE, lwd=1.8, frame.plot = F)#, ylim=c(-.2,.6))
# drawing the colors
rect(1-bpw/2, bp$stats[2, 1], 1+bpw/2, bp$stats[4, 1],col=colorU);
lines(c(.6,1.4),c(L1med,L1med),col="black")

# drawing the colors
rect(2-bpw/2, bp$stats[2, 2], 2+bpw/2, bp$stats[4, 2],col=colorD);
lines(c(1.6,2.4),c(L4med, L4med),col="black")
#dev.off()

barplot()