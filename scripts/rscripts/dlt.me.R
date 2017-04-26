pdf(file=paste(sep="","~/Desktop/mtls_and_histone_marks_",date.is,".pdf"))
par(mfrow=c(2,2))
main="MTL classification (proximal or distal)"
ylab="#"
mp <- barplot(t(x),beside=T,axes = FALSE,axisnames = FALSE,main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
n=3
text(seq(2,n*length(rownames(x)),by=n), par("usr")[3], labels = rownames(x), srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.5)
legend("topright",colnames(x),fill=c("darkgray","white"),cex=.7)

main="MTL classification (normalized)"
ylab="proportion"
mp <- barplot(t(x/apply(x,1,sum)),beside=T,axes = FALSE,axisnames = FALSE,ylim=c(0,1),main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
n=3
text(seq(2,n*length(rownames(x)),by=n), par("usr")[3], labels = rownames(x), srt = 45, adj = c(1,2), xpd = TRUE,cex=.5,cex.main=.9)

main="tss.enrich.total"
ylab="#"
mp <- barplot(tss.enrich.total,beside=T,axes = FALSE,axisnames = FALSE,main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
n=3
text(mp, par("usr")[3], labels = mtl.types, srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.3)

# main="tss.enrich.proximal"
# ylab="#"
# mp <- barplot(tss.enrich.prox,beside=T,axes = FALSE,axisnames = FALSE,main=main,ylab=ylab,cex.main=.9)
# axis(2)
# axis(4)
# n=3
# text(mp, par("usr")[3], labels = mtl.types, srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.3)

l <- sapply(tss.density.list,length)
main="tss enrichment"
ylab="z-score"
boxplot(tss.density.list,axes = FALSE,outline=T,main=main,ylab=ylab,cex.main=.9,cex=.4)
axis(2)
axis(4)
text(1:length(sam.list), par("usr")[3], labels = paste(sep="",names(sam.list)," (",l,")"), srt = 45, adj = c(1,-1), xpd = TRUE,cex=.5)

main="Th17 vs. Th0 differential expression"
ylab="z-score"
boxplot(sam.list,axes = FALSE,outline=F,main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
text(1:length(sam.list), par("usr")[3], labels = names(sam.list), srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.3)

main="Th17 vs. Th0 differential expression\n(MTLs affecting more than 1000 genes)"
ylab="zscore"
l <- x[,"proximal"]
l2 <- sapply(sam.list,function(i) length(i))
ix <- which(l>1000)
boxplot(sam.list[ix],axes = FALSE,outline=F,main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
text(1:length(sam.list[ix]), par("usr")[3], labels = paste(sep="",names(sam.list)[ix]," (",l2[ix],")"), srt = 45, adj = c(1,-1), xpd = TRUE,cex=.5)

par(mfrow=c(2,1))
main="Th17 RPKM (no outliers)"
ylab="RPKM"
boxplot(rpkm.th17.list,axes = FALSE,outline=F,cex=.1,main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
text(1:length(rpkm.th17.list), par("usr")[3], labels = names(rpkm.th17.list), srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.5)

main="Th17 RPKM (with outliers)"
ylab="RPKM"
boxplot(rpkm.th17.list,axes = FALSE,outline=T,cex=.1,main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
text(1:length(rpkm.th17.list), par("usr")[3], labels = names(rpkm.th17.list), srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.5)

main="Th0 RPKM (no outliers)"
ylab="RPKM"
boxplot(rpkm.th0.list,axes = FALSE,outline=F,cex=.1,main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
text(1:length(rpkm.th0.list), par("usr")[3], labels = names(rpkm.th0.list), srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.5)

main="Th0 RPKM (with outliers)"
ylab="RPKM"
boxplot(rpkm.th0.list,axes = FALSE,outline=T,cex=.1,main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
text(1:length(rpkm.th0.list), par("usr")[3], labels = names(rpkm.th0.list), srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.5)

par(mfrow=c(1,1))
main="Th17 vs. Th0 Fold Change (no outliers)"
ylab="log2(FC)"
par(mfrow=c(1,1))
boxplot(fc.list,axes = FALSE,outline=F,cex=.1,main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
text(1:length(fc.list), par("usr")[3], labels = names(fc.list), srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.5)

main="Th17 vs. Th0 Fold Change (with outliers)"
ylab="log2(FC)"
par(mfrow=c(1,1))
boxplot(fc.list,axes = FALSE,outline=T,cex=.1,main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
text(1:length(fc.list), par("usr")[3], labels = names(fc.list), srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.5)
dev.off()
