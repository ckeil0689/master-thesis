##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Nov 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

library(gregmisc) # to use permutations function
# set paths
path.mtl.list <- "input/th17/used_for_paper/chipClusterAnalysis/"

# file names
f.nm.ll.th17 <- paste(sep="",path.mtl.list,"tf_clusters_th17_100_Dec_6_2011.RData")
f.nm.ll.th0 <- paste(sep="",path.mtl.list,"tf_clusters_th0_100_Dec_6_2011.RData")
f.nm.sam <- paste(sep="",path.mtl.list,"../samTh17VsTh0Zscores.xls")
f.nm.specificity <- paste(sep="",path.mtl.list,"../specificity_mat_Sep_21_2011.xls")
f.nm.rpkm <- paste(sep="",path.mtl.list,"../ranseqDatasetNoQuartNorm.RData")

# load mtls lists
load(f.nm.ll.th17)
load(f.nm.ll.th0)
load(f.nm.rpkm)
sam.scores <- read.delim(file=f.nm.sam)
spec.scores <- read.delim(file=f.nm.specificity)


core.tfs <- c("BATF","IRF4","MAF","STAT3","RORC")
# create a matrix each row is an MTL (of up to five core TFs) and each column is a TF
# each pixel can take a value of zero if TF is not part of MTL, or 1:5 depending of the location of the TF in the MTL 1=most 5' 5=most 3'
mtl.th17 <- pval.th17 <- matrix(0,nr=length(ll.th17),nc=5)
mtl.th0 <- pval.th0 <- matrix(0,nr=length(ll.th0),nc=5)
sam.th17 <- spec.th17 <- dtss.th17 <- numeric(length=length(ll.th17))
sam.th0 <- spec.th0 <- dtss.th0 <- numeric(length=length(ll.th0))
colnames(mtl.th17) <- colnames(mtl.th0) <- colnames(pval.th17) <- colnames(pval.th0) <- core.tfs


for(i in 1:length(ll.th17)){
	mtl.i <- ll.th17[[i]]
	tfs <- mtl.i$tf.to.s
	n.tfs <- length(tfs)
	mtl.th17[i,tfs] <- 1:length(tfs)	
	pval.th17[i,tfs] <- mtl.i$pval
	tmp <- as.numeric(unlist(strsplit(mtl.i$d.tss,"_")))
	dtss.th17[i] <- min(tmp[which(tmp!="")])
	trgt <- mtl.i$trgts[which(mtl.i$trgts!="")][1]
	if(!is.na(trgt)){
		sam.th17[i] <- sam.scores[trgt,"Score_d"]
		spec.th17[i] <- spec.scores[trgt,"min.specificity"]
	}
}

# for(i in 1:length(ll.th17)){
# # for(i in 1:100){
# 	mtl.i <- ll.th17[[i]]
# 	tfs <- mtl.i$tf.to.s
# 	n.tfs <- length(tfs)
# 	if(length(tfs)==1) {
# 		x <- 1
# 	} else if (length(tfs)==2) {
# 		x <- 1:2
# 	} else if (length(tfs)==3) {
# 		x <- c(1,2,1)
# 	} else if (length(tfs)==4) {
# 		x <- c(1,2,2,1)
# 	} else if (length(tfs)==5) {
# 		x <- c(1,2,3,2,1)
# 	}
# 	mtl.th17[i,tfs] <- x
# 	pval.th17[i,tfs] <- mtl.i$pval
# 	tmp <- as.numeric(unlist(strsplit(mtl.i$d.tss,"_")))
# 	dtss.th17[i] <- min(tmp[which(tmp!="")])
# 	trgt <- mtl.i$trgts[which(mtl.i$trgts!="")][1]
# 	if(!is.na(trgt)){
# 		sam.th17[i] <- sam.scores[trgt,"Score_d"]
# 		spec.th17[i] <- spec.scores[trgt,"min.specificity"]
# 	}
# 	if( (i%%100) == 0) {cat(".")}
# }

sam.th17[which(is.na(sam.th17))] <- 0
spec.th17[which(is.na(spec.th17))] <- 0

for(i in 1:length(ll.th0)){
	mtl.i <- ll.th0[[i]]
	tfs <- mtl.i$tf.to.s
	n.tfs <- length(tfs)
	mtl.th0[i,tfs] <- 1:length(tfs)	
	pval.th0[i,tfs] <- mtl.i$pval
}

mtl.length.th17 <- apply(mtl.th17,1,function(i) length(which(i!=0))) 
ix <- which(mtl.length.th17==5)
plot.l.th17 <- list()
for(i in 1:length(core.tfs)){
	tf <- core.tfs[i]
	plot.l.th17[[tf]] <- list()
	plot.l.th17[[tf]]$pos <- mtl.th17[ix,tf]
	plot.l.th17[[tf]]$pval <- pval.th17[ix,tf]
	plot.l.th17[[tf]]$mean.pval.per.pos <- numeric(length=length(core.tfs))
	plot.l.th17[[tf]]$median.pval.per.pos <- numeric(length=length(core.tfs))
	plot.l.th17[[tf]]$sd.pval.per.pos <- numeric(length=length(core.tfs))
	names(plot.l.th17[[tf]]$mean.pval.per.pos) <- names(plot.l.th17[[tf]]$sd.pval.per.pos) <- 1:length(core.tfs)
	for(j in 1:length(core.tfs)){
		ixx <- which(plot.l.th17[[tf]]$pos == j)
		plot.l.th17[[tf]]$mean.pval.per.pos[j] <- mean(plot.l.th17[[tf]]$pval[ixx])
		plot.l.th17[[tf]]$sd.pval.per.pos[j] <- sd(plot.l.th17[[tf]]$pval[ixx])
		plot.l.th17[[tf]]$median.pval.per.pos[j] <- median(plot.l.th17[[tf]]$pval[ixx])
	}
}

x <- lapply(plot.l.th17,function(i) table(i$pos))
# y <- lapply(plot.l.th17,function(i) mean(i$pval))
# y <- lapply(ll.th17[ix],function(i) i$d.tss)
# y2 <- lapply(y,function(i) unlist(strsplit(i,"_")))
# y3 <- unlist(y2)
# y4 <- as.numeric(y3[which(y3!="")])

pdf(file="~/Desktop/5_way_MTL_tf_positions_newd.pdf")
par(mfrow=c(2,3))
for(i in 1:length(core.tfs)){
	tf <- core.tfs[i]
	barplot(x[[tf]],xlab="position in MTL (1=most 5', ... , 5=most 3')",ylab="occurrence number",main=paste(sep="",tf," (5 way MTLs)"))
}
dev.off()

pdf(file="~/Desktop/5_way_MTL_tf_pval_per_position_newd.pdf")
par(mfrow=c(2,3))
for(i in 1:length(core.tfs)){
	tf <- core.tfs[i]
	names(plot.l.th17[[tf]]$mean.pval.per.pos) <- 1:length(core.tfs)
	barplot(plot.l.th17[[tf]]$median.pval.per.pos,xlab="position in MTL (1=most 5', ... , 5=most 3')",ylab="median pval",main=paste(sep="",tf," (5 way MTLs)"))
}
dev.off()

ix <- which(mtl.length.th17==5)
x <- apply(mtl.th17[ix,],1, function(i) paste(i,collapse="_"))
y <- sam.th17[ix]
z <- spec.th17[ix]
w <- sort(table(x),decreasing=T)
w.sam <- w.spec <- vector("list",length(w))
names(w.sam) <- names(w.spec) <- names(w)

for(i in 1:length(w)){
	type <- names(w)[i]
	ixx <- which(x==type)
	w.sam[[i]] <- y[ixx][which(y[ixx]!=0)]
	w.spec[[i]] <- z[ixx][which(z[ixx]!=0)]
}

boxplot(w.sam,axes = FALSE,xlab=paste(core.tfs,collapse="_"))
axis(2)
text(1:length(w.sam), par("usr")[3], labels = names(w.sam), srt = 45, adj = 1, xpd = TRUE,cex=.5)

boxplot(w.spec,axes = FALSE,xlab=paste(core.tfs,collapse="_"))
axis(2)
text(1:length(w.sam), par("usr")[3], labels = names(w.sam), srt = 45, adj = 1, xpd = TRUE,cex=.5)

#############
d <- rnaseq.complete.matrix
th17.wt.ix <- grep("th17_ss_48hr_.{1,6}_wt",colnames(d),value=T,perl=T)
th0.wt.ix <-  grep("th0_ss_48hr_.{1,6}_wt",colnames(d),value=T,perl=T)
th17.wt.ix <- c(th17.wt.ix,"SL1858_th17_ts_48hr#1")
th0.wt.ix <- c(th0.wt.ix,"SL2673_th0_ts_48hr#1")
rpkm.th17 <- apply(d[,th17.wt.ix],1,median)
rpkm.th0 <- apply(d[,th0.wt.ix],1,median)

clust.types.th17 <- sapply(ll.th17,function(i) paste( sort(unique(i[["tf.to.s"]])),collapse="_") )
c.table.th17 <- table(clust.types.th17)
mtl.types <- character()
for(i in 1:length(core.tfs)){
	mtl.types <- c(mtl.types,unique(apply(permutations(5,i,v=colnames(mtl.th17)),1,function(i) paste(sort(i),collapse="_"))))
}

# get ix for down-regulated vs. up-regulated genes
rpkm.cut <- GLOBAL[["median.abs.cut"]] 
# which genes are included? based on median rpkm > median.abs.cut
ix.med <- which(sam.scores[,"median_rpkm"] > rpkm.cut )
# which genes are included? based on zscores cutoff
z.cut <- GLOBAL[["z.abs.cut"]]
w.z <- sam.scores[,"Score_d"]
w.r <- sam.scores[,"median_rpkm"]
names(w.z) <- names(w.r) <- rownames(sam.scores)
ix.z.up <- which((w.z > z.cut) &  (w.r > rpkm.cut) )
ix.z.down <- which((w.z < -z.cut)  & (w.r > rpkm.cut) )
ix.z.same <- which( (w.z < z.cut) & (w.z > -z.cut) & (w.r > rpkm.cut) )

sam.list <- vector("list",length=length(c.table.th17))
names(sam.list) <- mtl.types
rpkm.th17.list <- rpkm.th0.list <- fc.list <- tss.density.list <- vector("list",length=length(c.table.th17))
names(rpkm.th17.list) <- names(rpkm.th0.list) <- names(fc.list)  <- names(tss.density.list) <- mtl.types
mtl.proximal.distal <- numeric(length=length(clust.types.th17)) # 1 is proximal,2 is distal
x <- matrix(0,nr=length(c.table.th17),nc=2,dimnames=list(mtl.types,c("proximal","distal")))
for(i in 1:length(mtl.types)){
	ix <- which(clust.types.th17==mtl.types[i])
	y <- lapply(ll.th17[ix],function(i) i$trgts)
	y2 <- lapply(y,function(i) unlist(strsplit(i,"_")))
	mtl.proximal.distal[ix] <- sapply(y2,function(i) if(length(i)>0){1}else{2})
	y3 <- sapply(y2,function(i) length(i))
	tss.density.list[[mtl.types[i]]] <- y3[which(y3>0)]
	w2=unique(unlist(y2[which(mtl.proximal.distal[ix]==1)]))
	ixx <- which(w2 %in% rownames(sam.scores))
	sam.list[[mtl.types[i]]] <- sam.scores[w2[ixx],"Score_d"]
	names(sam.list[[mtl.types[i]]]) <- rownames(sam.scores[w2[ixx],])
	ixx <- which(w2 %in% names(rpkm.th17))
	rpkm.th17.list[[mtl.types[i]]] <- rpkm.th17[w2[ixx]]
	ixx <- which(w2 %in% names(rpkm.th17))
	rpkm.th0.list[[mtl.types[i]]] <- rpkm.th0[w2[ixx]]
	fc.list[[mtl.types[i]]] <- log2((rpkm.th17.list[[mtl.types[i]]]+1)/(rpkm.th0.list[[mtl.types[i]]]+1))
	x[mtl.types[i],"proximal"] <- length(which(mtl.proximal.distal[ix]==1))
	x[mtl.types[i],"distal"] <- length(which(mtl.proximal.distal[ix]==2))
}
# xp <- matrix(0,nr=length(c.table.th17),nc=6,dimnames=list(mtl.types,
# 	c("proximal.up","proximal.sim","proximal.down","distal.up","distal.sim","distal.down")))

# calc features
mtl.num.prox <- x[,"proximal"]
mtl.num.total <- apply(x,1,sum)
mtl.num.gene.affected <- sapply(sam.list,function(i) length(i))


tss.enrich.prox <- tss.enrich.total <- matrix(0,nr=length(mtl.types),nc=3,dimnames=list(mtl.types,c("up","down","similar")))
for(i in 1:ncol(tss.enrich.prox)){
	# calc tss enrichment (#tss/#mtl_total)
	tss.enrich.total[,i] <- mtl.num.gene.affected/mtl.num.total
	# calc tss enrichment (#tss/#mtl_proximal)
	tss.enrich.prox[,i] <- mtl.num.gene.affected/mtl.num.prox
	
}

pdf("~/Desktop/mtls_and_histone_marks_newd.pdf")
par(mfrow=c(2,2))
main="MTL classification (proximal or distal)"
ylab="#"
mp <- barplot(t(x),beside=T,axes = FALSE,axisnames = FALSE,main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
n=3
text(seq(2,n*length(mtl.types),by=n), par("usr")[3], labels = mtl.types, srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.3)
legend("topright",colnames(x),fill=c("darkgray","white"),cex=.7)

main="MTL classification (normalized)"
ylab="proportion"
mp <- barplot(t(x/apply(x,1,sum)),beside=T,axes = FALSE,axisnames = FALSE,ylim=c(0,1),main=main,ylab=ylab,cex.main=.9)
axis(2)
axis(4)
n=3
text(seq(2,n*length(mtl.types),by=n), par("usr")[3], labels = mtl.types, srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.3,cex.main=.9)

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







