run.it <- 7
if(run.it==1){
	# plot kc vs. sam scores
	date.combine.data.run <- "Apr_22_2012"
	path.input.combined <- paste(sep="","results/combinedAnalysis/",date.combine.data.run,"/")
	infl.nm <- paste(path.input.combined,"/results.combine.data.Rdata",sep="")
	# load input files
	load(infl.nm)
	# get deseq values for one TF (it has the same rownames for all tfs)
	gns.all <- names(res[[1]][[1]][[1]][,"DEseq"])
	cases <- c("activation")
	for(i in 1:length(cases)){
		cs <- cases[i]
		if (i==1) {
			# get node annotation data (stays the same for each case)
			m <- matrix(0,nr=nrow(res[[cs]][[tfs[1]]][[combine.type.1]]),nc=(length(tfs)+5))
			colnames(m) <- c(tfs,"SAM","th0_rpkm","th17_rpkm","highlight.gns","gn.function")
			rownames(m) <- rownames(res[[cs]][[tfs[1]]][[combine.type.1]])
			m[,c("SAM","th0_rpkm","th17_rpkm")] <- res[[cs]][[tfs[1]]][[combine.type.1]][,c("SAM","th0_rpkm","th17_rpkm")]
			m[which(rownames(m) %in% highlight.gns),"highlight.gns"] <- 1
			m[gns,"gn.function"] <- gn.function
			for(j in 1:length(tfs)){
				m[,tfs[j]] <- -1*res[[cs]][[tfs[j]]][[combine.type.1]][,"DEseq"]
			}
			gns <- gns[which(gns %in% rownames(m))]
			m.annot <- m[gns,]
		}

		m <- matrix(0,nr=length(gns),nc=length(tfs))
		colnames(m) <- core.2.tfs
		rownames(m) <- gns
		for(j in 1:length(tfs)){
			if(cs=="absolute"){
				m[,j] <- res[[cs]][[tfs[j]]][[combine.type.1]][gns,combine.type.2]*sign(res[[cs]][[tfs[j]]][[combine.type.1]][gns,"DEseq"])
			} else {
				m[,j] <- res[[cs]][[tfs[j]]][[combine.type.1]][gns,combine.type.2]
			}
		}
		# remove self regulation (we can't make it out with KO it will look like they all autoregulate)
		for(j in 1:length(tfs)){
			m[colnames(m)[j],j] <- 0
		}
		## define ouput file names
	}
	
} else if(run.it==2) {
load("/Users/aviv/Documents/nyu/littmanLab/th17_used_for_paper/input/th17/used_for_paper/ranseqDatasetNoQuartNorm.RData")
d=rnaseq.complete.matrix
ix.ts=c("SL1859_th0_exvivo_ss_#1",grep("ts",colnames(d),value=T)[1:8])
d=log2(d+1)
tfs <- c("BATF","IRF4","STAT3","RORC","MAF","NFE2L2")
max.y <- max(d[tfs,ix.ts])
min.y <- min(d[tfs,ix.ts])


pdf("ts_BIR.pdf")
plot(d["BATF",ix.ts],type="b",pch="B",ylim=c(min.y,max.y),xaxt="n",yaxt="n",lwd=5,font=2,col="orange",xlab="",ylab="",cex=2)
# ,ylab="mRNA level [RPKM]",xlab="time after stimulation [hr]",cex.lab=2)
axis(1,at=1:9,label=c(0,1,3,6,9,12,16,24,48),cex.axis=2)
# axis(2,at=seq(from=0,to=800,by=100))
axis(2,at=0:10,cex.axis=2)
lines(d["IRF4",ix.ts],type="b",pch="I",lwd=5,font=2,col="purple",cex=2)
lines(d["RORC",ix.ts],type="b",pch="R",lwd=5,font=2,col="darkgreen",cex=2)
legend("bottomright",c("Batf", "Irf4", "Rorc"), col = c("orange","purple","darkgreen"),fill=c("orange","purple","darkgreen"),cex=2)
dev.off()

lines(d["IL17A",ix.ts],type="b",pch=1,col="darkred")
lines(d["MAF",ix.ts],type="b",pch="M")
lines(c(0,d["IL17A",ix.ts[-1]]),type="b",pch=1,col="darkred")

lines(d["MAF",ix.ts],type="b",pch="M")
lines(d["STAT3",ix.ts],type="b",pch="S")

lines(d["NFE2L2",ix.ts],type="b",pch="N")
lines(d["STAT3",ix.ts],type="b",pch="S")

plot(d["BATF",ix.ts],type="b",pch="B",ylim=c(min.y,max.y),xaxt="n",yaxt="n",ylab="mRNA level [RPKM]",xlab="time after stimulation [hr]")
lines(d["NFE2L2",ix.ts],type="b",pch="N")
lines(d["HIF1A",ix.ts],type="b",pch="H")
lines(d["FOSL2",ix.ts],type="b",pch="F")
lines(d["KDM6B",ix.ts],type="b",pch="K")

lines(d["RORC",ix.ts],type="b",pch="R")
lines(d["MAF",ix.ts],type="b",pch="M")
lines(d["NFE2L2",ix.ts],type="b",pch="N")
} else if(run.it == 3){
	# w=read.table("input/th17/used_for_paper/DEseq_log2fc_May_21_2012.xls")
	w=read.table("input/th17/used_for_paper/DEseq_pval_signed_May_21_2012.xls")
	# plot pr comp for different lineages (DEseq log 2 FC)
	expt <- colnames(w)[1:4]
	sum.rows <- rowSums(w[,expt])
	ix.good <- which(sum.rows!=0)
	x <- w[ix.good,expt]
	expt.lables <- expt
	pr.lin=prcomp(t(w[,expt]))
	plot(pr.lin$x[,1],pr.lin$x[,2],cex=.2,pch=20,
	xlab=paste(sep="","PC1 (",round(pr.lin$sdev[1]/sum(pr.lin$sdev),2)*100,"%)"),
	ylab=paste(sep="","PC2 (",round(pr.lin$sdev[2]/sum(pr.lin$sdev),2)*100,"%)"),
	main="DEseq (log2FC)")
	text(pr.lin$x[,1],pr.lin$x[,2],labels=expt.lables)
	
	# plot pr comp for core tfs
	expt.ko <- c("Th17.batf.wt.vs.Th17.batf.ko",
	"Th17.maf.wt.vs.Th17.maf.ko",
	"Th17.irf4.wt.vs.Th17.irf4.ko",
	"Th17.stat3.wt.vs.Th17.stat3.ko",
	"Th17.rorc.wt.vs.Th17.rorc.ko",
	"Th17.nfe2l2.wt.vs.Th17.nfe2l2.ko",
	"Th17.fosl2.wt.vs.Th17.fosl2.ko",
	"Th17.hif1a.wt.vs.Th17.hif1a.ko",
	"Th17.ikzf3wt.vs.Th17.ikzf3ko")
	expt.kd <- c(
	"Th17.siNonTar.24hIMDM.vs.Th17.siRorc.24hIMDM",
	"Th17.siNonTar.24hIMDM.vs.Th17.kdm6b.24hIMDM",
	"Th17.siNonTar.24hRPMI.vs.Th17.siEtv6.24hRPMI",
	"Th17.siNonTar.24hRPMI.vs.Th17.sikdm6b.24hRPMI",
	"Th17.siNonTar.24hRPMI.vs.Th17.siJmjd6.24hRPMI",
	"Th17.siNonTar.24hRPMI.vs.Th17.siRorc.24hRPMI",
	"Th17.siNonTar.24hScreen.vs.Th17.siRorc.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siCcr6.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siId2.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siSki.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siLef1.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siAtf6.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siTrib3.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siSirt2.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siPrkrir.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siSatb1.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siEtv6.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siNfatc2.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siFosl2.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siCrem.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siAes.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siInhba.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siEgr2.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siBcl11b.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siJmjd3.24hScreen",
	"Th17.siNonTar.24hScreen.vs.Th17.siJmjd6.24hScreen"
	)
	expt <- c(expt.ko,expt.kd)
	# expt <- expt[1:9]
	expt <- expt.kd
	sum.rows <- rowSums(w[,expt])
	ix.good <- which(sum.rows!=0)
	x <- w[ix.good,expt]
	expt.lables.ko <- sapply(expt.ko,function(i) paste(sep="",strsplit(i,"\\.")[[1]][c(2)],collapse="_"))
	expt.lables.kd <- sapply(expt.kd,function(i) paste(sep="",strsplit(i,"\\.")[[1]][c(6)],collapse="_"))
	expt.labels <- expt.lables.kd
	# expt.labels <- expt.labels[1:9]
	expt.labels <- expt.labels[10:length(expt.labels)]
	pr.lin=prcomp(t(x),scale=T)
	pdf("~/Desktop/deseq_kd_pc_2d.pdf")
	plot(pr.lin$x[,1],pr.lin$x[,2],cex=.2,pch=20,
	xlab=paste(sep="","PC1 (",round(pr.lin$sdev[1]/sum(pr.lin$sdev),2)*100,"%)"),
	ylab=paste(sep="","PC2 (",round(pr.lin$sdev[2]/sum(pr.lin$sdev),2)*100,"%)"),main="DEseq (signed_pvals)")
	text(pr.lin$x[,1]-1,pr.lin$x[,2]-1,labels=expt.labels,cex=1)
	dev.off()
	
	# library(scatterplot3d)
	# library(rgl)
	# library(Rcmdr)
	scatter3d(pr.lin$x[,1],pr.lin$x[,2],pr.lin$x[,3],type="s",radius=1.5,
		xlab=paste(sep="","PC1 (",round(pr.lin$sdev[1]/sum(pr.lin$sdev),2)*100,"%)"),
		ylab=paste(sep="","PC2 (",round(pr.lin$sdev[2]/sum(pr.lin$sdev),2)*100,"%)"),
		zlab=paste(sep="","PC3 (",round(pr.lin$sdev[3]/sum(pr.lin$sdev),2)*100,"%)"),
		main="DEseq (log2FC)",col="orange",box=F, scale=F)
	# postscript("~/Desktop/deseq_kd_pc_2d.ps")
	# plot3d(pr.lin$x[,1],pr.lin$x[,2],pr.lin$x[,3],type="s",radius=1.5,
	# 	xlab=paste(sep="","PC1 (",round(pr.lin$sdev[1]/sum(pr.lin$sdev),2)*100,"%)"),
	# 	ylab=paste(sep="","PC2 (",round(pr.lin$sdev[2]/sum(pr.lin$sdev),2)*100,"%)"),
	# 	zlab=paste(sep="","PC3 (",round(pr.lin$sdev[3]/sum(pr.lin$sdev),2)*100,"%)"),
	# 	main="DEseq (log2FC)",col="orange",box=F, scale=F)
	# mv=5
	# text3d(pr.lin$x[,1]-mv,pr.lin$x[,2]-mv,pr.lin$x[,3]-mv,text=expt.labels)				
	# dev.off()
	load("input/th17/used_for_paper/ranseqDatasetNoQuartNorm.RData")
	# for core tfs
	w <- rnaseq.complete.matrix
	expt <- c("SL6496_th17_ss_48hr_ikzf3_ko#1",
	"SL3303_th17_ss_48hr_batf_ko#1","SL4224_th17_ss_48hr_batf_ko#2",
	"SL4220_th17_ss_48hr_cmaf_cre_ko#2","SL2765_th17_ss_48hr_cmaf_cre_ko#1",
	"SL2777_th17_ss_48hr_irf4_ko#1","SL3310_th17_ss_48hr_irf4_ko#2",
	"SL1846_th17_ss_48hr_rorc_ko#1","SL2682_th17_ss_48hr_rorc_ko#2",
	"SL1850_th17_ss_48hr_stat3_cre_ko#1","SL2675_th17_ss_48hr_stat3_cre_ko#2")
	sum.rows <- rowSums(w[,expt])
	ix.good <- which(sum.rows!=0)
	x <- w[ix.good,expt]
	expt.lables <- sapply(expt,function(i) paste(sep="",strsplit(i,"_")[[1]][c(1,5)],collapse="_"))
	pr.lin=prcomp(t(x),scale=T)
	plot(pr.lin$x[,1],pr.lin$x[,2],cex=.2,pch=20,
	xlab=paste(sep="","PC1 (",round(pr.lin$sdev[1]/sum(pr.lin$sdev),2)*100,"%)"),
	ylab=paste(sep="","PC2 (",round(pr.lin$sdev[2]/sum(pr.lin$sdev),2)*100,"%)"),main="RPKM")
	text(pr.lin$x[,1],pr.lin$x[,2],labels=expt.lables)
} else if (run.it==4){
	## create heatmaps for knock down experiments
	library("ggplot2")
	library(pheatmap)
	
	x <- unlist(strsplit(date()," +",perl=TRUE))
	date.is <- paste(x[2],x[3],x[5],sep="_")
	
	prcnt.chng <- read.table("input/th17/used_for_paper/DEseq/May_29_2012/DEseq_prcnt_chng_May_29_2012.xls")
	fc <- read.table("input/th17/used_for_paper/DEseq/May_29_2012/DEseq_log2fc_May_29_2012.xls")
	pval <- read.table("input/th17/used_for_paper/DEseq/May_29_2012/DEseq_pval_signed_May_29_2012.xls")
	conds.all <- grep("24hScreen",colnames(fc),value=T)
	conds.str <- c("Sirt2","Etv6","Nfatc2","Lef1","Ddit3","Crem","Aes","Prkrir","Atf6","Ski","Satb1","Bcl11b","Jmjd3","Rorc")

	
	conds.final <- character()
	for(i in 1:length(conds.str)){
		conds.final <- c(conds.final,grep(conds.str[i],conds.all,value=T,ignore.case=T))
	}
	ix.etv6 <- grep("etv6",conds.final,ignore.case=T)
	conds.final[ix.etv6] <- "Th17.siNonTar.24hRPMI.vs.Th17.siEtv6.24hRPMI"
	fc <- as.matrix(fc[,conds.final])
	pval <- as.matrix(pval[,conds.final])
	prcnt.chng <- as.matrix(prcnt.chng[,conds.final])
	colnames(fc) <- colnames(pval) <- colnames(prcnt.chng) <- sapply(sapply(colnames(fc),function(i) strsplit(i,"\\.")[[1]][6]), function(i) strsplit(i,"si")[[1]][2])
	x <- seq(0.0, 1.0, by=1/10) 
	orange.hsv = 0.075
	orange.2.black = hsv(orange.hsv,1.0,rev(x))

	blue.hsv = 0.65
	black.2.blue = hsv(blue.hsv,1.0,x)
	# col <- c(orange.2.black,black.2.blue)
	orange.2.black.staurated.top = c(orange.2.black,rep(orange.2.black[length(x)],length(x)/3))
	black.2.blue.staurated.top = c(black.2.blue,rep(black.2.blue[length(x)],length(x)/3))
	col <- c(orange.2.black.staurated.top,black.2.blue.staurated.top)

	# x <- read.table("~/Desktop/tmp/Signature_list_for_heatmap.txt",sep="\t")
	# gns <- unique(x[,"gene"])

	# m <- fc.gns
	# m.pval <- pval.gns
	m <- fc
	m.pval <- pval
	path.output <- "~/Desktop/tmp/"	
	# vec <- c(0.5,1,1.5,1.6,1.7,1.8,1.9,2)
	vec <- c(0.5,1,1.5,1.6,1.7,1.8,1.9,2)
	pval.cut <- 1
	gns.list <- list()
	gns.f.nm <- paste(sep="",path.output,"unsup_kds_gns_",date.is,".txt")
	for(i in 1:length(vec)) {
		fc.cut <- vec[i]
		fc.max.row <- apply(m,1,max)
		ix.1 <- which(abs(fc.max.row)>fc.cut)
		pval.max.row <- apply(m.pval,1,max)
		ix.2 <- which(abs(pval.max.row)>pval.cut)
		ix.final <- intersect(ix.1,ix.2)
		if(length(ix.final)>2){
			mdf <- as.data.frame(m[ix.final,])
			mdf$name <- as.factor(rownames(mdf))
			mgg <- melt(mdf)
			q <- ggplot(mgg, aes(variable, name)) + geom_tile(aes(fill = value)) + 
			scale_fill_gradient2(low = hsv(orange.hsv,1.0,1.0), mid="white",high = hsv(blue.hsv,1.0,1.0))
			q <- q + opts(axis.text.x=theme_text(angle=-90))
			f.nm <- paste(sep="",path.output,"unsup_kds_","fc_cut=",fc.cut,"_pval_cut=",pval.cut,".pdf")
			ggsave(f.nm,plot=q)
			# gns.list[[f.nm]] <- as.character(mdf$name)
			if(i==1){
				cat(file=gns.f.nm,sep="","> ",paste(sep="","fc_cut=",fc.cut,"_pval_cut=",pval.cut),"\n")
			} else {
				cat(file=gns.f.nm,sep="","> ",paste(sep="","fc_cut=",fc.cut,"_pval_cut=",pval.cut),"\n",append=T)
			}
			cat(file=gns.f.nm,sep="\n",as.character(mdf$name),append=T)
		}
	}


################ AM need to work this part out
# 	gns <- x.n[,1]
# 	fc.gns <- fc[gns,]
# 	pval.gns <- pval[gns,]
# 	m <- fc.gns
# 	m.pval <- pval.gns
# 	path.output <- "~/Desktop/tmp/"	
# 	# vec <- c(0.5,1,1.5,1.6,1.7,1.8,1.9,2)
# 	vec <- c(0.5,1,1.5,1.6,1.7,1.8,1.9,2)
# 	pval.cut <- 2
# 	for(i in 1:length(vec)) {
# 		fc.cut <- vec[i]
# 		fc.max.row <- apply(m,1,max)
# 		ix.1 <- which(abs(fc.max.row)>fc.cut)
# 		pval.max.row <- apply(m.pval,1,max)
# 		ix.2 <- which(abs(pval.max.row)>pval.cut)
# 		ix.final <- intersect(ix.1,ix.2)
# 		if(length(ix.final)>2){
# 			mdf <- as.data.frame(m[ix.final,])
# 			mdf$name <- as.factor(rownames(mdf))
# 			mgg <- melt(mdf)
# 			q <- ggplot(mgg, aes(variable, name)) + geom_tile(aes(fill = value)) + 
# 			scale_fill_gradient2(low = hsv(orange.hsv,1.0,1.0), mid="white",high = hsv(blue.hsv,1.0,1.0))
# 			q <- q + opts(axis.text.x=theme_text(angle=-90))
# 			ggsave(paste(sep="",path.output,"sup_kds_","fc_cut=",fc.cut,"_pval_cut=",pval.cut,".pdf"),plot=q)
# 		}
# 	}
# ################ AM need to work this part out
path.input <- "input/th17/used_for_paper/"
f.nm <- paste(sep="",path.input,"Signature_list_for_heatmap.txt")
x <- read.table(f.nm,sep="\t",header=T,as.is=T)

x.fc <- fc

# eff.max.val <- min( abs( c( max(x.fc), min(x.fc) ) ) )
eff.max.val <- 2
ix <- which(x.fc < -eff.max.val)
if(length(ix)){x.fc[ix] <- -eff.max.val}
ix <- which(x.fc > eff.max.val)
if(length(ix)){x.fc[ix] <- eff.max.val}

# set color scheme
n=10
tot <- n*4
or.col <- c(255,102,0)
bl.col <- c(0,51,255)
blue.col <- rgb(bl.col[1],bl.col[2],bl.col[3],max=255) 
orange.col <- rgb(or.col[1],or.col[2],or.col[3],max=255) 

###############
# j=4
# d <- dist(as.matrix(t(m[gns,])))
# hc <- hclust(d)
# col.order <- hc$labels[hc$order]
###############



fl.nm <- paste(sep="",path.output,"kd_heatmaps/selected_gns_",date.is,".pdf")
pdf(fl.nm)
m <- x.fc
for(j in 2:ncol(x)){
	nm <- colnames(x)[j]
	gns <- toupper(x$gene[which(x[,j]!=0)])
	d <- dist(as.matrix(m[gns,]))
	hc <- hclust(d)
	gn.order <- hc$labels[hc$order]

	ix <- which(x[,j]!=0)
	types <- sort(unique(x[ix,j]))
	types.all <- x[ix,j]
	if(max(types)>1){
		if(max(types)==4){
			labels <- c("Th1", "Th2","Th17","iTreg")
		} else if (max(types)==5) {
			labels <- c("Th1", "Th2","Th17","iTreg","Thelper")
		}
		ant <- data.frame(Type = factor(types.all, labels = labels))
		rownames(ant) = gns
		
		pheatmap(t(m[gns,]),scale="none",color=colorRampPalette(c(orange.col,"white",blue.col))(tot),
		breaks=seq(-2,2,by=0.1),cellwidth = 4, cellheight = 12,
		cex=0.35,cluster_cols=T,cluster_rows=F,annotation=ant,annotation_legend=F)
	} else {
		pheatmap(m[gn.order,],scale="none",color=colorRampPalette(c(orange.col,"white",blue.col))(tot),
		breaks=seq(-2,2,by=0.1),cellwidth = 12, cellheight = 4,treeheight_col=0,
		cex=0.35,cluster_cols=F,cluster_rows=F)			
	}
}
dev.off()


x <- read.table("~/Desktop/tmp/Signature_list_for_heatmap_aviv_2.txt",sep="\t",header=T,as.is=T)
x.n <- list()
for(i in 1:nrow(x)){
	x.n$gene[i] <- x$gene[i]
	if(x$Th1[i]==1) {
		x.n$pheno[i] <- "Th1"
	} else if(x$Th2[i]==1) {
		x.n$pheno[i] <- "Th2"
	} else if(x$Th17[i]==1) {
		x.n$pheno[i] <- "Th17"
	} else if(x$Tfh[i]==1) {
		x.n$pheno[i] <- "Tfh"
	} else if(x$iTreg[i]==1) {
		x.n$pheno[i] <- "iTreg"
	} else if(x$T_helper[i]==1) {
		x.n$pheno[i] <- "T_helper"
	} else if(x$T_cell_activator[i]==1) {
		x.n$pheno[i] <- "T_cell_activator"
	} else if(x$hotspot[i]==1) {
		x.n$pheno[i] <- "hotspot"
	}
}
x.n <- as.data.frame(x.n)
x.n$gene <- toupper(as.character(x.n$gene))


	gns <- x.n[,1]
	fc.gns <- fc[gns,]
	pval.gns <- pval[gns,]
	m <- fc.gns
	m.pval <- pval.gns
	for(i in 1:length(vec)) {
		fc.cut <- vec[i]
		fc.max.row <- apply(m,1,max)
		ix.1 <- which(abs(fc.max.row)>fc.cut)
		pval.max.row <- apply(m.pval,1,max)
		ix.2 <- which(abs(pval.max.row)>pval.cut)
		ix.final <- intersect(ix.1,ix.2)
		if(length(ix.final)>2){
			mdf <- as.data.frame(m[ix.final,])
			mdf <- mdf[gns[ix.final],]
			mdf[["gene.name"]] <- gns[ix.final]
			mdf[["pheno"]] <- x.n$pheno[ix.final]
			# mdf[["gene.name"]] <- as.character(rownames(mdf))
			mgg <- melt(mdf)
			mgg[["log2.FC"]] <- mgg$value
			mgg <- mgg[-which(names(mgg)=="value")]
			
		# 	base_size <- 9
		# 	ggplot(mgg, aes(variable, gene.name)) + geom_tile(aes(fill = log2.FC),colour = "white") + 
		# 	scale_fill_gradient2(low = hsv(orange.hsv,1.0,1.0), mid="white",high = hsv(blue.hsv,1.0,1.0))+
		# 	opts(axis.text.x=theme_text(angle=-90),panel.background = theme_blank(),
		# 	strip.text.y= theme_text(col="red"))
		# 	q + opts(axis.ticks = theme_blank())
		# 	
		# theme_blank(), axis.text.x = theme_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))
		# 	q + theme_grey(base_size = base_size) + labs(x = "",y = "") + scale_x_discrete(expand = c(0, 0)) +
		# 	     scale_y_discrete(expand = c(0, 0)) + opts(legend.position = "none", axis.ticks = 
		# 	theme_blank(), axis.text.x = theme_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))
			
			
			q <- ggplot(mgg, aes(variable, gene.name)) + geom_tile(aes(fill = log2.FC)) + 
			scale_fill_gradient2(low = hsv(orange.hsv,1.0,1.0), mid="white",high = hsv(blue.hsv,1.0,1.0))
			q <- q + opts(axis.text.x=theme_text(angle=-90))
			ggsave(paste(sep="",path.output,"sup_kds_","fc_cut=",fc.cut,"_pval_cut=",pval.cut,".pdf"),plot=q)
		}
	}		
		
	aviv.melt <- function(df,gns.list){
		o <- list()
		for(i in 1:length(gns.list)){
			e.nm <- names(gns.list)[i]
			o[["gene.name"]] <- df[!is.na(df[,e.nm])]
			o[["value"]] <- df[!is.na(df[,e.nm])]
			o[["variable"]] <- df[!is.na(df[,e.nm])]
			
		}
	}
		# m.pos <- m[ix.final,]
		# m.neg <- m[ix.final,]
		# m.pos[which(m.pos<=0)] <- 0 
		# m.neg[which(m.neg>=0)] <- 0 
		# if(length(ix.final)>=2){
		# 	pdf(paste(sep="",path.output,"kds_","fc_cut=",fc.cut,"_pval_cut=",pval.cut,".pdf"))
		# 		heatmap(m.pos,col=black.2.blue,cexRow=.5,scale="none",Rowv=NA,Colv=NA)
		# 		par(new=TRUE)
		# 		heatmap(m.neg,col=orange.2.black,cexRow=.5,scale="none",Rowv=NA,Colv=NA)
		# 	dev.off()
		# 	# pdf(paste(sep="",path.output,"kds_","fc_cut=",fc.cut,"_pval_cut=",pval.cut,".pdf"))
		# 	heatmap(m[rev(ix.final),],col=col,cexRow=.5,scale="none",Rowv=NA,Colv=NA)
		# 	heatmap(m[rev(ix.final),],col=col,cexRow=.5,scale="none",Rowv=NA,Colv=NA)
		# 	# dev.off()
		# }
	# }
} else if (run.it==6){
	x=as.matrix(read.table("~/Downloads/GWAS_Jun_22_2012_processed.txt",sep="\t",header=T,as.is=T))
	gns <- x[,"gene"]
	x <- as.matrix(x[,-which(colnames(x)=="gene")])
	m <- matrix(0,nr=nrow(x),nc=ncol(x),dimnames=dimnames(x))
	for(j in 1:ncol(x)){
		m[,j] <- as.numeric(x[,j])
	}
	tbl.gns <- sort(table(gns),decreasing=T)
	gns.non.uniq <- names(which(tbl.gns>1))
	gns.uniq <- names(which(tbl.gns==1))
	m.u <- matrix(0,nr=length(tbl.gns),nc=ncol(x))
	colnames(m.u) <- colnames(m)
	rownames(m.u) <- names(tbl.gns)
	for(i in 1:length(gns.uniq)){
		gn <- gns.uniq[i]
		ix <- which(gns==gn)
		m.u[gn,] <- m[ix,]
	}

	
	for(i in 1:length(gns.non.uniq)){
		gn <- gns.non.uniq[i]
		ix <- which(gns==gn)
		tmp <- numeric()
		for(j in 1:length(ix)){
			if(m[ix[j],"distance"]<200){
				tmp <- rbind(tmp,m[ix[j],])
			}	
		}
		if(length(tmp)>1){
			m.u[gn,] <- apply(tmp,2,max)
			m.u[gn,"distance"] <- min(tmp[,"distance"])
		}
	}
	rownames(m.u) <- toupper(rownames(m.u))
	write.table(m.u,"~/Downloads/GWAS_Jun_22_2012_processed_final.txt",sep="\t")
} else if (run.it==7){
	path.input <- "results/validation/Jun_23_2012/"
	fl.nm <- paste(sep="",path.input,"vldtn_perTF_AllInOne2_activation_th17_cut_prcnt_0_num_tfs_5_sam_0_deseq_cut_1_gs_data_COMBINED_Jun_22_2012.xls")
	x <- read.delim(file=fl.nm,sep="\t",header=T,as.is=T)
	r.names <- x[,"comb.type"]
	ix.rm <- which( (colnames(x)=="comb.type") | (colnames(x)=="X"))
	combo <- as.matrix(x[,-ix.rm])
	rownames(combo) <- r.names

	fl.nm <- paste(sep="",path.input,"vldtn_perTF_AllInOne2_activation_th17_cut_prcnt_0_num_tfs_5_sam_0_deseq_cut_1_gs_data_SINGLE_KC_Jun_22_2012.xls")
	x <- read.delim(file=fl.nm,sep="\t",header=T,as.is=T)
	r.names <- x[,"comb.type"]
	ix.rm <- which( (colnames(x)=="comb.type") | (colnames(x)=="X"))
	single <- as.matrix(x[,-ix.rm])
	rownames(single) <- r.names
	
	
	x <- combo["KC",grep("aucroc",colnames(combo),ignore.case=T)]
	names(x) <- sapply(names(x), function(i) strsplit(i,"_")[[1]][1])
	y <- single[grep("aucroc",rownames(single),ignore.case=T),]
	rownames(y) <- sapply(rownames(y), function(i) strsplit(i,"_")[[1]][2])
	
	fl.nm <- paste(sep="",path.input,"gwas_with_th17_cut_prcnt_0_num_tfs_5_sam_0_deseq_cut_1_Jun_22_2012.pdf")
	pdf(fl.nm)
	m <- nrow(y)
	ylim <- c(0,max(c(x,y))+0.1)
	col <- rainbow(6)[-2]
	barplot(x,axes=F,col="white",ylim=ylim,names.arg = NA,xlim=c(0,16),space=0.2)
	axis(2)
	for(j in 1:length(x)){
		w <- y[,j]
		if(j==1){
			s <- 0.7
			points(rep(s,length(w)),w,cex=0.5,pch=3,col=col)
			text(s, par("usr")[3], labels = names(x)[j], srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.8)
		}else {
			pl <- s+(j-1)*1.2
			points(rep(pl,length(w)),w,cex=0.5,pch=3,col=col)
			text(pl, par("usr")[3], labels = names(x)[j], srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.8)
		}
		
	}
	legend("topleft",legend=rownames(y),pch=3,col=col)
	
	x <- combo["KC",grep("aucpr",colnames(combo),ignore.case=T)]
	names(x) <- sapply(names(x), function(i) strsplit(i,"_")[[1]][1])
	y <- single[grep("aucpr",rownames(single),ignore.case=T),]
	rownames(y) <- sapply(rownames(y), function(i) strsplit(i,"_")[[1]][2])
	
	m <- nrow(y)
	ylim <- c(0,0.06)
	barplot(x,axes=F,col="white",ylim=ylim,names.arg = NA,xlim=c(0,16),space=0.2)
	axis(2)
	for(j in 1:length(x)){
		w <- y[,j]
		if(j==1){
			s <- 0.7
			points(rep(s,length(w)),w,cex=0.5,pch=3,col=col)
			text(s, par("usr")[3], labels = names(x)[j], srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.8)
		}else {
			pl <- s+(j-1)*1.2
			points(rep(pl,length(w)),w,cex=0.5,pch=3,col=col)
			text(pl, par("usr")[3], labels = names(x)[j], srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.8)
		}
	}
	legend("topleft",legend=rownames(y),pch=3,col=col)
	dev.off()
	
	
	x <- combo["KC",grep("aucroc",colnames(combo),ignore.case=T)]
	names(x) <- sapply(names(x), function(i) strsplit(i,"_")[[1]][1])
	y <- single[grep("aucroc",rownames(single),ignore.case=T),]
	rownames(y) <- sapply(rownames(y), function(i) strsplit(i,"_")[[1]][2])
	x <- x[-length(x)]
	y <- y[,-ncol(y)]
	
	fl.nm <- paste(sep="",path.input,"gwas_no_th17_cut_prcnt_0_num_tfs_5_sam_0_deseq_cut_1_Jun_22_2012.pdf")
	pdf(fl.nm)
	m <- nrow(y)
	ylim <- c(0,max(c(x,y))+0.1)
	col <- rainbow(6)[-2]
	barplot(x,axes=F,col="white",ylim=ylim,names.arg = NA,xlim=c(0,16),space=0.2)
	axis(2)
	for(j in 1:length(x)){
		w <- y[,j]
		if(j==1){
			s <- 0.7
			points(rep(s,length(w)),w,cex=0.5,pch=3,col=col)
			text(s, par("usr")[3], labels = names(x)[j], srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.8)
		}else {
			pl <- s+(j-1)*1.2
			points(rep(pl,length(w)),w,cex=0.5,pch=3,col=col)
			text(pl, par("usr")[3], labels = names(x)[j], srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.8)
		}
		
	}
	legend("topleft",legend=rownames(y),pch=3,col=col)
	
	x <- combo["KC",grep("aucpr",colnames(combo),ignore.case=T)]
	names(x) <- sapply(names(x), function(i) strsplit(i,"_")[[1]][1])
	y <- single[grep("aucpr",rownames(single),ignore.case=T),]
	rownames(y) <- sapply(rownames(y), function(i) strsplit(i,"_")[[1]][2])
	x <- x[-length(x)]
	y <- y[,-ncol(y)]
		
	m <- nrow(y)
	ylim <- c(0,0.06)
	barplot(x,axes=F,col="white",ylim=ylim,names.arg = NA,xlim=c(0,16),space=0.2)
	axis(2)
	for(j in 1:length(x)){
		w <- y[,j]
		if(j==1){
			s <- 0.7
			points(rep(s,length(w)),w,cex=0.5,pch=3,col=col)
			text(s, par("usr")[3], labels = names(x)[j], srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.8)
		}else {
			pl <- s+(j-1)*1.2
			points(rep(pl,length(w)),w,cex=0.5,pch=3,col=col)
			text(pl, par("usr")[3], labels = names(x)[j], srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.8)
		}
	}
	legend("topleft",legend=rownames(y),pch=3,col=col)
	dev.off()
}



































