##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Apr 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# R code to create quality control plots for entire network pipeline
library(caTools) # to use integration under curve (Trapezoid Rule Numerical Integration)

path.input <- "input/th17/used_for_paper/DEseq/"
## path.output <- paste(sep="","results/validation_",date.is,"/")
path.output <- paste(sep="","results/validation/",date.is,"/")
system(paste(sep="","mkdir ", path.output))			
path.input.sam <- "input/th17/used_for_paper/"
path.input.combined <- paste(sep="","results/combinedAnalysis/",date.combine.data.run,"/")

# make output directory
if(filter.by.sam==TRUE){
   add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",z.abs.cut,"_deseq_cut_",deseq.pval.cut)
} else {
   add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",0,"_deseq_cut_",deseq.pval.cut)
}


# make output directory
# make output directory
#if(any("P300" %in% tfs)){
#	num.tfs.tmp <- num.tfs-1
#} else {
#	num.tfs.tmp <- num.tfs
#}
#if(filter.by.sam==TRUE){
#	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs.tmp,"_sam_",z.abs.cut)
#} else {
#	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs.tmp,"_sam_",0)
#}
#add.str.2 <- paste(sep="","_add_p300_",any("P300" %in% tfs))
infl.nm <- paste(path.input.combined,"/results_combine_data_exprsn",add.str,".Rdata",sep="")
# load input files
load(infl.nm)
# keepers <- colnames(res[[1]][[1]][[1]])[1:16]
# ## gold standard files
# th17.gns.f.nm <- paste(sep="",path.input,"gold_standard_th17_genes_",gold.stdrd.date,".txt")
# gs.gns <- read.delim(sep="\t",header=T,th17.gns.f.nm, as.is=T)
# gs.gns$Gene <- toupper(gs.gns$Gene)

# expt.invivo <- "IL17a.gfp.plusve.SI.vs.IL17a.gfp.minusve.SI"
# expt.invitro <- "Th17.vs.Th0"
expt.invivo <- "IL17a.gfp.plusve.SI.vs.naive"
expt.invitro <- "Th17.vs.naive"
# load fc data for all experiments
fc <- read.table(paste(sep="",path.input,date.deseq.run,"/DEseq_log2fc_",date.deseq.run,".xls"))
# load pval data for all experiments
pval <- read.table(paste(sep="",path.input,date.deseq.run,"/DEseq_pval_signed_",date.deseq.run,".xls"))



## get different categoires of the predicted genes:
# positive impact on th17
# ix.pos <- grep("positive",gs.gns$Effect,ignore.case=T)
# ix.neg <- grep("negative",gs.gns$Effect,ignore.case=T)
# ## get different categories of gene function (if has more than 5 genes in category)
# gn.type <- table(gs.gns$Gene.type)
# run.types <- list()
# run.types[["pos"]] <- gs.gns$Gene[ix.pos]
# run.types[["neg"]] <- gs.gns$Gene[ix.neg]
# for(j in 1:length(gn.type)){
# 	if(gn.type[j]>=5){
# 		ix <- which(gs.gns$Gene.type == names(gn.type)[j])
# 		run.types[[names(gn.type)[j]]] <- gs.gns$Gene[ix]
# 	}
# }

#################################
for(i in 1:length(res)){ # go over activation,repression,or absolute
	tp <- names(res)[i]
	for(j in 1:length(tfs)){ # go over tfs
		tf <- tfs[j]
		for(k in 1:length(res[[tp]][[tf]])){ # go over chip types
			c.tp <- names(res[[tp]][[tf]])[k]
			if(tp=="repression"){
				res[[tp]][[tf]][[c.tp]][,"SAM"] <- -1*res[[tp]][[tf]][[c.tp]][,"SAM"]
			} else if (tp=="absolute"){
				res[[tp]][[tf]][[c.tp]][,"SAM"] <- abs(res[[tp]][[tf]][[c.tp]][,"SAM"])
			}
		}		
	}
}

# combine scores for multiple tfs
ix.add <- 1:15 # only add for K,C,R,I or combination (not for sam, deseq, etc)
o.sum.list <- list() # here we store the summed scores over all tfs
for(i in 1:length(res)){ # go over activation,repression,or absolute
	tp <- names(res)[i]
	o.sum.list[[tp]] <- list()
	for(k in 1:length(res[[tp]][[1]])){ # go over chip types (this is same for any tf)
		c.tp <- names(res[[tp]][[1]])[k]
		# o.sum.rel.list <- list()
		for(j in 1:length(tfs)){ # go over tfs
			tf <- tfs[j]
			if(j==1){ # if first tf start the summation
				o.sum.list[[tp]][[c.tp]] <- res[[tp]][[tf]][[c.tp]]
			} else {
				o.sum.list[[tp]][[c.tp]][,ix.add] <- o.sum.list[[tp]][[c.tp]][,ix.add] + res[[tp]][[tf]][[c.tp]][,ix.add]
			}			
		}
	}
}

ix.fc.vivo <- rownames(fc)[which(abs(fc[,expt.invivo])>1)]
ix.pval.vivo <- rownames(fc)[which(abs(pval[,expt.invivo])>1)]
ix.vivo <- intersect(ix.fc.vivo,ix.pval.vivo)

ix.fc.invitro <- rownames(fc)[which(abs(fc[,expt.invitro])>1)]
ix.pval.invitro <- rownames(fc)[which(abs(pval[,expt.invitro])>-log10(0.1))]
ix.invitro <- intersect(ix.fc.invitro,ix.pval.invitro)

ix.fc.not.diff.invitro <- rownames(fc)[which(abs(fc[,expt.invitro])<0.5)]
ix.pval.not.diff.invitro <- rownames(fc)[which(abs(pval[,expt.invitro]) < (-log10(0.3)) )]
ix.vivo.not.invitro <- intersect(ix.fc.not.diff.invitro,ix.pval.not.diff.invitro)

ix.in.vivo.uniq <- setdiff(ix.vivo,ix.invitro)
ix.in.vivo.uniq.strict <- intersect(ix.in.vivo.uniq,ix.vivo.not.invitro)
ix.in.vivo.common <- intersect(ix.vivo,ix.invitro)

tp <- comb.case
c.tp <- c.type.all.in.one.plot 
ix <- sam.vs.tp

f.nm.res.pdf <- paste(sep="",path.output, "sumscore_vs_diffexp_",c.tp,"_",ix,"_",tp,"_",expt.invitro,add.str,".pdf")
pdf(f.nm.res.pdf )
# show that there is a correlation between sum scores and diff expression
######### plot for in vitro ############

keepers <- intersect(which(abs(fc[,expt.invitro])>1),which(abs(o.sum.list[[tp]][[c.tp]][,ix])>5))
x <- o.sum.list[[tp]][[c.tp]][,ix]
y <- fc[,expt.invitro]
names(y) <- names(x) <- rownames(o.sum.list[[tp]][[c.tp]])

lm.all <- lm(y~x)
lm.diff <- lm(y[keepers]~x[keepers])
r.sq.all <- summary(lm.all)$r.squared
r.sq.diff <- summary(lm.diff)$r.squared  

ylim <- c(-5,5)
main="In vitro"
ylab <- "Fold Change Th17 vs. Th0 in-vitro [log2]"
xlab <- paste(sep="",ix," combined scores")
labels.pch <- c("non diff. genes","Sig. diff. in-vitro")
labels.lty <- c(paste(sep="","linear regression ALL genes (R^2=",round(r.sq.all,2),")"),
				paste(sep="","linear regression DIFF genes (R^2=",round(r.sq.diff,2),")"))


plot(x,y,ylim=ylim,cex=.5,pch=20,col="gray",main=main,ylab=ylab,xlab=xlab)
keepers.up <- ix.invitro[which(y[ix.invitro]>0)]
keepers.down <- ix.invitro[which(y[ix.invitro]<0)]
keepers <- union(keepers.up,keepers.down)
points(x[keepers.up],y[keepers.up],col="blue",pch=20,cex=0.5)
points(x[keepers.down],y[keepers.down],col="blue",pch=20,cex=0.5)
abline(lm.all,col="black",lwd=2,lty=2)
abline(lm.diff,col="blue",lwd=2,lty=2)
legend("topleft",labels.pch,col=c("gray","blue"),pch=20,cex=0.5)
legend("bottomrigh",labels.lty,col=c("gray","blue"),lty=2,lwd=2,cex=0.5)
######### ######### ######### ######### 

######### plot for in vivo ############
keepers <- intersect(which(abs(fc[,expt.invivo])>1),which(abs(o.sum.list[[tp]][[c.tp]][,ix])>5))
x <- o.sum.list[[tp]][[c.tp]][,ix]
y <- fc[,expt.invivo]
names(y) <- names(x) <- rownames(o.sum.list[[tp]][[c.tp]])

lm.all <- lm(y~x)
lm.diff <- lm(y[keepers]~x[keepers])
r.sq.all <- summary(lm.all)$r.squared
r.sq.diff <- summary(lm.diff)$r.squared

ylim <- c(-5,5)
main="In vivo"
ylab <- "Fold Change Th17 vs. Th0 in-vivo [log2]"
xlab <- paste(sep="",ix," combined scores")
labels.pch <- c("non diff. genes","Sig. diff. in-vitro")
labels.lty <- c(paste(sep="","linear regression ALL genes (R^2=",round(r.sq.all,2),")"),
				paste(sep="","linear regression DIFF genes (R^2=",round(r.sq.diff,2),")"))

plot(x,y,ylim=ylim,cex=.5,pch=20,main=main,col="gray",ylab=ylab,xlab=xlab)
keepers.up <- ix.vivo[which(y[ix.vivo]>0)]
keepers.down <- ix.vivo[which(y[ix.vivo]<0)]
keepers <- union(keepers.up,keepers.down)
points(x[keepers.up],y[keepers.up],col="blue",pch=20,cex=0.5)
points(x[keepers.down],y[keepers.down],col="blue",pch=20,cex=0.5)
abline(lm.all,col="black",lwd=2,lty=2)
abline(lm.diff,col="blue",lwd=2,lty=2)
legend("topleft",labels.pch,col=c("gray","blue"),pch=20,cex=0.5)
legend("bottomrigh",labels.lty,col=c("gray","blue"),lty=2,lwd=2,cex=0.5)
######### ######### ######### ######### 
dev.off()

f.nm.res.pdf <- paste(sep="",path.output, "vldtn_invivo_",c.tp,"_",ix,"_",tp,add.str,".pdf")
pdf(f.nm.res.pdf)
	w <- sort(o.sum.list[[tp]][[c.tp]][,ix],decreasing=T)
	xlab <- paste(sep="","TF sum ",ix," scores [summed over: Batf,Irf4,Maf,Stat3,and Rorc]")
	main <- paste(sep="","TF sum ",ix," scores vs For Subsets of Differentially expressed genes")
	cls <- c("darkgrey","blue","red","darkgreen","darkgreen","darkgreen")
	plot(density(w),lty=1,lwd=2,col=cls[1],xlab=xlab,
			main=main,cex.main=0.9)
	lines(density(w[ix.invitro]),lty=1,lwd=2,col=cls[2])
	lines(density(w[ix.vivo]),lty=1,lwd=2,col=cls[3])
	labels=c(paste(sep="","All genes (n = ",length(w),")"),
			paste(sep="","Th17 vs Th0 invitro (n = ",length(ix.invitro),")"),
			paste(sep="","Th17 vs Th0 exvivo (n = ",length(ix.vivo),")"))
			legend("topright",labels,lty=2,cex=0.8,lwd=2,col=cls)

	cls <- c("darkgrey","red","red","red")
	plot(density(w),lty=1,lwd=2,col=cls[1],xlab=xlab,
			main=main,cex.main=0.9)
	lines(density(w[ix.in.vivo.common]),lty=2,lwd=2,col=cls[2])
	lines(density(w[ix.in.vivo.uniq]),lty=2,lwd=2,col=cls[3])
	lines(density(w[ix.in.vivo.uniq.strict]),lty=4,lwd=2,col=cls[4])

	# add numbering for different exvivo densities
	pch <- c("0","1","2","3")
	x <- density(w)
	points(x=x$x[which.max(x$y)],y=x$y[which.max(x$y)],lty=2,lwd=2,col="black",type="o",pch=pch[1],cex=1)
	x <- density(w[ix.in.vivo.common])
	points(x=x$x[which.max(x$y)],y=x$y[which.max(x$y)],lty=2,lwd=2,col="black",type="o",pch=pch[2],cex=1)
	x <- density(w[ix.in.vivo.uniq])
	points(x=x$x[which.max(x$y)],y=x$y[which.max(x$y)],lty=2,lwd=2,col="black",type="o",pch=pch[3],cex=1)
	x <- density(w[ix.in.vivo.uniq.strict])
	points(x=x$x[which.max(x$y)],y=x$y[which.max(x$y)],lty=2,lwd=2,col="black",type="o",pch=pch[4],cex=1)


	labels=c(paste(sep="","All genes (n = ",length(w),")"),
			paste(sep="","Th17 vs Th0 intersect {exvivo,invitro} (n = ",length(ix.in.vivo.common),")"),
			paste(sep="","Th17 vs Th0 intersect {exvivo,NOT invitro} (n = ",length(ix.in.vivo.uniq),")"),
			paste(sep="","Th17 vs Th0 intersect {exvivo,NOT invitro*} (n = ",length(ix.in.vivo.uniq.strict),")"))

	legend("topright",labels,cex=0.8,col="black",pch=pch)
dev.off()

## plot correlations between chip scores schemes (Th0, Th17, and Th17 minus th0) with SAM
# idea is we're looking to see which one is more correlated with expression changes
## plot chip score distributions after normalizing by ranks
f.nm.res.pdf <- paste(sep="",path.output, "vldtn_chip_vs_sam",c.tp,"_",ix,"_",tp,add.str,".pdf")
pdf(f.nm.res.pdf)

c.tp.all <- c("th0","th17","th17_minus_th0")
# tfs <- c(tfs,"P300")
ix <- "C"
sam.scores <- res[[tp]][[tf]][[c.tp]][,"SAM"]
# c.cor.with.sam <- list()
c.cor.with.sam <- matrix(0,nc=length(tfs),nr=length(c.tp.all))
rownames(c.cor.with.sam) <- c.tp.all
colnames(c.cor.with.sam) <- tfs
for(j in 1:length(tfs)){ # go over tfs
	tf <- tfs[j]
	for(k in 1:length(c.tp.all)){ # go over chip types
		c.tp <- c.tp.all[k]
		# get chip scores for c.tp k and tf j
		x <- res[[tp]][[tf]][[c.tp]][,ix]
		c.cor.with.sam[c.tp,tf] <- cor(x,sam.scores)
	}
}

barplot(c.cor.with.sam, main=paste(sep="","cor with differential expression"),ylab="correlation", col=rainbow(3),ylim=c(min(c.cor.with.sam),0.4),
  legend = rownames(c.cor.with.sam), beside=TRUE)

for(k in 1:length(c.tp.all)){ # go over chip types
	c.tp <- c.tp.all[k]
	m <- matrix(0,nc=length(tfs),nr=length(sam.scores))
	rownames(m) <- rownames(sam.scores)
	colnames(m) <- tfs
	for(j in 1:length(tfs)){ # go over tfs
		tf <- tfs[j]
		# get chip scores for c.tp k and tf j
		m[,tf] <- res[[tp]][[tf]][[c.tp]][,ix]
	}
	boxplot(m,main=paste(sep="",c.tp," chip scores"))
	
	plot(y=c(0,0,1,1),x=c(0,5,0,5),main=paste(sep="",c.tp," chip scores of genes with peaks only"),xlab="",ylab="normalized score",type="n",axes=FALSE)
	for(i in 1:dim(m)[2]){
		boxplot(m[which(m[,i]!=0),i],add=T,at=i-0.5)
	}
	axis(1,at=1:i-0.5,labels=tfs)		
	
}
dev.off()










