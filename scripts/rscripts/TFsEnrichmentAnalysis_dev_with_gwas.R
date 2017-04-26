##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# here we evaluate predictions made by R, I, or RI vs. KC (for the five core tfs) and 
# predict which TFs are also part of the core TFs using KC sum scores ranked list as gold standard

# Np <- 1000 # number of permutations to estimate background distribution
N.core <- 7
# tfs <- toupper(c("Irf4","Stat3","Rorc","Batf","Maf"))
# type <- "k_r_i_activator"
type <- type.2
# gs <- "sum.KC"
# gs <- "rorc.KC"
# dataset <- "I"
# # dataset <- "R"
# gs <- "SAMdownreg"
# gs <- "SAMdownreg"
# z.abs.cut <- 2.5
library(multicore)
# source required functions
source("r_scripts/th17/used_for_paper/util.R")

# annotate output directory
if(filter.by.sam==TRUE){
	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",z.abs.cut,"_deseq_cut_",deseq.pval.cut)
} else {
	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",0,"_deseq_cut_",deseq.pval.cut)
}

# setup paths
# path.input <- "input/th17/used_for_paper/"
path.input <- "results/combinedAnalysis/"
path.input.sam <- "input/th17/used_for_paper/"
path.output <- paste(sep="","results/enrichmentAnalysis/",date.is,"/")
system(paste(sep="","mkdir ", path.output))
path.output <- paste(sep="",path.output, type,add.str,"/")
system(paste(sep="","mkdir ", path.output))


# set full file names (inputs)
c.f.nm <- paste(sep="",path.input,date.combine,"/",type.1,"_",type.2,"_C",add.str,"_",date.combine,".xls")
k.f.nm <- paste(sep="",path.input,date.combine,"/",type.1,"_",type.2,"_K",add.str,"_",date.combine,".xls")
r.f.nm <- paste(sep="",path.input,date.combine,"/",type.1,"_",type.2,"_R",add.str,"_",date.combine,".xls")
i.f.nm <- paste(sep="",path.input,date.combine,"/",type.1,"_",type.2,"_I",add.str,"_",date.combine,".xls")
sam.f.nm <- paste(sep="",path.input.sam,"samTh17VsTh0Zscores.xls")
gwas.f.nm <- paste(sep="",path.input.sam,"GWAS_processed_final_Jun_22_2012.xls")

# set full file names (outputs)
tf.enrichment.pvals.f.nm <- paste(sep="",path.output,"tf_enrichment_pvals_Np=",Np,"_Ntrgts_",N,"_gs=",gs,"_pl=",dataset,"_type=",type,"_date=",date.combine,".xls")


# read files
k <- as.matrix(read.table(k.f.nm,header=T,sep="\t"))
c <- as.matrix(read.table(c.f.nm,header=T,sep="\t"))
i <- as.matrix(read.table(i.f.nm,header=T,sep="\t"))
r <- as.matrix(read.table(r.f.nm,header=T,sep="\t"))
kc <- k[,tfs]+c[,tfs]

ix <- which(! (colnames(r) %in% colnames(i))) # remove tfs not found in immgen data
if(length(ix)>0){
	r <- r[,-ix]
}
ri <-  i
ri[,colnames(r)] <- ri[,colnames(r)]+i[,colnames(r)]
# get sam scores
scores.sam <- read.table(sam.f.nm,sep="\t",colClasses = "character")
c.nms.sam <- as.character(scores.sam[1,1:(dim(scores.sam)[2]-1)])
scores.sam <- scores.sam[-1,]
r.nms.sam <- scores.sam[,1]
scores.sam <- scores.sam[,-1]
m <- matrix(0,nr=nrow(scores.sam),nc=ncol(scores.sam))
rownames(m) <- r.nms.sam
colnames(m) <- c.nms.sam
for(j in 1:ncol(m)){
	m[,j] <- as.numeric(scores.sam[,j])
}
s <- m
# deal with SAM data (have same genes as other datatypes: genes do not apear if not diff exprsd)
gns <- rownames(k)
m <- matrix(0,nr=nrow(k),nc=ncol(s))
rownames(m) <- gns
colnames(m) <- colnames(s)
ix <- which(gns %in% rownames(s))
m[gns[ix],] <- s[gns[ix],]
s <- m

##################################################
# get rpkm values
f.nm.rpkm <- paste(sep="",path.input.sam,"ranseqDatasetNoQuartNorm.RData") # which experiment to compare

## read rpkm data
load(f.nm.rpkm)
d <- rnaseq.complete.matrix
expt.names <-colnames(d)
th17.wt.ix <- grep("th17_ss_48hr_.{1,6}_wt",colnames(d),value=T,perl=T)
th17.wt.ix <- c(th17.wt.ix,"SL1858_th17_ts_48hr#1")
th17.mean <- apply(d[,th17.wt.ix],1,mean)
th0.wt.ix <- grep("th0_ss_48hr_.{1,6}_wt",colnames(d),value=T,perl=T)
th0.wt.ix <- c(th0.wt.ix,"SL2673_th0_ts_48hr#1")
th0.mean <- apply(d[,th0.wt.ix],1,mean)
rpkm.mean <- pmax(th17.mean,th0.mean)

# th17
m <- numeric(length=nrow(k))
names(m) <- gns
ix <- which(gns %in% names(th17.mean))
m[gns[ix]] <- th17.mean[gns[ix]]
th17.rpkm <- m
# th0
m <- numeric(length=nrow(k))
names(m) <- gns
ix <- which(gns %in% names(th0.mean))
m[gns[ix]] <- th0.mean[gns[ix]]
th0.rpkm <- m
##################################################
## get gwas file
gwas=as.matrix(read.table(gwas.f.nm,sep="\t",header=T,as.is=T))


kc.sum <- apply(kc,1,sum)
ix.sam <- names(which(abs(s[,1])>z.abs.cut))
ix.kc <- names(which(apply(kc,1,sum)!=0))
ix.keep <- intersect(ix.sam,ix.kc)
# ix.rm <- which(apply(kc,1,sum)==0)
# ix.rm <- c(ix.rm,which(! colnames(r) %in% rownames(sam)))
kc <- kc[ix.keep,]
ri <- ri[ix.keep,]
i <- i[ix.keep,]
r <- r[ix.keep,]
s <- s[ix.keep,]


if(gs == "sum.KC"){
	ref.list <- apply(kc[,tfs],1,sum)
} else if(gs == "sum.KC.up"){
	if(!(dataset=="gwas")){
		w <- apply(kc[,tfs],1,sum)
		x <- w+abs(min(w))
		ref.list <- x		
	} else {
		w <- kc.sum
		x <- w+abs(min(w))
		ref.list <- x
	}
} else if(gs == "sum.KC.down"){
	w <- apply(kc[,tfs],1,sum)
	w <- -1*w
	x <- w+abs(min(w))
	ref.list <- x
} else if(gs == "sum.KC.abs"){
	w <- apply(abs(kc[,tfs]),1,sum)
	ref.list <- w[which(w>0)]
	# ref.list <- -1*apply(kc[,tfs],1,sum)
} else if(gs == "SAMupreg"){
	x <- s[,1]+abs(min(s[,1]))
	ref.list <- x
} else if (gs == "SAMdownreg") {
	tmp <- -1*s[,1]
	x <- tmp+abs(min(tmp))
	ref.list <- x
} else if (gs == "batf.KC") {
	ref.list <- kc[,"BATF"]
} else if (gs == "irf4.KC") {
	ref.list <- kc[,"IRF4"]
} else if (gs == "irf4.KC") {
	ref.list <- kc[,"MAF"]
} else if (gs == "rorc.KC") {
	ref.list <- kc[,"RORC"]
}

############### START: get permutation data for GSEA, AUROC, and AUPR ############
l1 <- numeric(nrow(i))
l1[1:N] <- 1 # this will be permutated in each run
# kc.sum <- apply(kc[,tfs],1,sum) 
l2=names(sort(ref.list,decreasing=T)) # this is the 'big' list we compare to
l2.scores=ref.list[l2] # for GSEA this is the vector of 'signal to noise' but for us it's confidence scores
gsea.prm <- mclapply(1:Np,GSEA.EnrichmentScore, gene.list=l2,gene.set=1:N,weighted.score.type = 1,correl.vector=l2.scores,mc.cores=N.core)
gS.for.pr.roc <- numeric(length(ref.list)) # this will be randomized to create the 'gold standard' for AUROC and AUPR
names(gS.for.pr.roc) <- names(ref.list)
gS.for.pr.roc[1:N] <- 1
pr.roc.prm <- mclapply(1:Np,evaluate.two.pred.lists.perm, testMat=as.matrix(ref.list), goldStdrdMat=gS.for.pr.roc,mc.cores=N.core)
ES.prm.vec <- sapply(gsea.prm,function(i) i$ES) # get random enrichment scores from gsea
aupr.prm.vec <- sapply(pr.roc.prm,function(i) i$AUPR) # get random enrichment scores from aupr
auroc.prm.vec <- sapply(pr.roc.prm,function(i) i$AUROC) # get random enrichment scores from auroc
############### DONE: get permutation data for GSEA, AUROC, and AUPR ############

############### START: for each tf calc enrichment score and pval based on GSEA, AUROC, and AUPR ############
# put data set you want to use for enrichment in M
gS <- ref.list
###
if(dataset=="I"){
	x <- i
	# x <- as.matrix(read.table(i.f.nm,header=T,sep="\t"))
	# x <- x[-ix.rm,]
	# x[,"sam.score"] <- x[,"sam.score"] + abs(min(x[,"sam.score"]))
} else if(dataset=="R"){
	x <- r
	# x <- as.matrix(read.table(r.f.nm,header=T,sep="\t"))
	# x <- x[-ix.rm,]
	# x[,"sam.score"] <- x[,"sam.score"] + abs(min(x[,"sam.score"]))
} else if(dataset=="gwas"){
	x <- gwas
}
###
# fix potential problems where gn names from dataset are not found in gS
# stop("AM")
gns.data <- rownames(x)
ix.bad <- which(! (gns.data %in% names(gS)))
if(length(ix.bad)>0){
	x <- x[-ix.bad,]
}
if(dataset=="gwas"){
	gS <- gS[rownames(x)]
}
###

M <- x
res.x <- list()
for(k in 1:ncol(M)){
# for(k in 1:10){
	cat(".")
	tf <- colnames(M)[k]
	M.k <- as.matrix(M[,k]) # here are the predictions we made for tf k
	ix.keep <- sort(M.k,decreasing=T,index.return=T)$ix[1:N] # take first N predictions (the rest are zero to find AUPR and AUROC)
	ix.keep <- ix.keep[which(M.k[ix.keep]>0)]
	# goldStdrdMat.k <- M.k
	M.k[-ix.keep]<- 0
	M.k[ix.keep]<- 1
	res.x[[tf]]$tf <- tf
	# change kc.sum with something generic
	res.x[[tf]] <- evaluate.two.pred.lists(testMat=as.matrix(gS), goldStdrdMat=M.k)
	res.x[[tf]]$pval.aupr <- length(which(aupr.prm.vec>res.x[[tf]]$AUPR))/Np
	res.x[[tf]]$pval.auroc <- length(which(auroc.prm.vec>res.x[[tf]]$AUROC))/Np
	gene.list <- names(sort(gS,decreasing=T))
	gene.list.scores <- gS[gene.list]
	gene.set <- names(sort(M[,k],decreasing=T))[1:N]
	o <- GSEA.EnrichmentScore(gene.list=gene.list,gene.set=gene.set,weighted.score.type = 1,correl.vector=gene.list.scores)
	res.x[[tf]]$ES <- o$ES
	res.x[[tf]]$RES <- o$RES
	res.x[[tf]]$arg.ES <- o$arg.ES
	res.x[[tf]]$indicator <- o$indicator
	res.x[[tf]]$pval.ES <- length( which(ES.prm.vec > res.x[[tf]]$ES) )/Np
}
cat("\n")

############### END: for each tf calc enrichment score and pval based on GSEA, AUROC, and AUPR ############
aupr <- sapply(res.x,function(i) i$AUPR)
auroc <- sapply(res.x,function(i) i$AUROC)
ES <- sapply(res.x,function(i) i$ES)

aupr.pvals <- sapply(res.x,function(i) i$pval.aupr)
auroc.pvals <- sapply(res.x,function(i) i$pval.auroc)
ES.pvals <- sapply(res.x,function(i) i$pval.ES)

# geometric.mean.pval <- (aupr.pvals*auroc.pvals*ES.pvals)^(1/3)
# take geometric mean of three pvals (in log space to avoid numerical instability issues)
geometric.mean.pval <- 10^(1/3*( log10(aupr.pvals+1/Np) + log10(auroc.pvals+1/Np) + log10(ES.pvals+1/Np) ))
stop("AM")
if(dataset!="gwas"){
	sam.score <- 0 
	kc.sum <- kc.sum[names(ES.pvals)]
	th17.rpkm <- th17.rpkm[names(ES.pvals)]
	th0.rpkm <- th0.rpkm[names(ES.pvals)]
	x <- cbind(aupr,auroc,ES,ES.pvals,aupr.pvals,auroc.pvals,geometric.mean.pval,sam.score,kc.sum,th17.rpkm,th0.rpkm)
	sam.gns <- names(ES.pvals)[which( names(ES.pvals) %in% rownames(s))]
	x[sam.gns,"sam.score"] <- s[sam.gns,1]
	rpkm <- pmax(th17.rpkm,th0.rpkm)
	ix <- which(rpkm>3)
	cat(file=tf.enrichment.pvals.f.nm,sep="\t",c("gene_name",colnames(x)),"\n")
	write.table(x[ix,],file=tf.enrichment.pvals.f.nm,sep="\t",col.names=F,append=T)
} else {
	x <- cbind(aupr,auroc,ES,ES.pvals,aupr.pvals,auroc.pvals,geometric.mean.pval)
}

