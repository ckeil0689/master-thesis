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

Np <- 1000 # number of permutations to estimate background distribution
N <- 100 # number of top genes to consider from gold standard (to see how enriched they are in predicted list)
N.core <- 7
core.tfs <- toupper(c("Irf4","Stat3","Rorc","Batf","Maf"))
date.combine.data.run <- "Sep_21_2011"
type <- "k_r_i_activator"
gs <- "sum.KC"
# gs <- "rorc.KC"
dataset <- "I"
# dataset <- "R"
# gs <- "SAMupreg"
# gs <- "SAMdownreg"
z.abs.cut <- 2.5
library(multicore)
# source required functions
source("r_scripts/th17/used_for_paper/util.R")
# setup paths
path.input <- "input/th17/used_for_paper/"
path.output <- paste(sep="","results/validation/")

# set full file names (inputs)
kc.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_KC_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
ri.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_RI_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
i.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_I_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
r.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_R_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
s.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_S_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
# sam.f.nm <- paste(sep="",path.input,"samTh17VsTh0Zscores.xls")

# set full file names (outputs)
tf.enrichment.pvals.f.nm <- paste(sep="",path.output,"tf_enrichment_pvals_Np=",Np,"_gs=",gs,"_pl=",dataset,"_type=",type,"_date=",date.combine.data.run,".xls")


# read files
kc <- as.matrix(read.table(kc.f.nm,header=T,sep="\t"))
s <- as.matrix(kc[,"sam.score"])
ri <- as.matrix(read.table(ri.f.nm,header=T,sep="\t"))
i <- as.matrix(read.table(i.f.nm,header=T,sep="\t"))
r <- as.matrix(read.table(r.f.nm,header=T,sep="\t"))

# remove unnecessary columns
col.rm <- c("sam.score","sum.score")
kc <- kc[,-which(colnames(kc) %in% col.rm)]
ri <- ri[,-which(colnames(ri) %in% col.rm)]
i <- i[,-which(colnames(i) %in% col.rm)]
r <- r[,-which(colnames(r) %in% col.rm)]

# s <- as.matrix(read.table(s.f.nm,header=T,sep="\t"))

# sam <- as.matrix(read.table(sam.f.nm,header=T,sep="\t"))

# remove genes that had a kc score of zero from comparisons
ix.rm <- which(apply(kc,1,sum)==0)
# stop("AM")
# ix.rm <- c(ix.rm,which(! colnames(r) %in% rownames(sam)))
kc <- kc[-ix.rm,]
ri <- ri[-ix.rm,]
i <- i[-ix.rm,]
r <- r[-ix.rm,]
s <- s[-ix.rm,]
# sam <- sam[rownames(r),]

if(gs == "sum.KC"){
	ref.list <- apply(kc[,core.tfs],1,sum)
} else if(gs == "SAMupreg"){
	x <- s+abs(min(ref.list))
	ref.list <- x
} else if (gs == "SAMdownreg") {
	tmp <- s*-1
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


# # make s all positive
# s.orig <- s
# s <- s+abs(min(s))

############### START: get permutation data for GSEA, AUROC, and AUPR ############
l1 <- numeric(nrow(i))
l1[1:N] <- 1 # this will be permutated in each run
# kc.sum <- apply(kc[,core.tfs],1,sum) 
l2=names(sort(ref.list,decreasing=T)) # this is the 'big' list we compare to
l2.scores=ref.list[l2] # for GSEA this is the vector of 'signal to noise' but for us it's confidence scores
gsea.prm <- mclapply(1:Np,GSEA.EnrichmentScore, gene.list=l2,gene.set=1:N,weighted.score.type = 1,correl.vector=l2.scores,mc.cores=N.core)
gS.for.pr.roc <- numeric(length(ref.list)) # this will be randomized to create the 'gold standard' for AUROC and AUPR
names(gS.for.pr.roc) <- names(ref.list)
gS.for.pr.roc[1:N] <- 1
pr.roc.prm <- mclapply(1:Np,evaluate.two.pred.lists.perm, testMat=as.matrix(ref.list), goldStdrdMat=gS.for.pr.roc,mc.cores=N.core)
ES.prm.vec <- sapply(gsea.prm,function(i) i$ES) # get random enrichment scores from gsea
aupr.prm.vec <- sapply(pr.roc.prm,function(i) i$AUPR) # get random enrichment scores from gsea
auroc.prm.vec <- sapply(pr.roc.prm,function(i) i$AUROC) # get random enrichment scores from gsea
############### DONE: get permutation data for GSEA, AUROC, and AUPR ############

############### START: for each tf calc enrichment score and pval based on GSEA, AUROC, and AUPR ############
# put data set you want to use for enrichment in M
gS <- ref.list
###
if(dataset=="I"){
	x <- as.matrix(read.table(i.f.nm,header=T,sep="\t"))
	x <- x[-ix.rm,]
	x[,"sam.score"] <- x[,"sam.score"] + abs(min(x[,"sam.score"]))
} else if(dataset=="R"){
	x <- as.matrix(read.table(r.f.nm,header=T,sep="\t"))
	x <- x[-ix.rm,]
	x[,"sam.score"] <- x[,"sam.score"] + abs(min(x[,"sam.score"]))
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

sam.score <-s[names(ES.pvals)]
kc.sum.score <- apply(kc[,core.tfs],1,sum)[names(ES.pvals)]

x <- cbind(aupr,auroc,ES,ES.pvals,aupr.pvals,auroc.pvals,geometric.mean.pval,sam.score,kc.sum.score)
cat(file=tf.enrichment.pvals.f.nm,sep="\t",c("gene_name",colnames(x)),"\n")
write.table(x,file=tf.enrichment.pvals.f.nm,sep="\t",col.names=F,append=T)


sort(geometric.mean.pval,decreasing=F)[1:50]
sort(aupr.pvals,decreasing=F)[1:30]
sort(auroc.pvals,decreasing=F)[1:30]
sort(ES.pvals,decreasing=F)[1:30]


stop("AM")

ix.gs <- sort(i.sum.f,decreasing=T,index.return=T)$ix[1:N]
goldStdrdMat <- i.sum.f
goldStdrdMat[-ix.gs]<- 0
goldStdrdMat[ix.gs]<- 1
goldStdrdMat.i <- as.matrix(goldStdrdMat)


testMat.i <- kc.sum[-ix.zeros.i]
gene.list=names(sort(testMat.i,decreasing=T))
correl.vector=testMat.i[gene.list]
o.p.i <- mclapply(1:Np,GSEA.EnrichmentScore, gene.list=gene.list,gene.set=1:100,weighted.score.type = 1,correl.vector=correl.vector,mc.cores=5)

# get sum for five TFs
kc.sum <- apply(kc[,core.tfs],1,sum)
ri.sum <- apply(ri[,core.tfs],1,sum)
i.sum <- apply(i[,core.tfs],1,sum)
r.sum <- apply(r[,core.tfs],1,sum)

