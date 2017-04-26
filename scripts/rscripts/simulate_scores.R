##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Apr 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

library(caTools) # to use integration under curve (Trapezoid Rule Numerical Integration)
source("r_scripts/th17/used_for_paper/simulations_util.R")

############################################
# initialize
############################################
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
infl.nm <- paste(path.input.combined,"results_combine_data_exprsn",add.str,".Rdata",sep="")
# load input files
load(infl.nm)
keepers <- colnames(res[[1]][[1]][[1]])[1:16]
# remove chip combination types that are not used (only keep th17 chip scores)
for(i in 1:length(res)){ # go over activation,repression,or absolute, or whole
	tp <- names(res)[i]
	for(j in 1:length(res[[tp]])){ # for each tf
		tf <- names(res[[tp]])[j]
		res[[tp]][[tf]] <- res[[tp]][[tf]][[c.type.all.in.one.plot]]
	}
}

############################################
# create list to store simulated scores
############################################
# statistics we will follow
quants <- paste(sep="",as.character(c(0,1,5,25,50,75,95,99,100)),"%")
fdrs <- paste(sep="","FDR-",seq(0.5,1.9,by=0.1))
stats <- c("mean.r","sd.r","mean.p","sd.p")
c.nms <- c(fdrs,quants,stats)
# make place for all statistics we'll collect
rp <- list() # results permutations
ix.add <- 1:15 # only add for K,C,R,I or combination (not for sam, deseq, etc)
res.p <- list() # here we store the summed scores over all tfs
m <- matrix(0,nr=Np,nc=length(c.nms))
colnames(m) <- c.nms
ix.comb.conds <- colnames(res[[1]][[1]])[1:15]
for(i in 1:length(res)){ # go over activation,repression,or absolute, or whole
	tp <- names(res)[i]
	res.p[[tp]] <- list()
	for(j in 1:length(res[[tp]])){ # for each tf
		tf <- names(res[[tp]])[j]
		res.p[[tp]][[tf]] <- list()
		for(k in 1:length(ix.comb.conds)){
			cond <- ix.comb.conds[k]
			if(tp=="whole"){
				res.p[[tp]][[tf]][[cond]] <- list()
				res.p[[tp]][[tf]][[cond]]$pos <- m
				res.p[[tp]][[tf]][[cond]]$neg <- m
				res.p[[tp]][[tf]][[cond]]$all <- m
			} else {
				res.p[[tp]][[tf]][[cond]] <- m
			}
		}
	}
}
############################################
# get random performance based on permutations
############################################
m <- res[[tp]][[tf]] # a sample matrix of results
# this are the allowed genes for Chip (all genes)
ix.chip <- 1:nrow(m)
# this are the allowed genes for KO (genes with rpkm > 3 in either th17 or th0) this is more conservative than what we acutally did (or >3 in either ko or wt)
# ix.ko <- which( (abs(m[,"th17_rpkm"])>min.rpkm) | (abs(m[,"th0_rpkm"])>min.rpkm) ) 
ix.ko <- 1:nrow(m)
# this are the allowed genes for Immgen and RNAseq (genes with absolute sam th17 vs. th0 bigger than 2.5)
ix.immgen <- which(abs(m[,"SAM"])>z.abs.cut)
ix.rnaseq <- which(abs(m[,"SAM"])>z.abs.cut)
# only combine on 1:15, the rest are not network score (e.g. sam score, fpkms, etc.)
ix.keep <- colnames(res[[1]][[1]])[1:15]


# keep track of all kc whole scores
x.r.kc.whole <- numeric()
x.p.kc.whole <- numeric()
x.r.ri.whole <- numeric()
x.p.ri.whole <- numeric()
# for each random permutation
for(p in 1:Np){
	# create random scores for datasets and datasets combinations (in similar manner to what we did for real networks)
	for(j in 1:length(res[[tp]])){ # for each tf
		ix.sample.chip <- sample(ix.chip)
		ix.sample.ko <- sample(ix.ko)
		ix.sample.immgen <- sample(ix.immgen)
		ix.sample.rnaseq <- sample(ix.rnaseq)
		tf <- names(res[[tp]])[j] # tf we work on
		# random combo scores for activation
		tp <- "activation"
		m.r.act <- res[[tp]][[tf]] # real score combination matrix for tf
		m.p.act <- simulate.combine.data(m.r.act,ix.chip,ix.ko,ix.immgen,ix.rnaseq,ix.sample.chip,ix.sample.ko,ix.sample.immgen,ix.sample.rnaseq)	
		# random combo scores for repression
		tp <- "repression"
		m.r.rep <- res[[tp]][[tf]] # real score combination matrix for tf		
		m.p.rep <- simulate.combine.data(m.r.rep,ix.chip,ix.ko,ix.immgen,ix.rnaseq,ix.sample.chip,ix.sample.ko,ix.sample.immgen,ix.sample.rnaseq)		
		# random combo scores for absolute
		tp <- "absolute"
		m.r.abs <- res[[tp]][[tf]] # real score combination matrix for tf		
		m.p.abs <- simulate.combine.data(m.r.abs,ix.chip,ix.ko,ix.immgen,ix.rnaseq,ix.sample.chip,ix.sample.ko,ix.sample.immgen,ix.sample.rnaseq)		
		# random combo scores for whole (whole is generated from activation and repression just like in combine script)
		# first assign all values to be the values in activated
		tp <- "whole"
		m.p.whl <- m.p.act 
		m.r.whl <- m.r.act
		# now if repression scores are larger (in absolute sense) switch to repression scores
		for(l in 1:length(ix.keep)){
			combo <- ix.keep[l]			
			ix  <- which(m.p.act[,combo] < m.p.rep[,combo])
			m.p.whl[ix,combo] <- -1*m.p.rep[ix,combo]
			ix  <- which(m.r.act[,combo] < m.r.rep[,combo])
			m.r.whl[ix,combo] <- -1*m.r.rep[ix,combo]
		}
		############################################
		# keep for histogram later
		if(p==1){
			x.r.kc.whole <- cbind(x.r.kc.whole,m.r.whl[,"KC"])
			x.p.kc.whole <- cbind(x.p.kc.whole,m.p.whl[,"KC"])
			x.r.ri.whole <- cbind(x.r.ri.whole,m.r.whl[,"RI"])
			x.p.ri.whole <- cbind(x.p.ri.whole,m.p.whl[,"RI"])
		} else {
			x.p.kc.whole <- cbind(x.p.kc.whole,m.p.whl[,"KC"])
			x.p.ri.whole <- cbind(x.p.ri.whole,m.p.whl[,"RI"])
		}
		
		############################################
		# calc statistics for each of whole abs act and rep and store in res.p
		for(l in 1:length(ix.keep)){
			combo <- ix.keep[l]
			tp <- "activation"
			res.p[[tp]][[tf]][[combo]][p,] <- calc.statistics.rand.vs.true(x.r=m.r.whl[,combo],x.p=m.p.whl[,combo])
			tp <- "repression"
			res.p[[tp]][[tf]][[combo]][p,] <- calc.statistics.rand.vs.true(x.r=m.r.whl[,combo],x.p=m.p.whl[,combo])
			tp <- "absolute"
			res.p[[tp]][[tf]][[combo]][p,] <- calc.statistics.rand.vs.true(x.r=m.r.whl[,combo],x.p=m.p.whl[,combo])
			tp <- "whole"
			res.p[[tp]][[tf]][[combo]]$all[p,] <- calc.statistics.rand.vs.true(x.r=m.r.whl[,combo],x.p=m.p.whl[,combo])
			x.r <- m.r.whl[,combo]
			x.p <- m.p.whl[,combo]
			x.r[which(x.r<0)] <- 0
			x.p[which(x.p<0)] <- 0
			res.p[[tp]][[tf]][[combo]]$pos[p,] <- calc.statistics.rand.vs.true(x.r=x.r,x.p=x.p)
			x.r <- m.r.whl[,combo]
			x.p <- m.p.whl[,combo]
			x.r[which(x.r>0)] <- 0
			x.p[which(x.p>0)] <- 0
			res.p[[tp]][[tf]][[combo]]$neg[p,] <- calc.statistics.rand.vs.true(x.r=abs(x.r),x.p=abs(x.p))
		}
	}
}
	
f.nm <- paste(sep="",path.output,"fdrs.",add.str,"_Np_",Np,".RData")
save(res.p,file=f.nm)

f.nm <- paste(sep="",path.output,"kc.scores.real.vs.sim.",add.str,"_Np_",Np,".pdf")
pdf(f.nm)
x <- abs(as.numeric(x.r.kc.whole))
h.r <- hist(x[which(x>1)],plot = F)
x <- as.numeric(abs(x.p.kc.whole))
h.p <- hist(x[which(x>1)],plot = F)
h.p$counts <- h.p$counts/ncol(x.p.kc.whole)*length(tfs)
plot(h.r, col=rgb(0,0,1,1/4),main="KC scores vs. Simulated KC scores",cex.axis=2,xlab="KC score",ylab="# targets",cex.lab=1.5)  # first histogram
plot( h.p, col=rgb(1,0,0,1/4), add=T)  # second
legend("top",c("Real","Simulated"),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),pch=15,cex=2)
dev.off()
f.nm <- paste(sep="",path.output,"ri.scores.real.vs.sim.",add.str,"_Np_",Np,".pdf")
pdf(f.nm)
x <- abs(as.numeric(x.r.ri.whole))
h.r <- hist(x[which(x>1)],plot = F)
x <- as.numeric(abs(x.p.ri.whole))
h.p <- hist(x[which(x>1)],plot = F)
h.p$counts <- h.p$counts/ncol(x.p.ri.whole)*length(tfs)
plot(h.r, col=rgb(0,0,1,1/4),main="RI scores vs. Simulated RI scores",cex.axis=2,xlab="RI score",ylab="# targets",cex.lab=1.5)  # first histogram
plot( h.p, col=rgb(1,0,0,1/4), add=T)  # second
legend("top",c("Real","Simulated"),col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),pch=15,cex=2)
dev.off()
