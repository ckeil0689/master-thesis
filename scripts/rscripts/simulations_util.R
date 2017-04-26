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

############################################
# helper functions
############################################
# this function takes a real combine data matrix and returns a randomly combined data list
# the way it works: core.scores (K,C,R,I): are shuffled (the rest of score combinations are based on the same suffle)
# input:
# 		m.r is the real combine data matrix
#		ix.chip: which genes are allowed to get a chip scores (in case there is some restriction)
#		ix.ko: which genes are allowed to get a chip scores (in case there is some restriction)
#		ix.immgen: which genes are allowed to get a chip scores (in case there is some restriction)
#		ix.rnaseq: which genes are allowed to get a chip scores (in case there is some restriction)
#		ix.sample.chip: a permutation of ix.chip
#		ix.sample.ko: a permutation of ix.ko
#		ix.sample.immgen: a permutation of ix.immgen
#		ix.sample.rnaseq: a permutation of ix.rnaseq
simulate.combine.data <- function(m.r,ix.chip,ix.ko,ix.immgen,ix.rnaseq,ix.sample.chip,ix.sample.ko,ix.sample.immgen,ix.sample.rnaseq){
	# a matrix where we will store random results
	m.p <- m.r
	m.p[which(m.p[,1:15]!=0)] <- 0 
	# permute core scores
	m.p[ix.chip,"C"] <- m.r[ix.sample.chip,"C"]
	m.p[ix.ko,"K"] <- m.r[ix.sample.ko,"K"]
	m.p[ix.immgen,"I"] <- m.r[ix.sample.immgen,"I"]
	m.p[ix.rnaseq,"R"] <- m.r[ix.sample.rnaseq,"R"]
	# generate all combo scores from core scores
	m.p[,"KR"] <- m.p[,"K"]+m.p[,"R"]
	m.p[,"KI"] <- m.p[,"K"]+m.p[,"I"]
	m.p[,"RI"] <- m.p[,"R"]+m.p[,"I"]
	m.p[,"KC"] <- m.p[,"K"]+m.p[,"C"]			
	m.p[,"CR"] <- m.p[,"C"]+m.p[,"R"]	
	m.p[,"CI"] <- m.p[,"C"]+m.p[,"I"]	
	m.p[,"KRI"] <- m.p[,"K"]+m.p[,"R"]+m.p[,"I"]
	m.p[,"KCR"] <- m.p[,"K"]+m.p[,"C"]+m.p[,"R"]	
	m.p[,"KCI"] <- m.p[,"K"]+m.p[,"C"]+m.p[,"I"]		
	m.p[,"CRI"] <- m.p[,"C"]+m.p[,"R"]+m.p[,"I"]
	m.p[,"KCRI"] <- m.p[,"K"]+m.p[,"C"]+m.p[,"R"]+m.p[,"I"]
	return(m.p)
}
# this function takes a real combine data matrix and returns a randomly combined data list
# the way it works: core.scores (K,C,R,I): are shuffled (the rest of score combinations are based on the same suffle)
# input:
# 		x.r is the real data score vector
# 		x.p is the permuted data score vector
# output:
#		a vector with various statistics such as fdrs (at various cutoffs), mean permuted/real score, and sd permuted/real scores
calc.statistics.rand.vs.true <- function(x.r,x.p){
	quants <- c(0,1,5,25,50,75,95,99,100)
	names(quants) <- as.character(quants)
	sum.scores.to.test <- seq(0.5,1.9,by=0.1)
	fdr.nms <- paste(sep="","FDR-",sum.scores.to.test)
	fdrs <- numeric(length(fdr.nms))
	names(fdrs) <- fdr.nms
	stats <- numeric(4)
	names(stats) <- c("mean.r","sd.r","mean.p","sd.p")
	# get quants
	quants <- quantile(x.p,probs=quants/100)
	# get stats
	stats <- c(mean(x.r),sd(x.r),mean(x.p),sd(x.p))	
	# get fdrs
	for(k in 1:length(fdrs)){
		score.k <- sum.scores.to.test[k]
		n.p <- length(which(abs(x.p)>score.k)) # how many expected by chance
		n.r <- length(which(abs(x.r)>score.k)) # how many we seen in real data
		if(n.r==0){
			fdrs[k] <- 0
		} else {
			fdrs[k] <- (n.p)/(n.p+n.r)*100 # fdr is proportion of FP's out of all positives
		}
	}
	return(c(fdrs,quants,stats))
}
