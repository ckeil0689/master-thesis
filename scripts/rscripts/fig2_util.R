##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Nov 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

boxplots.wt.vs.ko.rorc <- function(ll,sl.wt="SL..",sl.ko="SL..",tf.to.s,main1="",main2="",ylab1="",ylab2="",read.total,plot.pval=T,add.pval=T){
	clust.type <- sapply(ll,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
	trgts <- sapply(ll,function(i) paste(sep="",sort(unique(i$trgts)),collapse="_"))
	# get reads per cluster per milion (to normalize for library size) +1 pseudo count to avoid division by zeros
	p.c <- 1 # pseudo count to prevent division by zero
	wt <- (sapply(ll,function(i) i[[sl.wt]]) + 1)/read.total[sl.wt]
	ko <- (sapply(ll,function(i) i[[sl.ko]]) + 1)/read.total[sl.ko]

	ix.stat3 <- which(clust.type=="STAT3")
	ix.irf4 <- which(clust.type=="IRF4")
	ix.rorc <- which(clust.type=="RORC")
	ix.rorc.stat3 <- which(clust.type=="RORC_STAT3")
	ix.rorc.irf4 <- which(clust.type=="IRF4_RORC")
	ix.rorc.irf4.stat3 <- which(clust.type=="IRF4_RORC_STAT3")

	fc.stat3 <- log2((wt[ix.stat3])/(ko[ix.stat3]))
	fc.irf4 <- log2((wt[ix.irf4])/(ko[ix.irf4]))
	fc.rorc <- log2((wt[ix.rorc])/(ko[ix.rorc]))
	fc.rorc.stat3 <- log2((wt[ix.rorc.stat3])/(ko[ix.rorc.stat3]))
	fc.rorc.irf4 <- log2((wt[ix.rorc.irf4])/(ko[ix.rorc.irf4]))
	fc.rorc.irf4.stat3 <- log2((wt[ix.rorc.irf4.stat3])/(ko[ix.rorc.irf4.stat3]))
	
	w.fc <- list(fc.stat3,fc.irf4,fc.rorc,fc.rorc.stat3,fc.rorc.irf4,fc.rorc.irf4.stat3)
	w.fc.l <- sapply(w.fc,length)

	pval.stat3 <- sapply(ll[ix.stat3],function(i) unlist(i$pval))
	pval.irf4 <- sapply(ll[ix.irf4],function(i) unlist(i$pval))
	pval.rorc <- sapply(ll[ix.rorc],function(i) unlist(i$pval))
	pval.stat3.in.rorc.stat3 <- sapply(ll[ix.rorc.stat3],function(i) i$pval[which(i$tf.to.s=="STAT3")])
	pval.rorc.in.rorc.stat3 <- sapply(ll[ix.rorc.stat3],function(i) i$pval[which(i$tf.to.s=="RORC")])
	pval.irf4.in.rorc.irf4 <- sapply(ll[ix.rorc.irf4],function(i) i$pval[which(i$tf.to.s=="IRF4")])
	pval.rorc.in.rorc.irf4 <- sapply(ll[ix.rorc.irf4],function(i) i$pval[which(i$tf.to.s=="RORC")])
	pval.stat3.in.rorc.irf4.stat3 <- sapply(ll[ix.rorc.irf4.stat3],function(i) i$pval[which(i$tf.to.s=="STAT3")])
	pval.irf4.in.rorc.irf4.stat3 <- sapply(ll[ix.rorc.irf4.stat3],function(i) i$pval[which(i$tf.to.s=="IRF4")])
	pval.rorc.in.rorc.irf4.stat3 <- sapply(ll[ix.rorc.irf4.stat3],function(i) i$pval[which(i$tf.to.s=="RORC")])
	w.pval <- list(pval.stat3,pval.irf4,pval.rorc,
		pval.stat3.in.rorc.stat3,pval.rorc.in.rorc.stat3,
		pval.irf4.in.rorc.irf4,pval.rorc.in.rorc.irf4,
		pval.stat3.in.rorc.irf4.stat3,pval.irf4.in.rorc.irf4.stat3,pval.rorc.in.rorc.irf4.stat3)
	w.l <- sapply(w.pval,length)			
	if(plot.pval==T){
		wilcox.pvals <- list()
		wilcox.pvals[["stat3[S]"]] <- wilcox.test(pval.stat3,pval.stat3,alternative="less")
		wilcox.pvals[["irf4[I]"]] <- wilcox.test(pval.irf4,pval.irf4,alternative="less",paired=T)
		wilcox.pvals[["rorc[R]"]] <- wilcox.test(pval.rorc,pval.rorc,alternative="less")
		wilcox.pvals[["stat3[SR]"]] <- wilcox.test(pval.stat3,pval.stat3.in.rorc.stat3,alternative="less")
		wilcox.pvals[["rorc[SR]"]] <- wilcox.test(pval.rorc,pval.rorc.in.rorc.stat3,alternative="less")
		wilcox.pvals[["irf4[IR]"]] <- wilcox.test(pval.irf4,pval.irf4.in.rorc.irf4,alternative="less")
		wilcox.pvals[["rorc[IR]"]] <- wilcox.test(pval.rorc,pval.rorc.in.rorc.irf4,alternative="less")
		wilcox.pvals[["stat3[ISR]"]] <- wilcox.test(pval.stat3,pval.stat3.in.rorc.irf4.stat3,alternative="less")
		wilcox.pvals[["irf4[ISR]"]] <- wilcox.test(pval.irf4,pval.irf4.in.rorc.irf4.stat3,alternative="less")
		wilcox.pvals[["rorc[ISR]"]] <- wilcox.test(pval.rorc,pval.rorc.in.rorc.irf4.stat3,alternative="less")
		pvals <- sapply(wilcox.pvals,function(i) i$p.value)	
		labels.pval <- paste(sep="",c("stat3\n[S]\n","irf4\n[I]\n","rorc\n[R]\n",
		"stat3\n[SR]\n","rorc\n[SR]\n",
		"irf4\n[IR]\n","rorc\n[IR]\n",
		"stat3\n[ISR]\n","irf4\n[ISR]\n","rorc\n[ISR]\n"),"(",w.l,")\n",format(pvals,scientific=T,digits=1))	
		
		wilcox.fc <- list()
		wilcox.fc[["[S]"]] <- wilcox.test(fc.irf4,fc.stat3,alternative="l")
		wilcox.fc[["[I]"]] <- wilcox.test(fc.irf4,fc.irf4,alternative="l",paired=T)
		wilcox.fc[["[R]"]] <- wilcox.test(fc.irf4,fc.rorc,alternative="l")
		wilcox.fc[["[SR]"]] <- wilcox.test(fc.irf4,fc.rorc.stat3,alternative="l")
		wilcox.fc[["[IR]"]] <- wilcox.test(fc.irf4,fc.rorc.irf4,alternative="l")
		wilcox.fc[["[ISR]"]] <- wilcox.test(fc.irf4,fc.rorc.irf4.stat3,alternative="l")
		pvals <- sapply(wilcox.fc,function(i) i$p.value)
		labels.fc <- paste(sep="",c("[S]","[I]","[R]",
		"[SR]","[IR]","[ISR]"),"\n(",w.fc.l,")\n",format(pvals,scientific=T,digits=1))
	} else {
		labels.pvals <- paste(sep="",c("stat3\n[S]\n","irf4\n[I]\n","rorc\n[R]\n",
		"stat3\n[SR]\n","rorc\n[SR]\n",
		"irf4\n[IR]\n","rorc\n[IR]\n",
		"stat3\n[ISR]\n","irf4\n[ISR]\n","rorc\n[ISR]\n"),"(",w.l,")")
		labels.fc <- paste(sep="",c("[S]","[I]","[R]",
		"[SR]","[IR]","[ISR]"),"\n(",w.fc.l,")")
	}

	if(add.pval==TRUE){
		ix <- which(sapply(w.pval,function(i) if(length(i)==0){1}else{0})==1)
		if(length(ix)>0){
			for(i in 1:length(ix)){
				w.pval[[ix[i]]] <- numeric()
			}
		}
		boxplot(w.pval,outline=FALSE,main=main1,ylab=ylab2,col="gray",xaxt="n")	
		axis(1,at=1:length(w.pval),labels=labels.pvals,padj=+.5,cex.axis=.75)
	}

	
	boxplot(w.fc,outline=T,main=main2,ylab=ylab1,col="gray",xaxt="n")	
	axis(1,at=1:length(w.fc),labels=labels.fc,padj=+.5,cex.axis=.75)

}

boxplots.wt.vs.ko <- function(ll,sl.wt="SL..",sl.ko="SL..",tf.to.s,main1="",main2="",ylab1="",ylab2="",read.total){
	clust.type <- sapply(ll,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
	wt <- sapply(ll,function(i) i[[sl.wt]])
	ko <- sapply(ll,function(i) i[[sl.ko]])
	ix.batf <- which(clust.type=="BATF")
	ix.irf4 <- which(clust.type=="IRF4")
	ix.batf.irf4 <- which(clust.type=="BATF_IRF4")
	const <- reads.total[sl.wt]/reads.total[sl.ko]# normalize to library size of wt
	p.c <- 1 # pseudo count to prevent division by zero
	fc.batf <- log2((wt[ix.batf]+p.c)/((ko[ix.batf]+p.c)*const))
	fc.irf4 <- log2((wt[ix.irf4]+p.c)/((ko[ix.irf4]+p.c)*const))
	fc.batf.irf4 <- log2((wt[ix.batf.irf4]+p.c)/((ko[ix.batf.irf4]+p.c)*const))
	
	pval.batf <- sapply(ll[ix.batf],function(i) unlist(i$pval))
	pval.irf4 <- sapply(ll[ix.irf4],function(i) unlist(i$pval))
	pval.irf4.in.batf.irf4 <- sapply(ll[ix.batf.irf4],function(i) i$pval[which(i$tf.to.s=="IRF4")])
	pval.batf.in.batf.irf4 <- sapply(ll[ix.batf.irf4],function(i) i$pval[which(i$tf.to.s=="BATF")])
	
	boxplot(list(pval.batf,pval.irf4,pval.irf4.in.batf.irf4,pval.batf.in.batf.irf4),outline=FALSE,
	main=main1,
	ylab=ylab2,col="gray",xaxt="n")	
	axis(1,at=1:4,labels=c("batf\n[batf]","irf4\n[irf4]","irf4\n[irf4,batf]","batf\n[irf4,batf]"),padj=+.5)
	
	ix.batf.coop <- which(pval.batf > 500)
	ix.irf4.coop <- which(pval.irf4 > 500)
	ix.irf4.in.batf.irf4.coop <- which(pval.irf4.in.batf.irf4  > 500)
	ix.batf.in.batf.irf4.coop <- which(pval.batf.in.batf.irf4  > 500)
	ix.batf.irf4.coop <- intersect(ix.irf4.in.batf.irf4.coop,ix.batf.in.batf.irf4.coop)
	
	ix.batf.non.coop <- which(pval.batf < 200)
	ix.irf4.non.coop <- which(pval.irf4 < 200)
	ix.irf4.in.batf.irf4.non.coop <- which(pval.irf4.in.batf.irf4  < 200)
	ix.batf.in.batf.irf4.non.coop <- which(pval.batf.in.batf.irf4  < 200)
	ix.batf.irf4.non.coop <- intersect(ix.irf4.in.batf.irf4.non.coop,ix.batf.in.batf.irf4.non.coop)
	
	tmp <- list(fc.batf,fc.batf[ix.batf.non.coop],fc.batf[ix.batf.coop],
		fc.irf4,fc.irf4[ix.irf4.non.coop],fc.irf4[ix.irf4.coop],
		fc.batf.irf4,fc.batf.irf4[ix.batf.irf4.non.coop],fc.batf.irf4[ix.batf.irf4.coop]
		)
	names(tmp) <- 	c("[batf]\n(all pvals)","[batf]\n(pval<200)","[batf]\n(pval>500)",
					"[irf4]\n(all pvals)","[irf4]\n(pval<200)","[irf4]\n(pval>500)",
					"[irf4,batf]\n(all pvals)","[irf4,batf]\n(pval<200)","[irf4,batf]\n(pval>500)")
	boxplot(tmp,outline=F,main=main2,ylab=ylab1,col="gray",xaxt="n")
	axis(1,at=1:length(tmp),labels=names(tmp),cex.axis=0.6)	
}

count.reads.to.ll.list <- function(ll,o){
	tmp <- list()
	for(i in 1:length(o)){
		cr <- names(o)[i]
		tmp[[cr]] <- o[cr]
	}
	o <- tmp
	rm(tmp)
	w <- mclapply(o,count.reads.to.ll.list.for.single.chr, ll,mc.cores=n.p)
	counts <- unlist(w,use.names=FALSE)
}

count.reads.to.ll.list.for.single.chr <- function(o.r,ll){
	cr <- names(o.r)
	o.r <- o.r[[1]]
	# init
	# l goes over MTBS's from ll
	tmp <- sapply(ll,function(i) i$chr[1])
	ix <- which(tmp == cr)
	# the output is the read count per multi TF binding sites (MTBS)
	read.count <- numeric(length=length(ix))
	l <- ix[1]
	l.l <- ll[[l]]$start # left side of MTBS
	r.l <- ll[[l]]$end # right side of MTBS
	# j goes over reads
	j <- 1
	# read.cnt goes over the read.count vector
	read.cnt <- 1
	stop.l <- ix[length(ix)] + 1
	stop.j <- length(o.r) + 1
	while(l < stop.l){ # the real stop point is inside the loop; this is for readability
		if(o.r[j] < l.l){
			j <- j+1
			if(j == stop.j){
				break
			}
		} else {
			j.start <- j
			while(o.r[j] < r.l){
				j <- j+1
				if(j == stop.j){
					break
				}
			}
			read.count[read.cnt] <- j-j.start
			read.cnt <- read.cnt + 1 # add pseudo count
			l <- l+1
			if(l == stop.l){
				break
			}
			l.l <- ll[[l]]$start # left side of MTBS
			r.l <- ll[[l]]$end # right side of MTBS
			# back track a few reads if needed
			while(o.r[j] > l.l){
				j <- j-1
				if(j==0){
					j <- 1
					break
				}
			}
		}
		if(j%%10000==0){
			cat(".")
		}
	}
	cat("\n")
	return(read.count)
}


# get.pileup.per.chr.parallel
get.pileup.per.chr.parallel <- function(x.i,r.u=r.u,N,n.p){
	w.pos <- mclapply(x.i,get.pileup.for.one.chr.pos.strand,N=N,mc.cores=n.p)
	w.neg <- mclapply(x.i,get.pileup.for.one.chr.neg.strand,N=N,mc.cores=n.p)
	pl <- list()
	pl.pos <- w.pos[[1]]
	for(l in 2:length(w.pos)){
		pl.pos <- pl.pos+w.pos[[l]]
	}
	pl.neg <- w.neg[[1]]
	for(l in 2:length(w.neg)){
		pl.neg <- pl.neg+w.neg[[l]]
	}
	pl <- list(pl.pos,pl.neg)
	return(pl)
}

# deal with positive strand
get.pileup.for.one.chr.pos.strand <- function(x.i.i,N){
	r <- names(x.i.i)
	x.i.i <- x.i.i[[r]]
	ix <- which(r.u[,"chrom"]==r & r.u[,"strand"]=="+")
	v <- r.u[ix,"txStart"]
	pl.pos <- numeric(length=2*N)
	pl <- numeric(length=2*N)		
	for(l in 1:length(v)){
		ix <- which(abs(x.i.i-v[l])<N )
		ix.up <- which(x.i.i[ix]<=v[l])
		ix.down <- which(x.i.i[ix]>v[l])
		ix.up.norm <- x.i.i[ix[ix.up]]-v[l]+N
		ix.down.norm <- x.i.i[ix[ix.down]]-v[l]+N
		pl[unique(c(ix.up.norm,ix.down.norm))] <- pl[unique(c(ix.up.norm,ix.down.norm))] + 1
		# cat(".")
	}
	return(pl)
}
# deal with netative strand
get.pileup.for.one.chr.neg.strand <- function(x.i.i,N){
	r <- names(x.i.i)
	x.i.i <- x.i.i[[r]]
	ix <- which(r.u[,"chrom"]==r & r.u[,"strand"]=="-")
	v <- r.u[ix,"txStart"]
	pl.pos <- numeric(length=2*N)
	pl <- numeric(length=2*N)		
	for(l in 1:length(v)){
		ix <- which(abs(x.i.i-v[l])<N )
		ix.up <- which(x.i.i[ix]<=v[l])
		ix.down <- which(x.i.i[ix]>v[l])
		ix.up.norm <- x.i.i[ix[ix.up]]-v[l]+N
		ix.down.norm <- x.i.i[ix[ix.down]]-v[l]+N
		pl[unique(c(ix.up.norm,ix.down.norm))] <- pl[unique(c(ix.up.norm,ix.down.norm))] + 1
		# cat(".")
	}
	return(pl)
}

## split.matrix.to.chrmosomes
# input:
#       - x: a list of matrices each having one column named "chr"
#		- chr: a vector with the names of the chromosomes we are interested in extracting info for
# output:
#      - l: a list with as many elements as chr specified by input.
split.matrix.to.chrmosomes <- function(x,chr="chr1",n.p=3){
  l <- list()
  # cat("wroking on expt", e,"\n")
  l <- mclapply(chr,split.matrix.to.chromosome,x=x,mc.cores=n.p)
  names(l) <- chr
  return(l)
}
# a helper function for split.matrix.to.chrmosomes, gives for one chromosome the matrix's rows for expt e
# input:
#     - r is chromosome
#     - x is a list with two sub list named chr and another named read_5_prime_bp_location
split.matrix.to.chromosome <- function(r,x){
  # ix.chr.i <- which(m[,"chr"]==r)
  ix.chr.i <- which(x[["chr"]]==r)
  o <- x[["read_5_prime_bp_location"]][ix.chr.i]
  o <- sort(o)
  return(o)
  # cat("working on",r,"\n")
  # o <- list()
  # o[[r]] <- x[["read_5_prime_bp_location"]][ix.chr.i]
  # ix <- sort(o[[r]],index.return=TRUE)$ix
  # return(ix.chr.i[ix])
}

## make.ruler.for.reads makes a ruler: a line per chromosome with the locations of all reads (5' base) (sorted from start of chromosome to end)
# input:
#     - chr: a list of chromosome names
#     - reads.list.per.chrom: a list of reads for each chromosome for all expts
# output:
#     - o: a list each chormosome ruler as an element
make.ruler.for.reads <- function(chr,x){
  o <- list()
  for(j in 1:length(chr)){
    r <- chr[j] # chrmosome we go over
    s <- numeric()
    expt <- character()
    genetic.bg <- names(x) # the genetic bg of chip
    for(i in 1:length(genetic.bg)){
      o[[r]] <- list()
      e <- genetic.bg[i] #experiment name
      s <- c(s,x[[e]][[r]])
      expt <- c(expt,rep(e,length(x[[e]][[r]])))
    }
    ix <- sort(s,index.return=TRUE)$ix
    o[[r]] <- list(s=s[ix],expt=expt[ix])
  }
  return(o)
}
