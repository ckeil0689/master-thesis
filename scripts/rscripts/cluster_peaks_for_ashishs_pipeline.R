##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Feb 2012 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

rm(list=ls())
debug=F
script.version=0.1
print.error <- function(){
	cat("
DESCRIPTIION:	
	cluster_peaks.R takes MACS tab delimited files as input and produces one tab delimeted file (named mtls.xls) where 
	each row corresponds to a Multi TF Loci (MTL) in which peaks from different experiments (input MACS files)
	fall within a certain distance between summits from eachother.

INPUT:
	1.path_input=path to MACS files '::' delim [path_input=f1::f2::f3::...::fk]
	2.path_output=path to save generated MTL cluster file (where to save mtls.xls)
	3.expt_names=user specified names for MACS files '::' delim [expt_names=n1::n2::n3::...::nk]
	4.dist.summits=maximum distance between summits belonging to the same MTL
	5.n.autosome.chr=19 for mouse, 22 for human
	
EXAMPLE RUN: 
	cluster_peaks.R
	input_macs_files=SL2870_SL2871_peak.xls::SL2872_SL2876_peak.xls::SL3032_SL2871_peak.xls::SL3037_SL3036_peak.xls::SL3315_SL3319_peak.xls
	path_output=~/Desktop/
	expt_names=RORC_Th17::IRF4_Th17::MAF_Th17::BATF_Th17::STAT3_Th17
	dist_summits=100
	n_autosome_chr=19
	
Please cite us if you used this script: 
	The transcription factor network regulating Th17 lineage specification and function.
	Maria Ciofani, Aviv Madar, Carolina Galan, Kieran Mace, Agarwal, Kim Newberry, Richard M. Myers, 
	Richard Bonneau and Dan R. Littman et. al. (in preperation)\n\n")
}

############# helper functions:
## split.macs.list.to.chrmosomes
# input:
#      - macs.list: a list of macs expts: here are a few lines of one expt
## chr	start	end	length	summit	tags	#NAME?	fold_enrichment	FDR(%)
## chr1	4322210	4323069	860	494	55	158.95	6.03	0.05
## chr1	4797749	4798368	620	211	29	119.82	3.47	0.09
## chr1	4848182	4849113	932	494	46	105.42	2.9	0.09
#      - expts: a list of the expts names from macs.list that you want to process
#      - chr: chrmosomes names
# output:
#      - x: a list with as many elements as chr specified by input.
#      -x[[i]]:macs list per chr, with peak.id column added (these are the row numbers of MACS)
split.macs.list.to.chrmosomes.no.pml <- function(macs.list,expts,chr="chr1"){
  x <- list()
  n <- length(expts)
  for(i in 1:n){
    e <- expts[i] #experiment name
    cat("wroking on expt", e,"\n")
    x[[e]] <- lapply(chr,split.one.macs.expt.by.chromosome.no.pml,m=macs.list[[e]])
    names(x[[e]]) <- chr
  }
  return(x)
}
# a helper function for spliat.macs.list.to.chrmosomes, gives for one chromosome the macs rows for expt MACS matrix m
# input:
#     - r is chromosome
#     - m is macs matrix for expt e from above function
split.one.macs.expt.by.chromosome.no.pml <- function(r,m){
  ix.chr.i <- which(m[,"chr"]==r)
  # cat("working on",r,"\n")
  o <- list()
  o[[r]] <- m[ix.chr.i,]
  o[[r]]$peak.id <- ix.chr.i
  return(o[[r]])
}

## make.ruler makes a ruler: a line per chromosome with the locations of all tf binding sites (sorted from start of chromosome to end)
# also match these summit locations with corresponding:
# pvals, tfs, peak start and peak end trgts
# input:
#     - chr: a list of chromosome names
#     - macs.list.per.chrom: a list of macs peaks for each chromosome
# output:
#     - o: a list each chormosome ruler as an element
make.ruler.no.pml <- function(chr,macs.list.per.chrom){
  x <- macs.list.per.chrom
  o <- list()
  for(j in 1:length(chr)){
    r <- chr[j] # chrmosome we go over
    s <- numeric()
    pval <- numeric()
    tf.to.s <- character()
    trgt.prox <- character()
    trgt.dist <- character()
    dtss.prox <- character() 
    dtss.dist <- character()
    start <- numeric()
    end <- numeric()
    trtmnts <- names(x) # the treatments name (expt names)
    ## debug parameters ###
    ## which experiment peaks come from
    expt <- character()
    ## what was the peak id in that experiment
    peak.ids <- numeric()
    ## this will allow us to always back track from a cluster to the actual peaks in it
    ## debug params end ###
    for(i in 1:length(trtmnts)){
      o[[r]] <- list()
      e <- trtmnts[i] #experiment name
      tf <- strsplit(e,"_")[[1]][1]
      s <- c(s,x[[e]][[r]][,"start"]+x[[e]][[r]][,"summit"])
      pval <- c(pval,x[[e]][[r]][,7]) # the name of 7th column is X.10.log10.pvalue. this is pval
      # tf.to.s <- c(tf.to.s,rep(tf,length(x[[e]][[r]][,7]))) # all summits belong to tf
      start <- c(start,x[[e]][[r]][,"start"])
      end <- c(end,x[[e]][[r]][,"end"])
      expt <- c(expt,rep(e,length(x[[e]][[r]][,"end"])))
      peak.ids <- c(peak.ids,x[[e]][[r]][,"peak.id"])    
      trgt.prox <- c(trgt.prox,x[[e]][[r]][,"trgt.prox"])
	  trgt.dist <- c(trgt.dist,x[[e]][[r]][,"trgt.dist"])
      dtss.prox <- c(dtss.prox,x[[e]][[r]][,"dtss.prox"])
      dtss.dist <- c(dtss.dist,x[[e]][[r]][,"dtss.dist"])
    }
    ix <- sort(s,index.return=TRUE)$ix
    o[[r]] <- list(s=s[ix],pval=pval[ix],expt=expt[ix],start=start[ix],end=end[ix],peak.ids=peak.ids[ix],
				   trgt.prox=trgt.prox[ix],trgt.dist=trgt.dist[ix],dtss.prox=dtss.prox[ix],dtss.dist=dtss.dist[ix])
  }
  return(o)
}

## add cluster memberships based on ruler
## require no more than d.cut distance between tf summits
## cur.l is the current loci number (or cluster number)
assign.clusters.ids.to.peaks <- function(ruler,d.cut){
 cur.l <- 0
 for(j in 1:length(ruler)){
   s <- ruler[[j]][["s"]]
   l <- numeric(length=length(s))
   if(length(l)>0){
	   cur.l <- cur.l+1
   }
   if(length(l)==1){
   	l[1] <- cur.l
   } else if(length(l)>1) {
	l[1] <- cur.l
	   for(i in 2:length(l)){
	     d <- s[i]-s[i-1] # assumes s is sorted increasingly
	     if(d>d.cut){
	       cur.l <- cur.l+1
	     }
	     l[i] <- cur.l    
	   }
	}
   ruler[[j]][["l"]] <- l
 }
 return(ruler)
}

## here we create a list of TF enriched loci (clusters)
## input:
##     - ruler: a line per chromosome with the locations of all tf binding sites (sorted from start of chromosome to end)
## output:
##     -ll: a list of clusters ll where each cluster (elem in ll holds):
##        - l: the current loci (elem in ll)
##        - s: summits vector as before
##        - pval: pvals vector matching the summits
##        - tfs: a vector of tfs matching s, the summits in loci l
##        - spans: a vector of spans matching s, the summits in loci l, where spans is the dist between start and end of a peak
##        - trgts: genes that are targeted by the cluster (10kb upstream, in gene, 10kb downstream)
create.cluster.list.no.pml <- function(ruler){
  tmp <- list()
  ll <- list()
  for(j in 1:length(ruler)){
    r <- names(ruler)[j]
    # cat("working on ",r,"\n")
    x <- ruler[[j]] # for short typing let x stand for ruler[[chr[j]]]
    l.vec <- unique(x[["l"]]) # the clusters ids on j chr (only unique)
    n <- length(l.vec) # iterate over n clusters
    tmp[[j]] <- lapply(1:n,get.cluster.params.no.pml,x=x,l.vec=l.vec,r=r)
  }
  ## concatenate tmp into one list
  command <- paste(sep="","c(",paste(sep="","tmp[[",1:length(tmp),"]]",collapse=","),")")
  ll=eval(parse(text=command))
  return(ll)
}
get.cluster.params.no.pml <- function(i,x,l.vec,r){
  ix <- which(x[["l"]]==l.vec[i])
  l <- l.vec[i]
  start <- min(x[["start"]][ix])
  end <- max(x[["end"]][ix])
  s <- x[["s"]][ix]
  # tf.to.s <- x[["tf.to.s"]][ix]
  pval <- x[["pval"]][ix]
  pval.mean <- mean(pval)
  span.tfs <- x[["end"]][ix]-x[["start"]][ix]
  span.l <- end-start
  peak.ids <- x[["peak.ids"]][ix]
  expt <- x[["expt"]][ix]
  expt.alphanum.sorted <- sort(x[["expt"]][ix])
  trgt.prox <- unique(x[["trgt.prox"]][ix])
  trgt.dist <- unique(x[["trgt.dist"]][ix])
  dtss.prox <- unique(x[["dtss.prox"]][ix])
  dtss.dist <- unique(x[["dtss.dist"]][ix])
  chr <- rep(r,length(ix))
  return(list(l=l,chr=chr,expt=expt,expt.alphanum.sorted=expt.alphanum.sorted,start=start,
			  end=end,s=s,pval=pval,pval.mean=pval.mean,span.tfs=span.tfs,
			  span.l=span.l,peak.ids=peak.ids,trgt.prox=trgt.prox,trgt.dist=trgt.dist,
			  dtss.prox=dtss.prox,dtss.dist=dtss.dist))
}

## here we create a list of TF enriched loci (clusters)
## input:
##     - ruler: a line per chromosome with the locations of all tf binding sites (sorted from start of chromosome to end)
## output:
##     -ll: a list of clusters ll where each cluster (elem in ll holds):
##        - l: the current loci (elem in ll)
##        - s: summits vector as before
##        - pval: pvals vector matching the summits
##        - tfs: a vector of tfs matching s, the summits in loci l
##        - spans: a vector of spans matching s, the summits in loci l, where spans is the dist between start and end of a peak
##        - trgts: genes that are targeted by the cluster (10kb upstream, in gene, 10kb downstream)
create.cluster.list.no.pml <- function(ruler){
  tmp <- list()
  ll <- list()
  for(j in 1:length(ruler)){
    r <- names(ruler)[j]
    # cat("working on ",r,"\n")
    x <- ruler[[j]] # for short typing let x stand for ruler[[chr[j]]]
    l.vec <- unique(x[["l"]]) # the clusters ids on j chr (only unique)
    n <- length(l.vec) # iterate over n clusters
    tmp[[j]] <- lapply(1:n,get.cluster.params.no.pml,x=x,l.vec=l.vec,r=r)
  }
  ## concatenate tmp into one list
  command <- paste(sep="","c(",paste(sep="","tmp[[",1:length(tmp),"]]",collapse=","),")")
  ll=eval(parse(text=command))
  return(ll)
}
get.cluster.params.no.pml <- function(i,x,l.vec,r){
  ix <- which(x[["l"]]==l.vec[i])
  l <- l.vec[i]
  start <- min(x[["start"]][ix])
  end <- max(x[["end"]][ix])
  s <- x[["s"]][ix]
  # tf.to.s <- x[["tf.to.s"]][ix]
  pval <- x[["pval"]][ix]
  pval.mean <- mean(pval)
  span.tfs <- x[["end"]][ix]-x[["start"]][ix]
  span.l <- end-start
  peak.ids <- x[["peak.ids"]][ix]
  expt <- x[["expt"]][ix]
  expt.alphanum.sorted <- sort(x[["expt"]][ix])
  trgt.prox <- unique(x[["trgt.prox"]][ix])
  trgt.dist <- unique(x[["trgt.dist"]][ix])
  dtss.prox <- unique(x[["dtss.prox"]][ix])
  dtss.dist <- unique(x[["dtss.dist"]][ix])
  chr <- rep(r,length(ix))
  return(list(l=l,chr=chr,expt=expt,expt.alphanum.sorted=expt.alphanum.sorted,start=start,
			  end=end,s=s,pval=pval,pval.mean=pval.mean,span.tfs=span.tfs,
			  span.l=span.l,peak.ids=peak.ids,trgt.prox=trgt.prox,trgt.dist=trgt.dist,
			  dtss.prox=dtss.prox,dtss.dist=dtss.dist))
}
## pretty print (tab deim) to file requested elements out of chip cluster list, ll.
## input:
##     - ll: a cluster list
##     - f.nm: a file name (include path) to where you want files to print
##     - tfs: a list of the tfs we want to print the file for (the same as the tfs used for the peak clustering)
## output
##     - a tab delim file with clusers as rows and elems tab delim for each cluster
print.ll.verbose.all <- function(ll,f.nm="ll.xls"){
  options(digits=5)
  cat(file=f.nm,names(ll[[1]]),sep="\t")
  cat(file=f.nm,"\n",append=TRUE)
  for(i in 1:length(ll)){
	line <- ll[[i]][[1]] # put MTL number
	line <- paste(line,ll[[i]][[2]][1],sep="\t")
	for(j in 3:length(ll[[i]])){
		val <- ll[[i]][[j]]
		if(length(val) == 1){
			line <- paste(line,val,sep="\t")
		} else {
			line <- paste(line,paste(sep="",unlist(ll[[i]][j]),collapse="_"),sep="\t")
		}
	}
  	cat(file=f.nm,line,"\n",append=TRUE)	
  }
}
############# Code:
# retrieve args
if(debug==T){
	cmd.args <- c(
		"input_macs_files=SL2870_SL2871_peak.xls::SL2872_SL2876_peak.xls::SL3032_SL2871_peak.xls::SL3037_SL3036_peak.xls::SL3315_SL3319_peak.xls",
		"path_output=~/Desktop/",
		"expt_names=RORC_Th17::IRF4_Th17::MAF_Th17::BATF_Th17::STAT3_Th17",
		"dist_summits=100",
		"n_autosome_chr=19"
	)
} else {
	cmd.args <- commandArgs();      
}

if(length(grep("--version",cmd.args))){
	cat("version",script.version,"\n")
	q()
}

arg.names.cmd.line <- sapply(strsplit(cmd.args,"="),function(i) i[1])
args.val.cmd.line <- sapply(strsplit(cmd.args,"="),function(i) i[2])

arg.nms <- c("input_macs_files","path_output","expt_names","dist_summits","n_autosome_chr")
arg.val <- character(length=length(arg.nms))
for(i in 1:length(arg.nms)){
	ix <- which(arg.names.cmd.line==arg.nms[i])
	if(length(ix)==1){
		arg.val[i] <- args.val.cmd.line[ix]
	} else {
		stop("######could not find ",arg.nms[i]," arg######\n\n",print.error())
		
	} 
}
if(debug==T){
	print(paste(arg.nms,"=",arg.val))
}
# the files here adhere to tab delim format
input.macs.files <- strsplit(arg.val[1],"::")[[1]]
if(length(input.macs.files)==1){
	cat("only provided one MACS file to cluster.")
	print.error()
}
path.output <- arg.val[2]
expt.names <- strsplit(arg.val[3],"::")[[1]]
dist.summits <- as.numeric(arg.val[4])
n.autosome.chr <- as.numeric(arg.val[5])

# source("~/Documents/nyu/littmanLab/th17_used_for_paper/r_scripts/th17/used_for_paper/cluster_peaks_util.R")
chr <- c(paste(sep="","chr",1:n.autosome.chr),paste(sep="","chr",c("X","Y")))

# read MACS files
macs.list <- list()
for(i in 1:length(input.macs.files)){
	e <- expt.names[i]
	macs.list[[ e ]] <- read.delim(file=input.macs.files[i])
	macs.list[[ e ]][,"chr"] <- as.character(macs.list[[ e ]][,"chr"])
	macs.list[[ e ]][,"trgt.prox"] <- as.character(macs.list[[ e ]][,"trgt.prox"])
	macs.list[[ e ]][,"trgt.dist"] <- as.character(macs.list[[ e ]][,"trgt.dist"])
}

# take all macs files and put peaks together on each chromosome 
# (as if each chromosome is a ruler and we specify where in the ruler each peak summit falls)
cat("splitting macs files per chromosome\n")
x <- split.macs.list.to.chrmosomes.no.pml(macs.list=macs.list,expts=expt.names,chr=chr)
cat("adding peaks from all macs files into chromosome rulers\n")
ruler <- make.ruler.no.pml(chr,macs.list.per.chrom=x)
cat("add MTL membership to the ruler\n")
ruler <- assign.clusters.ids.to.peaks(ruler,d.cut=dist.summits)
for(i in 1:length(ruler)){
	ix <- which(is.na(ruler[[i]][["dtss.prox"]]))
	ruler[[i]][["dtss.prox"]][ix] <- ""
	ix <- which(is.na(ruler[[i]][["dtss.dist"]]))
	ruler[[i]][["dtss.dist"]][ix] <- ""
}
cat("creating MTL list\n")
ll <- create.cluster.list.no.pml(ruler=ruler)
cat("writing MTL table\n")
# f.nm <- paste(sep="",path.output,paste(expt.names,collapse="_"),"_MTLs.xls")
f.nm <- paste(sep="","mtls",".xls")
print.ll.verbose.all(ll=ll,f.nm=f.nm)









