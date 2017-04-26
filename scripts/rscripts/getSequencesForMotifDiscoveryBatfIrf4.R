##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Apr 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

cat("getting batf/irf4 mtls DNA sequences for figure 2\n")

library(BSgenome.Mmusculus.UCSC.mm9) # to get sequence around each peak
source("rscripts/utilChipScores.R")
source("rscripts/cluster_peaks_util.R")


# set paths
path.input <- "input/th17/used_for_paper/fig2/"
path.results <- "results/fig2/"

if(lineage=="Th17"){
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
	ll.lin <- ll.th17
} else if(lineage=="Th0") {
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th0_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
	ll.lin <- ll.th0
}
max.pval.mtls <- sapply(ll.lin,function(i) max(i$pval))
ix.sig.mtls <- which(max.pval.mtls > mtl.min.pval)
ll.lin <- ll.lin[ix.sig.mtls]

mtl.type <- sapply(ll.lin,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_")) 


ix.list <- list()
ix.list[["batf"]] <- which(mtl.type == "BATF")
ix.list[["irf4"]] <- which(mtl.type == "IRF4")
ix.list[["batf_irf4"]] <- which(mtl.type == "BATF_IRF4")
ix.list[["batf_irf4_plus"]] <- setdiff(grep("BATF_IRF4",mtl.type),ix.list[["batf_irf4"]])

seq.list <- list()
for(t in 1:2){
	if(t==1){
		type <- "p"
	} else {
		type <- "d"
	}
	for( i  in 1:length(ix.list)){
		ll <- ll.lin
		mtl <- names(ix.list)[i]
		ix <- ix.list[[mtl]]
		pval <- sapply(ll[ix],function(i) median(i$pval))
		ix.p <- sort(pval,index.return=T,decreasing=T)$ix
		cnt <- 1
		j <- 1
		ll.tmp <- list()
		max.cnt <- min(seq.count,length(pval))
		while(cnt<=seq.count){
			if( ( (ll[[ix[ix.p][j]]]$trgt[1] == "") & (type == "d") ) 
			 	|
				( (ll[[ix[ix.p][j]]]$trgt[1] != "") & (type == "p"))
				){
				ll.tmp <- c(ll.tmp,ll[ix[ix.p][j]])
				cnt <- cnt+1
			}
			j <- j+1
		}
		ll <- ll.tmp # limit to seq.count mtls with highest median pval
		chr <- sapply(ll,function(i) i$chr[1])
		summit <- sapply(ll,function(i) median(i$s))
		pval <- sapply(ll,function(i) median(i$pval))
		ll.id <- sapply(ll,function(i) median(i$l))
		ll.trgt <- sapply(ll,function(i) paste(sep="",sort(unique(i$trgts)),collapse="_")) 
		seq.list[[paste(sep="",type,"_",names(ix.list)[i])]] <- data.frame(stringsAsFactors=F,
		                       chr=chr,
		                       start=summit-peak.region.seq.length,
		                       end=summit+peak.region.seq.length,
		                       strand=rep("+",length(summit)),
							   ll.id=paste(ll.id,ll.trgt,sep="_"),
							   mtl=rep(paste(sep="",type,"_",names(ix.list)[i]),length(summit)),
							   pval=pval
		                       )
	}
}
stop("AM")
myseqs <- seq.list[[1]]
for(i in 2:length(seq.list)){
	myseqs <- rbind(myseqs,seq.list[[i]])
}
seqs <- getSeq(Mmusculus, names=myseqs$chr, start=myseqs$start, end=myseqs$end,strand=myseqs$strand)	

# write sequences for each type of mtl seperating proximal from distal
for(i in 1:length(seq.list)){
	mtl <- names(seq.list)[i]
	fl.nm <- paste(path.results,lineage,"_",mtl,"_","DNA_seq_",date.is,".fasta",sep="")
	ix <- which(myseqs$mtl==mtl)
	## write fasta file for results
	cat("writing peak sequences for ", mtl,"\n")
	append.to <- FALSE
	x <- myseqs[ix,]
	sequences <- seqs[ix]
	for(j in 1:nrow(x)){
		if(j > 1){append.to <- TRUE}
		cat(sep="",">",x$ll.id[j],"|","pval=",x$pval[j],"\n",as.character(sequences[j]),"\n",file=fl.nm,append=append.to)
	}
}

