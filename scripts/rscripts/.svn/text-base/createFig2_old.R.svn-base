##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Nov 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# script to create figure 2 
# 1) we look at BATF ko and see what happens to IRF4 peaks
# 2) we look at IRF4 ko and see what happens to BATF peaks
# 3) we look at RORC ko and see what happens to IRF4 and STAT3 peaks

library(multicore)
source("r_scripts/th17/used_for_paper/cluster_peaks_util.R")

# set paths
path.sam.big <- "/Volumes/Drobo/aviv/th17/sam_files/"
path.sam.small <- "/Volumes/Drobo/aviv/th17/sam_files_small/"
path.macs.tab.delim <- paste(sep="","input/th17/used_for_paper/rawData/MACS_tab_delim_",date.macs.processed,"/")
path.input.refseq <- "input/th17/used_for_paper/rawData/"
path.sam.files <- "/Volumes/Drobo/aviv/th17/sam_files_small_Nov_7_2011/"
#~~~~~~~~~~~~~~~~~~~~~~~~1~~~~~~~~~~~~~~~~~~~~~~~~~~~~#	
# to make script run faster and more mem efficient first lets 'treat' sam files to remove unused columns and rows
if(pre.process == TRUE){
	all.sam.files <- list.files(path.sam.big)
	skip.files <- list.files(path.sam.small)
	all.sam.files <- setdiff(all.sam.files,skip.files )
	cat("preprocessing:\n")
	for (i in 1:length(all.sam.files)){
		fl.nm <- all.sam.files[i]
		cat("working on",fl.nm,"\n")
		system(paste(sep="","awk '{if ($2==0 || $2==16) print $3 \"\t\" $4}' ",
							path.sam.big,fl.nm," > ",path.sam.small,fl.nm))
	} 
}
#~~~~~~~~~~~~~~~~~~~~~~~~2~~~~~~~~~~~~~~~~~~~~~~~~~~~~#	
# read refseq gene annotation table
## refseq <- read.delim(paste(sep="",path.input.refseq,"20101005_UCSC_mm9_kgXref.csv"))
refseq <- read.delim(paste(sep="",path.input.refseq,"UCSC_mm9_refseq_genes_",date.macs.processed,".txt"))
cat("reading refseq table\n")
# refseq can have many transcripts for each gene
# here i make it have only one transcript for each gene (the longest one)
refseq$name2 <- as.character(refseq$name2)
refseq$chrom <- as.character(refseq$chrom)
refseq$strand <- as.character(refseq$strand)
gn.nms <- unique(refseq$name2)
# create refseq unique r.u
n <- length(gn.nms)
r.u <- data.frame(cbind(chrom=rep("NA",n),strand=rep("NA",n),txStart=rep(0,n),txEnd=rep(0,n)),stringsAsFactors=FALSE,row.names=gn.nms)
cat("making a unique (per gene name) refseq table (by taking the longest transcript of each gene)\n")
for(i in 1:n){
  ix <- which(refseq$name2==gn.nms[i])
  if(length(ix)==0) {
    error("could not find gene", ng.nms[i], "in refseq table.  Bailing out...\n")
  } else if (length(ix)>1){
    l <- apply(refseq[ix,c("txStart","txEnd")],1,function(i) abs(i[1]-i[2]) )
    l.max <- max(l)
    ix <- ix[which(l==l.max)[1]]
  }
  r.u[gn.nms[i],c("chrom","strand","txStart","txEnd")] <- refseq[ix,c("chrom","strand","txStart","txEnd")]
}
r.u[,"txStart"] <- as.numeric(r.u[,"txStart"])
r.u[,"txEnd"] <- as.numeric(r.u[,"txEnd"])

# switch TSS and TES if we have chr "-"
ix <- which(r.u$strand=="-")
tmp.tes <- r.u$txStart[ix]
tmp.tss <- r.u$txEnd[ix]
r.u[ix,c("txStart","txEnd")] <- cbind(tmp.tss,tmp.tes)
#~~~~~~~~~~~~~~~~~~~~~~~~3~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
chr <- c(paste(sep="","chr",1:19),paste(sep="","chr",c("X","Y")))

# read sam files for BATF ko irf4 chip
batf.ko.irf.chip.files <- c("SL6629.sam","SL6633.sam")
names(batf.ko.irf.chip.files) <- c("batf.ko.irf.chip.1","batf.ko.irf.chip.2")
batf.wt.irf.chip.files <- c("SL6628.sam","SL6632.sam")
names(batf.wt.irf.chip.files) <- c("batf.wt.irf.chip.1","batf.wt.irf.chip.2")
all.files <- c(batf.ko.irf.chip.files[1],batf.wt.irf.chip.files[1])
sam.list <- list()
cat("create sam files list for:", all.files,"\n")
for(i in 1:length(all.files)){
	nm <- names(all.files)[i]
	cat("reading:", nm,"\n")
	sam.list[[nm]] <- read.table(paste(sep="",path.sam.files,all.files[i]),sep="\t")
	colnames(sam.list[[nm]]) <- c("chr","read_5_prime_bp_location")
}
#~~~~~~~~~~~~~~~~~~~~~~~~4~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
reads.list.per.chrom <- split.matrix.to.chrmosomes(sam.list[[1]],expts=names(sam.list)[1],chr=chr)
rm(sam.list)
cat("create read-pileups for:", all.files,"\n")
x <- list()
for(i in 1:length(all.files)){
	e <- names(all.files)[i]
	x[[e]] <- list()
	for (j in 1:length(chr)){
		r <- chr[j]
		x[[e]][[r]] <- list()
		x[[e]][[r]][[r]] <- reads.list.per.chrom[[e]][[r]]
	}
}
#~~~~~~~~~~~~~~~~~~~~~~~~5~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot density around TSS for both experiments
# create read-pileup around TSS for each experiments
p <- list()
N=10*10^3
cat("create read-pileups for:", all.files,"\n")
for(i in 1:length(all.files)){
	e <- names(all.files)[i]
	cat("working on expt:", nm,"\n")
	p[[e]] <- get.pileup.per.chr.parallel(x.i=x[[e]],r.u=r.u,N=N,n.p=n.p)
}

cls=rainbow(length(p))
y.max <- 0.00014
pdf(paste(sep="","results/tag_density_",date.is,".pdf"))
for(i in 1:length(p)){
	e <- names(p)[i]
	pl=p[[i]]
	x <- unique(pl)
	w <- numeric()
	for(u in 1:length(x)){
		ix <- which(pl==x[u])
		w <- c(w,rep(ix,x[u]))
	}
	w.d <- density(w,from=0, to=2*N)
	if(i==1){
		plot(w.d,xlim=c(0,2*N),xlab="bp position relative to TSS [kb]",xaxt="n",main="tag density around TSS",col=cls[i],ylim=c(0,0.00015))
		axis(1,at=c(0,N,2*N),labels=c(-N/1000,0,N/1000))
	} else {
		lines(w.d,col=cls[i])
	}
}
legend("topright", names(p), col=cls,lty=rep("solid",length(p)))
dev.off()
