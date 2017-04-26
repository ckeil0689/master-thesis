##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Dec 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# create tag densities around TSS for ChIP vs. ctrl to see how specific our ChIPs are
library(multicore)
source("r_scripts/th17/used_for_paper/cluster_peaks_util.R")
source("r_scripts/th17/used_for_paper/fig2_util.R")

# set paths
path.sam.big <- "/Volumes/Drobo/aviv/th17/sam_files/"
path.sam.small <- "/Volumes/Drobo/aviv/th17/sam_files_small/"
path.macs.tab.delim <- paste(sep="","input/th17/used_for_paper/rawData/MACS_tab_delim_",date.macs.processed,"/")
path.input.refseq <- "input/th17/used_for_paper/rawData/"
path.sam.files <- "input/th17/used_for_paper/rawData/sam_files_small_Nov_7_2011/"
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
# get refseq annotations to find all TSS locations
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
# which sam seq files to compare with densities around TSS
sam.treatment <- c("SL3031.sam","SL3032.sam","SL4424.sam","SL3031.sam","SL6628.sam")
names(sam.treatment) <- c("Maf_Th17_48hr","Maf_Th17_48hr","Maf_Th0_48hr","Maf_Th17_48hr","batf_wt_irf_chip")
sam.ctrl <- c("SL2871.sam","SL2871.sam","SL4425.sam","SL4424.sam","SL6629.sam")
names(sam.ctrl) <- c("ctrl_Th17_48hr","ctrl_Th17_48hr","ctrl_Th0_48hr","Maf_Th0_48hr","batf_ko_irf_chip")
# all.files <- unique(c(sam.treatment,sam.ctrl))
names.comparisons <- c("maf_th17_rep1_vs_ctrl","maf_th17_rep2_vs_ctrl","maf_th0_rep1_vs_ctrl",
					   "maf_th17_rep1_vs_th0_rep1","batf_wt_irf_chip_vs_batf_ko_irf_chip")
#~~~~~~~~~~~~~~~~~~~~~~~~3~~~~~~~~~~~~~~~~~~~~~~~~~~~~#	
for(v in 1:length(sam.treatment)){
	# read data
	sl.num.trtmnt <- strsplit(sam.treatment[v],"\\.")[[1]][1]
	sl.num.ctrl <- strsplit(sam.ctrl[v],"\\.")[[1]][1]
	e.nms <- c(names(sam.treatment[v]),names(sam.ctrl[v]))
	e.nms <- paste(sep="_",e.nms,c(sl.num.trtmnt,sl.num.ctrl))
	names(e.nms) <- c("trtmnt","ctrl")
	sam.list <- list()
	
	f.nm.trtmnt <- paste(sep="",path.sam.files,sl.num.trtmnt,".sam")
	cat("reading reads for sl number",sl.num.trtmnt,"\n")
	n.lines <- as.numeric(strsplit(system(paste(sep="","wc -l ",f.nm.trtmnt),intern=T),split=" ")[[1]][2])
	sam.list[[e.nms["trtmnt"]]] <- scan(file=f.nm.trtmnt,what <- list(chr=character(),read_5_prime_bp_location=numeric()),sep="\t",nlines=n.lines)
	
	f.nm.ctrl <- paste(sep="",path.sam.files,sl.num.ctrl,".sam")
	cat("reading reads for sl number",sl.num.ctrl,"\n")
	n.lines <- as.numeric(strsplit(system(paste(sep="","wc -l ",f.nm.ctrl),intern=T),split=" ")[[1]][2])
	sam.list[[e.nms["ctrl"]]] <- scan(file=f.nm.ctrl,what <- list(chr=character(),read_5_prime_bp_location=numeric()),sep="\t",nlines=n.lines)

	# pile up data per chr
	reads.list.per.chrom <- list()
	chr <- c(paste(sep="","chr",1:19),paste(sep="","chr",c("X","Y")))
	reads.list.per.chrom[[e.nms["trtmnt"]]] <- split.matrix.to.chrmosomes(sam.list[[e.nms["trtmnt"]]],chr=chr)
	reads.list.per.chrom[[e.nms["ctrl"]]] <- split.matrix.to.chrmosomes(sam.list[[e.nms["ctrl"]]],chr=chr)
	# rm(sam.list)
	gc()
	x <- list()
	for(i in 1:length(e.nms)){
		e <- e.nms[i]
		x[[e]] <- list()
		for (j in 1:length(chr)){
			r <- chr[j]
			x[[e]][[r]] <- list()
			x[[e]][[r]][[r]] <- reads.list.per.chrom[[e]][[r]]
		}
	}
	
	# plot density around TSS for both trtmnt and ctrl
	# create read-pileup around TSS for each experiments
	p <- list()
	N=10*10^3
	cat("create read-pileups for:", e.nms,"\n")
	for(i in 1:length(e.nms)){
		e <- e.nms[i]
		cat("working on expt:", e,"\n")
		p[[e]] <- get.pileup.per.chr.parallel(x.i=x[[e]],r.u=r.u,N=N,n.p=n.p)
	}
	
	cls=rainbow(length(p))
	pdf(paste(sep="","results/tag_density_",names.comparisons[v],"_",date.is,".pdf"))
	for(i in 1:length(p)){
		e <- names(p)[i]
		pl=p[[i]]
		w <- numeric()
		x <- unique(c(pl[[1]],pl[[2]]))
		for(u in 1:length(x)){
			ix <- which(pl[[1]]==x[u])
			w <- c(w,rep(ix,x[u]))
		}
		w.d <- density(w,from=0, to=2*N)
		y.max <- 0.00014
		if(i==1){
			plot(w.d,xlim=c(0,2*N),xlab="bp position relative to TSS [kb]",xaxt="n",main="tag density around TSS",col=cls[i],ylim=c(0,y.max))
			axis(1,at=c(0,N,2*N),labels=c(-N/1000,0,N/1000))
		} else {
			lines(w.d,col=cls[i])
		}
	}
	legend("topright", names(p), col=cls,lty=rep("solid",length(p)))
	dev.off()
	
}






