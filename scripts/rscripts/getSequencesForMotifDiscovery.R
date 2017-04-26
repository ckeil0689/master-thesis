##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Apr 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

library(BSgenome.Mmusculus.UCSC.mm9) # to get sequence around each peak
source("r_scripts/th17/used_for_paper/utilChipScores.R")

path.input.macs <- paste(sep="","input/th17/used_for_paper/rawData/MACS_tab_delim_",date.macs.data.compiled,"/") # the files here not adhere to tab delim format
path.input <- "input/th17/used_for_paper/"
path.input.clusters <- "input/th17/used_for_paper/chipClusterAnalysis/"
path.output.sequences <- "results/sequences_for_motif_analysis/"


cat(sep="","getting ",peak.region.seq.length," bps around each of the ",num.genes.with.highest.pval," top MACS picks (based on pval)\n")
# keepers <- c("P300_Th0_SL1948_SL1947","P300_Th17rorcwt_SL3594_SL3592",
#              "IRF4_Th0_SL1235_SL1234","IRF4_Th17rorcwt_SL2872_SL2876",
#              "BATF_Th0_SL3192_SL3190","BATF_Th17_SL3037_SL3036",
#              "MAF_Th0_SL4424_SL4425","MAF_Th17_SL3032_SL2871",
#              "STAT3_Th0_SL3780_SL3778","STAT3_Th17rorcwt_SL3315_SL3319",
#              "RORC_Th0_SL3779_SL3778","RORC_Th17_SL2870_SL2871",
# 			"FOSL2_Th0_SL6500_SL6499","FOSL2_Th17_SL6498_SL6497")
			
keepers <- c("FOSL2_Th17_SL6498_SL6497")

load(paste(sep="",path.input,"peaks_mat_List.RData"))
load(paste(sep="",path.input.clusters,"tf_clusters_th17_",d.cut,"_",date.cluster.peaks.run,".RData"))
load(paste(sep="",path.input.clusters,"tf_clusters_th0_",d.cut,"_",date.cluster.peaks.run,".RData"))

find.singletons <- function(ll,factor="IRF4"){
  ix.factor <- which( sapply(ll,function(i) i[["tf.to.s"]]) == factor)
  p <- sapply(ll[ix.factor],function(i) i[["peak.ids"]])
  return(p)
}
stop("AM")
## peak mat list with sequences
pmls <- peaks.mat.List[keepers]
## get top peaks based on num.genes.with.highest.pval
for(i in 1:length(pmls)){
  ## resort pmls based on pval
  ix <- sort(pmls[[i]][,"pval"],decreasing=TRUE,index.return=TRUE)$ix
  pmls[[i]] <- pmls[[i]][ix,]
  e <- strsplit(names(pmls)[i],"_")[[1]][1:2]

  ## if irf4 than find best singleton peaks
  ## if(e[1]=="IRF4"){
  ##   cat("for", e[1], "only consider singleton peaks\n")
  ##   if(e[2]=="Th0"){
  ##     peak.ids <- find.singletons(ll=ll.th0,factor="IRF4")
  ##   } else {
  ##      peak.ids <- find.singletons(ll=ll.th17,factor="IRF4")
  ##   }
  ##   ix <- which(pmls[[i]][,"peak_id"] %in% peak.ids)
  ##   pmls[[i]] <- pmls[[i]][ix,]
  ## }

  ## remove redundent peaks from pmls (i.e. remove peaks that hit more than one gene)
  ## and keep no more than num.genes.with.highest.pval peaks
  unique.peaks <- unique(pmls[[i]][,"peak_id"])
  if(length(unique.peaks)>num.genes.with.highest.pval){
    unique.peaks <- unique.peaks[1:num.genes.with.highest.pval]
  }
  ix <- numeric(length=length(unique.peaks))
  for(j in 1:length(unique.peaks)){
    peak <- unique.peaks[j]
    ix[j] <- which(pmls[[i]][,"peak_id"]==peak)[1]
  }
  pmls[[i]] <- pmls[[i]][ix,]
}
## add the chromosome associated with each peak to peaks.mat.List
for(i in 1:length(pmls)){
  cat("adding chromosome number for", names(pmls)[i],"to chip summary table\n")
  macs.file.nm <- paste(sep="",paste(sep="",strsplit(names(pmls)[i],"_")[[1]][3],"_",
                                     strsplit(names(pmls)[i],"_")[[1]][4],"_peaks.xls"))
  x <- read.table(paste(path.input.macs,macs.file.nm,sep=""),header=T)
                                        # add chr column
  chr <- strsplit(as.character(x$chr),"r")
  chrNumSexNamesMod <- sapply(chr,function(i) i[2])
  chrNumSexNamesMod[which(chrNumSexNamesMod=="X")] <- -1 # chrX == -1
  chrNumSexNamesMod[which(chrNumSexNamesMod=="Y")] <- -2 # chrY == -2
  chrNumSexNamesMod <- as.numeric(chrNumSexNamesMod) # make chr into a numeric variable
  names(chrNumSexNamesMod) <- rownames(x) # associate each chrNum with it's respective peak id
                                        # peak mat list is not ordered by peak number so get the order of peaks from peaks mat list
  ix.reordered <- pmls[[i]][,"peak_id"]
  pmls[[i]] <- cbind(pmls[[i]],chrNumSexNamesMod[as.character(ix.reordered)])
  colnames(pmls[[i]])[dim(pmls[[i]])[2]] <- "chrNum_X_-1_Y_-2"
}

## add binding sequence column
## n.p <- 3 # number of processors to use for get.sequences.parallel more procss would make memory footprint too large
for(i in 1:length(pmls)){
  chr <- pmls[[i]][,"chrNum_X_-1_Y_-2"]
  chr[which(chr==-1)] <- "X"
  chr[which(chr==-2)] <- "Y"
  chr <- paste("chr",chr,sep="")
  myseqs <- data.frame(
                       chr=chr,
                       start=pmls[[i]][,"Summit"]-peak.region.seq.length,
                       end=pmls[[i]][,"Summit"]+peak.region.seq.length,
                       strand=rep("+",dim(pmls[[i]])[1])
                       )
  nr <- dim(myseqs)[1]
  unit.skip = floor(nr/n.p)
  ## so we can run in parallel hard coded to four processors
  ix.list = list(seg1=1:unit.skip,
    seg2=(unit.skip+1):(2*unit.skip),
    seg3=(2*unit.skip+1):nr
    )
  x <- mclapply(ix.list, get.sequences.parallel, myseqs,mc.cores=n.p)
  sequences <- c(x[["seg1"]],x[["seg2"]],x[["seg3"]])
  gn_nm <- rownames(pmls[[i]])
  fl.nm <- paste(path.output.sequences,names(pmls)[i],"_",date.is,".fasta",sep="")
  ## write fasta file for results
  cat("writing peak sequences for ", names(pmls[i]),"\n")
  append.to <- FALSE
  for(j in 1:dim(pmls[[i]])[1]){
    if(j > 1){append.to <- TRUE}
    cat(sep="",">",gn_nm[j],"|",pmls[[i]][j,"peak_id"]," ",pmls[[i]][j,"pval"],"\n",sequences[j],"\n",file=fl.nm,append=append.to)
  }
}
