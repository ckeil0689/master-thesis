##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Apr 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# this code calculates clusters of peaks and their association with TFs

library(multicore)
source("r_scripts/th17/used_for_paper/cluster_peaks_util.R")
########## read data ##########
cat("reading data\n")
# set paths
path.input <- "input/th17/used_for_paper/"
path.input.save.results <- "input/th17/used_for_paper/chipClusterAnalysis/"
path.output <- "results/cluster_analysis/"
path.output.sungear.files <- "results/sungear/"
path.input.macs <-  paste(sep="","input/th17/used_for_paper/rawData/MACS_tab_delim_",date.macs.data.compiled,"/")
path.input.refseq <- "input/th17/used_for_paper/rawData/"

# get MACS file names
MACS.files <- list.files(path.input.macs)
# get mapping of MACS file names to meaningful (Th17 Th0..) names
macs.id.to.expt.name <- read.delim(paste(sep="",path.input,"rawData/macs_id_to_expt_map.txt"))

# read MACS output files
macs.list <- list()
n <- length(MACS.files)
# AM some experiments are not relevant for th17 project! remove them Oct 5
mistake.id.full <- c("SL2874_SL2875","SL2874_SL3315","SL2875_SL3316","SL3315_SL3316","SL1172_SL1171","SL6140_SL6142","SL6141_SL6142")
id.full <- character(length=n)
for(i in 1:n){
  id.full[i] <- gsub(".*(SL\\d+_SL\\d+).*","\\1",MACS.files[i],perl=T) # e.g. "SL1040_SL972" (expt.id_ctrl.id)
  if(!any(mistake.id.full == id.full[i])){
    expt.nm <- as.character(macs.id.to.expt.name$expt[which(macs.id.to.expt.name$id == id.full[i])]) # map id to expt name
    expt.nm <- gsub("(.*_SL\\d+)_\\d$.*","\\1",expt.nm,perl=T)
    macs.list[[ expt.nm ]] <- read.delim(file=paste(sep="",path.input.macs,MACS.files[i]))
    macs.list[[ expt.nm ]][,"chr"] <- as.character(macs.list[[ expt.nm ]][,"chr"])
  } else {
    cat("found a wrong MACS experiment",id.full[i],"\n")
  }
}

########## preprocess data ##########
cat("preprocessing data\n")
## focus on key experiments
keepers <- c("IRF4_Th0_SL1235_SL1234","IRF4_Th1_SL3501_SL3500","IRF4_Th2_SL3194_SL3193",
			"IRF4_Th17rorcwt_SL2872_SL2876","IRF4_iTreg_SL3196_SL3195")
## for each peak/TF we will map coressponding gene trgts if available, thus we will need peaks.mat.list
load(paste(sep="",path.input,"peaks_mat_List_pioneer.RData"))
pml <- peaks.mat.List # short hand
rm(peaks.mat.List)
## take macs.list and split each table into it's respective chr
## also add trgts for each peak as "trgt1_trgt2_trgt3_..."
# This step takes time so if you have done it before, you can just load the data
chr <- c(paste(sep="","chr",1:19),paste(sep="","chr",c("X","Y")))

if(load.macs.per.chrom.RData==FALSE){
  x <- split.macs.list.to.chrmosomes(macs.list=macs.list,expts=keepers,chr=chr,pml=pml,n.p=n.p,
	mtbs.tss.dist=mtbs.tss.dist,match.targets=match.targets)
  # save macs.per.chrom
  save(x,file=paste(sep="",path.input.save.results,"macs.per.chrom.th17.pioneer.RData"))
} else {
  # the name of macs.per.chrom.RData variable is x
  load(paste(sep="",path.input,"chipClusterAnalysis/macs.per.chrom.th17.RData"))
}
########## fill up ruler for each chromosome with peaks of all tfs ##########
cat("th17-creating combined-peaks (for all tfs) ruler for each chromosome\n")
# ruler is a line per chromosome with the locations of all tf binding sites (sorted from start of chromosome to end)
# also match these summit locations with corresponding:
# pvals, tfs, peak start and peak end
ruler.th17 <- make.ruler(chr,macs.list.per.chrom=x)

cat("th17-add cluster/loci membership to the ruler\n")
ruler.th17 <- assign.clusters.ids.to.peaks(ruler.th17,d.cut)
########## for each loci (cluster) get a bunch of params ###########
cat("th17-create the loci list ll.th17\n")
ll.thx<- create.cluster.list(ruler.th17)
cond.names <- sapply(keepers,function(i) paste(sep="",strsplit(i,"_")[[1]][1:2],collapse="_"))
cond.names[grep("Th17",cond.names,ignore.case=T)] <- "IRF4_Th17"
cat("writing cluster for th17\n")
f.nm <- paste(sep="",path.output,"tf_clusters_irf4_thx_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_pioneer_",date.is,".xls")
print.ll.verbose.pioneer(ll=ll.thx,f.nm=f.nm,tfs=cond.names)
save(ll.thx,file=paste(sep="",path.input.save.results,"tf_clusters_irf4_thx_",d.cut,"_d.tss_",mtbs.tss.dist,"_pioneer_",date.is,".RData"))

















