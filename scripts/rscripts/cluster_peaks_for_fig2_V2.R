##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2012 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# this code calculates clusters of peaks and their association with TFs

library(multicore)
source("r_scripts/th17/used_for_paper/cluster_peaks_util.R")
source("r_scripts/th17/used_for_paper/fig2_util.R")
########## read data ##########
cat("reading data\n")
# set paths
path.input <- "input/th17/used_for_paper/"
path.input.macs <-  paste(sep="","input/th17/used_for_paper/rawData/MACS_tab_delim_",date.macs.data.compiled,"/")


path.input.save.results <- "input/th17/used_for_paper/fig2/"
path.output <- "results/fig2/"

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
keepers.th17 <- c("IRF4_Th17rorcwt_SL2872_SL2876",
             "BATF_Th17_SL3037_SL3036",
             "MAF_Th17_SL3032_SL2871",
             "STAT3_Th17rorcwt_SL3315_SL3319",
             "RORC_Th17_SL2870_SL2871",
             "P300_Th17rorcwt_SL3594_SL3592",
             "FOSL2_Th17_SL6498_SL6497")
keepers.th0 <- c("IRF4_Th0_SL1235_SL1234",
             "BATF_Th0_SL3192_SL3190",
             "MAF_Th0_SL4424_SL4425",
             "STAT3_Th0_SL3780_SL3778",
             "RORC_Th0_SL3779_SL3778",
             "P300_Th0_SL1948_SL1947",
             "FOSL2_Th0_SL6500_SL6499")
tfs.th17 <- sapply(strsplit(keepers.th17,"_"),function(i) i[[1]])
tfs.th0 <- sapply(strsplit(keepers.th0,"_"),function(i) i[[1]])

tmp.th17 <- character()
tmp.th0 <- character()

for(i in 1:length(tfs.to.cluster.t)){
	tmp.th17 <- c(tmp.th17,keepers.th17[grep(tfs.to.cluster.t[i],tfs.th17,ignore.case=T)])
	tmp.th0 <- c(tmp.th0,keepers.th0[grep(tfs.to.cluster.t[i],tfs.th0,ignore.case=T)])
}
keepers.th17 <- tmp.th17 
keepers.th0 <- tmp.th0

tfs <- sapply(keepers.th17,function(i) strsplit(i,"_")[[1]][1])
## for each peak/TF we will map coressponding gene trgts if available, thus we will need peaks.mat.list
load(paste(sep="",path.input,"peaks_mat_List.RData"))
pml <- peaks.mat.List # short hand
## take macs.list and split each table into it's respective chr

chr <- c(paste(sep="","chr",1:19),paste(sep="","chr",c("X","Y")))
x <- split.macs.list.to.chrmosomes(macs.list=macs.list,expts=keepers.th17,chr=chr,pml=pml,n.p=n.p,
	mtbs.tss.dist=mtbs.tss.dist,match.targets=match.targets)
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
ll.th17 <- create.cluster.list(ruler.th17)
# cat("writing cluster for th17\n")
# f.nm <- paste(sep="",path.input.save.results,"tf_clusters_",paste(tfs,collapse="_"),"_th17_d.tfs_",d.cut,"_",date.is,".xls")
# # print.ll.verbose(ll=ll.th17,f.nm=f.nm,  tfs=tfs)
f.nm <- paste(sep="",path.input.save.results,"tf_clusters_",paste(tfs,collapse="_"),"_th17_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.is,".RData")
save(ll.th17,file=f.nm)

################################
## done with th17 now go to th0
################################
## take macs.list and split each table into it's respective chr
x <- split.macs.list.to.chrmosomes(macs.list=macs.list,expts=keepers.th0,chr=chr,pml=pml,n.p=n.p,mtbs.tss.dist=mtbs.tss.dist,match.targets=match.targets)
########## fill up ruler for each chromosome with peaks of all tfs ##########
cat("th0-creating combined-peaks (for all tfs) ruler for each chromosome\n")
ruler.th0 <- make.ruler(chr,macs.list.per.chrom=x)

cat("th0-add cluster/loci membership to the ruler.th0\n")
ruler.th0 <- assign.clusters.ids.to.peaks(ruler.th0,d.cut)
########## for each loci (cluster) get a bunch of params ###########
cat("th0-create the loci list ll.th0\n")
ll.th0 <- create.cluster.list(ruler.th0)
# cat("writing cluster for th0\n")
# f.nm <- paste(sep="",path.input.save.results,"tf_clusters_",paste(tfs,collapse="_"),"_th0_d.tfs_",d.cut,"_",date.is,".xls")
# # print.ll.verbose(ll=ll.th0,f.nm=f.nm,tfs=tfs)
f.nm <-  paste(sep="",path.input.save.results,"tf_clusters_",paste(tfs,collapse="_"),"_th0_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.is,".RData")
save(ll.th0,file=f.nm)




