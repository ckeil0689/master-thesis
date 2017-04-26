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
keepers.th17 <- c("IRF4_Th17rorcwt_SL2872_SL2876",
             "BATF_Th17_SL3037_SL3036",
             "MAF_Th17_SL3032_SL2871",
             "STAT3_Th17rorcwt_SL3315_SL3319",
             "RORC_Th17_SL2870_SL2871")
             # "P300_Th17rorcwt_SL3594_SL3592",
             # "FOSL2_Th17_SL6498_SL6497")
keepers.th0 <- c("IRF4_Th0_SL1235_SL1234",
             "BATF_Th0_SL3192_SL3190",
             "MAF_Th0_SL4424_SL4425",
             "STAT3_Th0_SL3780_SL3778",
             "RORC_Th0_SL3779_SL3778")
             # "P300_Th0_SL1948_SL1947",
             # "FOSL2_Th0_SL6500_SL6499")
tfs <- sapply(keepers.th17,function(i) strsplit(i,"_")[[1]][1])
## for each peak/TF we will map coressponding gene trgts if available, thus we will need peaks.mat.list
load(paste(sep="",path.input,"peaks_mat_List.RData"))
pml <- peaks.mat.List # short hand
rm(peaks.mat.List)
## take macs.list and split each table into it's respective chr
## also add trgts for each peak as "trgt1_trgt2_trgt3_..."
# This step takes time so if you have done it before, you can just load the data
chr <- c(paste(sep="","chr",1:19),paste(sep="","chr",c("X","Y")))

if(load.macs.per.chrom.RData==FALSE){
  x <- split.macs.list.to.chrmosomes(macs.list=macs.list,expts=keepers.th17,chr=chr,pml=pml,n.p=n.p,
	mtbs.tss.dist=mtbs.tss.dist,match.targets=match.targets)
  # save macs.per.chrom
  save(x,file=paste(sep="",path.input.save.results,"macs.per.chrom.th17.RData"))
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
ll.th17 <- create.cluster.list(ruler.th17)
cat("writing cluster for th17\n")
f.nm <- paste(sep="",path.output,"tf_clusters_th17_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.is,".xls")
print.ll.verbose(ll=ll.th17,f.nm=f.nm,  tfs=tfs)
save(ll.th17,file=paste(sep="",path.input.save.results,"tf_clusters_th17_",d.cut,"_",date.is,".RData"))
################################
## done with th17 now go to th0
################################
## take macs.list and split each table into it's respective chr
## also add trgts for each peak as "trgt1_trgt2_trgt3_..."
# This step takes time so if you have done it before, you can just load the data
if(load.macs.per.chrom.RData==FALSE){
  x <- split.macs.list.to.chrmosomes(macs.list=macs.list,expts=keepers.th0,chr=chr,pml=pml,n.p=n.p,mtbs.tss.dist=mtbs.tss.dist,match.targets=match.targets)
  # save macs.per.chrom
  save(x,file=paste(sep="",path.input.save.results,"macs.per.chrom.th0.RData"))
} else {
  # the name of macs.per.chrom.RData variable is x
  load(paste(sep="",path.input,"chipClusterAnalysis/macs.per.chrom.th0.RData"))
}
########## fill up ruler for each chromosome with peaks of all tfs ##########
cat("th0-creating combined-peaks (for all tfs) ruler for each chromosome\n")
ruler.th0 <- make.ruler(chr,macs.list.per.chrom=x)

cat("th0-add cluster/loci membership to the ruler.th0\n")
ruler.th0 <- assign.clusters.ids.to.peaks(ruler.th0,d.cut)
########## for each loci (cluster) get a bunch of params ###########
cat("th0-create the loci list ll.th17\n")
ll.th0 <- create.cluster.list(ruler.th0)
cat("writing cluster for th0\n")
f.nm <- paste(sep="",path.output,"tf_clusters_th0_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.is,".xls")
print.ll.verbose(ll=ll.th0,f.nm=f.nm,tfs=tfs)
save(ll.th0,file=paste(sep="",path.input.save.results,"tf_clusters_th0_",d.cut,"_",date.is,".RData"))


## add to ll.th0 and ll.th17 info about histone marks
if(add.histone.scores==T){
  source("r_scripts/th17/used_for_paper/cluster_peaks_with_histone_marks.R")

######### QC ##########
## sam <- as.matrix(read.table(paste(path.input,"samTh17VsTh0Zscores.xls",sep=""),header=T,sep="\t"))
## get list of tfs memberships in clusters from ll.th17 e.g. a clust with BATF and IRF4 will be "BATF_IRF4"
cat("calculating and ploting quality controls\n")
# th17
clust.types.th17 <- sapply(ll.th17,function(i) paste( sort(unique(i[["tf.to.s"]])),collapse="_") )
c.table.th17 <- table(clust.types.th17)
pvals.th17 <- list()
pvals.th17.per.tf <- list()
for(i in 1:length(tfs)){
	pvals.th17.per.tf[[tfs[i]]] <- list()
}
h.pvals.th17 <- list()
h.fc.th17 <- list()
for(i in 1:length(c.table.th17)){
  type <- names(c.table.th17)[i]
  ix <- which(clust.types.th17==type)
  pvals.th17[[type]] <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["pval"]] )))
  for(j in 1:length(tfs)){
	tf <- tfs[j]
  	pvals.th17.per.tf[[tf]][[type]] <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["pval"]][which(i$tf.to.s==tf)] )))
  }
  ## add histone marks pvals
  h.pvals.th17[[type]] <- list()
  h.pvals.th17[[type]]$P300 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.pvals"]]["SL1041_P300_Th17"] )))
  h.pvals.th17[[type]]$P300.rorcwt <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.pvals"]]["SL3594_P300_Th17rorcwt"] )))
  h.pvals.th17[[type]]$P300.rorcko <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.pvals"]]["SL3595_P300_Th17rorcko"] )))
  h.pvals.th17[[type]]$POL2 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.pvals"]]["SL1946_POL2_Th17"] )))
  h.pvals.th17[[type]]$CTCF <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.pvals"]]["SL3034_CTCF_Th17"] )))
  h.pvals.th17[[type]]$H3K4me1 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.pvals"]]["SL4030_H3K4me1_Th17rorcwt"] ))) ## activation
  h.pvals.th17[[type]]$H3K4me3 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.pvals"]]["SL4034_H3K4me3_Th17rorcwt"] ))) ## tss /activation
  h.pvals.th17[[type]]$H3K4me2 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.pvals"]]["SL2993_H3K4me2_Th17rorcwt"] ))) ##?
  h.pvals.th17[[type]]$H3AcK914 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.pvals"]]["SL4038_H3AcK914_Th17rorcwt"] )))
  h.pvals.th17[[type]]$H3K27me3 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.pvals"]]["SL2997_H3K27me3_Th17rorcwt"] ))) ## repressive

  ## add histone marks fc
  h.fc.th17[[type]] <- list()
  h.fc.th17[[type]]$P300 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.fc"]]["SL1041_P300_Th17"] )))
  h.fc.th17[[type]]$P300.rorcwt <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.fc"]]["SL3594_P300_Th17rorcwt"] )))
  h.fc.th17[[type]]$P300.rorcko <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.fc"]]["SL3595_P300_Th17rorcko"] )))
  h.fc.th17[[type]]$POL2 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.fc"]]["SL1946_POL2_Th17"] )))
  h.fc.th17[[type]]$CTCF <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.fc"]]["SL3034_CTCF_Th17"] )))
  h.fc.th17[[type]]$H3K4me1 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.fc"]]["SL4030_H3K4me1_Th17rorcwt"] ))) ## activation
  h.fc.th17[[type]]$H3K4me3 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.fc"]]["SL4034_H3K4me3_Th17rorcwt"] ))) ## tss /activation
  h.fc.th17[[type]]$H3K4me2 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.fc"]]["SL2993_H3K4me2_Th17rorcwt"] ))) ##?
  h.fc.th17[[type]]$H3AcK914 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.fc"]]["SL4038_H3AcK914_Th17rorcwt"] )))
  h.fc.th17[[type]]$H3K27me3 <- as.numeric(unlist(sapply(ll.th17[ix],function(i) i[["h.fc"]]["SL2997_H3K27me3_Th17rorcwt"] ))) ## repressive
}

# th0
clust.types.th0 <- sapply(ll.th0,function(i) paste( sort(unique(i[["tf.to.s"]])),collapse="_") )
c.table.th0 <- table(clust.types.th0)
pvals.th0 <- list()
pvals.th0.per.tf <- list()
for(i in 1:length(tfs)){
	pvals.th0.per.tf[[tfs[i]]] <- list()
}
h.pvals.th0 <- list()
h.fc.th0 <- list()
for(i in 1:length(c.table.th0)){
  type <- names(c.table.th0)[i]
  ix <- which(clust.types.th0==type)
  pvals.th0[[type]] <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["pval"]] )))
  pvals.th0[[type]] <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["pval"]] )))
  for(j in 1:length(tfs)){
	tf <- tfs[j]
  	pvals.th0.per.tf[[tf]][[type]] <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["pval"]][which(i$tf.to.s==tf)] )))
  }
  ## add histone marks pvals
  h.pvals.th0[[type]] <- list()
  h.pvals.th0[[type]]$P300 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.pvals"]]["SL1948_P300_Th0"] )))
  h.pvals.th0[[type]]$POL2 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.pvals"]]["SL1952_POL2_Th0"] )))
  h.pvals.th0[[type]]$CTCF <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.pvals"]]["SL3191_CTCF_Th0"] )))
  h.pvals.th0[[type]]$P300.rorcwt <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.pvals"]]["SL3594_P300_Th17rorcwt"] )))
  h.pvals.th0[[type]]$P300.rorcko <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.pvals"]]["SL3595_P300_Th17rorcko"] )))
  h.pvals.th0[[type]]$H3K4me1 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.pvals"]]["SL4031_H3K4me1_Th17rorcko"] ))) ## activation
  h.pvals.th0[[type]]$H3K4me3 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.pvals"]]["SL4035_H3K4me3_Th17rorcko"] ))) ## tss /activation
  h.pvals.th0[[type]]$H3K4me2 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.pvals"]]["SL2994_H3K4me2_Th17rorcko"] ))) ##?
  h.pvals.th0[[type]]$H3AcK914 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.pvals"]]["SL4039_H3AcK914_Th17rorcko"] )))
  h.pvals.th0[[type]]$H3K27me3 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.pvals"]]["SL2998_H3K27me3_Th17rorcko"] ))) ## repressive
  ## add histone marks fold.changes  
  h.fc.th0[[type]] <- list()
  h.fc.th0[[type]]$P300 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.fc"]]["SL1948_P300_Th0"] )))
  h.fc.th0[[type]]$POL2 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.fc"]]["SL1952_POL2_Th0"] )))
  h.fc.th0[[type]]$CTCF <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.fc"]]["SL3191_CTCF_Th0"] )))
  h.fc.th0[[type]]$P300.rorcwt <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.fc"]]["SL3594_P300_Th17rorcwt"] )))
  h.fc.th0[[type]]$P300.rorcko <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.fc"]]["SL3595_P300_Th17rorcko"] )))
  h.fc.th0[[type]]$H3K4me1 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.fc"]]["SL4031_H3K4me1_Th17rorcko"] ))) ## activation
  h.fc.th0[[type]]$H3K4me3 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.fc"]]["SL4035_H3K4me3_Th17rorcko"] ))) ## tss /activation
  h.fc.th0[[type]]$H3K4me2 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.fc"]]["SL2994_H3K4me2_Th17rorcko"] ))) ##?
  h.fc.th0[[type]]$H3AcK914 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.fc"]]["SL4039_H3AcK914_Th17rorcko"] )))
  h.fc.th0[[type]]$H3K27me3 <- as.numeric(unlist(sapply(ll.th0[ix],function(i) i[["h.fc"]]["SL2998_H3K27me3_Th17rorcko"] ))) ## repressive
}

## get distances from start site for all peaks that are in [10k downstream, 10k upstream] of TSS
## th17
dtss.th17 <- list()
for(i in 1:length(c.table.th17)){
  type <- names(c.table.th17)[i]
  ix <- which(clust.types.th17==type)
  x <- sapply(ll.th17[ix],function(i) i[["d.tss"]] )
  x2 <- lapply(x,function(i) as.numeric(strsplit(i,"_")[[1]]))
  x3 <- sapply(x2, function(i) if(length(i)>0) {if(abs(min(i))<10000){min(i)}else{NA}}else{NA})
  dtss.th17[[type]] <- x3[which(!is.na(x3))]
}
## get distances from start site for all peaks that are in [10k downstream, 10k upstream] of TSS
## th0
dtss.th0 <- list()
for(i in 1:length(c.table.th0)){
  type <- names(c.table.th0)[i]
  ix <- which(clust.types.th0==type)
  x <- sapply(ll.th0[ix],function(i) i[["d.tss"]] )
  x2 <- lapply(x,function(i) as.numeric(strsplit(i,"_")[[1]]))
  x3 <- sapply(x2, function(i) if(length(i)>0) {if(abs(min(i))<10000){min(i)}else{NA}}else{NA})
  dtss.th0[[type]] <- x3[which(!is.na(x3))]
}

## give clusters that are only recovered in th17 a 0 value in th0
## i want both c.table vectors to be of same length
c.table.th0.tmp <- numeric(length=length(c.table.th17))
names(c.table.th0.tmp) <- names(c.table.th17)
for(i in 1:length(c.table.th17)){
  e <- names(c.table.th17)[i]
  if(! e %in% names(c.table.th0) ){
    c.table.th0.tmp[e] <- 0
  } else {
    c.table.th0.tmp[e] <- c.table.th0[e]
  }
}
c.table.th0 <- c.table.th0.tmp


######## plot sam values of target genes of each clustesr with more than 500 counts ##########
## for each cluster going from "big" cluster (5 tfs) to smaller clusters
######## plot occurances of clustesr with more than 500 counts ##########
cut.th17 <- 500 ##(only plot barplots and boxplots for clusters with more than cut.th0 occurrances)
cut.th0 <- 500 ##(for th0 plot barplots as for th17 but only plot boxplots for clusters with more than cut.th0 occurrances)
ylim.pval <- c(0,2000)
ylim.pval.histones=c(0,3100)
pdf(file=paste(sep="",path.output,"tf_cluster_peaks_",d.cut,"_",date.is,".pdf"))
par(mfrow=c(2,2))

plot.qc.1(c.table.th17,c.table.th0,
          cut.th17=cut.th17,
          cut.th0=cut.th0,
          ylim.pval=ylim.pval,
          ylim.pval.histones=ylim.pval.histones,
          plot.histone.marks=add.histone.scores)
dev.off()

######## plot all occurances of clustesr ##########
# plot th17
pdf(file=paste(sep="",path.output,"tf_cluster_peaks_all_",d.cut,"_",date.is,".pdf"))
par(mfrow=c(2,1))
cut.th17 <- 0
plot.qc.2(c.table.th17,c.table.th0,cut.th17=cut.th17, ylim.pval=ylim.pval)
dev.off()


## c.table.perm.th0 <- matrix(0,nr=n.b,nc=length(c.table.th17))
## c.table.perm.th17 <- matrix(0,nr=n.b,nc=length(c.table.th17))
## colnames(c.table.perm.th0) <- colnames(c.table.perm.th17) <- names(c.table.th17)

if(calc.simulations==TRUE){
## get statistics on clusrters by permutation analysis
## for each cluster type from th17 perform permutatois to see significance
## store results in

# get min/max genome position with peak on each chromosome
min.chr.loc <- numeric(length=length(chr))
max.chr.loc <- numeric(length=length(chr))
len.th17.per.chr <- numeric(length=length(chr))
len.th0.per.chr <- numeric(length=length(chr))
tf.th17.per.chr <- list()
tf.th0.per.chr <- list()
pval.th17.per.chr <- list()
pval.th0.per.chr <- list()
names(len.th0.per.chr) <- names(len.th17.per.chr) <- names(min.chr.loc) <- names(max.chr.loc) <- chr
for(i in 1:length(chr)){
  len.th17.per.chr[i] <- length(ruler.th17[[i]]$start)
  tf.th17.per.chr[[i]] <- ruler.th17[[i]]$tf.to.s
  pval.th17.per.chr[[i]] <- ruler.th17[[i]]$pval
  tf.th0.per.chr[[i]] <- ruler.th0[[i]]$tf.to.s
  len.th0.per.chr[i] <- length(ruler.th0[[i]]$start)
  pval.th0.per.chr[[i]] <- ruler.th0[[i]]$pval
  min.chr.loc[i] <- ruler.th17[[i]]$start[1]
  max.chr.loc[i] <- ruler.th17[[i]]$start[len.th17.per.chr[i]]
}

  nms.th17.table <- names(c.table.th17)

  x <- run.perms(n.b=n.b,
                 len.th17.per.chr=len.th17.per.chr,
                 tf.th17.per.chr=tf.th17.per.chr,
                 pval.th17.per.chr=pval.th17.per.chr,
                 tf.th0.per.chr=tf.th0.per.chr,
                 len.th0.per.chr=len.th0.per.chr,
                 pval.th0.per.chr=pval.th0.per.chr,
                 min.chr.loc=min.chr.loc,
                 max.chr.loc=max.chr.loc,
                 nms.th17.table=nms.th17.table,
                 n.p=n.p)
  save(x,file=paste(sep="",path.output,"perm_runs_chip_clusters_nb_",n.b,".RData"))
  c.table.perm.th17 <- x$th17.perm.table
  c.table.perm.th0 <- x$th0.perm.table
  p.table.perm.th17 <- x$th17.pval.table
  p.table.perm.th0 <- x$th0.pval.table
  ## get means counts
  m.th17 <- apply(c.table.perm.th17,2,mean)
  m.th0 <-  apply(c.table.perm.th0,2,mean)

  ## get median counts
  median.th17.counts <- apply(c.table.perm.th17,2,median)+1
  median.th0.counts <-  apply(c.table.perm.th0,2,median)+1

  ## get sd counts
  sd.th17.counts <- apply(c.table.perm.th17,2,sd)+1
  sd.th0.counts <- apply(c.table.perm.th0,2,sd)+1

  ## get zscores counts
  zscores.th17.counts <- (c.table.th17-m.th17)/sd.th17.counts
  zscores.th0.counts <- (c.table.th0-m.th0)/sd.th0.counts

  ## get median pval
  median.pval.th17 <- apply(p.table.perm.th17,2,median)+1
  median.pval.th0 <- apply(p.table.perm.th0,2,median)+1

  ## get sd counts
  sd.pval.th17 <- apply(p.table.perm.th17,2,sd)+1
  sd.th0.counts <- apply(p.table.perm.th0,2,sd)+1

  ## get zscores counts
  zscores.th17.pvals <- (p.table.perm.th17-median.pval.th17)/sd.pval.th17
  zscores.th0.pvals <- (p.table.perm.th17-median.pval.th0)/sd.th0.counts 
}
}



















