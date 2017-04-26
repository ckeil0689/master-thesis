##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## May 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

library(BSgenome.Mmusculus.UCSC.mm9) # to get sequence around each peak
source("r_scripts/th17/used_for_paper/utilChipScores.R")
source("r_scripts/th17/used_for_paper/cluster_peaks_util.R")

x <- unlist(strsplit(date()," +",perl=TRUE))
date.is <- paste(x[2],x[3],x[5],sep="_")

path.input <- "input/th17/used_for_paper/"
path.input.clusters <- "input/th17/used_for_paper/chipClusterAnalysis/"

load(paste(sep="",path.input.clusters,"tf_clusters_th0.RData"))
load(paste(sep="",path.input.clusters,"tf_clusters_th17.RData"))

# get IRF4 peaks alone
find.mtl.exact <- function(ll,mtl="IRF4"){
  ix.factor <- which( sapply(ll,function(i) paste(sep="",sort(unique(i[["tf.to.s"]])),collapse="_") ) == mtl)
  return(ix.factor)
}

find.mtl.grep <- function(ll,mtl="IRF4"){
  ix.factor <- which( sapply(ll,function(i) length(grep(perl=TRUE,pattern=mtl,x=paste(sep="",sort(unique(i[["tf.to.s"]])),collapse="_")))==1 ) )
  return(ix.factor)
}

th0 <- list()
ix <- find.mtl.exact(ll.th0,mtl="IRF4")
th0[["IRF4"]] <- ll.th0[ix]
ix <- find.mtl.exact(ll.th0,mtl="BATF_IRF4")
th0[["BATF_IRF4"]] <- ll.th0[ix]
ix <- find.mtl.grep(ll.th0,mtl="^BATF_IRF4_.*")
th0[["BATF_IRF4_more"]] <- ll.th0[ix]


th17 <- list()
ix <- find.mtl.exact(ll.th17,mtl="IRF4")
th17[["IRF4"]] <- ll.th17[ix]
ix <- find.mtl.exact(ll.th17,mtl="BATF_IRF4")
th17[["BATF_IRF4"]] <- ll.th17[ix]
ix <- find.mtl.grep(ll.th17,mtl="^BATF_IRF4_.*")
th17[["BATF_IRF4_more"]] <- ll.th17[ix]

## print data
for(i in 1:length(th0)){
  f.nm <- paste(sep="","results/cluster_analysis/tf_clusters_th0_",names(th0)[i],"_",date.is,".xls")
  print.ll.verbose(th0[[i]],f.nm)
}

## print data
for(i in 1:length(th17)){
  f.nm <- paste(sep="","results/cluster_analysis/tf_clusters_th17_",names(th17)[i],"_",date.is,".xls")
  print.ll.verbose(th17[[i]],f.nm)
}







