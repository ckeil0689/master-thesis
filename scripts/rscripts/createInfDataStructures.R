##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## January 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# set paths
path.input <- "input/th17/used_for_paper/"
# source required functions
source("r_scripts/th17/used_for_paper/util.R")

file.nm.map <- paste(sep="",path.input,"combinedAnalysis/",combine.file.map)
map <- read.table(file.nm.map,sep="\t",colClasses = "character",header=T)
additional.regulators <- toupper(map$tf)

# create path for output
add.str.1 <- paste(sep="","input/th17/used_for_paper/infDataStructures/",date.is,"/")
system(paste(sep="","mkdir ",add.str.1))
if(GLOBAL[["type"]]=="rnaseq"){
	add.str.2 <- paste(sep="",add.str.1,"rnaseq/")
} else {
	add.str.2 <- paste(sep="",add.str.1,"immgen/")
}
system(paste(sep="","mkdir ",add.str.2))
path.output <- paste(add.str.2,"z_abs_cut_",GLOBAL[["z.abs.cut"]] ,"/",sep="")
system(paste(sep="","mkdir ",path.output))


if(GLOBAL[["type"]]=="rnaseq"){
  load(paste(path.input,"ranseqDatasetNoQuartNorm_",date.cufflinks.data.compiled,".RData",sep=""))
  # put time series at front
  ts.nms <- c("SL2654_th0_ts_1hr#1","SL2655_th0_ts_3hr#1","SL2656_th0_ts_6hr#1","SL2778_th0_ts_16hr#1","SL2673_th0_ts_48hr#1",
              "SL1851_th17_ts_1hr#1","SL1852_th17_ts_3hr#1","SL1853_th17_ts_6hr#1",
              "SL1854_th17_ts_9hr#1","SL1855_th17_ts_12hr#1","SL1856_th17_ts_16hr#1",
              "SL1857_th17_ts_24hr#1","SL1858_th17_ts_48hr#1")
  # steady state conditions indices
  ix.ss <- which(!colnames(rnaseq.complete.matrix) %in% ts.nms)
  # change rnaseq long dataset name to just dataset
  # put time series experiments first
  dataset <- cbind(rnaseq.complete.matrix[,ts.nms],rnaseq.complete.matrix[,ix.ss])
  rownames(dataset) = toupper(rownames(dataset))
  rm(rnaseq.complete.matrix)
} else {
  dataset <- as.matrix(read.table(paste(path.input,"immgenRatiosTable.txt",sep="")))
  rownames(dataset) = toupper(rownames(dataset))
}

##################################################
# get rpkm values
f.nm.rpkm <- paste(path.input,"ranseqDatasetNoQuartNorm_",date.cufflinks.data.compiled,".RData",sep="") # which experiment to compare

## read rpkm data
load(f.nm.rpkm)
d <- rnaseq.complete.matrix
expt.names <-colnames(d)
th17.wt.ix <- grep("th17_ss_48hr_.{1,6}_wt",colnames(d),value=T,perl=T)
th17.wt.ix <- c(th17.wt.ix,"SL1858_th17_ts_48hr#1")
th17.mean <- apply(d[,th17.wt.ix],1,mean)
th0.wt.ix <- grep("th0_ss_48hr_.{1,6}_wt",colnames(d),value=T,perl=T)
th0.wt.ix <- c(th0.wt.ix,"SL2673_th0_ts_48hr#1")
th0.mean <- apply(d[,th0.wt.ix],1,mean)
rpkm.mean <- pmax(th17.mean,th0.mean)
# 
# # th17
# m <- numeric(length=nrow(scores.ko))
# names(m) <- gns
# ix <- which(gns %in% names(th17.mean))
# m[gns[ix]] <- th17.mean[gns[ix]]
# th17.rpkm <- m
# # th0
# m <- numeric(length=nrow(scores.ko))
# names(m) <- gns
# ix <- which(gns %in% names(th0.mean))
# m[gns[ix]] <- th0.mean[gns[ix]]
# th0.rpkm <- m
##################################################
##################################################
# get SAM differential expression z scores
x <- as.matrix(read.table(paste(path.input,"samTh17VsTh0Zscores.xls",sep="")))
x.sam <- x[,"Score_d"]
ix.sam <- which(abs(x.sam) > GLOBAL[["z.abs.cut"]] )
gn.nms <- names(ix.sam)
##################################################
# keep genes with z score above cutoff and genes that are expressed in either th17 or th0

# AND expressed in either th17 or th0
# ix.th0 <- which(th0.mean > GLOBAL[["min.th17.or.th0.rpkm"]] )
# ix.th17 <- which(th17.mean > GLOBAL[["min.th17.or.th0.rpkm"]] )
# gn.nms <- intersect(names(ix.sam),names(ix.th0))
# gn.nms <- intersect(gn.nms,names(ix.th17))

# which genes are included? based on zscore > z.abs.cut (z.abs.cut defined in main)
# ix.z <- which(abs(x.sam) > GLOBAL[["z.abs.cut"]] )
# x.med <- x[,"median_rpkm"]
# # which genes are included? based on median rpkm > median.abs.cut (median.abs.cut defined in main)
# ix.med <- which(abs(x.med) > GLOBAL[["median.abs.cut"]])
# gn.nms <- intersect(names(ix.z),names(ix.med))
# if(GLOBAL[["type"]]=="rnaseq"){
#   x.med <- x[,"median_rpkm"]
#   # which genes are included? based on median rpkm > median.abs.cut (median.abs.cut defined in main)
#   ix.med <- which(abs(x.med) > GLOBAL[["median.abs.cut"]])
#   gn.nms <- intersect(names(ix.z),names(ix.med))
# } else {
#   gn.nms <- names(ix.z)
# }
gn.nms <- unique(c(gn.nms,GLOBAL[["known.tfs"]],GLOBAL[["known.gns"]],additional.regulators))

# get dataset for only required genes
cat("creating ratios\n")
bad.ix <- which(!(gn.nms %in% rownames(dataset)))
if(length(bad.ix)>0) {
  gn.nms <- gn.nms[-bad.ix]
}


if(GLOBAL[["type"]]=="rnaseq"){
  # log2 transform the data
  ratios <- log2(dataset[gn.nms,]+1)
} else {
  # already log2 transformed by RMA normalization
  ratios <- dataset[gn.nms,]
}

if(GLOBAL[["type"]]=="immgen") {
  # immgen has 508 cel files but only ~167 unique experiments.  Some experiments have many repeats.
  # To allow for more fair CV by inferelator and cut run time I take median over repeats.
  # This is not done for th17 b/c there the number of experiments is to small (our power loss will be to great).
  ratios <- avg.over.reps.immgen(ratios)
}


# create cluster stack
cat("creating clusterStack\n")
clusterStack <- createClusterStack(ratios)

# create cluster stack
cat("creating colMap\n")
is.ts <- logical(dim(ratios)[2])
del.t <- numeric(dim(ratios)[2])
is.ts[grep("ts",colnames(ratios),perl=TRUE)] <- TRUE
########################
# AM kludge
# switch 13 with 12
# AM this is hardwired using ts names
if(GLOBAL[["type"]]=="rnaseq"){
	del.t[1:5] <- c(0,120,180,600,1920) # for th0 time series
	del.t[6:13] <- c(0,120,180,180,180,240,480,1440) # for th17 time series
}
########################
colMap <- createColMap(ratios=ratios,is.ts=is.ts,del.t=del.t)

# create tfnames
cat("creating tfNames (human tfs found in dataset)\n")
load(paste(path.input,"humanTFs.RData",sep=""))
tfNames <- rownames(ratios)[which(rownames(ratios) %in% humanTFNames)]
additional.regulators <- additional.regulators[which(additional.regulators %in% rownames(ratios))]
# remove a falsely annotated histone from the list of possible tfs
tfNames  <- tfNames[-grep("HIST",tfNames)]
tfNames <- unique(c(tfNames,additional.regulators))

# create knockOutConfList from rnaseq data
cat("create knockOutConfList from rnaseq data")
knockOutConfList <- list()
ko.tf.nms <- GLOBAL[["known.tfs"]]
x <- as.matrix(read.table(paste(path.input,"rnaseqDiffExpMedianZscoresPerTf.xls",sep=""),sep="\t"))
for(i in 1:length(ko.tf.nms)){
  ix <- grep(ko.tf.nms[i],colnames(x),ignore.case=TRUE)
  knockOutConfList[[ko.tf.nms[i]]] <- sort(abs(x[,ix]),decreasing=TRUE)
}

cat("saving data structures (ratios,clusterStack,colMap,tfNames,knockOutConfList) in:\n",path.output,"\n")
save(ratios,file=paste(path.output,"ratios.RData",sep=""))
save(clusterStack,file=paste(path.output,"clusterStack.RData",sep=""))
save(colMap,file=paste(path.output,"colMap.RData",sep=""))
save(tfNames,file=paste(path.output,"tfNames.RData",sep=""))
save(knockOutConfList,file=paste(path.output,"knockOutConfList.RData",sep=""))




















