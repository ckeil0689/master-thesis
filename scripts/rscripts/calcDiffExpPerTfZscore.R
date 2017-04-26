##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## January 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

path.input <- "input/th17/used_for_paper/"
path.output <- "input/th17/used_for_paper/"
# get required functions
source("r_scripts/th17/used_for_paper/util.R")
# get dataset
load(paste(path.input,"ranseqDatasetNoQuartNorm.RData",sep=""))
d.raw <- rnaseq.complete.matrix
d <- log2(rnaseq.complete.matrix+pseudo.count)
# get expt namesd

expt.names <-colnames(d)
ko.tfs <- c("batf","maf","irf4","stat3","rorc")

# put results in this matrix
median.zscores.per.tf <- matrix(0,nr=dim(d)[1],nc=length(ko.tfs))
rownames(median.zscores.per.tf) <- rownames(d)
colnames(median.zscores.per.tf) <- paste(ko.tfs,"ko","median","zscore",sep="_")

# get all Th17 wt mean and sd
th17.wt.ix <- grep("th17_ss_48hr_.{1,6}_wt",colnames(d),value=T,perl=T)
th17.wt.ix <- c(th17.wt.ix,"SL1858_th17_ts_48hr#1")
wt.mean <- apply(d[,th17.wt.ix],1,mean)
wt.sd <- apply(d[,th17.wt.ix],1,sd)

# for each of these tfs calculate median diff expression zscores
ko.tfs <- c("batf","maf","irf4","stat3","rorc")
for(i in 1:length(ko.tfs)) {
  ## ix.th17.wt <- expt.names[grep(paste("th17","48hr",ko.tfs[i],"wt",sep=".*"),expt.names,ignore.case=TRUE)]
  th17.ko.ix <- expt.names[grep(paste("th17","48hr",ko.tfs[i],"ko",sep=".*"),expt.names,ignore.case=TRUE)]
  ko.mean <-  apply(d[,th17.ko.ix],1,mean)
  ix.good <- which( (ko.mean+wt.mean) > log2(rpkm.cut+pseudo.count))
  zscores <- (wt.mean[ix.good]-ko.mean[ix.good])/wt.sd[ix.good]
  median.zscores.per.tf[names(zscores),i] <- zscores
}

# get rid of na's that were put in by median.zscores.per.tf[names(median.zscores),i] <- median.zscores for genes not in names(median.zscores)
## median.zscores.per.tf[which( is.na(median.zscores.per.tf))] = 0
# get rid of genes that did not have a zscore in any of the ko tfs
rm.ix <-which(apply(abs(median.zscores.per.tf),1,sum)==0)
median.zscores.per.tf <- median.zscores.per.tf[-rm.ix,] 
write.table(median.zscores.per.tf,file=paste(sep="",path.output,"rnaseq_diff_exp_zscores_",  date.cufflinks.data.compiled,".xls"),sep="\t")
