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
# get expt names
expt.names <-colnames(rnaseq.complete.matrix)

# for each of these tfs calculate median diff expression zscores
## ko.tfs <- c("batf","maf","irf4","stat3","rorc","rora_and_rorc")
ko.tfs <- c("batf","maf","irf4","stat3","rorc")
# put results in this matrix
median.zscores.per.tf <- matrix(,nr=dim(rnaseq.complete.matrix)[1],nc=length(ko.tfs))
rownames(median.zscores.per.tf) <- rownames(rnaseq.complete.matrix)
colnames(median.zscores.per.tf) <- paste(ko.tfs,"ko","median","zscore",sep="_")
for(i in 1:length(ko.tfs)) {
  ix.th17.wt <- expt.names[grep(paste("th17","48hr",ko.tfs[i],"wt",sep=".*"),expt.names,ignore.case=TRUE)]
  ix.th17.ko <- expt.names[grep(paste("th17","48hr",ko.tfs[i],"ko",sep=".*"),expt.names,ignore.case=TRUE)]
  zscores <- numeric()
  for(j in 1:length(ix.th17.wt)){
    # calc differential expression based on diff in rpkm and diff in FC btwn wt and ctrl. c1 is cutoff such that only genes
    # with summed rpkm > c1 get zscores.  ps.cnt is a pseudo count so we can take log2(wt/ko) (if ko has 0...)
    zscores <- c(zscores,diff.exprs.analysis(wt=rnaseq.complete.matrix[,ix.th17.wt[j]],ko=rnaseq.complete.matrix[,ix.th17.ko[j]],c1=5,ps.cnt=1))
  }
  median.zscores <- get.median.zscore(zscores)
  median.zscores.per.tf[names(median.zscores),i] <- median.zscores
}
# get rid of na's that were put in by median.zscores.per.tf[names(median.zscores),i] <- median.zscores for genes not in names(median.zscores)
median.zscores.per.tf[which( is.na(median.zscores.per.tf))] = 0
# get rid of genes that did not have a zscore in any of the ko tfs
rm.ix <-which(apply(abs(median.zscores.per.tf),1,sum)==0)
median.zscores.per.tf <- median.zscores.per.tf[-rm.ix,] 
## write.table(median.zscores.per.tf,file=paste(sep="",path.output,"rnaseqDiffExpMedianZscoresPerTf.xls"),sep="\t")
write.table(median.zscores.per.tf,file=paste(sep="",path.output,"rnaseq_diff_exp_zscores_",  date.cufflinks.data.compiled,".xls"),sep="\t")
