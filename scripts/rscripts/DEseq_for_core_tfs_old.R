##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## May 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

## add 1 read psudo count
## take pvals instead of adjusted pvals


## helper function
plotDE <- function( res ) {
 plot( 
 res$baseMean, 
 res$log2FoldChange, 
 log="x", pch=20, cex=.1, 
 col = ifelse( res$padj < .1, "red", "black" ) ) 
}
########## read data ##########
library(DESeq)
cat("reading data\n")
# set paths
path.input <- "input/th17/used_for_paper/"
path.input.counts <- paste(sep="","input/th17/used_for_paper/rawData/rnaseq_ht_counts_",date.htseq.data.compiled,"/")
path.results <- "results/diff_expression/"

# set file names
f.nm.info <- paste(sep="",path.input,"diffExpressionExpts_",date.htseq.data.compiled,".txt") # which experiment to compare
f.nm.rpkm <- paste(sep="",path.input,"ranseqDatasetNoQuartNorm.RData") # which experiment to compare

## get count files (counting for each gene how many tags hit it)
count.files <- grep("counts",list.files(path.input.counts),value=T)

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
gns <- names(rpkm.mean)
gns.low.rpkm <- names(rpkm.mean)[which(rpkm.mean < rpkm.cut)]


## read tab delim file with experiments we want to compare (format below):
## name	expts
## Th17.rorawt.rorcwt.vs.Th17.rorako.rorcko	SL5874_SL5878::SL5876_SL5880
## Th17.dmso.vs.Th0.dmso.human	SL5252_SL5255::SL5251_SL5254
cat("read expt mat\n")
expt.mat <- read.delim(file=f.nm.info,colClasses = "character")

for(i in 1:nrow(expt.mat)){
  e.wt <- strsplit(strsplit(expt.mat[i,2],"::")[[1]],"_")[[1]]
  e.ko <- strsplit(strsplit(expt.mat[i,2],"::")[[1]],"_")[[2]]
  e.nm <- expt.mat[i,1]
  e.nm.wt <- strsplit(e.nm,".vs.")[[1]][1]
  e.nm.ko <- strsplit(e.nm,".vs.")[[1]][2]
  cat("working on",e.nm,"\n")
  f.nm.wt<- sapply(e.wt,function(i) grep(paste(sep="",i,"\\."),count.files,value=T))
  f.nm.ko<- sapply(e.ko,function(i) grep(paste(sep="",i,"\\."),count.files,value=T))
  if(i == 1){
    x <- read.delim(file=paste(sep="",path.input.counts,f.nm.wt[1]),header=T)
    n.row <- dim(x)[1]
    gn.nms <- toupper(as.character(x[,1]))
    s.mat.pval <- matrix(0,nc=nrow(expt.mat),nr=length(gn.nms),dimnames=list(gn.nms,expt.mat[,1]))
    s.mat.padj <- matrix(0,nc=nrow(expt.mat),nr=length(gn.nms),dimnames=list(gn.nms,expt.mat[,1]))
    s.mat.log2fc <- matrix(0,nc=nrow(expt.mat),nr=length(gn.nms),dimnames=list(gn.nms,expt.mat[,1]))
  }
  ## get ht experiments for e.wt into a matrix
  m.wt <- matrix(0,nr=n.row,nc=length(f.nm.wt))
  rownames(m.wt) <- gn.nms
  for(j in 1:length(f.nm.wt)){
    x <- read.delim(file=paste(sep="",path.input.counts,f.nm.wt[j]),header=T)
    m.wt[,j] <- as.numeric(x[,2])
  }
  ## get ht experiments for e.ko into a matrix
  m.ko <- matrix(0,nr=n.row,nc=length(f.nm.ko))
  rownames(m.ko) <- gn.nms
  for(j in 1:length(f.nm.ko)){
    x <- read.delim(file=paste(sep="",path.input.counts,f.nm.ko[j]),header=T)
    m.ko[,j] <- as.numeric(x[,2])
  }
  countData <- cbind(m.wt,m.ko) + pseudo.count
  colnames(countData) <- c(e.wt,e.ko)
  conditions <- factor(c(rep(e.nm.wt,length(e.wt)),rep(e.nm.ko,length(e.ko))),levels=c(e.nm.ko,e.nm.wt))
  cds <- newCountDataSet(countData, conditions, sizeFactors = NULL, phenoData = NULL, featureData=NULL)
  cds <- estimateSizeFactors( cds )
  ## if we have no replicates estimate gene's variance from both conditions as biological repeats
  if(dim(cds)[2]==2){
    cds <- estimateVarianceFunctions( cds,pool=T)
  } else {
    cds <- estimateVarianceFunctions( cds )
  }
  res <- nbinomTest( cds, levels(conditions)[1],levels(conditions)[2])
  colnames(res)[3:4] <- levels(conditions)
  ix.low.expressed <- which(res$id %in% gns.low.rpkm)
  res[ix.low.expressed,c("pval","padj")] <- 1
  gns.with.rpkm <- res$id[which(res$id %in% gns)]
  res$mean.rpkm.th17 <- res$mean.rpkm.th0 <- rep(0,length(res$id))
  res$mean.rpkm.th17[which(res$id %in% gns.with.rpkm)] <- th17.mean[gns.with.rpkm]
  res$mean.rpkm.th0[which(res$id %in% gns.with.rpkm)] <- th0.mean[gns.with.rpkm]
  s.mat.pval[,i] <- res$pval
  s.mat.padj[,i] <- res$padj
  s.mat.log2fc[,i] <- res$log2FoldChange
  write.table(res,file=paste(path.results,e.nm,"_",date.is,".xls",sep=""),sep="\t",row.names=FALSE)
}

## here we 'multiply' (element wise, A_i_j times B_i_j for all i and j) the pval matrix by the sign of the log2 fold change matrix
## this step allows us to view only the pvalue and say if we have down regulation or up regulation

# this matrix will be used for ChIP scores under combine data script
s.mat.pval.signed <- -log10(s.mat.pval) * sign(s.mat.log2fc)
write.table(s.mat.pval.signed,file=paste(path.input,"DEseq_pval_signed_",date.is,".xls",sep=""),sep="\t")
## other results...
write.table(s.mat.pval.signed,file=paste(path.results,"summary.pval.signed.xls",sep=""),sep="\t")
write.table(s.mat.pval,file=paste(path.results,"summary.pval.xls",sep=""),sep="\t")
write.table(s.mat.padj,file=paste(path.results,"summary.pval.adjusted.xls",sep=""),sep="\t")
write.table(s.mat.log2fc,file=paste(path.results,"summary.log2fc.xls",sep=""),sep="\t")

n <- ncol(s.mat.pval.signed)
ix.invivo <- grep("IL17a.gfp.plusve.SI.vs.IL17a.gfp.minusve.SI",colnames(s.mat.pval.signed))
fc.invivo <-  s.mat.log2fc[,ix.invivo]
pval.invivo <-  s.mat.pval.signed[,ix.invivo]
ix <- (n-2):n
s.mat.pval.signed.small <- s.mat.pval.signed[,ix]
s.mat.log2fc.small <- s.mat.log2fc[,ix]
# ix <- rev(colnames(s.mat.pval.signed.small))[1:3]
average.specificity <- apply(s.mat.pval.signed.small,1,mean)
sum.specificity <- apply(s.mat.pval.signed.small,1,sum)
max.specificity <- apply(abs(s.mat.pval.signed.small),1,max)*sign(apply(s.mat.pval.signed.small,1,max))
min.specificity <- apply(abs(s.mat.pval.signed.small),1,min)*sign(apply(s.mat.pval.signed.small,1,min))
s.mat.pval.signed.small <- cbind(s.mat.log2fc.small,s.mat.pval.signed.small,fc.invivo,pval.invivo,
	average.specificity,sum.specificity,max.specificity,min.specificity)
x <- colnames(s.mat.pval.signed.small)
colnames(s.mat.pval.signed.small)[1:3] <- paste(sep="","fc.",x[1:3])
colnames(s.mat.pval.signed.small)[4:6] <- paste(sep="","pval.",x[4:6])
s.mat.pval.signed.small <- s.mat.pval.signed.small[order(abs(s.mat.pval.signed.small[,"min.specificity"]),decreasing=T),]
write.table(s.mat.pval.signed.small,file=paste(path.input,"specificity_mat_",date.is,".xls",sep=""),sep="\t")
write.table(s.mat.pval.signed.small,file=paste(path.results,"specificity_mat_",date.is,".xls",sep=""),sep="\t")

