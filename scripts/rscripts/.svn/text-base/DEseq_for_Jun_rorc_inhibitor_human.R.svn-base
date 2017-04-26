l##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## May 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

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
path.output <- "results/"
path.input.counts <- paste(sep="","input/th17/used_for_paper/rawData/rnaseq_ht_counts_",date.htseq.data.compiled,"/")
path.results <- "results/diff_expression/"

date.htseq.data.compiled
## get count files (counting for each gene how many tags hit it)
count.files <- grep("counts",list.files(path.input.counts),value=T)

## read tab delim file with experiments we want to compare (format below):
## name	expts
## Th17.rorawt.rorcwt.vs.Th17.rorako.rorcko	SL5874_SL5878::SL5876_SL5880
## Th17.dmso.vs.Th0.dmso.human	SL5252_SL5255::SL5251_SL5254
expt.mat <- read.delim(file=paste(sep="",path.input,"diffExpressionExptsHuman_",date.htseq.data.compiled,".txt"),colClasses = "character")

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
    ## s.mat.pval[,1] <- gn.nms
    ## s.mat.padj[,1] <- gn.nms
    ## s.mat.log2fc[,1] <- gn.nms
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
  conditions <- as.factor(c(rep(e.nm.wt,length(e.wt)),rep(e.nm.ko,length(e.wt))))
  cds <- newCountDataSet(countData, conditions, sizeFactors = NULL, phenoData = NULL, featureData=NULL)
  cds <- estimateSizeFactors( cds )
  ## if we have no replicates estimate gene's variance from both conditions as biological repeats
  if(dim(cds)[2]==2){
    cds <- estimateVarianceFunctions( cds,pool=T)
  } else {
    cds <- estimateVarianceFunctions( cds )
  }
  res <- nbinomTest( cds, levels(conditions)[1],levels(conditions)[2])
  colnames(res)[3:4] <- c(e.nm.ko,e.nm.wt)
  s.mat.pval[,i] <- res$pval
  s.mat.padj[,i] <- res$padj
  s.mat.log2fc[,i] <- res$log2FoldChange
  write.table(res,file=paste(path.results,e.nm,".xls",sep=""),sep="\t",row.names=FALSE)
}
