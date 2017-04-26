##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 th17 project
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

##AM check gene order and make sure to use gene names when creating big rnaseq file ##
source("r_scripts/th17/used_for_paper/util.R")

path.input <- paste(sep="","input/th17/used_for_paper/rawData/cufflinks_no_quartnorm_",date.cufflinks.data.compiled,"/")
path.output <- "input/th17/used_for_paper/"

# get all rnaseq experiments with their sl number (e.g. SL1040)
rnaseq.expts <- list.files(path.input)

# currently datasets from cufflinks was run on datasets from chip expt. here i choose only rnaseq expts
## full.expt.names = as.character(read.delim(file=paste(sep="",path.input,"../rnaseq_colnames_Apr_18_2011.txt"),header=F,sep="")$V1)
full.expt.names <- as.character(read.delim(file=paste(sep="",path.input,"../rnaseq_colnames_",date.cufflinks.data.compiled,".txt"),header=F,sep="")$V1)


# no need in removing expts as they were removed in the annotations step
# # remove experiments with N2 inhibitor as they are not published yet by Jun
# full.expt.names <- full.expt.names[-grep("n2",full.expt.names,ignore.case=T)]
# # remove human experiments from dataset (we only treat mouse here)
# full.expt.names <- full.expt.names[-grep("human",full.expt.names,ignore.case=T)]

# get all rnaseq experiments with their sl number (e.g. SL1040)
good.expts <- sapply(strsplit(full.expt.names,"_"),function(i) i[1])
good.expts <- paste(sep="",good.expts,"_genes.expr")
ix.keep <- which(good.expts %in% rnaseq.expts)
good.expts <- good.expts[ix.keep]
full.expt.names <- full.expt.names[ix.keep]

# how many expts do we have
n <- length(good.expts)
cat("reading cufflinks rpkms\n")
for (i in 1:n) {
  # file name
  expt.file <- paste(path.input,good.expts[i],sep="")
  # get table
  x <- read.delim(expt.file)
  # if first iteration
  if(i == 1){
    # for first iteration get gene names for rownames later on
    gn.nms <- as.character(x[,"gene_id"])
    # variable for new complete rnaseq matrix
    rnaseq.complete.matrix <- matrix(,nr=length(gn.nms),nc=0)
    # variable name for column names of rnaseq.complete.matrix
    c.nms <- character()    
  }
  # get expt names into colnames
  c.nms <- c(c.nms,full.expt.names[i])
  # combine to growing matrix (based on gene names)
  tmp <- x[,"FPKM"]
  names(tmp) <- x[,"gene_id"]
  rnaseq.complete.matrix <- cbind(rnaseq.complete.matrix,tmp[gn.nms])
}
# add colnames to rnaseq.complete.matrix
colnames(rnaseq.complete.matrix) <- c.nms
# add rownames to rnaseq.complete.matrix
rownames(rnaseq.complete.matrix) <- toupper(gn.nms)
# take median gene of redundent gene names (about 230 genes with more than one obesrvation (row)) so we can have unique row names
rnaseq.complete.matrix <- take.median.gene.of.non.unique.genes(rnaseq.complete.matrix)$d
f.nm.1 <- paste(path.output,"ranseqDatasetNoQuartNorm_",date.is,".RData", sep="")
f.nm.2 <- paste("results/ranseq_dataset_no_quartnorm_",date.is,".xls",sep="")
cat("saving ranseqDatasetNoQuartNorm.RData\n")
save(rnaseq.complete.matrix,file=f.nm.1)
cat(file=f.nm.2,sep="\t",c("gene_name",colnames(rnaseq.complete.matrix)),"\n")
write.table(rnaseq.complete.matrix,file=f.nm.2,sep="\t",append=T,col.names=F)

