##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

rm(list=ls())
source("r_scripts/th17/used_for_paper/util.R")
path.input <- "input/th17/V2_used_for_paper/"

######### Combine chipseq with rnaSeq per tf####################
cat("creating summary table for chipSeq AND rnaSeq (ko data) per tf\n")

# getting rnaseq dataset
m.rna <- as.matrix(read.table(paste(path.input,"rnaseqDiffExpMedianZscoresPerTf.xls",sep=""),sep="\t"))

# getting chipseq dataset
m.chip <- as.matrix(read.table(paste(path.input,"chipScoresPerTF.xls",sep="")))
# leave only columns that have final score (in the same column order as m.rna)
tfs <- sapply(strsplit(colnames(m.rna),"_"),function(i) i[1])
ix.keep <- numeric()
for(i in 1:length(tfs)) {
  ix.keep <- c(ix.keep,grep( paste(tfs[i],"p300",sep=".*"), colnames(m.chip), ignore.case=TRUE ) )
}
m.chip <- m.chip[,ix.keep]

# create a matrix with rows for each gene that apears in any analysis (chip-seq or rna-seq)
gn.nms.rnaseq <- rownames(m.rna)
gn.nms.chipseq <- rownames(m.chip)
gn.nms.all <- unique(c(gn.nms.rnaseq,gn.nms.chipseq))

# create a final matrix for all chip-seq and rna-seq results
m <- matrix(0, nr=length(gn.nms.all),nc=(dim(m.rna)[2]+dim(m.chip)[2]))
rownames(m) <- gn.nms.all
colnames(m) <- c(paste("rna_seq_",tfs,"_zscore",sep=""),
                 paste("chip_seq_",tfs,"_score",sep=""))

# fill in rnaseq ko data
m[gn.nms.rnaseq,1:dim(m.rna)[2]] <- m.rna[gn.nms.rnaseq,]

# fill in chipseq data
m[gn.nms.chipseq,(dim(m.chip)[2]+1):(dim(m.chip)[2]+dim(m.chip)[2])] <- m.chip[gn.nms.chipseq,]

# add zscores to p300 rank based on chipseq
dist.zscores.chip <- sort(m[,1],decreasing =T)
for(i in 1:3){
   ix <- names(which(m.chip[,i]>0))
   ix2 <- names(m.chip[ix,i])[sort(m.chip[ix,i],index.return=T,decreasing=T)$ix]
   m[ix2,i+3] <- dist.zscores.chip[1:length(ix2)]
}

# remove negative values
m[which(m<0)] <- 0
rownames(m) <- toupper(rownames(m))
#write.table(m,file=paste(path.input,"/rnaseq_chipseq_scores.xls",sep="\t")

###############################################
cat("adding network inference scores for batf,maf,irf4,stat3, and rorc\n")
x <- as.matrix(read.table(paste(path.input,"infResults/th17_1000_mix_clr_inf_median_scores_2.xls",sep="")))
rownames(x) <- toupper(rownames(x))
ix <- rownames(x)[which(rownames(x) %in% rownames(m))]
m.colnames <- colnames(m)
m <- cbind(m,matrix(0,nr=dim(m)[1],nc=length(tfs)))
m.colnames <- c(m.colnames,paste(toupper(tfs),"Inf_mixCLR_zscores",sep="_"))
colnames(m) <- m.colnames
run1 <- toupper(tfs)
for(i in 1:length(run1)) {
  ix.m <- grep(paste(run1[i],"_Inf",sep=""),colnames(m))
  m[ix,ix.m] <- x[ix,run1[i]]
}
m <- m[ix,]

write.table(m,file=paste(path.input,dim(x)[1],"_gene_rnaseq_chipseq_InfMixCLR_zscores.xls",sep=""),sep="\t")

###############################################
cat("combining results for inf with chipseq")
ix.rnaseq <- grep("rna",colnames(m))
ix.chipseq <- grep("chip",colnames(m))
ix.inf <- grep("Inf",colnames(m))
m.colnames <- colnames(m)
m <- cbind(m,combine_mtrcs(m[,ix.chipseq],m[,ix.inf]))
m.colnames <- c(m.colnames,paste(toupper(tfs),"_Inf_Chip_pseudo_zscore",sep=""))
ix.chipseqInf <- grep("pseudo",m.colnames)

cat("combining results for infchipseq with rnaseq ko")
m <- cbind(m,combine_mtrcs(m[,ix.rnaseq],m[,ix.chipseqInf]))
m.colnames <- c(m.colnames,paste(toupper(tfs),"_Inf_Chip_Rnaseq_pseudo_zscore",sep=""))
colnames(m) <- m.colnames
write.table(m,file=paste(path.input,dim(x)[1],"_gene_rnaseq_chipseq_InfMixCLR_zscores.xls",sep=""),sep="\t")


































