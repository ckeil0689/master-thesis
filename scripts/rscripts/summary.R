##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Nov 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# create a summary table for results

# set paths
path.input <- "input/th17/used_for_paper/"
path.input.cytoscape <- paste(sep="","results/usedForPaper/cytoscape_",date.cytoscape.run,"/")
path.output <- "results/"

# set file names
f.nm.rpkm <- paste(sep="",path.input,"ranseqDatasetNoQuartNorm.RData") # which experiment to compare
f.nm.specificity <- paste(sep="",path.input.cytoscape,"specificity_mat_zcut_",z.abs.cut,"_",date.cytoscape.run,".txt") 
f.nm.sum.scores <- paste(sep="",path.input.cytoscape,"sum_scores_",date.cytoscape.run,".txt") 
f.nm.diff.exrpsn <- paste(sep="",path.input.cytoscape,"diff_exprssn_zcut_",z.abs.cut,"_",date.cytoscape.run,".txt") 

# get specificity data
specificity.mat <- read.delim(file=f.nm.specificity,sep="\t")
rownames(specificity.mat) <- specificity.mat[,"gene_id"]
specificity.mat <- specificity.mat[,-1]
# get sum scores data
sum.scores.mat <- read.delim(file=f.nm.sum.scores,sep="\t")
rownames(sum.scores.mat) <- sum.scores.mat[,"gene_id"]
sum.scores.mat <- sum.scores.mat[,-1]


# get differential expression data
diff.exrpsn.mat <- read.delim(file=f.nm.diff.exrpsn,sep="\t")
rownames(diff.exrpsn.mat) <- diff.exrpsn.mat[,"gene_id"]
diff.exrpsn.mat <- diff.exrpsn.mat[,-1]


# get rpkm data
load(f.nm.rpkm)
d <- rnaseq.complete.matrix
expt.names <-colnames(d)
th17.wt.ix <- grep("th17_ss_48hr_.{1,6}_wt",colnames(d),value=T,perl=T)
th17.wt.ix <- c(th17.wt.ix,"SL1858_th17_ts_48hr#1")
th17.mean <- apply(d[,th17.wt.ix],1,mean)
th0.wt.ix <- grep("th0_ss_48hr_.{1,6}_wt",colnames(d),value=T,perl=T)
th0.wt.ix <- c(th0.wt.ix,"SL2673_th0_ts_48hr#1")
th0.mean <- apply(d[,th0.wt.ix],1,mean)
max.mean.rpkm <- pmax(th17.mean,th0.mean)

# get gene names that are in all the four datasets
gn.nms.1 <- rownames(sum.scores.mat)
gn.nms.2 <- rownames(specificity.mat)
gn.nms.3 <- rownames(diff.exrpsn.mat)
gn.nms.4 <- names(th0.mean)
gene_id <- intersect(intersect(gn.nms.1,gn.nms.2),intersect(gn.nms.3,gn.nms.4))

# create summary table
x <- cbind(gene_id,specificity.mat[gene_id,],sum.scores.mat[gene_id,],
	diff.exrpsn.mat[gene_id,],th17.mean[gene_id],th0.mean[gene_id],max.mean.rpkm[gene_id])
colnames(x)[ncol(x)] <- "th0.th17.rpkm.max"
colnames(x)[ncol(x)-1] <- "th0.mean"
colnames(x)[ncol(x)-2] <- "th17.mean"

# write summary table to results
fl.nm <- paste(path.output,"summary_table_",z.abs.cut,"_",date.cytoscape.run,".xls",sep="")
cat(colnames(x),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(sep="\t",x,quote=FALSE,append=TRUE,col.names=FALSE,row.names=FALSE,file=fl.nm)	
	















