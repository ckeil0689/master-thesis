##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

source("r_scripts/th17/used_for_paper/util.R")
# setup paths
path.input <- "input/th17/used_for_paper/combinedAnalysis/"
path.output.in.input <- "input/th17/used_for_paper/enrichmentAnalysis/"
path.output.in.results <- "results/enrichmentAnalysis/"

# setup input file names
file.nm.combined.scores.ri <- paste(sep="","score_combine_RI_matrix_zcut_",z.abs.cut,"_",date.is,".xls")
file.nm.combined.scores.kc <- paste(sep="","score_combine_KC_matrix_zcut_",z.abs.cut,"_",date.is,".xls")
file.nm.combined.scores.i <- paste(sep="","score_combine_I_matrix_zcut_",z.abs.cut,"_",date.is,".xls")
file.nm.sam.scores <- "samTh17VsTh0Zscores.xls"
file.nm.d.r <- "rnaseqDataset.xls"
file.nm.d.i <- "immgenRatiosTable.txt"
# setup output file names
file.nm.qc.ri <- paste(sep="","qc_enrichment_scores_zcut_ri_",z.abs.cut,"_",date.is,".pdf")
file.nm.qc.i <- paste(sep="","qc_enrichment_scores_zcut_i_",z.abs.cut,"_",date.is,".pdf")
file.nm.enrchmnt.scores <- paste("enrichment_scores_zcut_",z.abs.cut,"_",date.is,".xls",sep="")
# load scores combined matrices
s.c.m.ri <- as.matrix(read.table(paste(sep="",path.input,file.nm.combined.scores.ri),sep="\t"))
s.c.m.kc <- as.matrix(read.table(paste(sep="",path.input,file.nm.combined.scores.kc),sep="\t"))
s.c.m.i <- as.matrix(read.table(paste(sep="",path.input,file.nm.combined.scores.i),sep="\t"))
# load SAM scores
s.c.m.sam <- as.matrix(read.table(paste(sep="",path.input,"../",file.nm.sam.scores),sep="\t"))

# load rnaseq data and immgen data
d.r <- as.matrix(read.table(paste(sep="",path.input,"../",file.nm.d.r),sep="\t"))
rownames(d.r) <- toupper(rownames(d.r))
d.i <- as.matrix(read.table(paste(sep="",path.input,"../",file.nm.d.i),sep="\t"))
rownames(d.i) <- toupper(rownames(d.i))

# get ranked list of core th17 target genes by summing over the scores each gene has with the core tfs
r.core.kc <- sort(s.c.m.kc[,"sum.score"],decreasing=T)
r.core.sam <- s.c.m.sam[,"Score_d"][sort(abs(s.c.m.sam[,"Score_d"]),decreasing=T,index.return=T)$ix]

# get a ranked list of targets for each tf (only non zero targets are considered) immgen and rnaseq data
r.new <- list()
new.tfs <- colnames(s.c.m.ri)
tfs.with.no.trgts.ix <- numeric()
for(i in 1:length(new.tfs)){
  tf <- new.tfs[i]
  # get list of ranked trgt genes for new tf i
  r.new[[ tf ]] <- sort(s.c.m.ri[ ,tf ],decreasing=T)
  ## if no targets found for tf remove it from enrichment analysis
  if (r.new[[tf]][1]==0){
    ## AM a kludge for tfs that have no targets
    tfs.with.no.trgts.ix <- c(tfs.with.no.trgts.ix,i)
  }
}
## remove tfs with no targets
if(length(tfs.with.no.trgts.ix)>0){
  new.tfs <- new.tfs[-tfs.with.no.trgts.ix]
  r.new <- r.new[-tfs.with.no.trgts.ix]
}
# s.ri: here i keep the statistic scores
s.ri <- list()
# l1 the list of core th17 genes with their scores
l1.th17 <- r.core.kc
pdf(file=paste(sep="","results/",file.nm.qc.ri))
for(i in 1:length(new.tfs)) {
  tf <- new.tfs[i]
  cat(tf," z.auc = ",sep="")
  # l2: list of trgt genes genes with their scores for tf i
  l2 <- r.new[[ tf ]]
  # perform enrichment analysis
  s.ri[[tf]] <- get.running.sum.statistic(l1.th17,l2,N=N,plot.it=TRUE,plot.name=tf)
  cat(s.ri[[tf]][["z.auc"]],"\n")
}
dev.off()
cat("\n")

# get a ranked list of targets for each tf (only non zero targets are considered) only immgen data
r.new <- list()
new.tfs <- colnames(s.c.m.i)
tfs.with.no.trgts.ix <- numeric()
for(i in 1:length(new.tfs)){
  tf <- new.tfs[i]
  # get list of ranked trgt genes for new tf i
  r.new[[ tf ]] <- sort(s.c.m.i[ ,tf ],decreasing=T)
  ## if no targets found for tf remove it from enrichment analysis
  if (r.new[[tf]][1]==0){
    ## AM a kludge for tfs that have no targets
    ## if they did not have trgts in RI they surely don't have in i alone so removed tfs from before are still valid
    tfs.with.no.trgts.ix <- c(tfs.with.no.trgts.ix,i)
  }
}
## remove tfs with no targets
if(length(tfs.with.no.trgts.ix)>0){
  new.tfs <- new.tfs[-tfs.with.no.trgts.ix]
  r.new <- r.new[-tfs.with.no.trgts.ix]
}
# s.i: here i keep the statistic scores
s.i <- list()
# l1 the list of core th17 genes with their scores
l1.th17 <- r.core.kc
pdf(file=paste(sep="","results/",file.nm.qc.i))
for(i in 1:length(new.tfs)) {
  tf <- new.tfs[i]
  cat(tf," z.auc = ",sep="")
  # l2: list of trgt genes genes with their scores for tf i
  l2 <- r.new[[ tf ]]
  # perform enrichment analysis
  s.i[[tf]] <- get.running.sum.statistic(l1.th17,l2,N=N,plot.it=TRUE,plot.name=tf)
  cat(s.i[[tf]][["z.auc"]],"\n")
}
dev.off()
cat("\n")

#~~~~~~~~~~~~~~~~~#
# find enrichment #
#~~~~~~~~~~~~~~~~~#
# define cutoffs:
# total genes ~2100
cut.KC.qunt <- .75 # score 64 or above, results in 508 genes
cut.SAM.qunt <- .94  # abs(z score) 5.21 or above, results in 575 genes
cut.RI.qunt <- 0.9
cut.I.qunt <- .9

  
# KC core with trgts from RI
m.core <- r.core.kc
m.tf <- s.c.m.ri
cutoff.tf <- quantile(m.tf[which(m.tf!=0)],cut.RI.qunt) # cutoff.tf >= 13 length = average tf trgt list = 162
cutoff.core <- quantile(abs(m.core[which(m.core!=0)]),cut.KC.qunt) # cutoff.KO > 64 length = 508
m.tf[which(m.tf<cutoff.tf)] <- 0
trgt.gn.lst <- names(m.core)[which(abs(m.core)>cutoff.core)]
x.KC.RI <- FisherEnrichedTFsToTargetGenes(m.tf, trgt.gn.lst)
# get correlations btwn tf and it's trgts list
x.r <- sapply(colnames(m.tf),getCorTfWithTrgts, d.r, m.tf, trgt.gn.lst)
x.KC.RI$cor.med.r <- sapply(x.r,median)
x.i <- sapply(colnames(m.tf),getCorTfWithTrgts, d.i, m.tf, trgt.gn.lst)
x.KC.RI$cor.med.i <- sapply(x.i,median)


# SAM core with trgts from RI
m.core <- r.core.sam
m.tf <- s.c.m.ri
cutoff.tf <- quantile(m.tf[which(m.tf!=0)],cut.RI.qunt) # cutoff.tf >= 13 length = average tf trgt list = 162
cutoff.core <- quantile(abs(m.core[which(m.core!=0)]),cut.SAM.qunt) # abs(cutoff.SAM) > 5.1 length = 575
m.tf[which(m.tf<cutoff.tf)] <- 0
trgt.gn.lst <- names(m.core)[which(abs(m.core)>cutoff.core)]
x.SAM.RI <- FisherEnrichedTFsToTargetGenes(m.tf, trgt.gn.lst)
# get correlations btwn tf and it's trgts list
x.r <- sapply(colnames(m.tf),getCorTfWithTrgts, d.r, m.tf, trgt.gn.lst)
x.SAM.RI$cor.med.r <- sapply(x.r,median)
x.i <- sapply(colnames(m.tf),getCorTfWithTrgts, d.i, m.tf, trgt.gn.lst)
x.SAM.RI$cor.med.i <- sapply(x.i,median)

# KC core with trgts from I
m.core <- r.core.kc
m.tf <- s.c.m.i
cutoff.tf <- quantile(m.tf[which(m.tf!=0)],cut.I.qunt) # cutoff.tf >= 9.6 length = average tf trgt list = 121
cutoff.core <- quantile(abs(m.core[which(m.core!=0)]),cut.KC.qunt) # cutoff.core <- 64 length = 508
m.tf[which(m.tf<cutoff.tf)] <- 0
trgt.gn.lst <- names(m.core)[which(abs(m.core)>cutoff.core)]
x.KC.I <- FisherEnrichedTFsToTargetGenes(m.tf, trgt.gn.lst)
# get correlations btwn tf and it's trgts list
x.r <- sapply(colnames(m.tf),getCorTfWithTrgts, d.r, m.tf, trgt.gn.lst)
x.KC.I$cor.med.r <- sapply(x.r,median)
x.i <- sapply(colnames(m.tf),getCorTfWithTrgts, d.i, m.tf, trgt.gn.lst)
x.KC.I$cor.med.i <- sapply(x.i,median)

# SAM core with trgts from I
m.core <- r.core.sam
m.tf <- s.c.m.i
cutoff.tf <- quantile(m.tf[which(m.tf!=0)],cut.I.qunt) # cutoff.tf >= 9.6 length = average tf trgt list = 121
cutoff.core <- quantile(abs(m.core[which(m.core!=0)]),cut.SAM.qunt) # abs(cutoff.SAM) > 5.1 length = 575
m.tf[which(m.tf<cutoff.tf)] <- 0
trgt.gn.lst <- names(m.core)[which(abs(m.core)>cutoff.core)]
x.SAM.I <- FisherEnrichedTFsToTargetGenes(m.tf, trgt.gn.lst)
# get correlations btwn tf and it's trgts list
x.r <- sapply(colnames(m.tf),getCorTfWithTrgts, d.r, m.tf, trgt.gn.lst)
x.SAM.I$cor.med.r <- sapply(x.r,median)
x.i <- sapply(colnames(m.tf),getCorTfWithTrgts, d.i, m.tf, trgt.gn.lst)
x.SAM.I$cor.med.i <- sapply(x.i,median)

# positive SAM core with trgts from RI
m.core <- r.core.sam
m.tf <- s.c.m.ri
cutoff.tf <- quantile(m.tf[which(m.tf!=0)],cut.RI.qunt) # cutoff.tf >= 13 length = average tf trgt list = 162
cutoff.core <- quantile(m.core[which(m.core!=0)],cut.SAM.qunt) # abs(cutoff.SAM) > 3.48 length = 575
m.tf[which(m.tf<cutoff.tf)] <- 0
trgt.gn.lst <- names(m.core)[which(m.core>cutoff.core)]
x.SAM.RI.pos <- FisherEnrichedTFsToTargetGenes(m.tf, trgt.gn.lst)
# get correlations btwn tf and it's trgts list
x.r <- sapply(colnames(m.tf),getCorTfWithTrgts, d.r, m.tf, trgt.gn.lst)
x.SAM.RI.pos$cor.med.r <- sapply(x.r,median)
x.i <- sapply(colnames(m.tf),getCorTfWithTrgts, d.i, m.tf, trgt.gn.lst)
x.SAM.RI.pos$cor.med.i <- sapply(x.i,median)

# negative SAM core with trgts from RI
m.core <- r.core.sam
m.tf <- s.c.m.ri
cutoff.tf <- quantile(m.tf[which(m.tf!=0)],cut.RI.qunt) # cutoff.tf >= 13 length = average tf trgt list = 162
cutoff.core <- quantile(m.core[which(m.core!=0)],1-cut.SAM.qunt) # abs(cutoff.SAM) > -3.7 length = 575
m.tf[which(m.tf<cutoff.tf)] <- 0
trgt.gn.lst <- names(m.core)[which(m.core<cutoff.core)]
x.SAM.RI.neg <- FisherEnrichedTFsToTargetGenes(m.tf, trgt.gn.lst)
# get correlations btwn tf and it's trgts list
x.r <- sapply(colnames(m.tf),getCorTfWithTrgts, d.r, m.tf, trgt.gn.lst)
x.SAM.RI.neg$cor.med.r <- sapply(x.r,median)
x.i <- sapply(colnames(m.tf),getCorTfWithTrgts, d.i, m.tf, trgt.gn.lst)
x.SAM.RI.neg$cor.med.i <- sapply(x.i,median)

# positive SAM core with trgts from I
m.core <- r.core.sam
m.tf <- s.c.m.i
cutoff.tf <- quantile(m.tf[which(m.tf!=0)],cut.I.qunt) # cutoff.tf >= 9.6 length = average tf trgt list = 121
cutoff.core <- quantile(m.core[which(m.core!=0)],cut.SAM.qunt) # abs(cutoff.SAM) > 3.48 length = 575
m.tf[which(m.tf<cutoff.tf)] <- 0
trgt.gn.lst <- names(m.core)[which(m.core>cutoff.core)]
x.SAM.I.pos <- FisherEnrichedTFsToTargetGenes(m.tf, trgt.gn.lst)
# get correlations btwn tf and it's trgts list
x.r <- sapply(colnames(m.tf),getCorTfWithTrgts, d.r, m.tf, trgt.gn.lst)
x.SAM.I.pos$cor.med.r <- sapply(x.r,median)
x.i <- sapply(colnames(m.tf),getCorTfWithTrgts, d.i, m.tf, trgt.gn.lst)
x.SAM.I.pos$cor.med.i <- sapply(x.i,median)

# negative SAM core with trgts from I
m.core <- r.core.sam
m.tf <- s.c.m.i
cutoff.tf <- quantile(m.tf[which(m.tf!=0)],cut.I.qunt) # cutoff.tf >= 9.6 length = average tf trgt list = 121
cutoff.core <- quantile(m.core[which(m.core!=0)],1-cut.SAM.qunt) # abs(cutoff.SAM) > -3.7 length = 575
m.tf[which(m.tf<cutoff.tf)] <- 0
trgt.gn.lst <- names(m.core)[which(m.core<cutoff.core)]
x.SAM.I.neg <- FisherEnrichedTFsToTargetGenes(m.tf, trgt.gn.lst)
# get correlations btwn tf and it's trgts list
x.r <- sapply(colnames(m.tf),getCorTfWithTrgts, d.r, m.tf, trgt.gn.lst)
x.SAM.I.neg$cor.med.r <- sapply(x.r,median)
x.i <- sapply(colnames(m.tf),getCorTfWithTrgts, d.i, m.tf, trgt.gn.lst)
x.SAM.I.neg$cor.med.i <- sapply(x.i,median)

# write enrichment output
z.auc.ri <- sapply(s.ri,function(i) i$z.auc)
s.auc.ri <- sapply(s.ri,function(i) i$s.auc)
z.auc.i <- sapply(s.i,function(i) i$z.auc)
s.auc.i <- sapply(s.i,function(i) i$s.auc)
## score.kc <- r.core.kc#[new.tfs]
score.core.tfs <- s.c.m.kc[,core.tfs]#[new.tfs,core.tfs]
score.sam <-  s.c.m.sam[,"Score_d"]
## ix.missing.in.sam <- which(! new.tfs %in% names(score.sam))
## score.sam <- c(score.sam,0)
## names(score.sam)[length(score.sam)] <- new.tfs[ix.missing.in.sam]
gn.nms <- names(sort(x.KC.RI$p.vals))
## score.core.tfs does not have a score for all the tfs that are diff expressed
## here I add them as to score.core.tfs with a zero score
miss.gns <- gn.nms[which(!gn.nms %in% rownames(score.core.tfs))]
m <- matrix(0,nr=length(miss.gns),nc=dim(score.core.tfs)[2])
rownames(m) <- miss.gns
score.core.tfs <- rbind(score.core.tfs,m)
rm(miss.gns,m)


e.mat <- cbind(score.sam[gn.nms],apply(score.core.tfs[gn.nms,],1,sum),
               x.KC.RI$p.vals[gn.nms],x.KC.I$p.vals[gn.nms],x.SAM.RI$p.vals[gn.nms],x.SAM.I$p.vals[gn.nms],
               x.SAM.RI.pos$p.vals[gn.nms],x.SAM.I.pos$p.vals[gn.nms],x.SAM.RI.neg$p.vals[gn.nms],x.SAM.I.neg$p.vals[gn.nms],
               s.auc.ri[gn.nms],z.auc.ri[gn.nms],s.auc.i[gn.nms],z.auc.i[gn.nms],
               score.core.tfs[gn.nms,],
               x.KC.RI$cor.med.r[gn.nms],x.KC.RI$cor.med.i[gn.nms],x.KC.I$cor.med.r[gn.nms],x.KC.I$cor.med.i[gn.nms],
               x.SAM.RI$cor.med.r[gn.nms], x.SAM.RI.pos$cor.med.r[gn.nms], x.SAM.RI.neg$cor.med.r[gn.nms], 
               x.SAM.RI$cor.med.i[gn.nms],x.SAM.RI.pos$cor.med.i[gn.nms],x.SAM.RI.neg$cor.med.i[gn.nms],
               x.SAM.I$cor.med.r[gn.nms],x.SAM.I.pos$cor.med.r[gn.nms],x.SAM.I.neg$cor.med.r[gn.nms],
               x.SAM.I$cor.med.i[gn.nms],x.SAM.I.pos$cor.med.i[gn.nms],x.SAM.I.neg$cor.med.i[gn.nms])
colnames(e.mat) <- c("score_sam","score_kc",
                     "pval_KC_RI","pval_KC_I","pval_SAM_RI","pval_SAM_I",
                     "pval_SAM_RI_pos","pval_SAM_I_pos","pval_SAM_RI_neg","pval_SAM_I_neg",
                     "s_auc_RI","z_auc_RI","s_auc_I","z_auc_I",
                     colnames(score.core.tfs),
                     "cor_KC_RI_r","cor_KC_RI_i","cor_KC_I_r","cor_KCI_i",
                     "cor_SAM_RI_r","pos_cor_SAM_RI_r","neg_cor_SAM_RI_r",
                     "cor_SAM_RI_i","pos_cor_SAM_RI_i","neg_cor_SAM_RI_i",
                     "cor_SAM_I_r","pos_cor_SAM_I_r","neg_cor_SAM_I_r",
                     "cor_SAM_I_i","pos_cor_SAM_I_i","neg_cor_SAM_I_i")
write.table(e.mat,file=paste(path.output.in.results,file.nm.enrchmnt.scores,sep=""),sep="\t")
## write.table(e.mat,file=paste(path.output.in.input,file.nm.enrchmnt.scores,sep=""),sep="\t")











