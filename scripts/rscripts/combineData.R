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
path.input.immgen <- "input/th17/used_for_paper/infResults/immgen/"
path.input.rnaseq <- "input/th17/used_for_paper/infResults/rnaseq/"
path.input.ko <- "input/th17/used_for_paper/"
path.input.chip <- "input/th17/used_for_paper/"
path.output.in.input <- "input/th17/used_for_paper/combinedAnalysis/"
path.output.in.results <- "results/combinedAnalysis/"
path.input.immgen.inf.data <-  paste(sep="","input/th17/used_for_paper/infDataStructures/immgen/z_abs_cut_",z.abs.cut,"/")
path.input.rnaseq.inf.data <-  paste(sep="","input/th17/used_for_paper/infDataStructures/rnaseq/z_abs_cut_",z.abs.cut,"/")
path.input.sam <- "input/th17/used_for_paper/"

# setup file names
file.nm.immgen <- paste(sep="","clrinf_cnfdnc_zcut_",z.abs.cut,"_Nb_",num.boots.immgen,"_",  date.inf.run,".xls")
file.nm.rnaseq <- paste(sep="","clrinf_cnfdnc_zcut_",z.abs.cut,"_Nb_",num.boots.rnaseq,"_",date.inf.run,".xls")
file.nm.ko <- paste(sep="","DEseq_pval_signed_",date.deseq.run,".xls")
file.nm.chip <- paste(sep="","chip_scores_",date.chip.run,".xls")
file.nm.sam <- "samTh17VsTh0Zscores.xls"

# load scores
#load sam diff expression
f.nm <- paste(sep="",path.input.sam,file.nm.sam)
cat("loading SAM scores file: ",f.nm,"\n")
scores.sam <- as.matrix(read.table(f.nm,sep="\t"))
# here we keep focusing on genes with abs diff expression zscore greater than z.abs.cut
gn.nms.diff <- rownames(scores.sam)[ which( abs(scores.sam[,"Score_d"]) > z.abs.cut ) ]
scores.sam <- scores.sam[gn.nms.diff,"Score_d"]

#immgen
f.nm <- paste(sep="",path.input.immgen,file.nm.immgen)
cat("loading immgen scores file: ",f.nm,"\n")
scores.immgen <- as.matrix(read.table(f.nm,sep="\t"))
#rnaseq
f.nm <- paste(sep="",path.input.rnaseq,file.nm.rnaseq)
cat("loading rnaseq scores file: ",f.nm,"\n")
scores.rnaseq <- as.matrix(read.table(f.nm,sep="\t"))
#KO
f.nm <- paste(sep="",path.input.ko,file.nm.ko)
cat("loading KO scores file: ",f.nm,"\n")
to.load <- c("Th17.batf.wt.vs.Th17.batf.ko","Th17.maf.wt.vs.Th17.maf.ko","Th17.irf4.wt.vs.Th17.irf4.ko","Th17.stat3.wt.vs.Th17.stat3.ko","Th17.rorc.wt.vs.Th17.rorc.ko")
scores.ko <- as.matrix(read.table(f.nm,sep="\t"))[,to.load]
#chip
# chip scores to load for each tf
f.nm <- paste(sep="",path.input.chip,file.nm.chip)
cat("loading ChIP scores file: ",f.nm,"\n")
if(type=="k_r_i_repressor" | type=="abs"){
  # load scores with no P300 boost
  to.load <- c("BATF_th17_minus_th0","IRF4_th17_minus_th0", # have th0 ctrl
               "MAF_th17_minus_th0","STAT3_th17_minus_th0","RORC_th17_minus_th0") # don't have th0 ctrl
  scores.chip <- as.matrix(read.table(f.nm,sep="\t"))[,to.load]
} else {
  # load scores with P300 boost
  to.load <- c("BATF_th17_minus_th0_plus_p300","IRF4_th17_minus_th0_plus_p300", # have th0 ctrl
               "MAF_th17_minus_th0_plus_p300","STAT3_th17_minus_th0_plus_p300","RORC_th17_minus_th0_plus_p300") # don't have th0 ctrl
  scores.chip <- as.matrix(read.table(f.nm,sep="\t"))[,to.load]
}

# simplify column names for ko and rnaseq (make into "BATF",...)
colnames(scores.chip) <- toupper(sapply(strsplit(colnames(scores.chip),"_"),function(i) i[1]))
colnames(scores.ko) <- toupper(sapply(strsplit(colnames(scores.ko),"\\."),function(i) i[2]))

# load datasets (rnaseq and immgen)
load(paste(sep="",path.input.immgen.inf.data,"ratios.RData"))
tfs <- colnames(scores.immgen)
sign.cor.mat.immgen <- sign(cor(t(ratios))[,tfs])
load(paste(sep="",path.input.rnaseq.inf.data,"ratios.RData"))
tfs <- colnames(scores.rnaseq)
sign.cor.mat.rnaseq <- sign(cor(t(ratios))[,tfs])
rm(tfs,ratios)

## put sign on confidence scores from immgen and ranseq
scores.immgen <- scores.immgen*sign(sign.cor.mat.immgen)
scores.rnaseq <- scores.rnaseq*sign(sign.cor.mat.rnaseq)

# convert scores to relative ranks: 1,...,0 where high relative rank means high confidence
# when possible we discriminate if regulatory role is activation or repression (or we discard such notation and just concentrate on the absolute confidence score).
# for KO's it is simple: if a tf KO leads to down-regulation of target then tf is an activator (direct of indirect), and vice versa for up-regulation/repressor
# for expression data we use the sign of the correlation between tf and target over the dataset
# abs regulation
if(type=="abs"){
  rel.rank.immgen <- convert.scores.to.relative.ranks.pos(abs(scores.immgen))
  rel.rank.rnaseq <- convert.scores.to.relative.ranks.pos(abs(scores.rnaseq))
  rel.rank.ko <- convert.scores.to.relative.ranks.pos(abs(scores.ko))
  rel.rank.sam <- convert.scores.to.relative.ranks.pos(abs(scores.sam))
} else if (type=="k_activator"){
  rel.rank.immgen <- convert.scores.to.relative.ranks.pos(abs(scores.immgen))
  rel.rank.rnaseq <- convert.scores.to.relative.ranks.pos(abs(scores.rnaseq))
  rel.rank.ko <- convert.scores.to.relative.ranks.pos(scores.ko)
  rel.rank.sam <- convert.scores.to.relative.ranks.pos(scores.sam)
} else if (type=="k_r_i_activator"){
  rel.rank.immgen <- convert.scores.to.relative.ranks.pos(scores.immgen)
  rel.rank.rnaseq <- convert.scores.to.relative.ranks.pos(scores.rnaseq)
  rel.rank.ko <- convert.scores.to.relative.ranks.pos(scores.ko)
  rel.rank.sam <- convert.scores.to.relative.ranks.pos(scores.sam)
} else if (type=="k_r_i_repressor"){
  rel.rank.immgen <- convert.scores.to.relative.ranks.pos(-1*scores.immgen)
  rel.rank.rnaseq <- convert.scores.to.relative.ranks.pos(-1*scores.rnaseq)
  rel.rank.ko <- convert.scores.to.relative.ranks.pos(-1*scores.ko)
  rel.rank.sam <- convert.scores.to.relative.ranks.pos(-1*scores.sam)
}
rel.rank.chip <- convert.scores.to.relative.ranks(scores.chip) # chip has no 'sign'

# each datatype gives you information for different genes (e.g. no probe for gene x in immgen may still be found in rnaseq).
# get union of all genes
gn.nms <- unique(c(rownames(rel.rank.immgen),
            rownames(rel.rank.rnaseq),
            rownames(rel.rank.ko),
            rownames(rel.rank.chip)))
tf.nms <- unique(c(colnames(rel.rank.immgen),
            colnames(rel.rank.rnaseq),
            colnames(rel.rank.ko),
            colnames(rel.rank.chip)))

## compile the list of genes we look at (differentially expressed OR core tfs OR known genes)
## very few genes don't survive this
gn.nms.keepers <- gn.nms.diff[which(gn.nms.diff %in% gn.nms)]
gn.nms.keepers <- unique(c(gn.nms.keepers,core.tfs))

## cat("combining SKCRI scores.\n")
## # combine SKCRI
## x <- combine_kcri(knockout=abs(rel.rank.ko),chip=abs(rel.rank.chip),rnaseq=abs(rel.rank.rnaseq),immgen=abs(rel.rank.immgen),
##                                       row.keepers=gn.nms.keepers,col.keepers=tf.nms)
## score.combined.matrix <- x[,which(apply(x,2,sum)>0)]
## sum.score.kcri <- apply(score.combined.matrix[,core.tfs],1,sum)
## sum.rel.ranks <- sum.score.kcri[order(sum.score.kcri,decreasing=T)]
## sam.score <- scores.sam[gn.nms.keepers]
## rel.rank.kcri <- convert.scores.to.relative.ranks(sum.rel.ranks)
## sum.score <- rel.rank.kcri[gn.nms.keepers] + rel.rank.sam[gn.nms.keepers]
## score.combined.matrix <- cbind(score.combined.matrix,sum.score.kcri,sam.score,sum.score)
## write.table(score.combined.matrix,file=paste(path.output.in.input,"score_combine_SKCRI_matrix_zcut_",z.abs.cut,"_",type,"_",date.is,".xls",sep=""),sep="\t")
## write.table(score.combined.matrix,file=paste(path.output.in.results,"score_combine_SKCRI_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),sep="\t")

cat("combining KCRI scores.\n")
# combine KCRI
score.combined.matrix <- combine_kcri(knockout=abs(rel.rank.ko),chip=abs(rel.rank.chip),rnaseq=abs(rel.rank.rnaseq),immgen=abs(rel.rank.immgen),
                                       row.keepers=gn.nms.keepers,col.keepers=tf.nms)
score.combined.matrix <- score.combined.matrix[,which(apply(score.combined.matrix,2,sum)>0)]
sum.score <- apply(score.combined.matrix[,core.tfs],1,sum)
sam.score <- scores.sam[gn.nms.keepers]
score.combined.matrix <- cbind(score.combined.matrix,sum.score,sam.score)
write.table(score.combined.matrix,file=paste(path.output.in.input,"score_combine_KCRI_matrix_zcut_",z.abs.cut,"_",type,"_",date.is,".xls",sep=""),sep="\t")
write.table(score.combined.matrix,file=paste(path.output.in.results,"score_combine_KCRI_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),sep="\t")

cat("combining KCR.\n")
# combine KCR
score.combined.matrix <- combine_kcri(knockout=abs(rel.rank.ko),chip=abs(rel.rank.chip),rnaseq=abs(rel.rank.rnaseq),immgen=0,
                                      row.keepers=gn.nms.keepers,col.keepers=tf.nms)
score.combined.matrix <- score.combined.matrix[,which(apply(score.combined.matrix,2,sum)>0)]
sum.score <- apply(score.combined.matrix[,core.tfs],1,sum)
sam.score <- scores.sam[gn.nms.keepers]
score.combined.matrix <- cbind(score.combined.matrix,sum.score,sam.score)
write.table(score.combined.matrix,file=paste(path.output.in.input,"score_combine_KCR_matrix_zcut_",z.abs.cut,"_",type,"_",date.is,".xls",sep=""),sep="\t")
write.table(score.combined.matrix,file=paste(path.output.in.results,"score_combine_KCR_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),sep="\t")

cat("combining KC.\n")
# combine KC
score.combined.matrix <- combine_kcri(knockout=abs(rel.rank.ko),chip=abs(rel.rank.chip),rnaseq=0,immgen=0,
                                      row.keepers=gn.nms.keepers,col.keepers=tf.nms)
score.combined.matrix <- score.combined.matrix[,which(apply(score.combined.matrix,2,sum)>0)]
sum.score <- apply(score.combined.matrix[,core.tfs],1,sum)
sam.score <- scores.sam[gn.nms.keepers]
score.combined.matrix <- cbind(score.combined.matrix,sum.score,sam.score)
write.table(score.combined.matrix,file=paste(path.output.in.input,"score_combine_KC_matrix_zcut_",z.abs.cut,"_",type,"_",date.is,".xls",sep=""),sep="\t")
write.table(score.combined.matrix,file=paste(path.output.in.results,"score_combine_KC_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),sep="\t")


cat("combining RI.\n")
# combine RI
score.combined.matrix <- combine_kcri(knockout=0,chip=0,rnaseq=abs(rel.rank.rnaseq),immgen=abs(rel.rank.immgen),
                                      row.keepers=gn.nms.keepers,col.keepers=tf.nms)
score.combined.matrix <- score.combined.matrix[,which(apply(score.combined.matrix,2,sum)>0)]
sum.score <- apply(score.combined.matrix[,core.tfs],1,sum)
sam.score <- scores.sam[gn.nms.keepers]
score.combined.matrix <- cbind(score.combined.matrix,sum.score,sam.score)
write.table(score.combined.matrix,file=paste(path.output.in.input,"score_combine_RI_matrix_zcut_",z.abs.cut,"_",type,"_",date.is,".xls",sep=""),sep="\t")
write.table(score.combined.matrix,file=paste(path.output.in.results,"score_combine_RI_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),sep="\t")

cat("output scores for each dataset alone: K,C,R,I.\n")
# knockout
score.combined.matrix <- combine_kcri(knockout=abs(rel.rank.ko),chip=0,rnaseq=0,immgen=0,
                                      row.keepers=gn.nms.keepers,col.keepers=tf.nms)
score.combined.matrix <- score.combined.matrix[,which(apply(score.combined.matrix,2,sum)>0)]
sum.score <- apply(score.combined.matrix[,core.tfs],1,sum)
sam.score <- scores.sam[gn.nms.keepers]
score.combined.matrix <- cbind(score.combined.matrix,sum.score,sam.score)
write.table(score.combined.matrix,file=paste(path.output.in.input,"score_combine_K_matrix_zcut_",z.abs.cut,"_",type,"_",date.is,".xls",sep=""),sep="\t")
write.table(score.combined.matrix,file=paste(path.output.in.results,"score_combine_K_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),sep="\t")
# chip
score.combined.matrix <- combine_kcri(knockout=0,chip=abs(rel.rank.chip),rnaseq=0,immgen=0,
                                      row.keepers=gn.nms.keepers,col.keepers=tf.nms)
score.combined.matrix <- score.combined.matrix[,which(apply(score.combined.matrix,2,sum)>0)]
sum.score <- apply(score.combined.matrix[,core.tfs],1,sum)
sam.score <- scores.sam[gn.nms.keepers]
score.combined.matrix <- cbind(score.combined.matrix,sum.score,sam.score)
write.table(score.combined.matrix,file=paste(path.output.in.input,"score_combine_C_matrix_zcut_",z.abs.cut,"_",type,"_",date.is,".xls",sep=""),sep="\t")
write.table(score.combined.matrix,file=paste(path.output.in.results,"score_combine_C_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),sep="\t")

# immgen
score.combined.matrix <- combine_kcri(knockout=0,chip=0,rnaseq=0,immgen=abs(rel.rank.immgen),
                                      row.keepers=gn.nms.keepers,col.keepers=tf.nms)
score.combined.matrix <- score.combined.matrix[,which(apply(score.combined.matrix,2,sum)>0)]
sum.score <- apply(score.combined.matrix[,core.tfs],1,sum)
sam.score <- scores.sam[gn.nms.keepers]
score.combined.matrix <- cbind(score.combined.matrix,sum.score,sam.score)
write.table(score.combined.matrix,file=paste(path.output.in.input,"score_combine_I_matrix_zcut_",z.abs.cut,"_",type,"_",date.is,".xls",sep=""),sep="\t")
write.table(score.combined.matrix,file=paste(path.output.in.results,"score_combine_I_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),sep="\t")
# rnaseq
score.combined.matrix <- combine_kcri(knockout=0,chip=0,rnaseq=abs(rel.rank.rnaseq),immgen=0,
                                      row.keepers=gn.nms.keepers,col.keepers=tf.nms)
score.combined.matrix <- score.combined.matrix[,which(apply(score.combined.matrix,2,sum)>0)]
sum.score <- apply(score.combined.matrix[,core.tfs],1,sum)
sam.score <- scores.sam[gn.nms.keepers]
score.combined.matrix <- cbind(score.combined.matrix,sum.score,sam.score)
write.table(score.combined.matrix,file=paste(path.output.in.input,"score_combine_R_matrix_zcut_",z.abs.cut,"_",type,"_",date.is,".xls",sep=""),sep="\t")
write.table(score.combined.matrix,file=paste(path.output.in.results,"score_combine_R_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),sep="\t")
# SAM
write.table(rel.rank.sam,file=paste(path.output.in.input,"score_combine_S_matrix_zcut_",z.abs.cut,"_",type,"_",date.is,".xls",sep=""),sep="\t")
write.table(rel.rank.sam,file=paste(path.output.in.results,"score_combine_S_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),sep="\t")
