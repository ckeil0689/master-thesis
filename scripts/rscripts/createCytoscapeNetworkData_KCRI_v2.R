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
# R code to create cyto-nets (after running combineData.R)
# output is:
#   core.net syf file
#   core.net eda file
#   cortex.net syf file
#   cortex.net eda file
#   dataset tab delim file
#   tf.scoresum.per.target: measure of how well supported was each gene by all tfs
#   tf.score.per.target.per.tf (from target.per.tf list): measure of how well supported was each gene by each tfs
#   isTf (tfs.in.cortex.net) : a list of genes with 1=tf 0=non.tf

## paths
path.input <- "input/th17/used_for_paper/"

## path for output
path.output <- paste(sep="","results/usedForPaper/cytoscape_",date.is,"/")
system(paste(sep="","mkdir ",path.output))

## define input file names
ri.f.nm <- paste(path.input,"combinedAnalysis/score_combine_RI_matrix_zcut_",z.abs.cut,"_","k_r_i_activator_",date.combine.data.run,".xls",sep="")
kc.f.nm <- paste(path.input,"combinedAnalysis/score_combine_KC_matrix_zcut_",z.abs.cut,"_","k_r_i_activator_",date.combine.data.run,".xls",sep="")
kcr.f.nm <- paste(path.input,"combinedAnalysis/score_combine_KCR_matrix_zcut_",z.abs.cut,"_","k_r_i_activator_",date.combine.data.run,".xls",sep="")
kcri.f.nm <- paste(path.input,"combinedAnalysis/score_combine_KCRI_matrix_zcut_",z.abs.cut,"_","k_r_i_activator_",date.combine.data.run,".xls",sep="")
ri.r.f.nm <- paste(path.input,"combinedAnalysis/score_combine_RI_matrix_zcut_",z.abs.cut,"_","k_r_i_repressor_",date.combine.data.run,".xls",sep="")
kc.r.f.nm <- paste(path.input,"combinedAnalysis/score_combine_KC_matrix_zcut_",z.abs.cut,"_","k_r_i_repressor_",date.combine.data.run,".xls",sep="")
kcr.r.f.nm <- paste(path.input,"combinedAnalysis/score_combine_KCR_matrix_zcut_",z.abs.cut,"_","k_r_i_repressor_",date.combine.data.run,".xls",sep="")
kcri.r.f.nm <- paste(path.input,"combinedAnalysis/score_combine_KCRI_matrix_zcut_",z.abs.cut,"_","k_r_i_repressor_",date.combine.data.run,".xls",sep="")
immgen.f.nm <- paste(sep="",path.input.immgen.inf.data,"ratios.RData")
rnaseq.f.nm <- paste(sep="",path.input.rnaseq.inf.data,"ratios.RData")
ko.f.nm <- paste(sep="","DEseq_pval_signed_",date.deseq.run,".xls")
specificity.f.nm <- paste(sep="",path.input,"specificity_mat_",date.deseq.run,".xls")

## define ouput file names
# for dataset KCRI
file.core.activate.kcri.eda <- paste(sep="",path.output,"kcri_core_net_activ_",date.is,".eda")
file.core.activate.kcri.sif <- paste(sep="",path.output,"kcri_core_net_activ_",date.is,".sif")
file.core.repress.kcri.eda <- paste(sep="",path.output,"kcri_core_net_repres_",date.is,".eda")
file.core.repress.kcri.sif <- paste(sep="",path.output,"kcri_core_net_repres_",date.is,".sif")
file.extend.activate.kcri.eda <- paste(sep="",path.output,"kcri_extend_net_activ_",date.is,".eda")
file.extend.activate.kcri.sif <- paste(sep="",path.output,"kcri_extend_net_activ_",date.is,".sif")
file.extend.repress.kcri.eda <- paste(sep="",path.output,"kcri_extend_net_repres_",date.is,".eda")
file.extend.repress.kcri.sif <- paste(sep="",path.output,"kcri_extend_net_repres_",date.is,".sif")
file.core.abs.kcri.eda <- paste(sep="",path.output,"kcri_core_net_abs_",date.is,".eda")
file.core.abs.kcri.sif <- paste(sep="",path.output,"kcri_core_net_abs_",date.is,".sif")
file.extend.abs.kcri.eda <- paste(sep="",path.output,"kcri_extend_net_abs_",date.is,".eda")
file.extend.abs.kcri.sif <- paste(sep="",path.output,"kcri_extend_net_abs_",date.is,".sif")
# for dataset KCR
file.core.activate.kcr.eda <- paste(sep="",path.output,"kcr_core_net_activ_",date.is,".eda")
file.core.activate.kcr.sif <- paste(sep="",path.output,"kcr_core_net_activ_",date.is,".sif")
file.core.repress.kcr.eda <- paste(sep="",path.output,"kcr_core_net_repres_",date.is,".eda")
file.core.repress.kcr.sif <- paste(sep="",path.output,"kcr_core_net_repres_",date.is,".sif")
file.extend.activate.kcr.eda <- paste(sep="",path.output,"kcr_extend_net_activ_",date.is,".eda")
file.extend.activate.kcr.sif <- paste(sep="",path.output,"kcr_extend_net_activ_",date.is,".sif")
file.extend.repress.kcr.eda <- paste(sep="",path.output,"kcr_extend_net_repres_",date.is,".eda")
file.extend.repress.kcr.sif <- paste(sep="",path.output,"kcr_extend_net_repres_",date.is,".sif")
file.core.abs.kcr.eda <- paste(sep="",path.output,"kcr_core_net_abs_",date.is,".eda")
file.core.abs.kcr.sif <- paste(sep="",path.output,"kcr_core_net_abs_",date.is,".sif")
file.extend.abs.kcr.eda <- paste(sep="",path.output,"kcr_extend_net_abs_",date.is,".eda")
file.extend.abs.kcr.sif <- paste(sep="",path.output,"kcr_extend_net_abs_",date.is,".sif")
# for dataset KC
file.core.activate.kc.eda <- paste(sep="",path.output,"kc_core_net_activ_",date.is,".eda")
file.core.activate.kc.sif <- paste(sep="",path.output,"kc_core_net_activ_",date.is,".sif")
file.core.repress.kc.eda <- paste(sep="",path.output,"kc_core_net_repres_",date.is,".eda")
file.core.repress.kc.sif <- paste(sep="",path.output,"kc_core_net_repres_",date.is,".sif")
file.extend.activate.kc.eda <- paste(sep="",path.output,"kc_extend_net_activ_",date.is,".eda")
file.extend.activate.kc.sif <- paste(sep="",path.output,"kc_extend_net_activ_",date.is,".sif")
file.extend.repress.kc.eda <- paste(sep="",path.output,"kc_extend_net_repres_",date.is,".eda")
file.extend.repress.kc.sif <- paste(sep="",path.output,"kc_extend_net_repres_",date.is,".sif")
file.core.abs.kc.eda <- paste(sep="",path.output,"kc_core_net_abs_",date.is,".eda")
file.core.abs.kc.sif <- paste(sep="",path.output,"kc_core_net_abs_",date.is,".sif")
file.extend.abs.kc.eda <- paste(sep="",path.output,"kc_extend_net_abs_",date.is,".eda")
file.extend.abs.kc.sif <- paste(sep="",path.output,"kc_extend_net_abs_",date.is,".sif")

file.node.annot.table <- paste(sep="","node_annot_table",date.is,".noa")

## get scores for interactions
ri <- as.matrix(read.table(ri.f.nm,header=T,sep="\t"))
kc <- as.matrix(read.table(kc.f.nm,header=T,sep="\t"))
kc.r <- as.matrix(read.table(kc.r.f.nm,header=T,sep="\t"))
kcr <- as.matrix(read.table(kcr.f.nm,header=T,sep="\t"))
kcri <- as.matrix(read.table(kcri.f.nm,header=T,sep="\t"))
kcr.r <- as.matrix(read.table(kcr.r.f.nm,header=T,sep="\t"))
kcri.r <- as.matrix(read.table(kcri.r.f.nm,header=T,sep="\t"))
# keep SAM scores
sam.score <- kcri[,"sam.score"]

# keep sum scores
kc.sum.score <- kc[,"sum.score"]
kc.r.sum.score <- kc.r[,"sum.score"]
kcr.sum.score <- kcr[,"sum.score"]
kcri.sum.score <- kcri[,"sum.score"]
kcr.r.sum.score <- kcr.r[,"sum.score"]
kcri.r.sum.score <- kcri.r[,"sum.score"]

# remove sum and SAM scores
n <- ncol(kc)
kc <- kc[,-c((n-1),n)]
n <- ncol(kc.r)
kc.r <- kc.r[,-c((n-1),n)]
n <- ncol(kcr)
kcr <- kcr[,-c((n-1),n)]
n <- ncol(kcr.r)
kcr.r <- kcr.r[,-c((n-1),n)]
n <- ncol(kcri)
kcri <- kcri[,-c((n-1),n)]
n <- ncol(kcri.r)
kcri.r <- kcri.r[,-c((n-1),n)]

#33333333333333333333333333333
#create sign matrix for kcri
#33333333333333333333333333333

# get gene and tf names
gns <- rownames(ri)
tfs.kc <- colnames(kc)
tfs.kcr <- colnames(kcr)
tfs.kcri <- colnames(kcri)

## get the sign (based on correlation between tfs.kcri and genes in both rnaseq and immgen dataset)
m.i <- matrix(0,nc=length(tfs.kcri),nr=nrow(kcri),dimnames=dimnames(kcri))
m.r <- matrix(0,nc=length(tfs.kcri),nr=nrow(kcri),dimnames=dimnames(kcri))
m.k <- matrix(0,nc=length(tfs.kcri),nr=nrow(kcri),dimnames=dimnames(kcri))
m.sign.kcri <- matrix(0,nc=length(tfs.kcri),nr=nrow(kcri),dimnames=dimnames(kcri))
# load datasets (immgen)
load(immgen.f.nm)
tf.immgen <- tfs.kcri[which(tfs.kcri %in% rownames(ratios))]
gn.immgen <- gns[which(gns %in% rownames(ratios))]
# get sign matrix for immgen
m.i[gn.immgen,tf.immgen] <- cor(t(ratios))[gn.immgen,tf.immgen]
# load datasets (rnaseq)
load(rnaseq.f.nm)
tf.rnaseq <- tfs.kcri[which(tfs.kcri %in% rownames(ratios))]
gn.rnaseq <- gns[which(gns %in% rownames(ratios))]
# get sign matrix for rnaseq
m.r[gn.rnaseq,tf.rnaseq] <- cor(t(ratios))[gn.rnaseq,tf.rnaseq]
# get sign matrix for five core tfs.kcri from the KO data
f.nm <- paste(sep="",path.input.ko,ko.f.nm)
cat("loading KO scores file: ",f.nm,"\n")
to.load <- c("Th17.batf.wt.vs.Th17.batf.ko","Th17.maf.wt.vs.Th17.maf.ko","Th17.irf4.wt.vs.Th17.irf4.ko","Th17.stat3.wt.vs.Th17.stat3.ko","Th17.rorc.wt.vs.Th17.rorc.ko")
scores.ko <- as.matrix(read.table(f.nm,sep="\t"))[,to.load]
colnames(scores.ko) <- toupper(sapply(strsplit(colnames(scores.ko),"\\."),function(i) i[2]))
gn.ko <- gns[which(gns %in% rownames(scores.ko))]
m.k <- scores.ko[gn.ko,core.tfs]

#33333333333333333333333333333
#create sign matrix for kcri
#33333333333333333333333333333
# make final sign matrix for each interactions as follows:
# for all genes other than core tf take the sign of the correlation over either rnaseq or immgen data (the sign of whichever has a larger abs(cor) value)
# for all core tfs take the sign from the KO data
# store final sign for each interactrion in m.sign.kcri (same dims as kcri)
for(i in 1:length(tfs.kcri)){
  tf <- tfs.kcri[i]
  v.i <- m.i[,tf]
  v.r <- m.r[,tf]
  ix.i <- which((abs(v.i)-abs(v.r))>=0)
  ix.r <- which((abs(v.i)-abs(v.r))<0)
  m.sign.kcri[ix.i,tf] <- v.i[ix.i]
  m.sign.kcri[ix.r,tf] <- v.r[ix.r]
}
m.sign.kcri[rownames(m.k),colnames(m.k)] <- m.k
m.sign.kcri <- sign(m.sign.kcri)
## some interaction may have no support by k r or i but may have support from Chip c
# set those interaction to have a default sign of 1
m.sign.kcri[which(m.sign.kcri==0)] <- 1
#33333333333333333333333333333
#create sign matrix for kcr
#33333333333333333333333333333
m.sign.kcr <- matrix(0,nc=length(tfs.kcr),nr=nrow(kcr),dimnames=dimnames(kcr))
m.sign.kcr <- m.r[,colnames(m.sign.kcr)]
m.sign.kcr[rownames(m.k),colnames(m.k)] <- m.k
m.sign.kcr <- sign(m.sign.kcr)
m.sign.kcr[which(m.sign.kcr==0)] <- 1
#33333333333333333333333333333
#create sign matrix for kc
#33333333333333333333333333333
m.sign.kc <- matrix(0,nc=length(tfs.kc),nr=nrow(kc),dimnames=dimnames(kc))
m.sign.kc[rownames(m.k),colnames(m.k)] <- m.k
m.sign.kc <- sign(m.sign.kc)
m.sign.kc[which(m.sign.kc==0)] <- 1


# get rnaseq mRNA levels for ko vs wt steady state experiments and th0 -> th17 time series experiment and put in d
load(paste(path.input,"ranseqDatasetNoQuartNorm.RData",sep=""))
rownames(rnaseq.complete.matrix) <- toupper(rownames(rnaseq.complete.matrix))
d <- rnaseq.complete.matrix
d.log <- log2(d+1)
rm(rnaseq.complete.matrix)


# add sign for activation/repression into score matrices
kcri.signed <- m.sign.kcri*kcri
kcri.r.signed <- -1*m.sign.kcri*kcri.r #(I want positive scores to mean repression here)
kcr.signed <- m.sign.kcr*kcr
kcr.r.signed <- -1*m.sign.kcr*kcr.r #(I want positive scores to mean repression here)
kc.signed <- m.sign.kc*kc
kc.r.signed <- -1*m.sign.kc*kc.r #(I want positive scores to mean repression here)

############ start writing sif and eda output files ####################
#33333333333333333333333333333
# FIRST BASED ON KCRI
#33333333333333333333333333333
m.cut.core.active <- quantile(kcri.signed[,core.tfs],probs=cut.kcri.core.activ.qunt)
m.cut.core.repress <- m.cut.core.active
m.cut.extend.active <- quantile(kcri.signed[,-which(colnames(kcri.signed) %in% core.tfs)],probs=cut.kcri.extend.activ.qunt)
m.cut.extend.repress <- m.cut.extend.active

pos.edge.core <- "kcri_core_active"
neg.edge.core <- "kcri_core_repress"
pos.edge.extend <- "kcri_extend_active"
neg.edge.extend <- "kcri_extend_repress"

##%%%%%%%% CORE NET ACTIVE %%%%%%%%##
f.eda <- file.core.activate.kcri.eda
f.sif <- file.core.activate.kcri.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
m <- kcri.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)

##%%%%%%%% CORE NET REPRESS %%%%%%%%##
f.eda <- file.core.repress.kcri.eda
f.sif <- file.core.repress.kcri.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
m <- kcri.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=F)

##%%%%%%%% EXTENDED NET ACTIVATE %%%%%%%%##
f.eda <- file.extend.activate.kcri.eda
f.sif <- file.extend.activate.kcri.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
## first build the core
m <- kcri.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)
## now extend for non core tfs based on RI scores (appending to file)
m <- kcri.signed[,-which(colnames(kcri.signed) %in% core.tfs)]
m.cut <- m.cut.extend.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.extend,neg.edge=neg.edge.extend,append=T)

##%%%%%%%% EXTENDED NET REPRESS %%%%%%%%##
f.eda <- file.extend.repress.kcri.eda
f.sif <- file.extend.repress.kcri.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
## first build the core
m <- kcri.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=F)
## now extend for now core tfs based on RI scores (appending to file)
m <- kcri.signed[,-which(colnames(kcri.signed) %in% core.tfs)]
m.cut <- m.cut.extend.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.extend,neg.edge=pos.edge.extend,append=T)

##%%%%%%%% COMBINED NET CORE ACTIVATE AND REPRESSED %%%%%%%%##
# first core active
f.eda <- file.core.abs.kcri.eda
f.sif <- file.core.abs.kcri.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
## build network for core tfs based on KCRI
m <- kcri.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)
# now append core repressed
m <- kcri.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=T)

##%%%%%%%% COMBINED NET EXTENDED ACTIVATE AND REPRESSED %%%%%%%%##
# first extended active
f.eda <- file.extend.abs.kcri.eda
f.sif <- file.extend.abs.kcri.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
m <- kcri.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)
## now extend for now core tfs based on RI scores (appending to file)
m <- kcri.signed[,-which(colnames(kcri.signed) %in% core.tfs)]
m.cut <- m.cut.extend.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.extend,neg.edge=neg.edge.extend,append=T)
# now extended repressed
## first build the core
m <- kcri.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=T)
## now extend for now core tfs based on RI scores (appending to file)
m <- kcri.signed[,-which(colnames(kcri.signed) %in% core.tfs)]
m.cut <- m.cut.extend.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.extend,neg.edge=pos.edge.extend,append=T)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############ writing sif and eda output files ####################
#33333333333333333333333333333
# BASED ON KCR
#33333333333333333333333333333
m.cut.core.active <- quantile(kcr.signed[,core.tfs],probs=cut.kcr.core.activ.qunt)
m.cut.core.repress <- m.cut.core.active
m.cut.extend.active <- quantile(kcr.signed[,-which(colnames(kcr.signed) %in% core.tfs)],probs=cut.kcr.extend.activ.qunt)
m.cut.extend.repress <- m.cut.extend.active

pos.edge.core <- "kcr_core_active"
neg.edge.core <- "kcr_core_repress"
pos.edge.extend <- "kcr_extend_active"
neg.edge.extend <- "kcr_extend_repress"

##%%%%%%%% CORE NET ACTIVE %%%%%%%%##
f.eda <- file.core.activate.kcr.eda
f.sif <- file.core.activate.kcr.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
m <- kcr.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)

##%%%%%%%% CORE NET REPRESS %%%%%%%%##
f.eda <- file.core.repress.kcr.eda
f.sif <- file.core.repress.kcr.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
m <- kcr.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=F)

##%%%%%%%% EXTENDED NET ACTIVATE %%%%%%%%##
f.eda <- file.extend.activate.kcr.eda
f.sif <- file.extend.activate.kcr.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
## first build the core
m <- kcr.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)
## now extend for non core tfs based on RI scores (appending to file)
m <- kcr.signed[,-which(colnames(kcr.signed) %in% core.tfs)]
m.cut <- m.cut.extend.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.extend,neg.edge=neg.edge.extend,append=T)

##%%%%%%%% EXTENDED NET REPRESS %%%%%%%%##
f.eda <- file.extend.repress.kcr.eda
f.sif <- file.extend.repress.kcr.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
## first build the core
m <- kcr.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=F)
## now extend for now core tfs based on RI scores (appending to file)
m <- kcr.signed[,-which(colnames(kcr.signed) %in% core.tfs)]
m.cut <- m.cut.extend.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.extend,neg.edge=pos.edge.extend,append=T)

##%%%%%%%% COMBINED NET CORE ACTIVATE AND REPRESSED %%%%%%%%##
# first core active
f.eda <- file.core.abs.kcr.eda
f.sif <- file.core.abs.kcr.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
## build network for core tfs based on kcr
m <- kcr.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)
# now append core repressed
m <- kcr.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=T)

##%%%%%%%% COMBINED NET EXTENDED ACTIVATE AND REPRESSED %%%%%%%%##
# first extended active
f.eda <- file.extend.abs.kcr.eda
f.sif <- file.extend.abs.kcr.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
m <- kcr.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)
## now extend for now core tfs based on RI scores (appending to file)
m <- kcr.signed[,-which(colnames(kcr.signed) %in% core.tfs)]
m.cut <- m.cut.extend.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.extend,neg.edge=neg.edge.extend,append=T)
# now extended repressed
## first build the core
m <- kcr.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=T)
## now extend for now core tfs based on RI scores (appending to file)
m <- kcr.signed[,-which(colnames(kcr.signed) %in% core.tfs)]
m.cut <- m.cut.extend.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.extend,neg.edge=pos.edge.extend,append=T)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#33333333333333333333333333333
# NOW BASED ON KC
#33333333333333333333333333333
m.cut.core.active <- quantile(kc.signed[,core.tfs],probs=cut.kc.core.activ.qunt)
m.cut.core.repress <- m.cut.core.active
m.cut.extend.active <- quantile(kc.signed[,-which(colnames(kc.signed) %in% core.tfs)],probs=cut.kc.extend.activ.qunt)
m.cut.extend.repress <- m.cut.extend.active

pos.edge.core <- "kc_core_active"
neg.edge.core <- "kc_core_repress"
pos.edge.extend <- "kc_extend_active"
neg.edge.extend <- "kc_extend_repress"

##%%%%%%%% CORE NET ACTIVE %%%%%%%%##
f.eda <- file.core.activate.kc.eda
f.sif <- file.core.activate.kc.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
m <- kc.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)

##%%%%%%%% CORE NET REPRESS %%%%%%%%##
f.eda <- file.core.repress.kc.eda
f.sif <- file.core.repress.kc.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
m <- kc.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=F)

##%%%%%%%% EXTENDED NET ACTIVATE %%%%%%%%##
f.eda <- file.extend.activate.kc.eda
f.sif <- file.extend.activate.kc.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
## first build the core
m <- kc.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)
## now extend for non core tfs based on RI scores (appending to file)
m <- kc.signed[,-which(colnames(kc.signed) %in% core.tfs)]
m.cut <- m.cut.extend.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.extend,neg.edge=neg.edge.extend,append=T)

##%%%%%%%% EXTENDED NET REPRESS %%%%%%%%##
f.eda <- file.extend.repress.kc.eda
f.sif <- file.extend.repress.kc.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
## first build the core
m <- kc.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=F)
## now extend for now core tfs based on RI scores (appending to file)
m <- kc.signed[,-which(colnames(kc.signed) %in% core.tfs)]
m.cut <- m.cut.extend.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.extend,neg.edge=pos.edge.extend,append=T)

##%%%%%%%% COMBINED NET CORE ACTIVATE AND REPRESSED %%%%%%%%##
# first core active
f.eda <- file.core.abs.kc.eda
f.sif <- file.core.abs.kc.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
## build network for core tfs based on kc
m <- kc.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)
# now append core repressed
m <- kc.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=T)

##%%%%%%%% COMBINED NET EXTENDED ACTIVATE AND REPRESSED %%%%%%%%##
# first extended active
f.eda <- file.extend.abs.kc.eda
f.sif <- file.extend.abs.kc.sif
system(paste(sep="","rm ",f.eda))  
system(paste(sep="","rm ",f.sif))
m <- kc.signed[,core.tfs]
m.cut <- m.cut.core.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.core,neg.edge=neg.edge.core,append=F)
## now extend for now core tfs based on RI scores (appending to file)
m <- kc.signed[,-which(colnames(kc.signed) %in% core.tfs)]
m.cut <- m.cut.extend.active
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=pos.edge.extend,neg.edge=neg.edge.extend,append=T)
# now extended repressed
## first build the core
m <- kc.r.signed[,core.tfs]
m.cut <- m.cut.core.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.core,neg.edge=pos.edge.core,append=T)
## now extend for now core tfs based on RI scores (appending to file)
m <- kc.signed[,-which(colnames(kc.signed) %in% core.tfs)]
m.cut <- m.cut.extend.repress
write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,pos.edge=neg.edge.extend,neg.edge=pos.edge.extend,append=T)


#33333333333333333333333333333333
# Now write additional info files
#33333333333333333333333333333333

## adding cell surface markers taken as a union of gene names in our dataset from three lists taken from human biomart, mouse biomart and GO mouse for
## external side of plasma membrane GO:0009897
srfc.gns.go <- toupper(unique(as.character(read.table(paste(sep="",path.input,"cytoscape/plsma_mmbrn_GO_0009897_Feb_21_2011.txt"),header=F,sep="\t")$V1)))
srfc.gns.biomart.mm <- toupper(unique(as.character(read.table(paste(sep="",path.input,"cytoscape/mart_export_plsma_mmbrn_mm_Feb_21_2011.txt"),header=F,sep="\t")$V1)))
srfc.gns.biomart.hg <- toupper(unique(as.character(read.table(paste(sep="",path.input,"cytoscape/mart_export_plsma_mmbrn_hg_Feb_21_2011.txt"),header=F,sep="\t")$V1)))
srfc.gns <- union(c(srfc.gns.biomart.mm,srfc.gns.biomart.hg),srfc.gns.go)
## catalytic activity GO:0003824
enzyme.gns <- unique(toupper(unique(as.character(read.table(paste(sep="",path.input,"cytoscape/catalytic_activity_GO_0003824_Mar_12_2011.txt"),header=F,sep="\t")$V1))))
## protein kinase activity GO:0004672
kinase.gns <- unique(toupper(unique(as.character(read.table(paste(sep="",path.input,"cytoscape/protein_kinase_activity_GO_0004672_Mar_12_2011.txt"),header=F,sep="\t")$V1))))

# tf=1, surface marker=2, kinase=3, enzyme=4, other=0 (if tf has priority than enzyme than srfc mrkr)
gn.function <- character(length=length(gns))
for (i in 1:length(gns)){
  if (gns[i] %in% tfs.kcri){
    gn.function[i] <- "tf"
  } else if (gns[i] %in% srfc.gns){
    gn.function[i] <- "srfc.mrkr"
  } else if (gns[i] %in% kinase.gns){
    gn.function[i] <- "kinase"
  } else if (gns[i] %in% enzyme.gns){
    gn.function[i] <- "enzyme"
  } else {
    gn.function[i] <- "other"
  }
}

## write kc.sum, show.gn.nm, and gn.function into sum scores
show.gn.nm <- numeric(length(gns))
names(show.gn.nm) <- gns
show.gn.nm[show.gn.names] <- 1
kc.sum.diff <- kc.sum.score-kc.r.sum.score
kcr.sum.diff <- kcr.sum.score-kcr.r.sum.score
kcri.sum.diff <- kcri.sum.score-kcri.r.sum.score
x <- cbind(kc.sum.score,kcr.sum.score,kcri.sum.score,
	kc.r.sum.score,kcr.r.sum.score,kcri.r.sum.score,
	kc.sum.diff,kcr.sum.diff,kcri.sum.diff,
	sam.score,gn.function,show.gn.nm)
fl.nm <- paste(path.output,"sum_scores_",date.is,".txt",sep="")
cat("gene_id",colnames(x),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(sep="\t",x,quote=FALSE,append=TRUE,col.names=FALSE,file=fl.nm)

## write differential expression levels
fl.nm <- paste(path.output,"diff_exprssn_zcut_",z.abs.cut,"_",date.is,".txt" ,sep="")
m.diff <- m.k
colnames(m.diff) <- paste(sep="","diff_exp_log10_pval_",colnames(m.diff))
ix.gns.not.in.m.diff <- which(!gns %in% rownames(m.diff)) # find gene names not in diff exp mat
tmp <- matrix(0,nr=length(ix.gns.not.in.m.diff),nc=dim(m.diff)[2])
rownames(tmp) <- gns[ix.gns.not.in.m.diff]
m.diff.inclusive <- rbind(m.diff,tmp) # add all gene names with 0 as diff exp values
cat("gene_id",colnames(m.diff),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(sep="\t",m.diff.inclusive[gns,],quote=FALSE,append=TRUE,col.names=FALSE,file=fl.nm)

## write specificity scores and save to variable specificity.mat
# load specificity matrix
f.nm <- specificity.f.nm
specificity.mat <- as.matrix(read.table(f.nm,sep="\t"))
fl.nm <- paste(path.output,"specificity_mat_zcut_",z.abs.cut,"_",date.is,".txt",sep="")
cat("gene_id",colnames(specificity.mat),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(specificity.mat,sep="\t",quote=FALSE,append=TRUE,col.names=FALSE,file=fl.nm)

## write all rnaseq data log2(rpkm)
fl.nm <- paste(path.output,"rnaseq_dataset_log2rpkm_zcut_",z.abs.cut,"_",date.is,".txt",sep="")
colnames(d.log) <- paste(sep="_",colnames(d.log),"log2_rpkm")
cat("gene_id",colnames(d.log),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(d.log[gns,],sep="\t",quote=FALSE,append=TRUE,col.names=FALSE,file=fl.nm)

## write all rnaseq data (rpkm)
fl.nm <- paste(path.output,"rnaseq_dataset_rpkm_zcut_",z.abs.cut,"_",date.is,".txt",sep="")
colnames(d) <- paste(sep="_",colnames(d),"rpkm")
cat("gene_id",colnames(d),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(d[gns,],sep="\t",quote=FALSE,append=TRUE,col.names=FALSE,file=fl.nm)
