##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

## R code to create specificity scores to suggest what regulatory events are specific to each lineage th17/th1/th2/itreg compared to th0
## output is:
##   specificity matrix with rows as genes and columns specificity scores for each lineage

## get diff expression function btwn two single replicate 
source("r_scripts/th17/used_for_paper/util.R")

## define file names
## file.specificity.cmprd.to.plate.noa <- paste(sep="","specificity_cmprd_to_plate_",date.is,".noa")
## file.specificity.cmprd.to.lineages.noa <- paste(sep="","specificity_cmprd_to_lineages_",date.is,".noa")

## get expt names
expt.names <- colnames(d.log)

## names of relevant condition in d
expt <- list()
expt[["th0.48hr"]] <- expt.names[grep(paste("th0","48hr","wt",sep=".*"),expt.names,ignore.case=T)] # our control for insilico conditions
expt[["th0.exvivo"]] <- "SL1859_th0_exvivo_ss_#1" # our control for exvivo condition of th17
expt[["th17.48hr"]] <- expt.names[grep(paste("th17","48hr","wt",sep=".*"),expt.names,ignore.case=TRUE)] # to find in silico th17 diff genes (11 reps)
expt[["th1.48hr"]] <- expt.names[grep(paste("th1_","48hr",sep=".*"),expt.names,ignore.case=TRUE)] # to find th1 diff genes and compare to th17 in silico diff genes (2 reps)
expt[["th2.48hr"]] <- expt.names[grep(paste("th2_","48hr",sep=".*"),expt.names,ignore.case=TRUE)] # to find th2 diff genes and compare to th17 in silico diff genes (2 reps)
expt[["itreg.48hr"]] <- expt.names[grep(paste("itreg","48hr",sep=".*"),expt.names,ignore.case=TRUE)] # to find itreg diff genes and compare to th17 in silico diff genes (2 reps)
expt[["th17.exvivo"]] <- "SL4083_IL17a_gfp_plus_exvivo_SI_ss_#1" # to find exvivo th17 diff genes (and compare to th0 exvivo diff genes)
expt[["itreg.exvivo"]] <- "SL4082_foxp3_gfp_plus_exvivo_spl_ln_ss_#1" # to find exvivo itreg diff genes (and compare to th0 exvivo diff genes)
expt[["gut.non.th17.exvivo"]] <- "SL3539_IL17a_gfp_minus_exvivo_SI_ss_#1" # to find exvivo gut (non th17) diff genes (and compare to th0 exvivo diff genes)


# put results in this matrix
x <- names(expt)[-c(1:2)]
m <- matrix(,nr=dim(d.log)[1],nc=length(x))
rownames(m) <- rownames(d.log)
colnames(m) <- x 
rm(x)
## get median and sd for each gene under th0.48hr condition
th0.48hr.mean <- apply(d.log[,expt[["th0.48hr"]] ],1,mean)
th0.48hr.sd <- apply(d.log[,expt[["th0.48hr"]] ],1,sd)

## diff expression with th0.exvivo or th0.48hr
for(i in 1:dim(m)[2]) {
  cond <- colnames(m)[i]
  zscores <- numeric()
  is.exvivo <- as.logical(length(grep("exvivo",cond)))
  if(is.exvivo==TRUE){
    for(j in 1:length(expt[[cond]])){
      zscores <- cbind(zscores, diff.exprs.analysis(wt=d[, expt[[cond]][j] ],ko=d[,expt[["th0.exvivo"]] ],c1=5,ps.cnt=1))
    }
  } else {
    if(cond=="th17.48hr"){
      zscores <- cbind(zscores, apply(d.log,1,function(i) t.test(i[expt[[cond]]],i[expt[["th0.48hr"]]])$statistic))
    } else {
      for(j in 1:length(expt[[cond]])){
        zscores <- cbind(zscores, (d.log[, expt[[cond]][j] ]-th0.48hr.mean)/th0.48hr.sd)
      }
    }
  }
  median.zscores <- apply(zscores,1,median)
  m[names(median.zscores),i] <- median.zscores
}
colnames(m) <- paste(sep=".",colnames(m),c("vs.th0[z]","vs.th0[z]","vs.th0[z]","vs.th0[z]","vs.th0.exvivo[pz]","vs.th0.exvivo[pz]","vs.th0.exvivo[pz]"))
## add differential expression analysis between th17.48hr invitro to th17.48hr exvivo
## ix.invitro <- which(apply(d.log[,expt[["th17.48hr"]] ],1,median)> 2.321928) ## consider genes with > 2.32 log2(rpkm) i.e. bigger than 5 rpkms.
## ix.exvivo <- which(d.log[,expt[["th17.exvivo"]] ]> 2.321928) ## consider genes with > 2.32 log2(rpkm) i.e. bigger than 5 rpkms.
## ix <- intersect(ix.invitro,ix.exvivo)
## th17.48hr.mean <- apply(d.log[,expt[["th17.48hr"]] ],1,mean)
## th17.48hr.sd <- apply(d.log[,expt[["th17.48hr"]] ],1,sd)
## zscores <- (d.log[ix, expt[["th17.exvivo"]][j] ]-th17.48hr.mean[ix])/th17.48hr.sd[ix]
m.add <- matrix(,nr=dim(m)[1],nc=3)
rownames(m.add) <- rownames(m)
colnames(m.add) <- c("th17.exvivo.vs.gut.exvivo[pz]","th17.exvivo.vs.itreg.exvivo[pz]","th17.exvivo.vs.th17.invitro[pz]")
## add diff expression for th17.exvivo vs. gut.non.th17.exvivo
x1 <- diff.exprs.analysis(wt=d[, expt[["th17.exvivo"]] ],ko=d[,expt[["gut.non.th17.exvivo"]] ],c1=5,ps.cnt=1)
## add diff expression for th17.exvivo vs. itreg.exvivo
x2 <- diff.exprs.analysis(wt=d[, expt[["th17.exvivo"]] ],ko=d[,expt[["itreg.exvivo"]] ],c1=5,ps.cnt=1)
th17.48hr.median <- apply(d[,expt[["th17.48hr"]] ],1,median)
x3 <- diff.exprs.analysis(wt=d[, expt[["th17.exvivo"]] ],ko=th17.48hr.median,c1=5,ps.cnt=1)
m.add[names(x1),1] <- x1
m.add[names(x2),2] <- x2
m.add[names(x3),3] <- x3
m <- cbind(m,m.add)

## get rid of na's that were put in by median.zscores.per.tf[names(median.zscores),i] <- median.zscores for genes not in names(median.zscores)
m[which( is.na(m))] = 0

specificity.mat <- m
rm(m)
