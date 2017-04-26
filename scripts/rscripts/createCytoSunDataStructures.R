##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

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
path.output <- paste(sep="","results/usedForPaper/cytoscape_",date.is,"/")
system(paste(sep="","mkdir ",path.output))

## define file names
file.all.eda <- paste(sep="",path.output,"all_net_",date.is,".eda")
file.all.sif <- paste(sep="",path.output,"all_net_",date.is,".sif")
file.extend.eda <- paste(sep="",path.output,"extend_net_",date.is,".eda")
file.extend.sif <- paste(sep="",path.output,"extend_net_",date.is,".sif")
file.core.eda <- paste(sep="",path.output,"core_net_",date.is,".eda")
file.core.sif <- paste(sep="",path.output,"core_net_",date.is,".sif")
file.all.tfs.only.eda <- paste(sep="",path.output,"all_net_tfs_only_",date.is,".eda")
file.all.tfs.only.sif <- paste(sep="",path.output,"all_net_tfs_only_",date.is,".sif")
file.core.tfs.only.eda <- paste(sep="",path.output,"core_net_tfs_only_",date.is,".eda")
file.core.tfs.only.sif <- paste(sep="",path.output,"core_net_tfs_only_",date.is,".sif")



file.node.annot.table <- paste(sep="","node_annot_table",date.is,".noa")

## get scores for interactions
RI <- as.matrix(read.table(paste(path.input,"combinedAnalysis/score_combine_RI_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),header=T,sep="\t"))
KC <- as.matrix(read.table(paste(path.input,"combinedAnalysis/score_combine_KC_matrix_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),header=T,sep="\t"))

gns <- rownames(RI)
# get core net tf names
core.tfs <- c("BATF","MAF","IRF4","STAT3","RORC")

# get inferelator betas for interactions (a part of the inferelator confidence score above)
beta.mat.i <- as.matrix(read.table(paste(path.input,"infResults/immgen/inf_betas_zcut_",z.abs.cut,"_Nb_",num.boot.i,"_",date.is,".xls",sep=""),header=T,sep="\t"))
beta.mat.r <- as.matrix(read.table(paste(path.input,"infResults/rnaseq/inf_betas_zcut_",z.abs.cut,"_Nb_",num.boot.r,"_",date.is,".xls",sep=""),header=T,sep="\t"))

# get diff expression zscores for each tf
tfs.diff.mat <- as.matrix(read.table(paste(path.input,"rnaseqDiffExpMedianZscoresPerTf.xls",sep=""),header=T,sep="\t"))
tfs.diff.mat <- tfs.diff.mat[,-grep("rora_and_rorc",colnames(tfs.diff.mat))]
colnames(tfs.diff.mat) <- core.tfs 

# get rnaseq mRNA levels for ko vs wt steady state experiments and th0 -> th17 time series experiment and put in d
load(paste(path.input,"ranseqDatasetNoQuartNorm.RData",sep=""))
rownames(rnaseq.complete.matrix) <- toupper(rownames(rnaseq.complete.matrix))
d <- rnaseq.complete.matrix
d.log <- log2(d+1)
rm(rnaseq.complete.matrix)

# get diff expression (SAM) between th17 and th0 (both at 48 hr)
sam <- as.matrix(read.table(paste(path.input,"samTh17VsTh0Zscores.xls",sep=""),header=T,sep="\t"))

tfs <- colnames(RI)
############ start writing output ####################
## write sif and eda files
system(paste(sep="","rm ",file.core.eda))  
system(paste(sep="","rm ",file.core.sif))
## build network for core tfs based on KC
## set treshold for significant interactions
KC.cut <- 17
## set treshold for assigning repression/activation edge
diff.z.cut.pos <- .5
diff.z.cut.neg <- -.5

trgt.gns <- core.tfs # keep track of genes that made it into network given cutoff
num.edges.in.net <- 0
num.edges.in.total <- 0
cat( "confidence_score\n", file = file.core.eda)
for(i in 1:length(core.tfs)){
  tf <- core.tfs[i]
  tf.targets.ix <- which( KC[,tf] > KC.cut )
  num.edges.in.net <- num.edges.in.net+length(tf.targets.ix)
  num.edges.in.total <- num.edges.in.total+length(which( KC[,tf] > 0 ))
  for(j in 1:length(tf.targets.ix)){
    trgt.gn <- names(tf.targets.ix)[j]
    trgt.gns <- c(trgt.gns,trgt.gn)
    if(any(rownames(tfs.diff.mat) == trgt.gn)) {
      if(tfs.diff.mat[trgt.gn,tf] > diff.z.cut.pos) {
        edge.type <- "KC_activator"
      } else if (tfs.diff.mat[trgt.gn,tf] < diff.z.cut.neg) {
        cat("in")
        edge.type <- "KC_repressor"
      } else {
        edge.type <- "KC_unsigned"
      }
    } else {
      edge.type <- "KC_unsigned"
    }
    eda.line <- paste(tf," (", edge.type,") ", trgt.gn, " = ", 
                      KC[tf.targets.ix[j],tf], sep = "")
    eda.line <- paste(eda.line,"\n",sep="")
    sif.line <- paste( tf, edge.type, trgt.gn, sep = " ")
    sif.line <- paste(sif.line,"\n",sep="")
    cat( eda.line, file = file.core.eda, append = TRUE)    
    cat( sif.line, file = file.core.sif, append = TRUE)
  }
}
cat(sep="","core net done! #nodes= ", length(unique(trgt.gns))," #edges= ", num.edges.in.net,
    " percent edges shown (top hits) out of total edges= ",round(num.edges.in.net/num.edges.in.total*100,2),"%\n")

system(paste(sep="","rm ",file.extend.eda))  
system(paste(sep="","rm ",file.extend.sif))

## build network for other tfs based on RI
## set treshold for significant interactions
RI.cut <- 17

other.tfs <- rownames(KC)[which(rownames(KC) %in% tfs)]
other.tfs <- names(which(apply(KC[other.tfs,],1,max)>KC.cut))
other.tfs <- other.tfs[which(! other.tfs %in% core.tfs)]
num.edges.in.net <- 0
num.edges.in.total <- 0
trgt.gns <- core.tfs # keep track of genes that made it into network given cutoff
cat( "confidence_score\n", file = file.extend.eda)
for(i in 1:length(other.tfs)){
  tf <- other.tfs[i]
  tf.targets.ix <- which( RI[,tf] > RI.cut )
  num.edges.in.net <- num.edges.in.net+length(tf.targets.ix)
  num.edges.in.total <- num.edges.in.total+length(which( RI[,tf] > 0 ))
  if(length(tf.targets.ix)>0){
    for(j in 1:length(tf.targets.ix)){
      trgt.gn <- names(tf.targets.ix)[j]
      trgt.gns <- c(trgt.gns,trgt.gn)
      ## check for inf edge with rnaseq data and immgen data
      if(any(rownames(beta.mat.r) == trgt.gn)) {    
        if(beta.mat.r[trgt.gn,tf]>0) {
          edge.type <- "RI_activator"
        } else if (beta.mat.r[trgt.gn,tf] < 0) {
          edge.type <- "RI_repressor"
        } else {
          edge.type <- "RI_unsigned"
        }
      } else if(any(rownames(beta.mat.i) == trgt.gn)) {
        if(beta.mat.i[trgt.gn,tf]>0) {
          edge.type <- "RI_activator"
        } else if (beta.mat.i[trgt.gn,tf] < 0) {
          edge.type <- "RI_repressor"
        } else {
          edge.type <- "RI_unsigned"
        }
      } else {
        edge.type <- "RI_unsigned"
      }
      ## eda.line <- paste(tf, edge.type, trgt.gn, "=", RI[tf.targets.ix[j],tf], sep = " ")
      eda.line <- paste(tf," (", edge.type,") ", trgt.gn, " = ", 
                        RI[tf.targets.ix[j],tf], sep = "")
      eda.line <- paste(eda.line,"\n",sep="")
      sif.line <- paste( tf, edge.type, trgt.gn, sep = " ")
      sif.line <- paste(sif.line,"\n",sep="")
      cat( eda.line, file = file.extend.eda, append = TRUE)
      cat( sif.line, file = file.extend.sif, append = TRUE)      
    }
  }
}
cat(sep="","extended net done! #nodes= ", length(unique(trgt.gns))," #edges= ", num.edges.in.net,
    " percent edges shown (top hits) out of total edges= ",round(num.edges.in.net/num.edges.in.total*100,2),"%\n")

## build network for all tfs based on RI
## set treshold for significant interactions
RI.cut <- 17
system(paste(sep="","rm ",file.all.eda))  
system(paste(sep="","rm ",file.all.sif))

other.tfs <- rownames(KC)[which(rownames(KC) %in% tfs)]
num.edges.in.net <- 0
num.edges.in.total <- 0
trgt.gns <- core.tfs # keep track of genes that made it into network given cutoff
cat( "confidence_score\n", file = file.all.eda)
for(i in 1:length(other.tfs)){
  tf <- other.tfs[i]
  tf.targets.ix <- which( RI[,tf] > RI.cut )
  num.edges.in.net <- num.edges.in.net+length(tf.targets.ix)
  num.edges.in.total <- num.edges.in.total+length(which( RI[,tf] > 0 ))
  if(length(tf.targets.ix)>0){
    for(j in 1:length(tf.targets.ix)){
      trgt.gn <- names(tf.targets.ix)[j]
      trgt.gns <- c(trgt.gns,trgt.gn)
      ## check for inf edge with rnaseq data and immgen data
      if(any(rownames(beta.mat.r) == trgt.gn)) {    
        if(beta.mat.r[trgt.gn,tf]>0) {
          edge.type <- "RI_activator"
        } else if (beta.mat.r[trgt.gn,tf] < 0) {
          edge.type <- "RI_repressor"
        } else {
          edge.type <- "RI_unsigned"
        }
      } else if(any(rownames(beta.mat.i) == trgt.gn)) {
        if(beta.mat.i[trgt.gn,tf]>0) {
          edge.type <- "RI_activator"
        } else if (beta.mat.i[trgt.gn,tf] < 0) {
          edge.type <- "RI_repressor"
        } else {
          edge.type <- "RI_unsigned"
        }
      } else {
        edge.type <- "RI_unsigned"
      }
      ## eda.line <- paste(tf, edge.type, trgt.gn, "=", RI[tf.targets.ix[j],tf], sep = " ")
      eda.line <- paste(tf," (", edge.type,") ", trgt.gn, " = ", 
                        RI[tf.targets.ix[j],tf], sep = "")
      eda.line <- paste(eda.line,"\n",sep="")
      sif.line <- paste( tf, edge.type, trgt.gn, sep = " ")
      sif.line <- paste(sif.line,"\n",sep="")
      cat( eda.line, file = file.all.eda, append = TRUE)
      cat( sif.line, file = file.all.sif, append = TRUE)      
    }
  }
}
cat(sep="","all net done! #nodes= ", length(unique(trgt.gns))," #edges= ", num.edges.in.net,
    " percent edges shown (top hits) out of total edges= ",round(num.edges.in.net/num.edges.in.total*100,2),"%\n")

## build network for tfs only based on RI
## set treshold for significant interactions
RI.cut <- 17
system(paste(sep="","rm ",file.all.tfs.only.eda))  
system(paste(sep="","rm ",file.all.tfs.only.sif))

other.tfs <- rownames(KC)[which(rownames(KC) %in% tfs)]
num.edges.in.net <- 0
num.edges.in.total <- 0
trgt.gns <- core.tfs # keep track of genes that made it into network given cutoff
cat( "confidence_score\n", file = file.all.tfs.only.eda)
for(i in 1:length(other.tfs)){
  tf <- other.tfs[i]
  tf.targets.ix <- which( RI[,tf] > RI.cut )
  tf.targets.ix <- tf.targets.ix[which( names(tf.targets.ix) %in% tfs )]
  num.edges.in.net <- num.edges.in.net+length(tf.targets.ix)
  num.edges.in.total <- num.edges.in.total+length(which( RI[tfs[which(tfs%in%rownames(RI))],tf] > 0 ))
  if(length(tf.targets.ix)>0){
    for(j in 1:length(tf.targets.ix)){
      trgt.gn <- names(tf.targets.ix)[j]
      trgt.gns <- c(trgt.gns,trgt.gn)
      ## check for inf edge with rnaseq data and immgen data
      if(any(rownames(beta.mat.r) == trgt.gn)) {    
        if(beta.mat.r[trgt.gn,tf]>0) {
          edge.type <- "RI_activator"
        } else if (beta.mat.r[trgt.gn,tf] < 0) {
          edge.type <- "RI_repressor"
        } else {
          edge.type <- "RI_unsigned"
        }
      } else if(any(rownames(beta.mat.i) == trgt.gn)) {
        if(beta.mat.i[trgt.gn,tf]>0) {
          edge.type <- "RI_activator"
        } else if (beta.mat.i[trgt.gn,tf] < 0) {
          edge.type <- "RI_repressor"
        } else {
          edge.type <- "RI_unsigned"
        }
      } else {
        edge.type <- "RI_unsigned"
      }
      ## eda.line <- paste(tf, edge.type, trgt.gn, "=", RI[tf.targets.ix[j],tf], sep = " ")
      eda.line <- paste(tf," (", edge.type,") ", trgt.gn, " = ", 
                        RI[tf.targets.ix[j],tf], sep = "")
      eda.line <- paste(eda.line,"\n",sep="")
      sif.line <- paste( tf, edge.type, trgt.gn, sep = " ")
      sif.line <- paste(sif.line,"\n",sep="")
      cat( eda.line, file = file.all.tfs.only.eda, append = TRUE)
      cat( sif.line, file = file.all.tfs.only.sif, append = TRUE)      
    }
  }
}
cat(sep="","all net tfs only done! #nodes= ", length(unique(trgt.gns))," #edges= ", num.edges.in.net,
    " percent edges shown (top hits) out of total edges= ",round(num.edges.in.net/num.edges.in.total*100,2),"%\n")

## write sif and eda files
system(paste(sep="","rm ",file.core.tfs.only.eda))  
system(paste(sep="","rm ",file.core.tfs.only.sif))
## build network for core tfs based on KC
## set treshold for significant interactions
KC.cut <- 17
## set treshold for assigning repression/activation edge
diff.z.cut.pos <- .5
diff.z.cut.neg <- -.5

trgt.gns <- core.tfs # keep track of genes that made it into network given cutoff
num.edges.in.net <- 0
num.edges.in.total <- 0

cat( "confidence_score\n", file = file.core.tfs.only.eda)
for(i in 1:length(core.tfs)){
  tf <- core.tfs[i]
  tf.targets.ix <- which( KC[,tf] > KC.cut )
  tf.targets.ix <- tf.targets.ix[which( names(tf.targets.ix) %in% tfs )]
  num.edges.in.net <- num.edges.in.net+length(tf.targets.ix)
  num.edges.in.total <- num.edges.in.total+length(which( KC[tfs[which(tfs%in%rownames(KC))],tf] > 0 ))
  for(j in 1:length(tf.targets.ix)){
    trgt.gn <- names(tf.targets.ix)[j]
    trgt.gns <- c(trgt.gns,trgt.gn)
    if(any(rownames(tfs.diff.mat) == trgt.gn)) {
      if(tfs.diff.mat[trgt.gn,tf] > diff.z.cut.pos) {
        edge.type <- "KC_activator"
      } else if (tfs.diff.mat[trgt.gn,tf] < diff.z.cut.neg) {
        cat("in")
        edge.type <- "KC_repressor"
      } else {
        edge.type <- "KC_unsigned"
      }
    } else {
      edge.type <- "KC_unsigned"
    }
    eda.line <- paste(tf," (", edge.type,") ", trgt.gn, " = ", 
                      KC[tf.targets.ix[j],tf], sep = "")
    eda.line <- paste(eda.line,"\n",sep="")
    sif.line <- paste( tf, edge.type, trgt.gn, sep = " ")
    sif.line <- paste(sif.line,"\n",sep="")
    cat( eda.line, file = file.core.tfs.only.eda, append = TRUE)    
    cat( sif.line, file = file.core.tfs.only.sif, append = TRUE)
  }
}
cat(sep="","core net tfs only done! #nodes= ", length(unique(trgt.gns))," #edges= ", num.edges.in.net,
    " percent edges shown (top hits) out of total edges= ",round(num.edges.in.net/num.edges.in.total*100,2),"%\n")

## get out degree for each tf (also weighted by RI score)
out.weighted.degree.sum <- apply(RI[,tfs],2,sum)
out.degree.sum <- apply(RI[,tfs],2, function(i) length(which(i>0)) )
degree.mat <- matrix(0,nr=length(gns),nc=2, dimnames=list(gns,c("out_weighted_degree_sum","out_degree_sum")))
good.tfs <- colnames(RI)[which(apply(RI,2,max)>10)]
degree.mat[good.tfs,"out_degree_sum"] <- out.degree.sum[good.tfs]
degree.mat[good.tfs,"out_weighted_degree_sum"] <- out.weighted.degree.sum[good.tfs]

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
  if (gns[i] %in% tfs){
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

## get node attributes for which gn names to show on network (for presentations)
show.gn.nm <- numeric(length=length(gns))
for(i in 1:length(gns)){
  if(gns[i] %in% core.tfs){
    show.gn.nm[i] <- 2
  } else if( gns[i] %in% show.gn.names){
    show.gn.nm[i] <- 1
  } else {
    show.gn.nm[i] <- 0
  }
}

## get KC.sum of each gene
kc.sum <- apply(KC[,core.tfs],1,sum)
## write kc.sum, show.gn.nm, and gn.function into other.mrna
fl.nm <- paste(path.output,"other_zcut_",z.abs.cut,"_",date.is,".txt",sep="")
x <- cbind(kc.sum,show.gn.nm,gn.function,degree.mat)
colnames(x) <- c("kc.sum","show.gn.nm","gn.function","out_weighted_degree_sum","out_degree_sum")
cat("gene_id",colnames(x),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(sep="\t",x,quote=FALSE,append=TRUE,col.names=FALSE,file=fl.nm)


## write differential expression levels as .mrna file
fl.nm <- paste(path.output,"file_diff_exprssn_zcut_",z.abs.cut,"_",date.is,".txt" ,sep="")
m.diff <- as.matrix(read.table(paste(path.input,"rnaseqDiffExpMedianZscoresPerTf.xls",sep=""),header=T,sep="\t"))
ix.gns.not.in.m.diff <- which(!gns %in% rownames(m.diff)) # find gene names not in diff exp mat
tmp <- matrix(0,nr=length(ix.gns.not.in.m.diff),nc=dim(m.diff)[2])
rownames(tmp) <- gns[ix.gns.not.in.m.diff]
m.diff.inclusive <- rbind(m.diff,tmp) # add all gene names with 0 as diff exp values
cat("gene_id",colnames(m.diff),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(sep="\t",m.diff.inclusive[gns,],quote=FALSE,append=TRUE,col.names=FALSE,file=fl.nm)

## calc specificity scores and save to variable specificity.mat
source("r_scripts/th17/used_for_paper/calcSpecificityScores.R")
fl.nm <- paste(path.output,"specificity_mat_zcut_",z.abs.cut,"_",date.is,".txt",sep="")
cat("gene_id",colnames(specificity.mat),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(specificity.mat,sep="\t",quote=FALSE,append=TRUE,col.names=FALSE,file=fl.nm)


## write all rnaseq data log2(rpkm)
fl.nm <- paste(path.output,"rnaseq_dataset_log2rpkm_zcut_",z.abs.cut,"_",date.is,".txt",sep="")
colnames(d.log) <- paste(sep="_",colnames(d.log),"log2_rpkm")
cat("gene_id",colnames(d.log),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(d.log[,],sep="\t",quote=FALSE,append=TRUE,col.names=FALSE,file=fl.nm)

## write all rnaseq data (rpkm)
fl.nm <- paste(path.output,"rnaseq_dataset_rpkm_zcut_",z.abs.cut,"_",date.is,".txt",sep="")
colnames(d) <- paste(sep="_",colnames(d),"rpkm")
cat("gene_id",colnames(d),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(d[,],sep="\t",quote=FALSE,append=TRUE,col.names=FALSE,file=fl.nm)

# write all enrichment results
e.mat <- as.matrix(read.table(paste("results/enrichmentAnalysis/enrichment_scores_zcut_",z.abs.cut,"_",date.is,".xls",sep=""),header=T,sep="\t"))
cat("\n",append=TRUE,file=fl.nm)
ix.na <- which(is.na(e.mat))
e.mat[ix.na] <- 0
c.names <- colnames(e.mat)[c(3:8,12,14)]
fl.nm <- paste(path.output,"enrichment_scores_zcut_",z.abs.cut,"_",date.is,".txt",sep="")
cat("gene_id",c.names,sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(e.mat[,c.names],sep="\t",quote=FALSE,append=TRUE,col.names=FALSE,file=fl.nm)

## ## create histogram of SAM diff expression scores for SKI and SMAD3 regulated genes
## tf="SKI"
## t.ski <- names(which(RI[,tf]>RI.cut))
## tf="SMAD3"
## t.smad3 <- names(which(RI[,tf]>RI.cut))
## t.ski.and.smad3 <- intersect(t.ski,t.smad3)
## tf="RORC"
## t.rorc <- names(which(RI[,tf]>RI.cut))
## tf="IRF4"
## t.irf4 <- names(which(RI[,tf]>RI.cut))
## tf="BATF"
## t.batf <- names(which(RI[,tf]>RI.cut))

## s.sam.rorc <- sam[t.rorc,1]
## s.sam.irf4 <- sam[t.irf4,1]
## s.sam.batf <- sam[t.batf,1]

## s.sam.ski <- sam[t.ski,1]
## s.sam.smad3 <- sam[t.smad3,1]
## s.sam.both <- sam[t.ski.and.smad3,1]

## s.kc.ski <- kc.sum[t.ski]
## s.kc.smad3 <- kc.sum[t.smad3]
## s.kc.both <- kc.sum[t.ski.and.smad3]








