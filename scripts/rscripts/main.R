##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# all that is needed for running network inference for th17 project
#----------------#
# 0-
#----------------#
# set global variables
cat("\n- setting global variables\n")
rm(list=ls())
GLOBAL <- list()
#GLOBAL[["run.these.steps"]] <- c("12.1.1","14.1.1","14.2.1") # for a complete run c(1:last_step)
GLOBAL[["run.these.steps"]] <- c(7)
# stamp the date on this run
x <- unlist(strsplit(date()," +",perl=TRUE))
GLOBAL[["date.is"]] <- paste(x[2],x[3],x[5],sep="_")
# data set size: take genes with abs(zscore) based on diff expression btwn th17_48hr and th0_48hr
GLOBAL[["z.abs.cut"]] <- 2.5 #zscore
# Only used to determine rnaseq data set size
# GLOBAL[["median.abs.cut"]] <- 3 #rpkm
GLOBAL[["min.th17.or.th0.rpkm"]] <- 3 #rpkm
# The number of bootstraps for inferelator rnaseq and immgen
GLOBAL[["num.boots.rnaseq"]] <- 200
GLOBAL[["num.boots.immgen"]] <- 200
GLOBAL[["num.perms.enrichment.analysis"]] <- 10
# do we want to output the sequences around the peaks?
GLOBAL[["get.sequence"]] <- TRUE
# what is the minimum distance between peaks to cluster them together
## GLOBAL[["tfs.min.dist"]] <- c(50,100,150,200,250,300,500,1000)
GLOBAL[["tfs.min.dist"]] <- c(100)
# multi TF binding sites distance (up/downstrem) from TSS for gene MTBS mapping
GLOBAL[["mtbs.tss.dist"]] <- 5000
GLOBAL[["mm9.effective.genome.size"]] <- 1.87e9

# genes that we want in as a rule even if they have low median expression or z.score that kicks them out
GLOBAL[["known.tfs"]] <- c("BATF","MAF","IRF4","STAT3","RORC")
GLOBAL[["known.gns"]] <- toupper(c("Il17a","Il17f","Il23r","Il22","Il21","Foxp3","Rora"))
GLOBAL[["use.multicore"]] <- TRUE
if(GLOBAL[["use.multicore"]]==TRUE){
  library(multicore)
}


#----------------#
# 1-
#----------------#
#Put all the rnaseq data that we have in different files into one data matrix.
#    -- script:
#       - r_scripts/th17/used_for_paper/createRnaseqCompleteDataset.R
#    -- input: 
#       - input/th17/rawData/no_quartnorm/SL* (one of many rnaseq experiment, e.g. SL1040)
#    -- output:
#       - input/th17/V2_used_for_paper/ranseqDatasetNoQuartNorm.RData (this matrix is not log2 transformed)
if(any(GLOBAL[["run.these.steps"]]==1)) {
  cat("\n- creating rnaseq dataset from raw data\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  date.is <- GLOBAL[["date.is"]]
  date.cufflinks.data.compiled <- "Jul_30_2012"
  source("r_scripts/th17/used_for_paper/createRnaseqCompleteDataset.R")
}

#----------------#
# 2-
#----------------#
#Calculate differential expression.  For each gene compare its expression levels between th17_48hr and th0_48hr using SAM to find differentially expressed genes that are specifically important for th17 phenotype (we have a wopping number of 9 replicates for each condition!!). rpkm are log2 transformed before differential expression by SAM
#    -- script:
#       - r_scripts/th17/used_for_paper/calcDiffExpSAM.R
#    -- input: 
#       - input/th17/used_for_paper/ranseqDatasetNoQuartNorm.RData
#    -- output:
#       - input/th17/used_for_paper/samTh17VsTh0Zscores.xls
#       - input/th17/used_for_paper/gnsMedian.xls
if(any(GLOBAL[["run.these.steps"]]==2)) {
  cat("\n- calculating SAM differential expression zscores between th17_48hr and th0_48hr\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  source("r_scripts/th17/used_for_paper/calcDiffExpSAM.R")
}

#----------------#
# 3-
#----------------#
#Calculate for differential expression for each tf between th17_48hr_tf_wt to th17_48hr_tf_ko to find potential targets for each gene.
#    -- script:
#       - r_scripts/th17/used_for_paper/DEseq_for_core_tfs.R
#    -- input: 
#       - input/th17/used_for_paper/
#    -- output:
#       - input/th17/used_for_paper/
if(any(GLOBAL[["run.these.steps"]]==3)) {
  cat("\n- using DEseq to calculate differential expression for tf-ko vs. tf-wt\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  date.htseq.data.compiled <- "May_27_2012"
  date.is <- GLOBAL[["date.is"]]
  fpkm.cut <- 3 # pvals for genes that have low mean rpkm (lower than rpkm.cut) in both WT and KO recieve pval of 1
  # fc.cut.ko <- 30 # pvals for genes that have low FC (lower than fc.cut % up or down) in between WT and KO recieve pval of 1
  # fc.cut.kd <- 20 # pvals for genes that have low FC (lower than fc.cut % up or down) in between WT and KO recieve pval of 1
  pseudo.count <- 1   ## we add 1 pseudo count to all pixles to resolve division by zero
  feature.rm <- c("NO_FEATURE","AMBIGUOUS","TOO_LOW_AQUAL","NOT_ALIGNED","ALIGNMENT_NOT_UNIQUE")
  source("r_scripts/th17/used_for_paper/DEseq_for_core_tfs_paired.R")
}

#----------------#
# 4-
#----------------#
## this code is now part of Ashish's pipeline and does not need to be run
# map chipseq MACS peaks to genes. This takes time (2-3 mins per chip)
#    -- script:
#       - r_scripts/th17/used_for_paper/createChipScores.R
#    -- input: 
#       - input/th17/input/th17/used_for_paper/rawData/MACS/*
#    -- output:
#       - input/used_for_paper/
if(any(GLOBAL[["run.these.steps"]]==4)) {
  cat("\n- mapping chip peaks to genes\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  date.is <- GLOBAL[["date.is"]]
  date.macs.data.compiled <- "Sep_15_2011"
  peak.dist.upstream <- 10000
  peak.dist.downstream <- 10000
  source("r_scripts/th17/used_for_paper/mapChipPeaksToGenes.R")
}

#----------------#
# 5-
#----------------#
#Calculate for each tf how good it binds to each gene.  Calculation is based on the ranks of three parameters of chip-seq:
#  number and significance of pvaluethat fall inside a gene region, and the number of randomly expected number of peaks of equal or higher mean_pvalue (of all peaks inside gene region)
#  in the gene region.
#    -- script:
#       - r_scripts/th17/used_for_paper/createChipScores.R
#    -- input: 
#       - input/th17/chipSeqV10/geneCentric/* (one of many chipseq experiments mapped to genes, e.g. BATF-Th0-117-anno-20101123-1444.txt)
#    -- output:
#       - input/V2_used_for_paper/chipScoresPerGenePerTfList.RData
if(any(GLOBAL[["run.these.steps"]]==5)) {
  cat("\n- calculating chipseq scores\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  peak.dist.upstream <- 10000
  peak.dist.downstream <- 10000
  core.tfs <- GLOBAL[["known.tfs"]]
  known.genes <- GLOBAL[["known.gns"]]
  use.multicore <- GLOBAL[["use.multicore"]]
  n.pr <- 7 # number of processors to use if use.multicore is TRUE
  mm9.effective.genome.size <- GLOBAL[["mm9.effective.genome.size"]]
  date.macs.data.compiled <- "Sep_15_2011"
  source("r_scripts/th17/used_for_paper/calcChipScores.R")
}

#----------------#
# 6-
#----------------#
# This code calculates clusters of peaks and their association with TFs.
#    -- script:
#       - r_scripts/th17/used_for_paper/cluster_peaks.R
#    -- input: 
#       - input/th17/used_for_paper/peaks_mat_List.RData
#    -- output:
#       - input/th17/used_for_paper/tf_clusters_per_loci_list.RData
if(any(GLOBAL[["run.these.steps"]]==6)) {
  cat("\n- calculating tfs chipseq clusters\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  ## core.tfs <- GLOBAL[["known.tfs"]]
  date.is <- GLOBAL[["date.is"]]
  date.macs.data.compiled <- "Sep_15_2011"
  mm9.genome.size <- GLOBAL[["mm9.effective.genome.size"]]
  # make mtbs.tss.dist Inf to use original association of chip peaks with trgt genes (ie. all gene region +/-10kb)
  mtbs.tss.dist <- GLOBAL[["mtbs.tss.dist"]]
  n.p <- 4 # number of processors to use
  n.b <- 4 # number of bootstraps
  load.macs.per.chrom.RData <- FALSE
  add.histone.scores <- FALSE
  match.targets <- TRUE # do you want to include the targets of each cluster?
  calc.simulations <- FALSE # do you want to run 10000 simulations to determine significance?
  for(iter in 1:length(GLOBAL[["tfs.min.dist"]])){
    d.cut <- GLOBAL[["tfs.min.dist"]][iter]
    cat("working on d.cut ",d.cut,"\n")
    source("r_scripts/th17/used_for_paper/cluster_peaks.R")
  }
}

if(any(GLOBAL[["run.these.steps"]]==6.1)) {
  cat("\n- calculating tfs chipseq clusters\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  ## core.tfs <- GLOBAL[["known.tfs"]]
  date.is <- GLOBAL[["date.is"]]
  date.macs.data.compiled <- "Sep_15_2011"
  mm9.genome.size <- GLOBAL[["mm9.effective.genome.size"]]
  mtbs.tss.dist <- GLOBAL[["mtbs.tss.dist"]] 
  n.p <- 4 # number of processors to use
  n.b <- 4 # number of bootstraps
  load.macs.per.chrom.RData <- FALSE
  add.histone.scores <- F
  match.targets <- TRUE # do you want to include the targets of each cluster?
  calc.simulations <- FALSE # do you want to run 10000 simulations to determine significance?
  for(iter in 1:length(GLOBAL[["tfs.min.dist"]])){
    d.cut <- GLOBAL[["tfs.min.dist"]][iter]
    cat("working on d.cut ",d.cut,"\n")
    source("r_scripts/th17/used_for_paper/cluster_peaks_with_histone_marks.R")
  }
}

#----------------#
# 7-
#----------------#
#Calculate for each tf how good it binds to each gene.  Calculation is based on the ranks of three parameters of chip-seq: pvalue, foldchange, and distance from transcrption start site. 
#    -- script:
#       - r_scripts/th17/used_for_paper/createChipScores.R
#    -- input: 
#       - input/th17/chipSeqV10/geneCentric/* (one of many chipseq experiments mapped to genes, e.g. BATF-Th0-117-anno-20101123-1444.txt)
#    -- output:
#       - input/V2_used_for_paper/chipScoresPerGenePerTfList.RData
if(any(GLOBAL[["run.these.steps"]]==7)) {
  cat("\n- getting sequences for motif discovery\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  date.is <- GLOBAL[["date.is"]]
  date.macs.data.compiled <- "Sep_15_2011"
  date.cluster.peaks.run <- "Sep_18_2011"
  # next value is hard wired atm, do not change!
  n.p <- 3 # number of processors to use to get sequences (mem intensive don't use more than 3 for 8gigs of mem)
  peak.region.seq.length <- 50
  num.genes.with.highest.pval <- 600
  d.cut <- 100 # only do it for d.cut = 100
  source("r_scripts/th17/used_for_paper/getSequencesForMotifDiscovery.R")
}

#----------------#
# 8-
#----------------#
#Create inferelator data structures.
#    -- script:
#       - r_scripts/th17/used_for_paper/createInfDataStructures.R
#    -- input: 
#       - input/th17/used_for_paper/ranseqDatasetNoQuartNorm.RData (also: humanTFs.RData,samTh17VsTh0Zscores.xls,rnaseqDiffExpMedianZscoresPerTf.xls)
#    -- output:
#       - input/th17/used_for_paper/infDataStructures/rnaseq/1000Genes/ratios.RData (also: clusterStack.RData,colMap.RData,tfNames.RData,knockOutConfList.RData)
if(any(GLOBAL[["run.these.steps"]]==8)) {
  cat("\n- creating inferelator data structures for rnaseq run\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  date.is <- GLOBAL[["date.is"]]
  # for rnaseq or immgen?
  GLOBAL[["type"]] <- "rnaseq"
  date.cufflinks.data.compiled <- "Jul_30_2012"
  combine.file.map <- paste(sep="","combine_kcri_core_",28,"_tfs.xls") # or combine_kcri_core_5tfs.xls
  # additional.regulators <- c("KDM6B","JMJD6","LEF1","INHBA","TRIB3","SIRT2","AES","EGR2")
  source("r_scripts/th17/used_for_paper/createInfDataStructures.R")
  GLOBAL[["dataset.path"]] <- path.output
}
#----------------#
# 9-
#----------------#
#Run mix-CLR->Inferelator.  Allow for calculation of an extended network (i.e. not just for tfs for which we have ko chipseq and rnaseq experiments)
#    -- script:
#       - r_scripts/th17/used_for_paper/runNetInference.R
#    -- input: 
#       - input/th17/used_for_paper/infDataStructures/rnaseq/1000Genes/* (or immgen/* for either 100gene or 1000gene)
#    -- output:
#       - input/th17/used_for_paper/infDataStructures/rnaseq/1000Genes/ratios.RData (also: clusterStack.RData,colMap.RData,tfNames.RData,knockOutConfList.RData)
if(any(GLOBAL[["run.these.steps"]]==9)) {
  cat("\n- running inferelator on rnaseq data\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  PARAMS = list()
  PARAMS[["general"]] = list()
  PARAMS[["clr"]] = list()
  PARAMS[["lars"]] = list()
  PARAMS[["output"]] = list()
  ######################### most useful params ##################################
  # what dataset do you want to run? 
  # options:
  # - th17_1000
  # - immgen_1000
  # (this must agree with n and type from step 8)!!!!!!!!!!
  PARAMS[["general"]][["data"]] = GLOBAL[["dataset.path"]]
  # how many predictors (tf's) for elastic net to choose from?
  PARAMS[["lars"]][["max_single_preds"]] = 30
  # what l2 norm weights to use in elastic net? (lambda =0 same as LARS)
  PARAMS[["lars"]][["lambda"]] = c(0) # (0,1,100)
  # how many bootstrap runs? (one for each bootstraped dataset)
  PARAMS[["general"]][["numBoots"]] = GLOBAL[["num.boots.rnaseq"]]
  # what is the maximum delta T to consider (in minutes)? (if delta T is bigger than it will be treated as steady state)
  PARAMS[["general"]][["delT_max"]] = 7*60
  # what is the time order (tau) of the reactinos (i.e. in what time do you expect reactions to happen in)?
  PARAMS[["general"]][["tau"]] = 3*60
  # for tfs for which we have ko data.  What zscore to consider for regulation (bigger than this number in abs terms)?
  PARAMS[["general"]][["rnaseq_ko_zscore_abs_filter"]] = 1
  # how many processors to use? if use more than one need to install the multicore package
  PARAMS[["general"]][["processorsNumber"]] = 6
  # What n fold cross validation for inferelator should we use?
  PARAMS[["lars"]][["nCv"]] = 10
  # How many bins to use to calculate mutual information?
  PARAMS[["clr"]][["n.bins"]] = 10
  # do you want to use priors from rnaseq ko data for each tf?
  PARAMS[["general"]][["use.priors"]] =FALSE
################################################################################
  source("r_scripts/th17/used_for_paper/runNetInference.R")
}
#----------------#
# 10-
#----------------#
#Create inferelator data structures.
#    -- script:
#       - r_scripts/th17/used_for_paper/createInfDataStructures.R
#    -- input: 
#       - input/th17/used_for_paper/ranseqDatasetNoQuartNorm.RData (also: humanTFs.RData,samTh17VsTh0Zscores.xls,rnaseqDiffExpMedianZscoresPerTf.xls)
#    -- output:
#       - input/th17/used_for_paper/infDataStructures/rnaseq/1000Genes/ratios.RData (also: clusterStack.RData,colMap.RData,tfNames.RData,knockOutConfList.RData)
if(any(GLOBAL[["run.these.steps"]]==10)) {
  cat("\n- creating inferelator data structures for immgen run\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  date.is <- GLOBAL[["date.is"]]
  # for rnaseq or immgen?
  GLOBAL[["type"]] <- "immgen"
  date.cufflinks.data.compiled <- "Jul_30_2012"
  combine.file.map <- paste(sep="","combine_kcri_core_",28,"_tfs.xls") # or combine_kcri_core_5tfs.xls
  # additional.regulators <- c("KDM6B","JMJD6","LEF1")
  source("r_scripts/th17/used_for_paper/createInfDataStructures.R")
  GLOBAL[["dataset.path"]] <- path.output
}
#----------------#
# 11-
#----------------#
#Run mix-CLR->Inferelator.  Allow for calculation of an extended network (i.e. not just for tfs for which we have ko chipseq and rnaseq experiments)
#    -- script:
#       - r_scripts/th17/used_for_paper/runNetInference.R
#    -- input: 
#       - input/th17/used_for_paper/infDataStructures/rnaseq/1000Genes/* (or immgen/* for either 100gene or 1000gene)
#    -- output:
#       - input/th17/used_for_paper/infDataStructures/rnaseq/1000Genes/ratios.RData (also: clusterStack.RData,colMap.RData,tfNames.RData,knockOutConfList.RData)
if(any(GLOBAL[["run.these.steps"]]==11)) {
  cat("\n- running inferelator on immgen data\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  PARAMS = list()
  PARAMS[["general"]] = list()
  PARAMS[["clr"]] = list()
  PARAMS[["lars"]] = list()
  PARAMS[["output"]] = list()
  ######################### most useful params ##################################
  # what dataset do you want to run? 
  # options:
  # - th17_1000
  # - immgen_1000
  # (this must agree with n and type from step 5)!!!!!!!!!!
  PARAMS[["general"]][["data"]] = GLOBAL[["dataset.path"]]
  # how many predictors (tf's) for elastic net to choose from?
  PARAMS[["lars"]][["max_single_preds"]] = 30
  # what l2 norm weights to use in elastic net? (lambda =0 same as LARS)
  PARAMS[["lars"]][["lambda"]] = c(0) # (0,1,100)
  # how many bootstrap runs? (one for each bootstraped dataset)
  PARAMS[["general"]][["numBoots"]] = GLOBAL[["num.boots.immgen"]]
  # what is the maximum delta T to consider (in minutes)? (if delta T is bigger than it will be treated as steady state)
  PARAMS[["general"]][["delT_max"]] = 7*60
  # what is the time order (tau) of the reactinos (i.e. in what time do you expect reactions to happen in)?
  PARAMS[["general"]][["tau"]] = 3*60
  # for tfs for which we have ko data.  What zscore to consider for regulation (bigger than this number in abs terms)?
  PARAMS[["general"]][["rnaseq_ko_zscore_abs_filter"]] = 1
  # how many processors to use? if use more than one need to install the multicore package
  PARAMS[["general"]][["processorsNumber"]] = 6
  # What n fold cross validation for inferelator should we use?
  PARAMS[["lars"]][["nCv"]] = 10
  # How many bins to use to calculate mutual information?
  PARAMS[["clr"]][["n.bins"]] = 10
  # do you want to use priors from rnaseq ko data for each tf?
  PARAMS[["general"]][["use.priors"]] =FALSE
  ################################################################################
  source("r_scripts/th17/used_for_paper/runNetInference.R")
}

#----------------#
# 12-
#----------------#
# Combine results from: (1) KO differential expression (for core tfs),(2) chipseq analysis (for core tfs),(3) rnaseq mixCLR->Inf analysis,(4) immgen mixCLR->Inf analysis.
#    -- script:
#       - r_scripts/th17/used_for_paper/combineData.R
#    -- input:
#       - input/th17/V2_used_for_paper/infResults/rnaseq/z_abs_cut_2.5_mix_clr_inf_median_scores_Nb_200.xls (for rnaseq results)
#       - input/th17/V2_used_for_paper/infResults/immgen/z_abs_cut_2.5_mix_clr_inf_median_scores_Nb_200.xls (for immgen results)
#       - input/th17/V2_used_for_paper/rnaseqDiffExpMedianZscoresPerTf.xls (for KO diff expression results)
#       - input/th17/V2_used_for_paper/chipScores.xls (for chipseq results)
#    -- output:
#       - results/th17/used_for_paper/combined_analysis_date/combinedScores.xls
if(any(GLOBAL[["run.these.steps"]]==12)) {
  cat("\n- combine scores from: ko,chip,rnaseq,and immgen\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
  date.is <- GLOBAL[["date.is"]]
  num.boots.immgen <- GLOBAL[["num.boots.immgen"]]
  num.boots.rnaseq <- GLOBAL[["num.boots.rnaseq"]]
  core.tfs <- GLOBAL[["known.tfs"]]
  core.gns <- GLOBAL[["known.gns"]]
  date.inf.run <- "Sep_16_2011"
  date.chip.run <- "Sep_15_2011"
  date.deseq.run <- "Sep_21_2011"
  types <- c("abs","k_activator","k_r_i_activator","k_r_i_repressor")
  for(t in 1:length(types)){
    type <- types[t]
    cat("working on type: ",type,"\n")
    source("r_scripts/th17/used_for_paper/combineData.R")
  }
}

if(any(GLOBAL[["run.these.steps"]]==12.1)) {
  cat("\n- combine scores from: ko,chip,rnaseq,and immgen\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  num.tfs <- 8 # 8 or 5 (5=core.tfs,8=core.tfs+etv6,fosl2,hif1a)
  combine.file.map <- paste(sep="","combine_kcri_core_",num.tfs,"_tfs.xls") # or combine_kcri_core_5tfs.xls
  date.is <- GLOBAL[["date.is"]]
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
  num.boots.immgen <- GLOBAL[["num.boots.immgen"]]
  num.boots.rnaseq <- GLOBAL[["num.boots.rnaseq"]]
  date.inf.run <- "Sep_16_2011"
  #date.deseq.run <- "May_29_2012"
  date.deseq.run <- "Jul_12_2012"
  chip.integration <- "genewide_pois_model_pval"
  # chip.integration <- "prox_pois_model_pval"
  prcnt.chng.cut <- 0 # filter out all scores K,C,R,I that have a KO support that doesn't lead to more than x% change of wt expression are set to zero!
  filter.by.sam <- F
  source("r_scripts/th17/used_for_paper/combine_kc_generic.R")
}

if(any(GLOBAL[["run.these.steps"]]==12.2)) {
	cat("\n- combine scores from: ko,chip,rnaseq,and immgen\n")
	for(iter in 1:4){
		rm(list=ls()[-which(ls() %in% c("GLOBAL","iter"))])
		date.is <- GLOBAL[["date.is"]]
		num.boots.immgen <- GLOBAL[["num.boots.immgen"]]
		num.boots.rnaseq <- GLOBAL[["num.boots.rnaseq"]]
		# date.inf.run <- "Sep_16_2011"
		date.inf.run <- "Aug_6_2012"		
		# date.deseq.run <- "Jun_18_2012"
		date.deseq.run <- "Aug_2_2012"
		chip.integration <- "genewide_pois_model_pval"
		# chip.integration <- "prox_pois_model_pval"
		z.abs.cut.old <- GLOBAL[["z.abs.cut"]]
		z.abs.cut <- 0
		if(iter==1){prcnt.chng.cut <- 0;filter.by.sam <- F;num.tfs <- 5;deseq.pval.cut <- 1;}
		else if(iter==2){prcnt.chng.cut <- 20;filter.by.sam <- F;num.tfs <- 5;deseq.pval.cut <- 0.25;}
		else if(iter==3){prcnt.chng.cut <- 20;filter.by.sam <- F;num.tfs <- 28;deseq.pval.cut <- 0.25;}
		else if(iter==4){prcnt.chng.cut <- 0;filter.by.sam <- F;num.tfs <- 28;deseq.pval.cut <- 1;}
		combine.file.map <- paste(sep="","combine_kcri_core_",num.tfs,"_tfs.xls") # or combine_kcri_core_5tfs.xls
		source("r_scripts/th17/used_for_paper/combine_kc_generic_no_p300_no_th0.R")
	}
}
#----------------#
# 13-
#----------------#
# Perform enrichemnt analysis.  Identify new tfs, (i.e. not core tfs we already no about), that may be involved in th17 differentiation.
# Only two of our four data types,rnaseq and immgen, contain information with respect to not core tfs.
# We use the rnaseq and immgen results to rank the targets of each putative tf and see if the list is enriched (similar) for the list of ranked targets that results from the core TF.
#    -- script:
#       - r_scripts/th17/used_for_paper/newTFsEnrichmentAnalysis.R
#    -- input:
#       - input/th17/V2_used_for_paper/infResults/rnaseq/z_abs_cut_2.5_mix_clr_inf_median_scores_Nb_200.xls (for rnaseq results)
#       - input/th17/V2_used_for_paper/infResults/immgen/z_abs_cut_2.5_mix_clr_inf_median_scores_Nb_200.xls (for immgen results)
#       - input/th17/V2_used_for_paper/rnaseqDiffExpMedianZscoresPerTf.xls (for KO diff expression results)
#       - input/th17/V2_used_for_paper/chipScores.xls (for chipseq results)
#    -- output:
#       - results/th17/used_for_paper/combined_analysis_date/combinedScores.xls
if(any(GLOBAL[["run.these.steps"]]==13)) {
  cat("\n- performing new TFs enrichment analysis\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
  date.is <- GLOBAL[["date.is"]]
  core.tfs <- GLOBAL[["known.tfs"]][-which(GLOBAL[["known.tfs"]]=="FOSL2")]
  N <- GLOBAL[["num.perms.enrichment.analysis"]]
  source("r_scripts/th17/used_for_paper/TFsEnrichmentAnalysis.R")
}

# development version: still need work
if(any(GLOBAL[["run.these.steps"]]==13.1)) {
  cat("\n- performing new TFs enrichment analysis\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
  prcnt.chng.cut <- 20
  filter.by.sam <- F
  # tfs = c("BATF","MAF","IRF4","STAT3","RORC")
  # OR
  # tfs = c("BATF","MAF","IRF4","STAT3","RORC","FOSL2","HIF1A","ETV6")
  tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
  num.tfs <- length(tfs) # 8 or 5 (5=core.tfs,8=core.tfs+etv6,fosl2,hif1a)
  date.is <- GLOBAL[["date.is"]]
  # core.tfs <- GLOBAL[["known.tfs"]][-which(GLOBAL[["known.tfs"]]=="FOSL2")]
  # core.tfs <- GLOBAL[["known.tfs"]]
  # Np <- GLOBAL[["num.perms.enrichment.analysis"]]
  #date.combine <- "Jun_21_2012"
  # date.combine <- "Jul_11_2012"
  date.combine <- "Aug_8_2012"
  type.1 <- "th17" # "th17_minus_th0"
  # gs <- "SAMdownreg" # what to create gold standard list from?
  dataset <- "I" # compare gold standardto what data?
  # number of top genes to consider from gold standard (to see how enriched they are in predicted list)
  N <- 100 
  Np <- 10000 # num permutation for pval calculations
  deseq.pval.cut <- 0.25
  # all.gs <- c("SAMdownreg","SAMupreg","sum.KC.up","sum.KC.down")
  all.gs <- c("sum.KC.up")
  # all.tps <- c("whole","activation" ,"repression", "absolute")
  all.tps <- c("activation")
  for(gs in all.gs){
	  for(tp in all.tps){
	  	type.2 <- tp # "repression" ,"activation", "absolute", "whole"
	  	source("r_scripts/th17/used_for_paper/TFsEnrichmentAnalysis_dev.R")
	  }
  }
}

#----------------#
# 14-
#----------------#
# Create datastructures for cytoscape and for sungear.
#    -- script:
#       - r_scripts/th17/used_for_paper/newTFsEnrichmentAnalysis.R
#    -- input:
#       - input/th17/used_for_paper/infResults/rnaseq/z_abs_cut_2.5_mix_clr_inf_median_scores_Nb_200.xls (for rnaseq results)
#       - input/th17/used_for_paper/infResults/immgen/z_abs_cut_2.5_mix_clr_inf_median_scores_Nb_200.xls (for immgen results)
#       - input/th17/used_for_paper/rnaseqDiffExpMedianZscoresPerTf.xls (for KO diff expression results)
#       - input/th17/used_for_paper/chipScores.xls (for chipseq results)
#    -- output:
#       - results/th17/used_for_paper/combined_analysis_date/combinedScores.xls
if(any(GLOBAL[["run.these.steps"]]==14)) {
  cat("\n- creating cytoscape and sungear data structinos\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  date.is <- GLOBAL[["date.is"]]
  date.combine.data.run <- "Sep_21_2011"
  date.inf.run <- "Sep_16_2011"
  date.deseq.run <- "Sep_21_2011"
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
  num.boot.i = GLOBAL[["num.boots.immgen"]]
  num.boot.r = GLOBAL[["num.boots.rnaseq"]]
  show.gn.names <- c("BATF","MAF","IRF4","STAT3","RORC")
  core.tfs <- GLOBAL[["known.tfs"]]
  cut.kcri.core.activ.qunt <- .8 # get network for 20% top edges for core
  cut.kcri.extend.activ.qunt <- .97
  cut.kcr.core.activ.qunt <- .8 # get network for 20% top edges for core
  cut.kcr.extend.activ.qunt <- .97
  cut.kc.core.activ.qunt <- .8 # get network for 20% top edges for core
  cut.kc.extend.activ.qunt <- .97
  # for repression use the exact cut values (in absolute term, not quontile, to make the networks comparable)
  source("r_scripts/th17/used_for_paper/createCytoscapeNetworkData_KCRI.R")
}

if(any(GLOBAL[["run.these.steps"]]==14.1)) {
  cat("\n- creating cytoscape data structinos\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  date.is <- GLOBAL[["date.is"]]
  # date.combine.data.run <- "Jun_6_2012"
  date.combine.data.run <- "Jun_17_2012"
  # tfs = c("BATF","MAF","IRF4","STAT3","RORC")
  # OR
  # tfs = c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")

  # filter out all scores K,C,R,I that have a KO support that doesn't lead to more than x% change of wt expression are set to zero!
  prcnt.chng.cut <- 0 
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
  filter.by.sam <- F
  combine.type.1 <- "th17" # "th17_minus_th0"

  combine.type.2 <- "KC"
  cut.abs <- 1.5 # get network with edges with score larger than this

  tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
  num.tfs <- length(tfs) # 8 or 5 (5=core.tfs,8=core.tfs+etv6,fosl2,hif1a)

  core.1.tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
  core.2.tfs <- sapply(tfs,function(i) strsplit(i,"_")[[1]][1])
  show.gn.names <- tfs
  core.tfs <- GLOBAL[["known.tfs"]]
  # cut.abs <- 1.5 # get network with edges with score larger than this
  source("r_scripts/th17/used_for_paper/createCytoscapeNetworkData.R")
}

if(any(GLOBAL[["run.these.steps"]]==14.2)) {
  cat("\n- creating cytoscape data structinos\n")
  for(iter in 1:4){
  rm(list=ls()[-which(ls() %in% c("GLOBAL","iter"))])
  date.is <- GLOBAL[["date.is"]]
  #date.combine.data.run <- "Jun_20_2012"
  date.combine.data.run <- "Jul_12_2012"
  core.1.tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
  core.tfs <- GLOBAL[["known.tfs"]]
  # filter out all scores K,C,R,I that have a KO support that doesn't lead to more than x% change of wt expression are set to zero!
  combine.type.1 <- "th17" # "th17_minus_th0"
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
	deseq.pval.cut <- 1 
	plot.qc=T
      if(iter == 1){
         prcnt.chng.cut <- 0;filter.by.sam <- F;combine.type.2 <- "KC";tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
      } else if(iter == 2){
         prcnt.chng.cut <- 0;filter.by.sam <- F;combine.type.2 <- "KC";tfs <- c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")
      } else if(iter == 3){
         prcnt.chng.cut <- 20;filter.by.sam <- F;combine.type.2 <- "KC";tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
      }  else if(iter == 4){
         prcnt.chng.cut <- 20;filter.by.sam <- F;combine.type.2 <- "KC";tfs <- c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")
      }

		cut.abs <- 1.5 # get network with edges with score larger than this
	  	core.2.tfs <- sapply(tfs,function(i) strsplit(i,"_")[[1]][1])
		num.tfs <- length(tfs) # 8 or 5 (5=core.tfs,8=core.tfs+etv6,fosl2,hif1a)
	  	show.gn.names <- tfs
		source("r_scripts/th17/used_for_paper/createCytoscapeNetworkData.R")
	}
}

if(any(GLOBAL[["run.these.steps"]]==14.3)) {
	 cat("\n- creating cytoscape data structinos\n")
	 rm(list=ls()[-which(ls()=="GLOBAL")])
	 date.is <- GLOBAL[["date.is"]]
	 # date.combine.data.run <- "Jun_6_2012"
	 date.combine.data.run <- "Aug_2_2012"
	  # filter out all scores K,C,R,I that have a KO support that doesn't lead to more than x% change of wt expression are set to zero!
	 prcnt.chng.cut <- 20
	 deseq.pval.cut <- 0.25
	 z.abs.cut <- GLOBAL[["z.abs.cut"]]
	 filter.by.sam <- F
	 prcnt.chng.cut <- 20
	 combine.type.1 <- "th17" # "th17_minus_th0"
	 cut.abs <- 1.5 # get network with edges with score larger than this
	 num.tfs <- 28
	 core.tfs <- GLOBAL[["known.tfs"]]
	 source("r_scripts/th17/used_for_paper/createCytoscapeNetworkData_fig7.R")
}


if(any(GLOBAL[["run.these.steps"]]==14.4)) {
  cat("\n- creating cytoscape data heatmaps\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  # date.combine.data.run <- "Jun_6_2012"
  date.combine.data.run <- "Jun_17_2012"
  date.is <- GLOBAL[["date.is"]]
  # tfs = c("BATF","MAF","IRF4","STAT3","RORC")
  # OR
  # tfs = c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")
  tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
  num.tfs <- length(tfs) # 8 or 5 (5=core.tfs,8=core.tfs+etv6,fosl2,hif1a)
  # cut.abs <- 1.5 # get network with edges with score larger than this
  combine.type.1 <- "th17"
  combine.type.2 <- "KC"
  core.2.tfs <- sapply(tfs,function(i) strsplit(i,"_")[[1]][1])
  show.gn.names <- tfs
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
  filter.by.sam <- F
  prcnt.chng.cut <- 0 
  make.binary.heatmap <- FALSE
  source("r_scripts/th17/used_for_paper/createCytoscapeNetworkData_heatmaps.R")
}
if(any(GLOBAL[["run.these.steps"]]==14.5)) {
  cat("\n- creating cytoscape data heatmaps\n")
  for(iter in 1:4){
  rm(list=ls()[-which(ls() %in% c("GLOBAL","iter"))])
  date.is <- GLOBAL[["date.is"]]
  #date.combine.data.run <- "Jun_20_2012"
  date.combine.data.run <- "Jul_12_2012"
  show.perm <- T
  # tfs = c("BATF","MAF","IRF4","STAT3","RORC")
  # OR
  # tfs = c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")
  core.1.tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
  core.tfs <- GLOBAL[["known.tfs"]]
  # filter out all scores K,C,R,I that have a KO support that doesn't lead to more than x% change of wt expression are set to zero!
  combine.type.1 <- "th17" # "th17_minus_th0"
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
	deseq.pval.cut <- 1
      if(iter == 1){
         prcnt.chng.cut <- 0;filter.by.sam <- F;combine.type.2 <- "KC";tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
      } else if(iter == 2){
         prcnt.chng.cut <- 0;filter.by.sam <- F;combine.type.2 <- "KC";tfs <- c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")
      } else if(iter == 3){
         prcnt.chng.cut <- 20;filter.by.sam <- F;combine.type.2 <- "KC";tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
      }  else if(iter == 4){
         prcnt.chng.cut <- 20;filter.by.sam <- F;combine.type.2 <- "KC";tfs <- c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")
      }
		cut.abs <- 1.5 # get network with edges with score larger than this
	  core.2.tfs <- sapply(tfs,function(i) strsplit(i,"_")[[1]][1])
		num.tfs <- length(tfs) # 8 or 5 (5=core.tfs,8=core.tfs+etv6,fosl2,hif1a)
	  show.gn.names <- tfs
  	make.binary.heatmap <- FALSE
  	source("r_scripts/th17/used_for_paper/createCytoscapeNetworkData_heatmaps.R")
	}
}

#----------------#
# 15-
#----------------#
# run net Validation script.
#    -- script:
#       - r_scripts/th17/used_for_paper/netValidation.R
#    -- input:
#       - input/th17/
#    -- output:
#       - results/th17/
if(any(GLOBAL[["run.these.steps"]]==15)) {
  cat("\n- calculating quality controls\n")
  types <- c("abs","k_activator","k_r_i_activator","k_r_i_repressor")
  for(t in 1:length(types)){
    rm(list=ls()[-which(ls() %in% c("GLOBAL","types","t"))])
    type <- types[t]  
    date.is <- GLOBAL[["date.is"]]
	gold.stdrd.date <- "May_31_2012"
    date.combine.data.run <- "Sep_21_2011"
    z.abs.cut <- GLOBAL[["z.abs.cut"]]
    core.tfs <- GLOBAL[["known.tfs"]]
    source("r_scripts/th17/used_for_paper/netValidation.R")
  }
}

#----------------#
# 15.1-
#----------------#
# create net validation plots:
# - datatype combos vs. performance/tf combos vs. performance
# - Diff expression vs. kc scores for both invivo and invitro samples
if(any(GLOBAL[["run.these.steps"]]==15.1)) {
	rm(list=ls()[-which(ls()=="GLOBAL")])
	cat("\n- calculating validation plots for targets recovery\n")
	date.is <- GLOBAL[["date.is"]]
    cut.dist <- 100 # only consider snps within 100kb
	gs.list <- "th17"
	#gs.list <- "type2diab"	
	gold.stdrd.date <- "Jun_22_2012"
	fl.nm.gs <- paste(sep="","gold_standard_",gs.list,"_genes_",gold.stdrd.date,".txt")
	#date.deseq.run <- "Jun_18_2012"
	# date.deseq.run <- "Aug_05_2012"
	# date.combine.data.run <- "Jun_3_2012"
	# date.combine.data.run <- "Jun_6_2012"	
	# date.combine.data.run <- "Jul_12_2012"	
	date.combine.data.run <- "Aug_5_2012"	
	core.tfs <- GLOBAL[["known.tfs"]]
	prcnt.chng.cut <- 0
	deseq.pval.cut <- 1
	# c.type.all.in.one.plot <- "th17_minus_th0"
	c.type.all.in.one.plot <- "th17"
	# tfs = c("BATF","MAF","IRF4","STAT3","RORC")
	# OR
	# tfs = c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")
	tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
	num.tfs <- length(tfs) # 8 or 5 (5=core.tfs,8=core.tfs+etv6,fosl2,hif1a)
	z.abs.cut <- GLOBAL[["z.abs.cut"]]
	filter.by.sam <- F
	comb.case <- "activation" # or activation or repression or absolute
	source("r_scripts/th17/used_for_paper/netValidation1_v2.R")
}

#----------------#
# 15.2-
#----------------#
# create net validation plots:
# - datatype combos vs. performance/tf combos vs. performance
# - Diff expression vs. kc scores for both invivo and invitro samples
if(any(GLOBAL[["run.these.steps"]]==15.2)) {
	gs.lists <- c("altz","sciz","asthma","celiac","type2diab","type1diab","SLE","psoriasis","CD","MS","RA","UC","th17")
	roc.sam.vec <- pr.sam.vec <- pr.rand.95.prcntile.vec <- roc.rand.95.prcntile.vec <- numeric(length(gs.lists))
	names(roc.sam.vec) <- names(pr.sam.vec) <- names(pr.rand.95.prcntile.vec) <- names(roc.rand.95.prcntile.vec) <- gs.lists
	m.final <- numeric() # here we will hold the results for all these gs lists (5 way combined results only)
	r.names <- c("AUCROC_BATF","AUCROC_MAF","AUCROC_IRF4","AUCROC_STAT3","AUCROC_RORC",
			 "AUCPR_BATF","AUCPR_MAF","AUCPR_IRF4","AUCPR_STAT3","AUCPR_RORC")
	m.final.2 <- matrix(0,nr=length(r.names),nc=length(gs.lists))
	colnames(m.final.2) <- gs.lists
	rownames(m.final.2) <- r.names
	for(iter in 1:length(gs.lists)){
	  	rm(list=ls()[-which(ls() %in% c("GLOBAL","iter","gs.lists","m.final","m.final.2",
			"pr.rand.95.prcntile.vec","roc.rand.95.prcntile.vec","roc.sam.vec","pr.sam.vec"))])
		cat("\n- calculating validation plots for targets recovery\n")
		date.is <- GLOBAL[["date.is"]]
		cut.dist <- 100 # only consider snps within 100kb
		# gs.list <- "th17"
		gs.list <- gs.lists[iter]	
		gold.stdrd.date <- "Jun_22_2012"
		fl.nm.gs <- paste(sep="","gold_standard_",gs.list,"_genes_",gold.stdrd.date,".txt")
		# date.combine.data.run <- "Jun_6_2012"	
		date.combine.data.run <- "Aug_8_2012"	
		core.tfs <- GLOBAL[["known.tfs"]]
		prcnt.chng.cut <- 0
		deseq.pval.cut <- 1
		# c.type.all.in.one.plot <- "th17_minus_th0"
		c.type.all.in.one.plot <- "th17"
		# tfs = c("BATF","MAF","IRF4","STAT3","RORC")
		# OR
		# tfs = c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")
		tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
		num.tfs <- length(tfs) # 8 or 5 (5=core.tfs,8=core.tfs+etv6,fosl2,hif1a)
		z.abs.cut <- GLOBAL[["z.abs.cut"]]
		filter.by.sam <- F
		comb.case <- "activation" # or activation or repression or absolute
		# source("r_scripts/th17/used_for_paper/netValidation1_v2.R")
		Np=15
		source("r_scripts/th17/used_for_paper/simulate_netValidation1_v2.R")
		roc.rand.95.prcntile.vec[gs.lists[iter]] <- roc.rand.95.prcntile
		pr.rand.95.prcntile.vec[gs.lists[iter]] <- pr.rand.95.prcntile
		pr.sam.vec[gs.lists[iter]] <- aucPR.sam
		roc.sam.vec[gs.lists[iter]] <- aucROC.sam
		
		#####################################
		# keep values for ploting later
		aucpr <- unlist(AUCPR.5.way)
		aucroc <- unlist(AUCROC.5.way)
		names(aucpr) <- names(aucroc) <- names(aucpr)
		m <- cbind(aucpr,aucroc)
		colnames(m) <- paste(sep="",gs.list ,"_",colnames(m))
		m.final <- cbind(m.final,m)
		tfs <- names(AUCROC)
		for(i in 1:length(tfs)){
			tf <- tfs[i]
			ix <- grep(paste(sep="","AUCPR_",tf),rownames(m.final.2)) 
			m.final.2[ix,gs.list] <- AUCPR[[tf]][["KCRI"]]
			ix <- grep(paste(sep="","AUCROC_",tf),rownames(m.final.2)) 
			m.final.2[ix,gs.list] <- AUCROC[[tf]][["KCRI"]]
		}
		#####################################	
	}
	stop("AM")
	rand.roc <- median(roc.rand.95.prcntile.vec)
	rand.pr <- median(pr.rand.95.prcntile.vec)
	sam.pr <- median(pr.sam.vec[which(names(pr.sam.vec)!="SAM")])
	sam.roc <- median(roc.sam.vec[which(names(roc.sam.vec)!="SAM")])
	source("r_scripts/th17/used_for_paper/make_gwas_validation_plot.R")	
}

#----------------#
# 16-
#----------------#
# create a tab delimited table summary of many smaller tables (taken the intersection of row names, i.e. for genes that appear in all tables)
#    -- script:
#       - r_scripts/th17/used_for_paper/summary.R
#    -- input:
#       - input/th17/
#    -- output:
#       - results/th17/

if(any(GLOBAL[["run.these.steps"]]==16)) {
  cat("\n- calculating quality controls\n")
  rm(list=ls()[-which(ls() %in% c("GLOBAL","types","t"))])
  date.is <- GLOBAL[["date.is"]]
  date.cytoscape.run <- "Oct_25_2011"
  date.diff.expsn.run <- "Sep_21_2011"
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
  date.chip.run <- "Sep_15_2011"
  stop("AM")
  source("r_scripts/th17/used_for_paper/summary.R")
  source("r_scripts/th17/used_for_paper/summary2.R")
}

#----------------#
# 17-
#----------------#
# create Figure 2 part 1
# Calculate clusters of peaks for tfs of interest (e.g. BATF and IRF4)
#    -- script:
#       - r_scripts/th17/used_for_paper/cluster_peaks)_for_fig2.R
#    -- input: 
#       - input/th17/used_for_paper/peaks_mat_List.RData
#    -- output:
#       - input/th17/used_for_paper/tf_clusters_per_loci_list.RData
if(any(GLOBAL[["run.these.steps"]]==17)) {
  cat("\n- calculating tfs chipseq clusters\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  tfs <-  GLOBAL[["known.tfs"]]
  tfs.to.cluster <- list()
  tfs.to.cluster[["BATF_IRF4"]] <- c("BATF","IRF4")
  tfs.to.cluster[["RORC_IRF4"]] <- c("RORC","IRF4")
  tfs.to.cluster[["RORC_STAT3"]] <- c("RORC","STAT3")
  tfs.to.cluster[["IRF4_STAT3"]] <- c("IRF4","STAT3")
  tfs.to.cluster[["BATF_STAT3"]] <- c("BATF","STAT3")
  tfs.to.cluster[["STAT3_IRF4"]] <- c("STAT3","IRF4")
  tfs.to.cluster[["STAT3_BATF"]] <- c("STAT3","BATF")
  tfs.to.cluster[[ paste(sep="",tfs,collapse="_") ]] <- tfs
  date.is <- GLOBAL[["date.is"]]
  # date.macs.data.compiled <- "Sep_15_2011"
  n.p <- 6 # number of processors to use
  match.targets <- TRUE # for figure 2 we don't care if peak is near a TSS or not
  mtbs.tss.dist <- GLOBAL[["mtbs.tss.dist"]]
  load.macs.per.chrom.RData <- FALSE
  d.cut <- 100
  mm9.effective.genome.size <- GLOBAL[["mm9.effective.genome.size"]]
  chr <- c(paste(sep="","chr",1:19),paste(sep="","chr",c("X","Y")))
  # for(type in 1:length(tfs.to.cluster)){
	type=1
	tfs.to.cluster.t <- tfs.to.cluster[[type]]
	# source("r_scripts/th17/used_for_paper/cluster_peaks_for_fig2_v2.R")
	# load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
	# load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th0_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_3_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th0_d.tfs_100_d.tss_5000_Apr_3_2012.RData")
	source("r_scripts/th17/used_for_paper/add_fc_and_pval_for_mtls_v2.R")
  # }
  # source("r_scripts/th17/used_for_paper/calc_pval_for_histone_marks_in_peak_clusters_v2.R")	
}

#----------------#
# 18-
#----------------#
# create Figure 2 part 2
# Calculate clusters of peaks for tfs of interest (e.g. BATF and IRF4)
#    -- script:
#       - r_scripts/th17/used_for_paper/cluster_peaks)_for_fig2.R
#    -- input: 
#       - input/th17/used_for_paper/peaks_mat_List.RData
#    -- output:
#       - input/th17/used_for_paper/tf_clusters_per_loci_list.RData
if(any(GLOBAL[["run.these.steps"]]==18)) {
  cat("\n- calculating tfs chipseq clusters\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  genetic.bg.trgt <- c("irf4_p300","stat3_p300","irf4_stat3","batf_stat3","stat3_irf4","stat3_batf","batf_irf4","irf4_batf","rorc_irf4","rorc_stat3","rorc_all") # small case letters pls :) first is bg other is associated TF to look at
  date.is <- GLOBAL[["date.is"]]
  n.p <- 6 # number of processors to use
  mtbs.tss.dist <- GLOBAL[["mtbs.tss.dist"]]
  mm9.effective.genome.size <- GLOBAL[["mm9.effective.genome.size"]]
  mtl.min.pval <- 300 # the minimum pval to include mtl in analysis (we only consider high significance mtls)
  core.tfs <- GLOBAL[["known.tfs"]]
  for(g in 1:length(genetic.bg.trgt)){
  	# source("r_scripts/th17/used_for_paper/calc_pval_for_histone_marks_in_peak_clusters_v3.R")
	source("r_scripts/th17/used_for_paper/calc_pval_for_histone_marks_in_peak_clusters_v4.R")
  }
}

#----------------#
# 19-
#----------------#
# create Figure 2 part 2
# Calculate clusters of peaks for tfs of interest (e.g. BATF and IRF4)
#    -- script:
#       - r_scripts/th17/used_for_paper/cluster_peaks)_for_fig2.R
#    -- input: 
#       - input/th17/used_for_paper/peaks_mat_List.RData
#    -- output:
#       - input/th17/used_for_paper/tf_clusters_per_loci_list.RData
if(any(GLOBAL[["run.these.steps"]]==19)) {
  cat("\n- calculating tfs chipseq clusters\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  tf.pair <- c("batf","irf4")
  date.is <- GLOBAL[["date.is"]]
  n.p <- 3 # number of processors to use
  mtbs.tss.dist <- GLOBAL[["mtbs.tss.dist"]]
  mm9.effective.genome.size <- GLOBAL[["mm9.effective.genome.size"]]
  mtl.min.pval <- 300 # the minimum pval to include mtl in analysis (we only consider high significance mtls)
  seq.count <- 500
  peak.region.seq.length <- 50
  lineage <- "Th17" # "Th0" or "Th17"
  source("r_scripts/th17/used_for_paper/getSequencesForMotifDiscoveryBatfIrf4.R")
}

#----------------#
# 20-
#----------------#
# create Figure 2 supplemental
#    -- script:
#       - r_scripts/th17/used_for_paper/createFig2.R
#    -- input:
#       - input/th17/
#    -- output:
#       - results/th17/
if(any(GLOBAL[["run.these.steps"]]==20)) {
  cat("\n- calculating quality controls\n")
  rm(list=ls()[-which(ls() %in% c("GLOBAL","types","t"))])
  date.is <- GLOBAL[["date.is"]]
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
  date.macs.processed <- "Sep_15_2011"
  pre.process <- FALSE
  n.p <- 4
  source("r_scripts/th17/used_for_paper/chipSpecificityPlots.R")
  # source("r_scripts/th17/used_for_paper/createFig2.R")
}

#----------------#
# 21-
#----------------#
# create Figure 2 supplemental
#    -- script:
#       - r_scripts/th17/used_for_paper/createFig2.R
#    -- input:
#       - input/th17/
#    -- output:
#       - results/th17/
if(any(GLOBAL[["run.these.steps"]]==21)) {
  cat("\n- calculating quality controls\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  date.combine.data.run <- "Apr_22_2012"
  date.is <- GLOBAL[["date.is"]]
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
  date.macs.processed <- "Sep_15_2011"
  pre.process <- FALSE
  n.p <- 4
  source("r_scripts/th17/used_for_paper/exploreMTLs.R")
}

#----------------#
# 22-
#----------------#
# run various tmp functions and scripts
#    -- script:
#       - r_scripts/th17/used_for_paper/tmp.R
#    -- input:
#       - input/th17/
#    -- output:
#       - results/th17/
if(any(GLOBAL[["run.these.steps"]]==22)) {
  cat("\n- calculating quality controls\n")
  rm(list=ls()[-which(ls() %in% c("GLOBAL","types","t"))])
  date.combine.data.run <- "May_10_2012"
  date.is <- GLOBAL[["date.is"]]
  cut.qunt <- .8 # get network for 20% top edges for core
  combine.type.1 <- "th17_w_p300_minus_th0_w_p300"
  combine.type.2 <- "KC"
  tfs <- c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","NFE2L2_ETOH","HIF1A","ETV6_KD_RPMI_24","KDM6B_KD_RPMI_24","JMJD6_KD_RPMI_24")
  core.1.tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
  core.2.tfs <- sapply(tfs,function(i) strsplit(i,"_")[[1]][1])
  show.gn.names <- tfs
  source("r_scripts/th17/used_for_paper/tmp1.R")
}

#----------------#
# 23-
#----------------#
# run various tmp functions and scripts
#    -- script:
#       - r_scripts/th17/used_for_paper/calc_tf_cor_with_diff_exp.R
#    -- input:
#       - input/th17/
#    -- output:
#       - results/th17/
if(any(GLOBAL[["run.these.steps"]]==23)) {
	rm(list=ls()[-which(ls()=="GLOBAL")])
	cat("\n- calculating validation plots for targets recovery\n")
	date.is <- GLOBAL[["date.is"]]
	# date.deseq.run <- "Jul_12_2012"
	date.combine.data.run <- "Jul_12_2012"	
	core.tfs <- GLOBAL[["known.tfs"]]
	z.abs.cut <- 0
	prcnt.chng.cut <- 20
	deseq.pval.cut <- 1
	c.type <- "th17"
	comb.case <- "whole" # or activation or repression or absolute
	comb.type <- "KC"
	comb.type.cut <- 1.5
	diff.exp.by.sam.or.fc <- "FC"
	# tfs = c("BATF","IRF4","STAT3","MAF","RORC")
	# OR
	# tfs = c("BATF","IRF4","STAT3","MAF","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")
	tfs <- c("BATF","IRF4","STAT3","MAF","RORC","FOSL2_TH17")
	num.tfs <- 8 # 8 or 5 (5=core.tfs,8=core.tfs+etv6,fosl2,hif1a)
    source("r_scripts/th17/used_for_paper/calc_tf_cor_with_diff_exp.R")
}

#----------------#
# 24-
#----------------#
# create diff chip vocanoplots looking at what happens to chip of x_i in ko of y_j
#    -- script:
#       - r_scripts/th17/used_for_paper/calc_tf_cor_with_diff_exp.R
#    -- input:
#       - input/th17/
#    -- output:
#       - results/th17/
if(any(GLOBAL[["run.these.steps"]]==24)) {
	cat("\n- calculating tfs chipseq clusters\n")
	rm(list=ls()[-which(ls()=="GLOBAL")])
	mm9.effective.genome.size <- GLOBAL[["mm9.effective.genome.size"]]
	date.is <- GLOBAL[["date.is"]]
	n.p <- 6 # number of processors to use
	mtbs.tss.dist <- GLOBAL[["mtbs.tss.dist"]]
	d.cut <- GLOBAL[["tfs.min.dist"]]
	mtl.min.pval <- 300 # the minimum pval to include mtl in analysis (we only consider high significance mtls)
	core.tfs <- GLOBAL[["known.tfs"]]
	core.tfs <- c("IRF4","STAT3","RORC","BATF","MAF","P300")
	count.reads.per.mtl <- F
	date.last.count <- "Jul_24_2012"
	expts.file.nm <- "input/th17/used_for_paper/diff_chip_volcanoplots/map_sam_files_to_expt_Jul_23_2012.txt"
	if(count.reads.per.mtl==T){
		# add counts for mtls (th17 and th0)
		source("r_scripts/th17/used_for_paper/add_fc_and_pval_for_expts_to_mtl_list.R")
		source("r_scripts/th17/used_for_paper/calc_mapped_reads_for_expts.R")
		date.last.count <- date.is
	}
	# which genes to highlight
	highlight.genes <- c("IL17A","IL17F","IL23R","FURIN","CXCL3","CCL20","IL1R1","IL12RB2","IL22")		
	# do you want to calc adjustable lambda (more strict)
	calc.local.lambda  <- T
	source("r_scripts/th17/used_for_paper/create_volcanoplots_diff_chip_v2.R")
}

#----------------#
# 25-
#----------------#
# simulate k,c,r,i networks to estimate FDRs and random scores
#    -- script:
#       - r_scripts/th17/used_for_paper/simulate_scores.R
#    -- input:
#       - input/th17/
#    -- output:
#       - results/th17/
if(any(GLOBAL[["run.these.steps"]]==25)) {
	rm(list=ls()[-which(ls()=="GLOBAL")])
	cat("\n- running simulated scores to get FDRs \n")
	date.is <- GLOBAL[["date.is"]]
	cut.dist <- 100 # only consider snps within 100kb
	gs.list <- "th17"
	gold.stdrd.date <- "Jun_22_2012"
	# date.deseq.run <- "Jul_11_2012"
	date.combine.data.run <- "Jul_12_2012"	
	core.tfs <- GLOBAL[["known.tfs"]]
	prcnt.chng.cut <- 0
	deseq.pval.cut <- 1
	c.type.all.in.one.plot <- "th17"
	# tfs = c("BATF","MAF","IRF4","STAT3","RORC")
	# OR
	# tfs = c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")
	tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
	num.tfs <- length(tfs) # 8 or 5 (5=core.tfs,8=core.tfs+etv6,fosl2,hif1a)
	z.abs.cut <- GLOBAL[["z.abs.cut"]]
	min.rpkm <- GLOBAL[["min.th17.or.th0.rpkm"]]
	prcnt.chng.cut <- 20
	deseq.pval.cut <- 0.25
	filter.by.sam <- F
	comb.case <- "activation" # or activation or repression or absolute
	Np <- 2 # number of random sets
	source("r_scripts/th17/used_for_paper/simulate_scores.R")
}

#----------------#
# 26-
#----------------#
# simulate validation plots to get random expected aucPR and aucROC
#    -- script:
#       - r_scripts/th17/used_for_paper/simulate_netValidation1_v2
#    -- input:
#       - input/th17/
#    -- output:
#       - results/th17/
if(any(GLOBAL[["run.these.steps"]]==26)) {
	rm(list=ls()[-which(ls()=="GLOBAL")])
	cat("\n- calculating validation plots for targets recovery\n")
	date.is <- GLOBAL[["date.is"]]
    cut.dist <- 100 # only consider snps within 100kb
	gs.list <- "th17"
	#gs.list <- "type2diab"	
	gold.stdrd.date <- "Jun_22_2012"
	fl.nm.gs <- paste(sep="","gold_standard_",gs.list,"_genes_",gold.stdrd.date,".txt")
	#date.deseq.run <- "Jun_18_2012"
	# date.deseq.run <- "Jul_12_2012"
	# date.deseq.run <- "Jul_10_2012"
	# date.combine.data.run <- "Jul_12_2012"	
	date.combine.data.run <- "Aug_8_2012"	
	core.tfs <- GLOBAL[["known.tfs"]]
	prcnt.chng.cut <- 0
	deseq.pval.cut <- 1
	c.type.all.in.one.plot <- "th17"
	# tfs = c("BATF","MAF","IRF4","STAT3","RORC")
	# OR
	# tfs = c("BATF","MAF","IRF4","STAT3","RORC","FOSL2_TH17","HIF1A","ETV6_KD_RPMI_24")
	tfs <- c("BATF","MAF","IRF4","STAT3","RORC")
	num.tfs <- length(tfs) # 8 or 5 (5=core.tfs,8=core.tfs+etv6,fosl2,hif1a)
	z.abs.cut <- GLOBAL[["z.abs.cut"]]
	filter.by.sam <- F
	comb.case <- "activation" # or activation or repression or absolute
	Np=15
	source("r_scripts/th17/used_for_paper/simulate_netValidation1_v2.R")
}

#----------------#
# 24-
#----------------#
# find essential targets for Rorc
#    -- script:
#       - r_scripts/th17/used_for_paper/find_rorc_essential_targets.R
#    -- input:
#       - input/th17/
#    -- output:
#       - results/th17/
if(any(GLOBAL[["run.these.steps"]]==27)) {
	cat("\n- calculating tfs chipseq clusters\n")
	rm(list=ls()[-which(ls()=="GLOBAL")])
	date.is <- GLOBAL[["date.is"]]
	date.deseq.run <- "Aug_2_2012"
	examined.tf <- c("RORC")
	expts.file.nm <- "results/diff_chip_volcanoplots/Jul_27_2012/th17_wt_vs_ko_min_mtl_pval_300_Jul_27_2012_rorc.xls"
	source("r_scripts/th17/used_for_paper/find_rorc_essential_targets.R")
}














