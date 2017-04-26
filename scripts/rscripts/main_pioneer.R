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
GLOBAL[["run.these.steps"]] <- 6 # for a complete run c(1:last_step)
# stamp the date on this run
x <- unlist(strsplit(date()," +",perl=TRUE))
GLOBAL[["date.is"]] <- paste(x[2],x[3],x[5],sep="_")
# data set size: take genes with abs(zscore) based on diff expression btwn th17_48hr and th0_48hr
GLOBAL[["z.abs.cut"]] <- 2.5 #zscore
# Only used to determine rnaseq data set size
GLOBAL[["median.abs.cut"]] <- 3 #rpkm
# The number of bootstraps for inferelator rnaseq and immgen
GLOBAL[["num.boots.rnaseq"]] <- 200
GLOBAL[["num.boots.immgen"]] <- 200
GLOBAL[["num.perms.enrichment.analysis"]] <- 100
# do we want to output the sequences around the peaks?
GLOBAL[["get.sequence"]] <- TRUE
# what is the minimum distance between peaks to cluster them together
## GLOBAL[["tfs.min.dist"]] <- c(50,100,150,200,250,300,500,1000)
GLOBAL[["tfs.min.dist"]] <- c(100)
# multi TF binding sites distance (up/downstrem) from TSS for gene MTBS mapping
GLOBAL[["mtbs.tss.dist"]] <- 5000
GLOBAL[["mm9.effective.genome.size"]] <- 1.87e9

# genes that we want in as a rule even if they have low median expression or z.score that kicks them out
GLOBAL[["known.tfs"]] <- toupper(c("Irf4","Stat3","Rorc","Batf","Maf"))
GLOBAL[["known.gns"]] <- toupper(c("Il17a","Il17f","Il23r","Il22","Il21","Foxp3","Rora"))
GLOBAL[["use.multicore"]] <- TRUE
if(GLOBAL[["use.multicore"]]==TRUE){
  library(multicore)
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
  source("r_scripts/th17/used_for_paper/calc_chip_scores_pioneer.R")
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
  mtbs.tss.dist <- GLOBAL[["mtbs.tss.dist"]] #<- Inf
  n.p <- 4 # number of processors to use
  n.b <- 4 # number of bootstraps
  load.macs.per.chrom.RData <- FALSE
  add.histone.scores <- FALSE
  match.targets <- TRUE # do you want to include the targets of each cluster?
  calc.simulations <- FALSE # do you want to run 10000 simulations to determine significance?
  for(iter in 1:length(GLOBAL[["tfs.min.dist"]])){
    d.cut <- GLOBAL[["tfs.min.dist"]][iter]
    cat("working on d.cut ",d.cut,"\n")
    source("r_scripts/th17/used_for_paper/cluster_peaks_pioneer.R")
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
  add.histone.scores <- TRUE
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
  # for rnaseq or immgen?
  GLOBAL[["type"]] <- "rnaseq"
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
  # for rnaseq or immgen?
  GLOBAL[["type"]] <- "immgen"
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
  source("r_scripts/th17/used_for_paper/newTFsEnrichmentAnalysis.R")
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
   ## show.gn.names <- toupper(c("BATF","MAF","IRF4","STAT3","RORC","Ahr",
  ##                            "Runx1","RORa","IL2RA","IL2RB","Got1",
  ##                            "Pank2","HK1","IL21","Il17f","Il17a","Il23R")) # this was decided by maria and I on March 16th via email
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
    date.combine.data.run <- "Sep_21_2011"
    z.abs.cut <- GLOBAL[["z.abs.cut"]]
    core.tfs <- GLOBAL[["known.tfs"]]
    source("r_scripts/th17/used_for_paper/netValidation.R")
  }
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
  tfs.to.cluster <- list()
  tfs.to.cluster[["BATF_IRF4"]] <- c("BATF","IRF4")
  tfs.to.cluster[["RORC_IRF4"]] <- c("RORC","IRF4")
  tfs.to.cluster[["RORC_STAT3"]] <- c("RORC","STAT3")
  tfs.to.cluster[[ paste(sep="",GLOBAL[["known.tfs"]],collapse="_") ]] <- GLOBAL[["known.tfs"]]
  date.is <- GLOBAL[["date.is"]]
  date.macs.data.compiled <- "Sep_15_2011"
  n.p <- 6 # number of processors to use
  match.targets <- TRUE # for figure 2 we don't care if peak is near a TSS or not
  mtbs.tss.dist <- GLOBAL[["mtbs.tss.dist"]]
  load.macs.per.chrom.RData <- FALSE
  d.cut <- 100
  mm9.effective.genome.size <- GLOBAL[["mm9.effective.genome.size"]]
  # for(type in 1:length(tfs.to.cluster)){
	type=1
	tfs.to.cluster.t <- tfs.to.cluster[[type]]
	source("r_scripts/th17/used_for_paper/cluster_peaks_for_fig2_v2.R")
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
  genetic.bg.trgt <- c("batf_irf4","irf4_batf","rorc_irf4","rorc_stat3","rorc_all") # small case letters pls :) first is bg other is associated TF to look at
  date.is <- GLOBAL[["date.is"]]
  n.p <- 6 # number of processors to use
  mtbs.tss.dist <- GLOBAL[["mtbs.tss.dist"]]
  mm9.effective.genome.size <- GLOBAL[["mm9.effective.genome.size"]]
  mtl.min.pval <- 300 # the minimum pval to include mtl in analysis (we only consider high significance mtls)
  core.tfs <- GLOBAL[["known.tfs"]]
  for(g in 1:length(genetic.bg.trgt)){
  	source("r_scripts/th17/used_for_paper/calc_pval_for_histone_marks_in_peak_clusters_v3.R")
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
  rm(list=ls()[-which(ls() %in% c("GLOBAL","types","t"))])
  date.is <- GLOBAL[["date.is"]]
  z.abs.cut <- GLOBAL[["z.abs.cut"]]
  date.macs.processed <- "Sep_15_2011"
  pre.process <- FALSE
  n.p <- 4
  source("r_scripts/th17/used_for_paper/exploreMTLs.R")
}














