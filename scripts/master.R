#!/usr/bin/env Rscript
# Master script for processing TF-target gene interaction data into network files for Cytoscape 
# (see Computational Methods paper, http://www.cell.com/action/showImagesData?pii=S0092-8674%2812%2901123-3 ) 
# Aviv Madar has provided R-scripts from the original project. Some procedures have been taken and/or adapted 
# from his code.
#
# Input: 
#       C = customized MACS output (ChIP-seq - provided on GEO) 
#       K = DEseq files for core TFs (RNA-seq KO - provided on GEO but should be exchangeable) 
#       R = Inferelator matrix (RNA-seq compendium - provided on GEO)
#       I = Inferelator matrix (2011 Immgen microarray data - provided on GEO)

# Add required libraries
source("addLibraries.R")
# Load global variables
source("setGlobalVars.R")
# Load utility functions
source("util.R")

# Input file directories
scriptdir <- getwd()
rnaseqfile <- paste0(getwd(), "/../suppl/data/inferelator/GSE40918_Inferelator_RNAseq.txt")
immgenfile <- paste0(getwd(), "/../suppl/data/inferelator/GSE40918_Inferelator_Immgen.txt")
zscores_filepath <- paste0(getwd(), "/../suppl/mmc5.xls")

# Debug output directory
outpath.debug <- paste0(getwd(), "/../suppl/data/analysis/debug/")
if(!dir.exists(outpath.debug)) {
  dir.create(outpath.debug)
}

# Cytoscape output directory
outpath.cyt <- paste0(getwd(), "/../suppl/data/analysis/cyt/")
if(!dir.exists(outpath.cyt)) {
  dir.create(outpath.cyt)
}

# Get command line arguments
# Possible options:
# noload - skip generating new ChIP-seq and DEseq (knockout) confidence score matrices and attempt to load existing files
args = commandArgs(trailingOnly=TRUE)
if(length(args) > 1) {
  stop("Only one argument allowed. Stopping.")
}

# ChIP and DEseq data is generated from GEO input, if this is true.
# If the user enters the 'noload' option, existing matrices will be loaded (if they exist)
shouldRegenerateKCData <- TRUE

# Can only be length 0 or 1 here
if(length(args) == 1) {
  if(args[[1]] == "noload") {
    shouldRegenerateKCData <- FALSE
    print("Attemtpting to use existing data for knockout and ChIP data.")
  } else {
    stop("Option not recognized. Available options: noload")
  }
} else {
  print("Newly generating the KC data.")
}

# START OF PROCEDURE
# --------------
# 1) Load data from each selected data type to create confidence score matrix S
# --------------
print("--------------------------------------------------------------------")
print(">>>>>>>>>>>>>>>>> 1) Loading all initial data <<<<<<<<<<<<<<<<<<<<<<")
print("--------------------------------------------------------------------")
# Load z-scores from SAM stored in mmc5 table of original authors (available on the Cell page, link on top of this script). 
# Z-scores provide a color scheme for nodes in Cytoscape which shows differential expression Th17 vs Th0 after 48h.
zscores.all <- load.zscores(zscores_filepath)

# Z-scores may also be used for filtering of included genes
genes.final <- filter.genes.by.zscore(zscores.all, GLOBAL[["z.abs.cut"]])
print(paste("Total number of genes with z-scores:", length(rownames(zscores.all)), "--- Genes with abs(zscore) > 2.50:", length(genes.final)))
  
ko.scores <- NULL
chip.scores.activator <- NULL
chip.scores.repressor <- NULL
rna.scores <- NULL
immgen.scores <- NULL

# ChIP always has positive values
k_sign_mat <- NULL
r_sign_mat <- NULL
i_sign_mat <- NULL

# Load RNA-seq knockout scores
loadKOData <- function() {
  setwd(scriptdir)
  source(paste0(getwd(), "/" , "deseqExtract-fun.R"))
  print("Generating knockout scores.")
  scores <- load.deseq()
  k_sign_mat <- sign(as.data.frame(scores))
  if(GLOBAL[["DEBUG"]]) write.mat(scores, outpath.debug, "K", "smat")
  return(scores)
}

if(shouldRegenerateKCData) {
  ko.scores <- loadKOData()
} else {
  print("Looking for existing knockout data matrix.")
  ko.file <- paste0(outpath.debug, "K_smat.txt")
  if(file.exists(ko.file)) {
    print("Found knockout matrix file.")
    ko.scores <- as.matrix(read.table(ko.file, sep="\t", header=TRUE, row.names = 1))
    k_sign_mat <- sign(ko.scores)
  } else {
    print("Could not find knockout data matrix. Generating new matrix.")
    ko.scores <- loadKOData()
  }
}
 
# Load ChIP-seq scores  
loadChIPData <- function(type) {
  setwd(scriptdir)
  source(paste0(getwd(), "/" , "chipExtract-fun.R"))
  print("Loading ChIP scores.")
  if(type == "activator") {
    scores <- load.chip(boost.p300 = TRUE)
    if(GLOBAL[["DEBUG"]]) write.mat(scores, outpath.debug, "C", "activator_smat")
  } else if(type == "repressor") {
    scores <- load.chip(boost.p300 = FALSE)
    if(GLOBAL[["DEBUG"]]) write.mat(scores, outpath.debug, "C", "repressor_smat")
  } else {
    stop("Unknown type when loading ChIP data.")
  }
  return(scores)
}

if(shouldRegenerateKCData) {
  chip.scores.activator <- loadChIPData("activator")
  chip.scores.repressor <- loadChIPData("repressor")
} else {
  # Activator
  print("Looking for existing ChIP activator data matrix.")
  chip.activator.file <- paste0(outpath.debug, "C_activator_smat.txt")
  if(file.exists(chip.activator.file)) {
    print("Found ChIP activator matrix file.")
    chip.scores.activator <- as.matrix(read.table(chip.activator.file, sep="\t", header=TRUE, row.names = 1))
  } else {
    print("Could not find ChIP activator data matrix. Generating new matrix.")
    chip.scores.activator <- loadChIPData("activator")
  }
  
  # Repressor
  print("Looking for existing ChIP repressor data matrix.")
  chip.repressor.file <- paste0(outpath.debug, "C_repressor_smat.txt")
  if(file.exists(chip.repressor.file)) {
    print("Found ChIP repressor matrix file.")
    chip.scores.repressor <- as.matrix(read.table(chip.repressor.file, sep="\t", header=TRUE, row.names = 1))
  } else {
    print("Could not find ChIP repressor data matrix. Generating new matrix.")
    chip.scores.repressor <- loadChIPData("repressor")
  }
}

# Load RNA-compendium Inferelator scores - directly provided on GEO
print("Loading RNA-seq inferelator scores.")
rna.scores <- as.matrix(read.table(rnaseqfile, sep="\t", header=TRUE, row.names = 1))
r_sign_mat <- sign(as.data.frame(rna.scores))

# Load Immgen microarray Inferelator scores - directly provided on GEO 
print("Loading Immgen microarray inferelator scores.")
immgen.scores <- as.matrix(read.table(immgenfile, sep="\t", header=TRUE, row.names = 1))
i_sign_mat <- sign(as.data.frame(immgen.scores))

# --------------
# 2) Perform ranking on each confidence score matrix S
# --------------
print("--------------------------------------------------------------------")
print(">>>>>>>>>>>>> 2) Perform quantile ranking on each data set <<<<<<<<<<<<<<<<<<")
print("--------------------------------------------------------------------")
# Reset because previous functions may globally change working directory and source() breaks
setwd(scriptdir)

# Wrapper for performing ranking. Write ranked matrix if desired.
do.quantile.rank <- function(mat, positiveOnly = FALSE, prefix) {
  if(!is.null(mat)) {
    qmat <- calc.quantile.ranks(mat, positiveOnly = positiveOnly)
    if(GLOBAL[["DEBUG"]]) write.mat(qmat, outpath.debug, prefix, "_qmat")
    return(qmat)
  }
  return(NULL)
}

if(GLOBAL[["use.nyu.rank"]]) {
  print(">>>>>>>>>>>>>>>>> DEBUG: Calculating Aviv Madar's rank-score algorithm (in place of own ranking) <<<<<<<<<<<<<<<<<<<<<<")
  # Utilizing ranking methods from AM for testing purposes
  source(paste0(getwd(), "/external/rscripts/rscripts/" , "util.R"))
  
  # KO scores
  print("Ranking knockout scores.")
  print(ko.scores[1:3,])
  ko_qmat.activator <- abs(convert.scores.to.relative.ranks.pos(ko.scores))
  ko_qmat.repressor <- abs(convert.scores.to.relative.ranks.pos(-1*ko.scores))
  write.mat(ko_qmat.activator, outpath.debug, "K", "activator_nyu_qmat")
  write.mat(ko_qmat.repressor, outpath.debug, "K", "repressor_nyu_qmat")
  
  # ChIP scores
  print("Ranking ChIP scores.")
  chip_qmat.activator <- abs(convert.scores.to.relative.ranks(chip.scores.activator)) # all values!
  chip_qmat.repressor <- abs(convert.scores.to.relative.ranks(chip.scores.repressor)) # all values!
  write.mat(chip_qmat.activator, outpath.debug, "C", "activator_nyu_qmat")
  write.mat(chip_qmat.repressor, outpath.debug, "C", "repressor_nyu_qmat")
  
  # RNA compendium scores
  print("Ranking RNA compendium scores.")
  rna_qmat.activator <- abs(convert.scores.to.relative.ranks.pos(rna.scores))
  rna_qmat.repressor <- abs(convert.scores.to.relative.ranks.pos(-1*rna.scores))
  write.mat(rna_qmat.activator, outpath.debug, "R", "activator_nyu_qmat")
  write.mat(rna_qmat.repressor, outpath.debug, "R", "repressor_nyu_qmat")
  
  # Immgen microarray scores
  print("Ranking Immgen microarray scores.")
  immgen_qmat.activator <- abs(convert.scores.to.relative.ranks.pos(immgen.scores))
  immgen_qmat.repressor <- abs(convert.scores.to.relative.ranks.pos(-1*immgen.scores))
  write.mat(immgen_qmat.activator, outpath.debug, "I", "activator_nyu_qmat")
  write.mat(immgen_qmat.repressor, outpath.debug, "I", "repressor_nyu_qmat")
  
} else {
  source(paste0(getwd(), "/" , "quantileRank-fun.R"))
  print("Creating quantil rank matrices (Q).")
  print("Ranking knockout scores.")
  ko_qmat.activator <- do.quantile.rank(ko.scores, positiveOnly = TRUE, "K_activator") # positive only!
  ko_qmat.repressor <-do.quantile.rank(-1*ko.scores, positiveOnly = TRUE, "K_repressor") # positive only!
  
  print("Ranking ChIP scores.")
  chip_qmat.activator <- do.quantile.rank(chip.scores.activator, positiveOnly = FALSE, "C_activator")
  chip_qmat.repressor <- do.quantile.rank(chip.scores.repressor, positiveOnly = FALSE, "C_repressor")
  
  print("Ranking RNA compendium scores.")
  rna_qmat.activator <- do.quantile.rank(rna.scores, positiveOnly = TRUE, "R_activator") # positive only!
  rna_qmat.repressor <- do.quantile.rank(-1*rna.scores, positiveOnly = TRUE, "R_repressor") # positive only!
  
  print("Ranking Immgen microarray scores.")
  immgen_qmat.activator <- do.quantile.rank(immgen.scores, positiveOnly = TRUE, "I_activator") # positive only!
  immgen_qmat.repressor <- do.quantile.rank(-1*immgen.scores, positiveOnly = TRUE, "I_repressor") # positive only!
}

# --------------
# 3) Combine data according to various data type combinations
# --------------
print("--------------------------------------------------------------------")
print(">>>>>>> Combining Q-matrices to a single interaction matrix <<<<<<<<")
print("--------------------------------------------------------------------")
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "createCombinedMat-fun.R"))

# KC
kc.activator <- createCombinedMat(combo = "kc", type = "activator", ko_qmat = ko_qmat.activator, chip_qmat = chip_qmat.activator, 
                                  rna_qmat = NULL, immgen_qmat = NULL, genes.final)
kc.repressor <- createCombinedMat(combo = "kc", type = "repressor", ko_qmat = ko_qmat.repressor, chip_qmat = chip_qmat.repressor, 
                                  rna_qmat = NULL, immgen_qmat = NULL, genes.final)

# KCRI
kcri.activator <- createCombinedMat(combo = "kcri", type = "activator", ko_qmat = ko_qmat.activator, chip_qmat = chip_qmat.activator, 
                                    rna_qmat = rna_qmat.activator, immgen_qmat = immgen_qmat.activator, genes.final)
kcri.repressor <- createCombinedMat(combo = "kcri", type = "repressor", ko_qmat = ko_qmat.repressor, chip_qmat = chip_qmat.repressor, 
                                    rna_qmat = rna_qmat.repressor, immgen_qmat = immgen_qmat.repressor, genes.final)

# --------------
# 4) Write a copy of mmc5 Th17 vs. Th0 (both at 48h) z-scores to a table, which should be loaded in Cytoscape as 'Node table' 
# --------------
print("--------------------------------------------------------------------")
print(">>>>>>>>>>>>>>>>>>>>> Writing z-score table <<<<<<<<<<<<<<<<<<<<<<<<")
print("--------------------------------------------------------------------")
print("Writing zscores from mmc5.")
write.mat(zscores.all, outpath.cyt, "", "zscores")

# --------------
# 5) From combine data matrix, create a list of node-node-value interactions for Cytoscape
# --------------
print("--------------------------------------------------------------------")
print(">>>>>>>>> Writing interaction lists for Cytoscape network <<<<<<<<<<")
print("--------------------------------------------------------------------")
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "createInteractions-fun.R"))

# Positive edges are surpassing the cut value for the edges. They are tagged as 'positive_KC' if >1.65 (or any other defined value) is about activation, and as 'negative_KC' if >1.50 is
# about repression. This has been imitated from Aviv Madar's original code in an effort to get the procedure right, but I find the variable handling and
# naming very inconvenient and confusing. If there is time, I will change this.

# KC
print("Writing KC interactions as single list...")
create.interactions(kc.activator, outpath.cyt, "kc", "single", pos.edge= "positive_KC", neg.edge = "negative_KC", append = FALSE)
create.interactions(kc.repressor, outpath.cyt, "kc", "single", pos.edge= "negative_KC", neg.edge = "positive_KC", append = TRUE)

print("Writing KC interactions as separate lists...")
create.interactions(kc.activator, outpath.cyt, "kc", "activator", pos.edge= "positive_KC", neg.edge = "negative_KC", append = FALSE)
create.interactions(kc.repressor, outpath.cyt, "kc", "repressor", pos.edge= "negative_KC", neg.edge = "positive_KC", append = FALSE)

# KCRI
print("Writing KCRI interactions as single list...")
create.interactions(kcri.activator, outpath.cyt, "kcri", "single", pos.edge= "positive_KCRI", neg.edge = "negative_KCRI", append = FALSE)
create.interactions(kcri.repressor, outpath.cyt, "kcri", "single", pos.edge= "negative_KCRI", neg.edge = "positive_KC", append = TRUE)

print("Writing KCRI interactions as separate lists...")
create.interactions(kcri.activator, outpath.cyt, "kcri", "activator", pos.edge= "positive_KCRI", neg.edge = "negative_KCRI", append = FALSE)
create.interactions(kcri.repressor, outpath.cyt, "kcri", "repressor", pos.edge= "negative_KCRI", neg.edge = "positive_KCRI", append = FALSE)

print("---------------------")
print("Done. Now files can be loaded into Cytoscape.")