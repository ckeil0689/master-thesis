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

# Perform some intial setups
list.of.packages <- c("data.table", "reshape2")#, "xlsx")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# install if missing
if(length(new.packages)) {
  print("Installing data.table package...")
  install.packages(new.packages, repos="http://cran.rstudio.com/")
}

library(data.table)
library(reshape2)
# library(xlsx) // causes issues with rJava package installs which rely on specific Java dev kit installs --> unreliable mess

# Set some global variables
GLOBAL <- list()
GLOBAL[["DEBUG"]] <- FALSE
GLOBAL[["z.abs.cut"]] <- 0.00 # carried over from original Aviv Madar code

# Input file directories
scriptdir <- getwd()
deseqdir <- paste0(getwd(), "/../suppl/data/deseq/")
chipdir <- paste0(getwd(), "/../suppl/data/chipseq/")
rnaseqfile <- paste0(getwd(), "/../suppl/data/inferelator/GSE40918_Inferelator_RNAseq.txt")
immgenfile <- paste0(getwd(), "/../suppl/data/inferelator/GSE40918_Inferelator_Immgen.txt")
ref_filepath <- paste0(getwd(), "/../suppl/mmc4.csv")
zscores_filepath <- paste0(getwd(), "/../suppl/mmc5.xls")

# Put all output into analysis folder
outpath.debug <- paste0(getwd(), "/../suppl/data/analysis/debug/")
if(!dir.exists(outpath.debug)) {
  dir.create(outpath.debug)
}

# What will eventually be used in Cytoscape
outpath.cyt <- paste0(getwd(), "/../suppl/data/analysis/cyt/")
if(!dir.exists(outpath.cyt)) {
  dir.create(outpath.cyt)
}

# Core target transcription factors
# CORE_TFS <- c("batf", "irf4", "stat3", "maf", "rorc")
CORE_TFS <- c("batf", "irf4", "stat3", "maf", "rorc", "fosl2", "hif1a")
   
# Load confidence score matrices by option
write.mat <- function(mat, outpath, prefix, suffix) {
  filename = paste0(outpath, prefix, suffix, ".txt")
  print(paste("Writing matrix to file:", filename))
  write.table(mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)
}

# --------------
# 1) Load data from each selected data type to create confidence score matrix S
# --------------
# Load z-scores from SAM stored in mmc5 table of original authors (available on the Cell page, link on top of this script). 
# Z-scores provide a color scheme for nodes in Cytoscape which shows differential expression Th17 vs Th0 after 48h.
zscore.table <- read.table(zscores_filepath, sep="\t", header=TRUE)
zscore_col <- "Th17.Th0.zscores"
zscores.all <- as.matrix(zscore.table[, zscore_col])
rownames(zscores.all) <- zscore.table[, "Gene_id"]
colnames(zscores.all) <- "Th17_vs_Th0_Zscores" # matches name used in Cytoscape Style for example KC.cys

# Z-scores may also be used for filtering of included genes
genes.final.idx <- which(abs(zscores.all) > GLOBAL[["z.abs.cut"]])
genes.final <- zscore.table[genes.final.idx, "Gene_id"]

print(paste("Original gene number:", length(zscore.table[, "Gene_id"]), "--- Genes with abs(zscore) > 2.50:", length(genes.final)))
print(paste())

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
setwd(scriptdir)
source(paste0(getwd(), "/" , "deseqExtract-fun.R"))
ko.scores <- load.deseq(dir = deseqdir, CORE_TFS)
k_sign_mat <- sign(as.data.frame(ko.scores))
if(GLOBAL[["DEBUG"]]) write.mat(ko.scores, outpath.debug, "K", "_smat")
 
# Load ChIP-seq scores  
setwd(scriptdir)
source(paste0(getwd(), "/" , "chipExtract-fun.R"))
chip.scores.activator <- load.chip(dir = chipdir, reflibfile = ref_filepath, boost.p300 = TRUE, CORE_TFS)
chip.scores.repressor <- load.chip(dir = chipdir, reflibfile = ref_filepath, boost.p300 = FALSE, CORE_TFS)
if(GLOBAL[["DEBUG"]]) write.mat(chip.scores.activator, outpath.debug, "C", "_activator_smat")
if(GLOBAL[["DEBUG"]]) write.mat(chip.scores.repressor, outpath.debug, "C", "_repressor_smat")

# Load RNA-compendium Inferelator scores - directly provided on GEO
rna.scores <- as.data.frame(read.table(rnaseqfile, sep="\t", header=TRUE))
# remove first column or sign()-function will explode
rownames(rna.scores) <- rna.scores[, 1]
rna.scores[, 1] <- NULL
r_sign_mat <- sign(as.data.frame(rna.scores))

# Load Immgen microarray Inferelator scores - directly provided on GEO 
immgen.scores <- as.data.frame(read.table(immgenfile, sep="\t", header=TRUE))
# remove first column or sign()-function will explode
rownames(immgen.scores) <- immgen.scores[, 1]
immgen.scores[, 1] <- NULL
i_sign_mat <- sign(as.data.frame(immgen.scores))

# --------------
# 2) Perform ranking on each confidence score matrix S
# --------------
# Reset because previous functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "rankSmat-fun.R"))

# Wrapper for performing ranking. Write ranked matrix if desired.
do.rank <- function(mat, prefix) {
  if(!is.null(mat)) {
    mat_ranked <- rank.smat(mat)
    if(GLOBAL[["DEBUG"]]) write.mat(mat_ranked, outpath.debug, prefix, "_ranked")
    return(mat_ranked)
  }
  return(NULL)
}

print("Creating rank matrices.")
ko.scores_ranked <- do.rank(ko.scores, "K")
chip.scores.activator_ranked <- do.rank(chip.scores.activator, "C_activator")
chip.scores.repressor_ranked <- do.rank(chip.scores.repressor, "C_repressor")
rna.scores_ranked <- do.rank(rna.scores, "R")
immgen.scores_ranked <- do.rank(immgen.scores, "I")

# --------------
# 3) Calculate quantiles for each ranked matrix to obtain the Q-matrix.
# --------------
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "qmat-fun.R"))

# Wrapper for calculating quantile scores from ranked matrices.
# @param mat - The original confidence score S-matrix for the data type
do.qcalc <- function(scores, scores_ranked, prefix) {
  if(!is.null(scores_ranked)) {
    qmat <- calc.qmat(scores, scores_ranked)
    if(GLOBAL[["DEBUG"]]) write.mat(qmat, outpath.debug, prefix, "_qmat")
    return(as.matrix(qmat))
  }
  return(NULL)
}

print("Calculating Q-matrices.")
# ko_qmat <- do.qcalc(ko.scores, ko.scores_ranked, "K")
# chip_qmat <- do.qcalc(chip.scores, chip.scores_ranked, "C")
# rna_qmat.activator <- do.qcalc(rna.scores, rna.scores_ranked, "R")
# immgen_qmat.activator <- do.qcalc(immgen.scores, immgen.scores_ranked, "I")

# Utilizing ranking methods from AM for testing purposes
source(paste0(getwd(), "/external/rscripts/rscripts/" , "util.R"))

# KO scores
print("Ranking knockout scores.")
ko_qmat.activator <- abs(convert.scores.to.relative.ranks.pos(ko.scores))
ko_qmat.repressor <- abs(convert.scores.to.relative.ranks.pos(-1*ko.scores))
write.mat(ko_qmat.activator, outpath.debug, "K", "_activator_nyu_qmat")
write.mat(ko_qmat.repressor, outpath.debug, "K", "_repressor_nyu_qmat")

# ChIP scores
print("Ranking ChIP scores.")
chip_qmat.activator <- abs(convert.scores.to.relative.ranks(chip.scores.activator))
chip_qmat.repressor <- abs(convert.scores.to.relative.ranks(chip.scores.repressor))
write.mat(chip_qmat.activator, outpath.debug, "C", "_activator_nyu_qmat")
write.mat(chip_qmat.repressor, outpath.debug, "C", "_repressor_nyu_qmat")

# RNA compendium scores
print("Ranking RNA compendium scores.")
rna_qmat.activator <- abs(convert.scores.to.relative.ranks(rna.scores))
rna_qmat.repressor <- abs(convert.scores.to.relative.ranks(-1*rna.scores))
write.mat(rna_qmat.activator, outpath.debug, "R", "_activator_nyu_qmat")
write.mat(rna_qmat.repressor, outpath.debug, "R", "_repressor_nyu_qmat")

# Immgen microarray scores
print("Ranking Immgen microarray scores.")
immgen_qmat.activator <- abs(convert.scores.to.relative.ranks(immgen.scores))
immgen_qmat.repressor <- abs(convert.scores.to.relative.ranks(-1*immgen.scores))
write.mat(immgen_qmat.activator, outpath.debug, "I", "_activator_nyu_qmat")
write.mat(immgen_qmat.repressor, outpath.debug, "I", "_repressor_nyu_qmat")

# --------------
# 4) Combine data according to various data type combinations
# --------------
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "createCombinedMat-fun.R"))

# KC
kc.activator <- createCombinedMat(combo = "kc", type = "activator", ko_qmat = ko_qmat.activator, chip_qmat = chip_qmat.activator, 
                                  rna_qmat = NULL, immgen_qmat = NULL, genes.final, CORE_TFS)
kc.repressor <- createCombinedMat(combo = "kc", type = "repressor", ko_qmat = ko_qmat.repressor, chip_qmat = chip_qmat.repressor, 
                                  rna_qmat = NULL, immgen_qmat = NULL, genes.final, CORE_TFS)

# KCRI
kcri.activator <- createCombinedMat(combo = "kcri", type = "activator", ko_qmat = ko_qmat.activator, chip_qmat = chip_qmat.activator, 
                                    rna_qmat = rna_qmat.activator, immgen_qmat = immgen_qmat.activator, genes.final, CORE_TFS)
kcri.repressor <- createCombinedMat(combo = "kcri", type = "repressor", ko_qmat = ko_qmat.repressor, chip_qmat = chip_qmat.repressor, 
                                    rna_qmat = rna_qmat.repressor, immgen_qmat = immgen_qmat.repressor, genes.final, CORE_TFS)

# --------------
# 6) Write a copy of mmc5 Th17 vs. Th0 (both at 48h) z-scores to a table, which should be loaded in Cytoscape as 'Node table' 
# --------------
print("Writing zscores from mmc5.")
write.mat(zscores.all, outpath.cyt, "zscores")

# --------------
# 7) From combine data matrix, create a list of node-node-value interactions for Cytoscape
# --------------
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