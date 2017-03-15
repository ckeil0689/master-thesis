#!/usr/bin/env Rscript
# Master script for processing TF-target gene interaction data into network files for Cytoscape 
# (see Computational Methods paper, http://www.cell.com/action/showImagesData?pii=S0092-8674%2812%2901123-3 ) 
# Aviv Madar has provided R-scripts from the original project. Some procedures have been taken and/or adapted 
# from his code.
#
# Input: 
#       C = customized MACS output (ChIP-seq) 
#       K = dESeq files for core TFs (RNA-seq KO) 
#       R = Inferelator matrix (RNA-seq compendium)
#       I = Inferelator matrix (Immgen microarray data)
list.of.packages <- c("data.table", "reshape2")#, "xlsx")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# install if missing
if(length(new.packages)) {
  print("Installing data.table package...")
  install.packages(new.packages, repos="http://cran.rstudio.com/")
}

library(data.table)
library(reshape2)
# library(xlsx)

# Global variables
GLOBAL <- list()
GLOBAL[["DEBUG"]] <- TRUE
GLOBAL[["z.abs.cut"]] <- 2.50 # carried over from original Aviv Madar code

# Input file directories
scriptdir <- getwd()
deseqdir <- paste0(getwd(), "/../suppl/data/deseq/")
chipdir <- paste0(getwd(), "/../suppl/data/chipseq/")
rnaseqfile <- paste0(getwd(), "/../suppl/data/inferelator/GSE40918_Inferelator_RNAseq.txt")
immgenfile <- paste0(getwd(), "/../suppl/data/inferelator/GSE40918_Inferelator_Immgen.txt")
ref_filepath <- paste0(getwd(), "/../suppl/mmc4.csv")
zscores_filepath <- paste0(getwd(), "/../suppl/mmc5.xls")

# Put all output into analysis folder
outpath <- paste0(getwd(), "/../suppl/data/analysis/")
if(!dir.exists(outpath)) {
  dir.create(outpath)
}

# For input checking
ALLOWED_COMBOS <- c("c", "k", "ri", "kc", "kcr", "kcri")

# Core target transcription factors
# CORE_TFS <- c("batf", "irf4", "stat3", "maf", "rorc")
CORE_TFS <- c("batf", "irf4", "stat3", "maf", "rorc", "fosl2", "hif1a")

# Get user input 
args <- commandArgs(trailingOnly=TRUE)
combo <- tolower(as.character(args[1]))

# Check user input for valiitiy
print_usage <- function() {
  print(paste("Usage: ./script.r <combo>", "(k = rnaseq-ko, c = chipseq, r = rna-compendium, i = immgen)"))
}

if(length(args) != 1){
  print_usage()
  stop("Incorrect number of arguments.")
}
   
if(combo == "" || !combo %in% ALLOWED_COMBOS) {
  err <- c("Problem with argument. Enter a valid combination: ", ALLOWED_COMBOS)
  print(paste(err, collapse = " "))
  print_usage()
  stop("invalid argument")
}

# Load confidence score matrices by option
write.mat <- function(mat, prefix, suffix) {
  filename = paste0(outpath, prefix, suffix, ".txt")
  print(paste("Writing matrix to file:", filename))
  write.table(mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)
}

# --------------
# 1) Load data from each selected data type to create confidence score matrix S
# --------------
zscore.table <- read.table(zscores_filepath, sep="\t", header=TRUE)
# Column is '8' because R has trouble reading "Th17/Th0 zscores" and doesn't accept the String it produces when printing colnames()...
zscore_col <- "Th17.Th0.zscores"
genes.final.idx <- which(abs(zscore.table[, zscore_col]) > GLOBAL[["z.abs.cut"]])
genes.final <- zscore.table[genes.final.idx, "Gene_id"]

print(paste("Original gene number:", length(zscore.table[, "Gene_id"])))
print(paste("Genes with abs(zscore) > 2.50:", length(genes.final)))

ko.scores <- NULL
chip.scores.activator <- NULL
chip.scores.repressor <- NULL
rna.scores <- NULL
immgen.scores <- NULL

# ChIP always has positive values
k_sign_mat <- NULL
r_sign_mat <- NULL
i_sign_mat <- NULL

# Split into separate characters
opts = unlist(strsplit(combo, ""))

for(opt in opts) {
  # reset because functions may globally change working directory and source() breaks
  setwd(scriptdir)
  
  if(opt == "k") {
    source(paste0(getwd(), "/" , "deseqExtract-fun.R"))
    ko.scores <- load.deseq(dir = deseqdir, CORE_TFS)
    k_sign_mat <- sign(as.data.frame(ko.scores))
    if(GLOBAL[["DEBUG"]]) write.mat(ko.scores, "K", "_smat")
    
  } else if(opt == "c") {
    source(paste0(getwd(), "/" , "chipExtract-fun.R"))
    chip.scores.activator <- load.chip(dir = chipdir, reflibfile = ref_filepath, boost.p300 = TRUE, CORE_TFS)
    chip.scores.repressor <- load.chip(dir = chipdir, reflibfile = ref_filepath, boost.p300 = FALSE, CORE_TFS)
    if(GLOBAL[["DEBUG"]]) write.mat(chip.scores.activator, "C", "_activator_smat")
    if(GLOBAL[["DEBUG"]]) write.mat(chip.scores.repressor, "C", "_repressor_smat")
    
  } else if(opt == "r") {
    rna.scores <- as.data.frame(read.table(rnaseqfile, sep="\t", header=TRUE))
    # remove first column or sign()-function will explode
    rownames(rna.scores) <- rna.scores[, 1]
    rna.scores[, 1] <- NULL
    r_sign_mat <- sign(as.data.frame(rna.scores))
    
  } else if(opt == "i") {
    immgen.scores <- as.data.frame(read.table(immgenfile, sep="\t", header=TRUE))
    # remove first column or sign()-function will explode
    rownames(immgen.scores) <- immgen.scores[, 1]
    immgen.scores[, 1] <- NULL
    i_sign_mat <- sign(as.data.frame(immgen.scores))
    
  } else {
    stop(paste("Option not recognized:", opt))
  }
}

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
    if(GLOBAL[["DEBUG"]]) write.mat(mat_ranked, prefix, "_ranked")
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
    if(GLOBAL[["DEBUG"]]) write.mat(qmat, prefix, "_qmat")
    return(as.matrix(qmat))
  }
  return(NULL)
}

print("Calculating Q-matrices.")
# ko_qmat <- do.qcalc(ko.scores, ko.scores_ranked, "K")
# chip_qmat <- do.qcalc(chip.scores, chip.scores_ranked, "C")
# rna_qmat <- do.qcalc(rna.scores, rna.scores_ranked, "R")
# immgen_qmat <- do.qcalc(immgen.scores, immgen.scores_ranked, "I")

# Utilizing ranking methods from AM for testing purposes
source(paste0(getwd(), "/external/rscripts/rscripts/" , "util.R"))

# KO scores
ko_qmat.activator <- abs(convert.scores.to.relative.ranks.pos(ko.scores))
ko_qmat.repressor <- abs(convert.scores.to.relative.ranks.pos(-1*ko.scores))
write.mat(ko_qmat.activator, "K", "_activator_nyu_qmat")
write.mat(ko_qmat.repressor, "K", "_repressor_nyu_qmat")

# ChIP scores
chip_qmat.activator <- abs(convert.scores.to.relative.ranks(chip.scores.activator))
chip_qmat.repressor <- abs(convert.scores.to.relative.ranks(chip.scores.repressor))
write.mat(chip_qmat.activator, "C", "_activator_nyu_qmat")
write.mat(chip_qmat.repressor, "C", "_repressor_nyu_qmat")

rna_qmat <- NULL
immgen_qmat <- NULL

# --------------
# 4) Combine data according to chosen data type combination
# --------------
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "combineQmats-fun.R"))

print("Combining Q-matrices to a single matrix.")
combined_mat.activator <- combine.qmats(ko_qmat.activator, chip_qmat.activator, rna_qmat, immgen_qmat, CORE_TFS)
combined_mat.repressor <- combine.qmats(ko_qmat.repressor, chip_qmat.repressor, rna_qmat, immgen_qmat, CORE_TFS)

# Filter the combined matrices by zscores from mmc5
filtered.genes.idx <- which(rownames(combined_mat.activator) %in% genes.final)
filtered.genes.activator <- rownames(combined_mat.activator)[filtered.genes.idx]
combined_mat.activator <- combined_mat.activator[filtered.genes.idx,]
rownames(combined_mat.activator) <- filtered.genes.activator

filtered.genes.idx <- which(rownames(combined_mat.repressor) %in% genes.final)
filtered.genes.repressor <- rownames(combined_mat.repressor)[filtered.genes.idx]
combined_mat.repressor <- combined_mat.repressor[filtered.genes.idx,]
rownames(combined_mat.repressor) <- filtered.genes.repressor

if(GLOBAL[["DEBUG"]]) write.mat(combined_mat.activator, combo, "_mat_activator")
if(GLOBAL[["DEBUG"]]) write.mat(combined_mat.repressor, combo, "_mat_repressor")

# --------------
# 5) Apply sign matrix
# --------------
print("Applying signs to activator matrix.")
#  Activator
# Empty matrix with same dimension as combined matrix
m.sign.kc <- matrix(0, nc=ncol(combined_mat.activator), nr=nrow(combined_mat.activator), dimnames=dimnames(combined_mat.activator))
# Only set values which also appear in KO matrix (TF-target gene pairs)
ko.genes <- rownames(ko.scores)[which(filtered.genes.activator %in% rownames(ko.scores))]
m.sign.kc[rownames(ko.genes), colnames(ko.scores)] <- ko.scores[ko.genes,]
# The knockout values will give us signs, everything else treated as positive (ChIP!)
m.sign.kc <- sign(m.sign.kc)
m.sign.kc[which(m.sign.kc==0)] <- 1
if(GLOBAL[["DEBUG"]]) write.mat(m.sign.kc, combo, "_activator_signmat")

print("Applying signs to repressor matrix.")
# Repressor
# Empty matrix with same dimension as combined matrix
m.sign.kc.r <- matrix(0, nc=ncol(combined_mat.repressor), nr=nrow(combined_mat.repressor), dimnames=dimnames(combined_mat.repressor))
# Only set values which also appear in KO matrix (TF-target gene pairs)
ko.genes <- rownames(ko.scores)[which(filtered.genes.repressor %in% rownames(ko.scores))]
m.sign.kc.r[rownames(ko.genes), colnames(ko.scores)] <- ko.scores[ko.genes,]
# The knockout values will give us signs, everything else treated as positive (ChIP!)
m.sign.kc.r <- sign(m.sign.kc.r)
m.sign.kc.r[which(m.sign.kc.r==0)] <- 1
if(GLOBAL[["DEBUG"]]) write.mat(m.sign.kc, combo, "_repressor_signmat")

print("Checking dimensions...")
if(!identical(dim(m.sign.kc), dim(combined_mat.activator))) {
  print(paste("Dimension sign_mat:", dim(m.sign.kc)))
  print(paste("Dimension combined_mat.activator:", dim(combined_mat.activator)))
  print(paste("Dimension ko_qmat.activator:", dim(ko_qmat.activator)))
  print(paste("Dimension chip_qmat:", dim(chip_qmat.activator)))
  stop("Sign matrix does not have the same dimension as combined matrix, things will break. Stopping.")
}

print("Applying sign matrix to combined activator matrix...")
# Element-wise multiplication with sign matrix
combined_mat.activator <- combined_mat.activator * as.vector(m.sign.kc)

print("Checking dimensions...")
if(!identical(dim(m.sign.kc.r), dim(combined_mat.repressor))) {
  print(paste("Dimension sign_mat:", dim(m.sign.kc.r)))
  print(paste("Dimension combined_mat.activator:", dim(combined_mat.repressor)))
  print(paste("Dimension ko_qmat.activator:", dim(ko_qmat.repressor)))
  print(paste("Dimension chip_qmat:", dim(chip_qmat.repressor)))
  stop("Sign matrix does not have the same dimension as combined matrix, things will break. Stopping.")
}

print("Applying sign matrix to combined repressor matrix...")
# Positive scores in repressor mean repression. Multiply by -1 so repressor edges >1.50 will be filtered as negative in createInteractions
combined_mat.repressor <- combined_mat.repressor * as.vector(m.sign.kc.r * -1) # element-wise multiplication

if(GLOBAL[["DEBUG"]])  {
  write.mat(combined_mat.activator, combo, "_signed_activator")
  write.mat(combined_mat.repressor, combo, "_signed_repressor")
}

# --------------
# 6) From combine data matrix, create a list of node-node-value interactions for Cytoscape
# --------------
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "createInteractions-fun.R"))

# Positive edges are surpassing the cut value for the edges. They are tagged as 'positive_KC' if >1.50 is about activation, and as 'negative_KC' if >1.50 is
# about repression. This has been imitated from Aviv Madar's original code in an effort to get the procedure right, but I find the variable handling and
# naming very inconvenient and confusing. If there is time, I will change this.
print("Writing interactions as single list...")
create.interactions(combined_mat.activator, outpath, combo, "single", pos.edge= "positive_KC", neg.edge = "negative_KC", append = FALSE)
create.interactions(combined_mat.repressor, outpath, combo, "single", pos.edge= "negative_KC", neg.edge = "positive_KC", append = TRUE)

print("Writing interactions as separate lists...")
create.interactions(combined_mat.activator, outpath, combo, "activator", pos.edge= "positive_KC", neg.edge = "negative_KC", append = FALSE)
create.interactions(combined_mat.repressor, outpath, combo, "repressor", pos.edge= "negative_KC", neg.edge = "positive_KC", append = FALSE)
print("Done.")