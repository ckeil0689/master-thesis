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

GLOBAL <- list()
GLOBAL[["DEBUG"]] <- TRUE
GLOBAL[["BOOST_P300"]] <- TRUE

# Input file directories
scriptdir <- getwd()
deseqdir <- paste0(getwd(), "/../suppl/data/deseq/")
chipdir <- paste0(getwd(), "/../suppl/data/chipseq/")
rnaseqfile <- paste0(getwd(), "/../suppl/data/inferelator/GSE40918_Inferelator_RNAseq.txt")
immgenfile <- paste0(getwd(), "/../suppl/data/inferelator/GSE40918_Inferelator_Immgen.txt")
ref_filepath <- paste0(getwd(), "/../suppl/mmc4.csv")

# Put all output into analysis folder
outpath <- paste0(getwd(), "/../suppl/data/analysis/")
if(!dir.exists(outpath)) {
  dir.create(outpath)
}

zz <- file(paste0(outpath,"/all.Rout"), open="wt")
sink(zz, type="message")

# For input checking
ALLOWED_COMBOS <- c("c", "k", "ri", "kc", "kcr", "kcri")

# Core target transcription factors
# CORE_TFS <- c("batf", "irf4", "stat3", "maf", "rorc")
CORE_TFS <- c("batf", "irf4", "stat3", "maf", "rorc", "fosl2")

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
k_mat <- NULL
c_mat <- NULL
r_mat <- NULL
i_mat <- NULL

# ChIP always has positive values
k_sign_mat <- NULL
r_sign_mat <- NULL
i_sign_mat <- NULL

# Split into separate characters
opts = unlist(strsplit(combo, ""))

if(length(opts) == 2) {
  GLOBAL[["abs.cut"]] <- 1.50
} else if(length(opts) == 3) {
  GLOBAL[["abs.cut"]] <- 2.00
} else if(length(opts) == 4) {
  GLOBAL[["abs.cut"]] <- 2.50
} else {
  GLOBAL[["abs.cut"]] <- 0.75
}

print(paste("Determined absolute cutoff from chosen data types:", GLOBAL[["abs.cut"]]))

for(opt in opts) {
  # reset because functions may globally change working directory and source() breaks
  setwd(scriptdir)
  
  if(opt == "k") {
    source(paste0(getwd(), "/" , "deseqExtract-fun.R"))
    k_mat <- load.deseq(dir = deseqdir, CORE_TFS)
    k_sign_mat <- sign(as.data.frame(k_mat))
    if(GLOBAL[["DEBUG"]]) write.mat(k_mat, "K", "_smat")
    
  } else if(opt == "c") {
    source(paste0(getwd(), "/" , "chipExtract-fun.R"))
    c_mat <- load.chip(dir = chipdir, reflibfile = ref_filepath, CORE_TFS)
    if(GLOBAL[["DEBUG"]]) write.mat(c_mat, "C", "_smat")
    
  } else if(opt == "r") {
    r_mat <- as.data.frame(read.table(rnaseqfile, sep="\t", header=TRUE))
    # remove first column or sign()-function will explode
    rownames(r_mat) <- r_mat[, 1]
    r_mat[, 1] <- NULL
    r_sign_mat <- sign(as.data.frame(r_mat))
    
  } else if(opt == "i") {
    i_mat <- as.data.frame(read.table(immgenfile, sep="\t", header=TRUE))
    # remove first column or sign()-function will explode
    rownames(i_mat) <- i_mat[, 1]
    i_mat[, 1] <- NULL
    i_sign_mat <- sign(as.data.frame(i_mat))
    
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
k_mat_ranked <- do.rank(k_mat, "K")
c_mat_ranked <- do.rank(c_mat, "C")
r_mat_ranked <- do.rank(r_mat, "R")
i_mat_ranked <- do.rank(i_mat, "I")

# --------------
# 3) Calculate quantiles for each ranked matrix to obtain the Q-matrix.
# --------------
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "qmat-fun.R"))

# Wrapper for calculating quantile scores from ranked matrices.
# @param mat - The original confidence score S-matrix for the data type
do.qcalc <- function(mat, mat_ranked, prefix) {
  if(!is.null(mat_ranked)) {
    qmat <- calc.qmat(mat, mat_ranked)
    if(GLOBAL[["DEBUG"]]) write.mat(qmat, prefix, "_qmat")
    return(as.matrix(qmat))
  }
  return(NULL)
}

print("Calculating Q-matrices.")
# k_qmat <- do.qcalc(k_mat, k_mat_ranked, "K")
# c_qmat <- do.qcalc(c_mat, c_mat_ranked, "C")
# r_qmat <- do.qcalc(r_mat, r_mat_ranked, "R")
# i_qmat <- do.qcalc(i_mat, i_mat_ranked, "I")

# Utilizing ranking methods from AM for testing purposes
source(paste0(getwd(), "/external/rscripts/rscripts/" , "util.R"))
k_qmat.activator <- abs(convert.scores.to.relative.ranks.pos(k_mat))
k_qmat.repressor <- abs(convert.scores.to.relative.ranks.pos(-1*k_mat))
write.mat(k_qmat.activator, "K", "_nyu_qmat_activator")
write.mat(k_qmat.repressor, "K", "_nyu_qmat_repressor")
c_qmat <- abs(convert.scores.to.relative.ranks(c_mat))
write.mat(c_qmat, "C", "_nyu_qmat")
r_qmat <- NULL
i_qmat <- NULL

# --------------
# 4) Combine data according to chosen data type combination
# --------------
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "combineQmats-fun.R"))

print("Combining Q-matrices to a single matrix.")
combined_mat.activator <- combine.qmats(k_qmat.activator, c_qmat, r_qmat, i_qmat, CORE_TFS)
combined_mat.repressor <- combine.qmats(k_qmat.repressor, c_qmat, r_qmat, i_qmat, CORE_TFS)
if(GLOBAL[["DEBUG"]]) write.mat(combined_mat.activator, combo, "_mat_activator")
if(GLOBAL[["DEBUG"]]) write.mat(combined_mat.repressor, combo, "_mat_repressor")

# --------------
# 5) Apply sign matrix
# --------------
# Empty matrix with same dimension as combined matrix
m.sign.kc <- matrix(0, nc=ncol(combined_mat.activator), nr=nrow(combined_mat.activator), dimnames=dimnames(combined_mat.activator))
# Only set values which also appear in KO matrix (TF-target gene pairs)
m.sign.kc[rownames(k_mat), colnames(k_mat)] <- k_mat
# The knockout values will give us signs, everything else treated as positive (ChIP!)
m.sign.kc <- sign(m.sign.kc)
m.sign.kc[which(m.sign.kc==0)] <- 1

if(GLOBAL[["DEBUG"]]) write.mat(m.sign.kc, combo, "_signmat")

print("Checking dimensions...")
if(!identical(dim(m.sign.kc), dim(combined_mat.activator))) {
  print(paste("Dimension sign_mat:", dim(sign_mat)))
  print(paste("Dimension combined_mat.activator:", dim(combined_mat.activator)))
  print(paste("Dimension k_qmat.activator:", dim(k_qmat.activator)))
  print(paste("Dimension c_qmat:", dim(c_qmat)))
  stop("Sign matrix does not have the same dimension as combined matrix, things will break. Stopping.")
}

print("Applying sign matrix to combined matrix...")
# Element-wise multiplication with sign matrix
combined_mat.activator <- combined_mat.activator * as.vector(m.sign.kc)
# Positive scores in repressor mean repression. Multiply by -1 so repressor edges >1.50 will be filtered as negative in createInteractions
combined_mat.repressor <- combined_mat.repressor * as.vector(m.sign.kc*-1) # element-wise multiplication

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
close(zz)
