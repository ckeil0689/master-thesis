#!/usr/bin/env Rscript
# Master script for processing TF-target gene interaction data into network files for Cytoscape 
# (see Computational Methods paper, http://www.cell.com/action/showImagesData?pii=S0092-8674%2812%2901123-3 ) 
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

print(paste("Chosen absolute cutoff:", GLOBAL[["abs.cut"]]))

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
k_qmat <- do.qcalc(k_mat, k_mat_ranked, "K")
c_qmat <- do.qcalc(c_mat, c_mat_ranked, "C")
r_qmat <- do.qcalc(r_mat, r_mat_ranked, "R")
i_qmat <- do.qcalc(i_mat, i_mat_ranked, "I")

# source(paste0(getwd(), "/external/rscripts/rscripts/" , "util.R"))
# k_qmat <- convert.scores.to.relative.ranks.pos(k_mat)
# write.mat(k_qmat, "K", "_nyu_qmat")
# c_qmat <- convert.scores.to.relative.ranks.pos(c_mat)
# write.mat(c_qmat, "C", "_nyu_qmat")
# r_qmat <- NULL
# i_qmat <- NULL
# # write.mat(c_nyu_qmat, "C", "_nyu_qmat")

# --------------
# 4) Combine data according to chosen data type combination
# --------------
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "combineQmats-fun.R"))

print("Combining Q-matrices to a single matrix.")
combined_mat <- combine.qmats(k_qmat, c_qmat, r_qmat, i_qmat, CORE_TFS)
if(GLOBAL[["DEBUG"]]) write.mat(combined_mat, combo, "_mat")

# --------------
# 5) Apply sign matrix 
# --------------
sign_mat <- k_sign_mat #temp
sign_mat[sign_mat == 0] <- 1
if(GLOBAL[["DEBUG"]]) write.mat(sign_mat, combo, "_signmat")

print("Checking dimensions...")

# if(!identical(dim(sign_mat), dim(combined_mat))) {
#   stop("Sign matrix does not have the same dimension as combined matrix, things will break. Stopping.")
# }

print("Applying sign matrix to combined matrix...")
combined_mat <- combined_mat[rownames(k_qmat),] * as.vector(sign_mat) # element-wise multiplication

if(GLOBAL[["DEBUG"]]) write.mat(combined_mat, combo, "_signed")

# --------------
# 6) From combine data matrix, create a list of node-node-value interactions for Cytoscape
# --------------
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "createInteractions-fun.R"))

print("Creating interactions...")
create.interactions(combined_mat, outpath, combo)
print("Done.")
close(zz)
