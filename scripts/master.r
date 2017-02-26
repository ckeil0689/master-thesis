#!/usr/bin/env Rscript
# Master script for processing TF-target gene interaction data into network files for Cytoscape 
# (see Computational Methods paper, http://www.cell.com/action/showImagesData?pii=S0092-8674%2812%2901123-3 ) 
# 
# Input: 
#       C = customized MACS output (ChIP-seq) 
#       K = dESeq files for core TFs (RNA-seq KO) 
#       R = Inferelator matrix (RNA-seq compendium)
#       I = Inferelator matrix (Immgen microarray data)
list.of.packages <- c("data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# install if missing
if(length(new.packages)) {
  print("Installing data.table package...")
  install.packages(new.packages, repos="http://cran.rstudio.com/")
}

library(data.table)

DEBUG = TRUE

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

ALLOWED_COMBOS <- c("c", "k", "ri", "kc", "kcr", "kcri")
ALLOWED_CELLS <- c("th17", "th0")

# Get user input 
args <- commandArgs(trailingOnly=TRUE)
combo <- tolower(as.character(args[1]))
thx <- tolower(as.character(args[2]))

# Check user input for valiitiy
print_usage <- function() {
  print(paste("Usage: ./script.r <combo> <cell-type>", "(k = rnaseq-ko, c = chipseq, r = rna-compendium, i = immgen) (Th17 | Th0)"))
}

if(length(args) != 2){
  print_usage()
  stop("Incorrect number of arguments.")
}
   
if(combo == "" || !combo %in% ALLOWED_COMBOS) {
  err <- c("Problem with argument. Enter a valid combination: ", ALLOWED_COMBOS)
  print(paste(err, collapse = " "))
  print_usage()
  stop("invalid argument")
}

if(thx == "" || !thx %in% ALLOWED_CELLS) {
  err <- c("Problem with argument. Enter a valid cell type: ", ALLOWED_CELLS)
  print(paste(err, collapse = " "))
  print_usage()
  stop("invalid argument")
}

# Laod confidence score matrices by option
write.mat <- function(mat, prefix, suffix) {
  filename = paste0(outpath, prefix, thx, suffix, ".txt")
  print(paste("Writing matrix to file:", filename))
  write.table(mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)
}

# Generate zero initialized matrix
tfs <- c()
genes <- c()

# Core target transcription factors
CORE_TFS <- c("batf", "irf4", "stat3", "hif1a", "maf", "fosl2", "rorc")

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

for(opt in opts) {
  # reset because functions may globally change working directory and source() breaks
  setwd(scriptdir)
  
  if(opt == "k") {
    source(paste0(getwd(), "/" , "deseqExtract-fun.R"))
    k_mat <- load.deseq(dir = deseqdir)
    k_sign_mat <- sign(as.data.frame(k_mat))
    if(DEBUG) write.mat(k_mat, "K_", "_smat")
    
  } else if(opt == "c") {
    source(paste0(getwd(), "/" , "chipExtract-fun.R"))
    c_mat <- load.chip(dir = chipdir, reflibfile = ref_filepath, thx = thx)
    if(DEBUG) write.mat(c_mat, "C_", "_smat")
    
  } else if(opt == "r") {
    r_mat <- as.data.frame(read.table(rnaseqfile, sep="\t", header=TRUE))
    r_sign_mat <- sign(as.data.frame(r_mat))
    
  } else if(opt == "i") {
    i_mat <- as.data.frame(read.table(immgenfile, sep="\t", header=TRUE))
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
do.rank <- function(mat, prefix, shouldWrite = FALSE) {
  if(!is.null(mat)) {
    mat_ranked <- rank.smat(mat)
    if(shouldWrite) write.mat(mat_ranked, prefix, "_ranked")
    return(mat_ranked)
  }
  return(NULL)
}

k_mat_ranked <- do.rank(k_mat, "K_")
c_mat_ranked <- do.rank(c_mat, "C_")
r_mat_ranked <- do.rank(r_mat, "R_")
i_mat_ranked <- do.rank(i_mat, "I_")

# --------------
# 3) Calculate quantiles for each ranked matrix to obtain the Q-matrix.
# --------------
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "qmat-fun.R"))

# Wrapper for calculating quantile scores from ranked matrices.
do.qcalc <- function(mat, mat_ranked, prefix, shouldWrite = FALSE) {
  if(!is.null(mat_ranked)) {
    qmat <- calc.qmat(k_mat, k_mat_ranked)
    if(shouldWrite) write.mat(qmat, prefix, "_qmat")
    return(qmat)
  }
  return(NULL)
}

k_qmat <- do.qcalc(k_mat, k_mat_ranked, "K_")
c_qmat <- do.qcalc(c_mat, c_mat_ranked, "C_")
r_qmat <- do.qcalc(r_mat, r_mat_ranked, "R_")
i_qmat <- do.qcalc(i_mat, i_mat_ranked, "I_")

# --------------
# 4) Combine data according to chosen data type combination
# --------------
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "combineQmats-fun.R"))

combined_mat <- combine.qmats(c(k_qmat, c_qmat, r_qmat, i_qmat))
if(DEBUG) write.mat(combined_mat, combo, "")

# --------------
# 5) From combine data matrix, create a list of node-node-value interactions for Cytoscape
# --------------
# Reset because functions may globally change working directory and source() breaks
setwd(scriptdir)
source(paste0(getwd(), "/" , "createInteractions-fun.R"))

edges <- create.interactions(combined_mat)
write.mat(edges, combo, "_edges")

print("Done.")
