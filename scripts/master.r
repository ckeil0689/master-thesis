#!/usr/bin/env Rscript
# Master script for processing TF-target gene interaction data into network files for Cytoscape 
# (see Computational Methods paper, http://www.cell.com/action/showImagesData?pii=S0092-8674%2812%2901123-3 ) 
# 
# Input: 
#       C = customized MACS output (ChIP-seq) 
#       K = dESeq files for core TFs (RNA-seq KO) 
#       R = Inferelator matrix (RNA-seq compendium)
#       I = Inferelator matrix (Immgen microarray data)

# Input file directories
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

# Generate zero initialized matrix
tfs <- c()
genes <- c()

# Core target transcription factors
CORE_TFS <- c("batf", "irf4", "stat3", "hif1a", "maf", "fosl2", "rorc")

# Laod confidence score matrices by option
k_mat <- NULL
c_mat <- NULL
r_mat <- NULL
i_mat <- NULL

# Split into separate characters
opts = unlist(strsplit(combo, ""))

for(opt in opts) {
  if(opt == "k") {
    k_mat <- loaddeseq()
    
  } else if(opt == "c") {
    source("chipExtract-fun.R")
    c_mat <- load.chip(dir = chipdir, reflibfile = ref_filepath, thx = thx)
    # write matrices to a tab-delimited file
    filename=paste0(outpath, "C_", thx, "_mat.txt")
    write.table(c_mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)
    print(paste("Wrote:", filename))
    
  } else if(opt == "r") {
    r_mat <- as.data.frame(read.table(rnaseqfile, sep="\t", header=TRUE))
    
  } else if(opt == "i") {
    i_mat <- as.data.frame(read.table(immgenfile, sep="\t", header=TRUE))
    
  } else {
    stop(paste("Option not recognized:", opt))
  }
}

loaddeseq <- function() {
  setwd(deseqdir)
}
