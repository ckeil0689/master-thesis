#!/usr/bin/env Rscript
# Takes a combined matrix file and transforms it to a list of node-node interactions that can be loaded into Cytoscape
setwd(paste0(getwd(), "/../suppl/data/analysis/"))

args <- commandArgs(trailingOnly=TRUE)

combo <- tolower(as.character(args[1]))

ALLOWED_COMBOS <- c("ri", "kc", "kcr", "kcri")

if(length(args) != 1 || combo == "" || !combo %in% ALLOWED_COMBOS) {
  err <- c("Problem with argument. Enter a valid combination: ", ALLOWED_COMBOS)
  print(paste(err, collapse = " "))
  print(paste("Usage: ./combine_data.r <combo>", "(k = rnaseq-ko, c = chipseq, r = rna-compendium, i = immgen)"))
  stop("invalid argument")
}

matpath <- paste0(getwd(), "/", combo, ".txt")

if(!file.exists(matpath)) {
  stop(paste("File does not exist:", matpath))
}

# read in the combined matrix file
cst <- as.data.frame(read.table(matpath, sep="\t", header=TRUE))
rownames(cst) <- cst[, 1]
cst[, 1] <- NULL

