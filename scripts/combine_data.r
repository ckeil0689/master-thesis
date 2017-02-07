#!/usr/bin/env Rscript
# Combine previously generated Q-matrices
setwd(paste0(getwd(), "/../suppl/data/analysis/"))

args <- commandArgs(trailingOnly=TRUE)
combo <- casefold(as.character(args[1]), upper = FALSE)

ALLOWED_COMBOS <- c("ri", "kc", "kcr", "kcri")

if(length(args) != 1 || combo == "" || !combo %in% ALLOWED_COMBOS) {
  err <- c("Problem with argument. Enter a valid combination: ", ALLOWED_COMBOS)
  print(paste(err, collapse = " "))
  print(paste("Usage: ./combine_data.r <combo>", "(k = rnaseq-ko, c = chipseq, r = rna-compendium, i = immgen)"))
  stop("invalid argument")
}

opts = strsplit(combo, "")

# generate zero initialized matrix
tfs <- c()
genes <- c()

# figure out maximal, unique set of genes and TFs from q-matrix files
for(opt in opts) {
  if(opt == "k") {
    
  }
}