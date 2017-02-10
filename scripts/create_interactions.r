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

edges <- NULL
for (i in 1:nrow(cst)) {
  for (j in 1:ncol(cst)) {
    edges <- rbind(edges, c(rownames(cst)[i], rownames(cst)[j], cst[i,j]))
  }
}

colnames(edges) <- c('node1', 'node2', 'value')
filename <- paste0(combo, "_edges.txt")
write.table(edges, file = filename, row.names=FALSE, quote=FALSE, sep='\t')