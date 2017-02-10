#!/usr/bin/env Rscript
# Takes a combined matrix file and transforms it to a list of node-node interactions that can be loaded into Cytoscape
print("Checking for data.table package.")
list.of.packages <- c("data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# install if missing
if(length(new.packages)) {
  print("Installing data.table package...")
  install.packages(new.packages, repos="http://cran.rstudio.com/")
}

library(data.table)

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
print(paste("Reading matrix:", matpath))
cst <- as.data.frame(read.table(matpath, sep="\t", header=TRUE))
rownames(cst) <- cst[, 1]
cst[, 1] <- NULL

print("Transforming matrix to node-node-value list.")
edge_num <- nrow(cst) * ncol(cst)
print(paste("Edges to write:", edge_num))

#pre-allocate data table since dimensions are known
edges <- data.table("node1"=rep(NA,nrow(cst)), "node2"=rep(NA,nrow(cst)), "value"==rep(0,nrow(cst)))
filename <- paste0(combo, "_edges.txt")
write.table(edges, file = filename, row.names=FALSE, quote=FALSE, sep='\t')
stop()

for (i in 1:nrow(cst)) {
  if(i%%100==0) cat("\r", paste0("Progress: ", round((i*100/nrow(cst)), digits = 0), "%"))
  for (j in 1:ncol(cst)) {
    edges <- rbind(edges, c(rownames(cst)[i], colnames(cst)[j], cst[i,j]))
  }
}
print("Done. Writing file...")
colnames(edges) <- c('node1', 'node2', 'value')
filename <- paste0(combo, "_edges.txt")
write.table(edges, file = filename, row.names=FALSE, quote=FALSE, sep='\t')