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

opts = unlist(strsplit(combo, ""))

# generate zero initialized matrix
tfs <- c()
genes <- c()

# figure out maximal, unique set of genes and TFs from q-matrix files
matfiles <- c()
for(opt in opts) {
  
  if(opt == "k") {
    filename <- paste0(getwd(), "/", "K_Th17_mat_ranked_q.txt") # unix only!
    
  } else if(opt == "c") {
    filename <- paste0(getwd(), "/", "C_Th17_mat_ranked_q.txt") # unix only!
    
  } else {
    stop(paste("Option not recognized:", opt))
  }
  
  matfiles <- c(matfiles, filename)
}

# Vectors for row and column names of final Thx (x=0/=17) matrix
all_genes_thx_wt <- c()
tfs_thx <- c()

print("Finding all genes and TFs to build a matrix.")
# Iteration 1: create unique, maximal list of TFs and genes
for(i in matfiles) {
  # read in the data
  cst <- as.data.frame(read.table(i, sep="\t", header=TRUE))
  rownames(cst) <- cst[, 1]
  cst[, 1] <- NULL
  
  tfs_thx <- c(tfs_thx, colnames(cst))
  all_genes_thx_wt <- append(all_genes_thx_wt, rownames(cst))
}

print("Generate zero-filled matrix skeleton.")
# generate the empty Th17 and Th0 matrices for ChIP-seq
all_genes_thx_unique <- toupper(sort(unique(all_genes_thx_wt)))
tfs_thx_unique <- toupper(sort(unique(tfs_thx)))

# 0-initialized matrix  
combined_mat <- matrix(0, nrow = length(all_genes_thx_unique), ncol = length(tfs_thx_unique))
# unique gene list makes up rows
rownames(combined_mat) <- all_genes_thx_unique
# unique transcription factor list makes up columns
colnames(combined_mat) <- tfs_thx_unique

# Iteration 2: filling data by adding values from each file at the 
# matching index in the combined matrix
for(i in matfiles) {
  cst <- as.data.frame(read.table(i, sep="\t", header=TRUE))
  rownames(cst) <- cst[, 1]
  cst[, 1] <- NULL
  
  cAB <- union(colnames(combined_mat), colnames(cst))
  rAB <- union(rownames(combined_mat), rownames(cst))
  
  A1 <- matrix(0, ncol=length(cAB), nrow=length(rAB), dimnames=list(rAB, cAB))
  B1 <- A1
  
  indxA <- outer(rAB, cAB, FUN=paste) %in% outer(rownames(combined_mat), colnames(combined_mat), FUN=paste) 
  indxB <- outer(rAB, cAB, FUN=paste) %in% outer(rownames(cst), colnames(cst), FUN=paste)
  A1[indxA] <- combined_mat
  B1[indxB] <- as.matrix(cst)
  
  combined_mat <- A1 + B1
}

print("Writing combined matrix to file.")
# write matrices to a tab-delimited file
filename=paste0(combo, ".txt")
write.table(combined_mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)