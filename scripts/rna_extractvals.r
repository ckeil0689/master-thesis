#!/usr/bin/env Rscript
# TODO: could be combined with generate-chipseq-mat.r to reduce code
# Read a DeSeq files from dir
setwd(paste0(getwd(), "/../suppl/data/deseq"))

args <- commandArgs(trailingOnly=TRUE)
thx <- toupper(as.character(args[1]))
thx <-'Th17'

# fixed target transcription factors to use
target_tfs <- c("batf", "irf4", "stat3", "hif1a", "ikzf3", "cmaf", "maf", "fosl2", "rorc", "rorg")

# ensure we have a usable argument to work with
if(thx != 'Th17' && thx != 'Th0') {
  print("Incorrect or missing argument.")
  print("Usage: script.R <(Th17|Th0)>")
  stop()
}

# DEseq results files
deseqfiles <- list.files(getwd())

# Vectors for row and column names of final Thx (x=0/=17) matrix
all_genes_thx_wt <- c()
tfs_thx <- c()

print("Reading DESeq files to create complete list of genes.")
# Iteration 1: get list of all tested genes and TFs so a matrix can be set up
for(i in deseqfiles) {
  
  # skip non-DEseq files, regex tests for format as found on GEO Series GSE40918
  if(!grepl('^GSE[0-9]+_(Th[0-9]{1,2})\\..+\\.wt\\.vs\\.(Th[0-9]{1,2})\\..+\\.ko_.+_(20\\d\\d)\\.txt$', i)) {
    next
  }
  
  # read in the data and extract the library name
  cst <- read.table(i, sep="\t", header=TRUE)
  
  # get TF name and add it to list
  tfcol <- colnames(cst)[3]
  tf <- strsplit(tfcol, "[.]")[[1]][2] # assumes convention column name (e.g. Th17.batf.wt -> batf)
  
  # only use target TFs
  if(!tf %in% target_tfs) {
    print(paste("Transcription factor not in target group:", tf, "(skipped)"))
    next
  }
  
  tfs_thx <- c(tfs_thx, tf)
  
  # create full gene list by adding all genes from this file
  genes_thx_wt <- as.character(cst$id)
  all_genes_thx_wt <- append(all_genes_thx_wt, genes_thx_wt)
}

print("Generate zero-filled matrix skeleton.")
# generate the empty Th17 and Th0 matrices for ChIP-seq
all_genes_thx_unique <- sort(unique(all_genes_thx_wt))

# 0-initialized matrix  
thx_mat <- matrix(0, nrow = length(all_genes_thx_unique), ncol = length(tfs_thx))
# unique gene list makes up rows
rownames(thx_mat) <- all_genes_thx_unique
# unique transcription factor list makes up columns
colnames(thx_mat) <- tfs_thx

print("Extract DESeq p-values (non-adjusted) and log2 foldchange from files.")
# iteration 2: extract p-values and log2 values from DESeq results files and 
# fill confidence score matrix according to formula in Computational Methods
for(i in deseqfiles) {
  
  # read in the data and extract the library name
  cst <- read.table(i, sep="\t", header=TRUE)
  
  # get TF name and add it to list
  tfcol <- colnames(cst)[3]
  tf <- strsplit(tfcol, "[.]")[[1]][2] # assumes convention column name (e.g. Th17.batf.wt -> batf)
  
  # get the DESeq p-values by iterating and accessing matrix via id and TF-name
  idx <- 1
  for(j in cst$id) {
    thx_mat[j, tf] <- -log10(cst$pval[idx]) * sign(cst$log2FoldChange[idx])
    idx <- idx + 1
  }
}

# replace NA and Inf values in matrix with 0s (for later ranking procedure)
thx_mat[thx_mat == Inf] <- 0
thx_mat[is.na(thx_mat)] <- 0

print("Writing K-matrix to file.")
# write matrices to a tab-delimited file
filename=paste0("K_", thx, "_mat.txt")
write.table(thx_mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)