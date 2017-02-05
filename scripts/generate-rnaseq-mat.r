#!/usr/bin/env Rscript
# TODO: could be combined with generate-chipseq-mat.r to reduce code
# Read a DeSeq files from dir
setwd("/home/chrisk/Documents/uni/thesis/suppl/data/deseq")

args <- commandArgs(trailingOnly=TRUE)
thx <- toupper(as.character(args[1]))
thx <-'Th17'

# fixed target transcription factors to use
target_tfs <- c("batf", "irf4", "stat3", "hif1a", "cmaf", "maf", "fosl2", "rorc", "rorg")

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
# Iteration 1: create unique, maximal list of tested TFs and genes with MACS peaks
for(i in deseqfiles) {
  
  # skip non-DEseq files, regex tests for format as found on GEO Series GSE40918
  if(!grepl('^GSE[0-9]+_(Th[0-9]{1,2})\\..+\\.wt\\.vs\\.(Th[0-9]{1,2})\\..+\\.ko_.+_(20\\d\\d)\\.txt$', i)) {
    next
  }
  
  # read in the data and extract the library name
  cst <- read.table(i, sep="\t", header=TRUE)
  
  # TODO get TF name and add it to list
  tfcol <- colnames(cst)[3]
  tf <- strsplit(tfcol, "[.]")[[1]][2] # assumes convention column name (e.g. Th17.batf.wt -> batf)
  
  # only use target TFs
  if(!tf %in% target_tfs) {
    print(paste("Transcription factor not in target group:", tf, "(skipped)"))
    next
  }
  
  # create full gene list by adding all genes from this file
  genes_thx_wt <- as.character(cst$id)
  all_genes_thx_wt <- append(all_genes_thx_wt, genes_thx_wt)
}

print("Generate dummy matrix skeleton.")
# generate the empty Th17 and Th0 matrices for ChIP-seq
all_genes_thx_unique <- sort(unique(all_genes_thx_wt))

# 0-initialized matrix  
thx_mat <- matrix(0, nrow = length(all_genes_thx_unique), ncol = length(tfs_thx))
# unique gene list makes up rows
rownames(thx_mat) <- all_genes_thx_unique
# unique transcription factor list makes up columns
colnames(thx_mat) <- tfs_thx

print("Extract DESeq p-values and log2 foldchange from files.")
#TODO implement