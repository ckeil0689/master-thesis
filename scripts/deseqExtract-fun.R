# ----------------
# Constants
# ----------------
# DEseq file directory relative to /scripts/
deseqdir <- paste0(getwd(), "/../suppl/data/deseq/")
# Ensure we are in correct directory
if(!dir.exists(deseqdir)) stop("Cannot load ChIP-files because the directory does not exist.")
setwd(deseqdir)
deseqfiles <- list.files(getwd())

extract.tf <- function(tfcol) {
  # get TF name and add it to list
  tf <- toupper(strsplit(tfcol, "[.]")[[1]][2]) # assumes convention column name (e.g. Th17.batf.wt -> batf)
  return(tf)
}

get.skel.mat <- function() {
  print(paste("Reading", length(deseqfiles),"DESeq files to create complete list of genes."))
  
  # Vectors for row and column names of final KO-matrix
  all.genes <- c()
  all.tfs <- c()
  
  # Iteration 1: get list of all tested genes and TFs so a matrix can be pre-allocated
  for(i in deseqfiles) {
    # read in the data and extract the library name
    cst <- read.table(i, sep="\t", header=TRUE)
    tf <- extract.tf(colnames(cst)[3])
    
    # only use target TFs
    if(!tolower(tf) %in% CORE_TFS) {
      print(paste("Transcription factor not in target group:", tf, "(skipped)"))
      next
    }
    
    all.tfs <- c(all.tfs, tf)
    
    # create full gene list by adding all genes from this file
    file.genes <- as.character(cst$id)
    all.genes <- append(all.genes, file.genes)
  }
  
  print("Generate zero-filled matrix skeleton.")
  # generate the empty Th17 and Th0 matrices for ChIP-seq
  all.genes.unique <- toupper(sort(unique(all.genes)))
  
  # 0-initialized matrix  
  scores.empty <- matrix(0, nrow = length(all.genes.unique), ncol = length(all.tfs))
  # unique gene list makes up rows
  rownames(scores.empty) <- all.genes.unique
  # unique transcription factor list makes up columns
  colnames(scores.empty) <- toupper(all.tfs)
  
  return(scores.empty)
}

# Extract p-values and log2(foldchange) values from DESeq results files and 
# fill pre-allocated confidence score matrix according to formula in Computational Methods:
# score = pval * sign(log2(foldchange))
populate.ko.scores <- function(ko.scores) {
  print("Extract DESeq non-adjusted p-values and log2(foldchange) from files.")
  
  for(i in deseqfiles) {
    # read in the data and extract the library name
    cst <- read.table(i, sep="\t", header=TRUE)
    tf <- extract.tf(colnames(cst)[3])
    
    # only use target TFs
    if(!tolower(tf) %in% CORE_TFS) {
      print(paste("Transcription factor not in target group:", tf, "(skipped)"))
      next
    }
    
    # get the DESeq p-values by iterating and accessing matrix via id and TF-name (genes are not ordered by name in DESeq files!)
    idx <- 1
    for(j in cst$id) {
      ko.scores[j, tf] <- -log10(cst$pval[idx]) * sign(cst$log2FoldChange[idx])
      idx <- idx + 1
    }
  }
  
  # replace NA and Inf values in matrix with 0s (for later ranking procedure)
  print("Replacing Inf and NA values with 0 KO-score.")
  ko.scores[ko.scores == Inf] <- 0
  ko.scores[is.na(ko.scores)] <- 0
  
  # drop 0-only-rows
  # print(dim(ko.scores))
  # ko.scores <- ko.scores[rowSums(abs(ko.scores[, -1]))>(1e-10),]
  # print("DROPPED---------------------------------------")
  # print(dim(ko.scores))
  
  return(ko.scores)
}

# Read a DEseq files from dir
load.deseq <- function(CORE_TFS) {
  empty.scores <- get.skel.mat()
  ko.scores <- populate.ko.scores(empty.scores)
  return(ko.scores)
}