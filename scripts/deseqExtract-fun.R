# ----------------
# Constants
# ----------------
# DEseq file directory relative to /scripts/
deseqdir <- paste0(getwd(), "/../suppl/data/deseq/")
# Ensure we are in correct directory
if(!dir.exists(deseqdir)) stop("Cannot load ChIP-files because the directory does not exist.")
setwd(deseqdir)

# We consider all DESeq result files (as opposed to selective ChIP-seq loading)
deseqfiles <- list.files(getwd())

# ----------------
# Functions
# ----------------
# Extracts the TF target from a DESeq file according to their format
# Assumes convention column name (e.g. Th17.batf.wt -> batf)
extract.tf <- function(tfcol) {
  splitcol <- strsplit(tfcol, "[.]")
  tf <- splitcol[[1]][2]
  # only use target TFs
  if(!tf %in% GLOBAL[["CORE_TFS"]]) {
    warning(paste("Could not load DESeq file for:", tf, "(skipped)"))
    return(NA)
  }
  return(toupper(tf))
}

# Scans DESeq files for genes and transcription factors (matching CORE_TFS) to generate a pre-allocated skeleton matrix
get.skel.mat <- function() {
  # Ensure we are in correct directory
  if(!dir.exists(deseqdir)) stop("Cannot load ChIP-files because the directory does not exist.")
  setwd(deseqdir)
  print(paste("Reading", length(deseqfiles),"DESeq files to create complete list of genes."))
  
  # Vectors for row and column names of final KO-matrix
  all.genes <- c()
  all.tfs <- c()
  
  for(i in deseqfiles) {
    # read in the data and extract the library name
    cst <- read.table(i, sep="\t", header=TRUE)
    tf <- extract.tf(colnames(cst)[3])
    if(is.na(tf)) next
    
    all.tfs <- c(all.tfs, tf)
    
    # Collect all genes from this file
    file.genes <- as.character(cst$id)
    all.genes <- append(all.genes, file.genes)
  }
  
  print("Generating zero-filled KO-matrix skeleton.")
  all.genes.unique <- toupper(sort(unique(all.genes)))
  all.tfs.unique <- toupper(sort(unique(all.tfs)))
  
  # 0-initialized matrix  
  scores.empty <- matrix(0, nrow = length(all.genes.unique), ncol = length(all.tfs.unique))
  rownames(scores.empty) <- all.genes.unique
  colnames(scores.empty) <- all.tfs.unique
  
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
    if(is.na(tf)) next
    
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

# ----------------
# Main function: load & process DEseq data and return the confidence score matrix S(KO)
# ----------------
load.deseq <- function() {
  empty.scores <- get.skel.mat()
  ko.scores <- populate.ko.scores(empty.scores)
  return(ko.scores)
}