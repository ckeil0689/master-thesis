# DEseq file directory relative to /scripts/
deseqdir <- paste0(getwd(), "/../suppl/data/deseq/")
# Ensure we are in correct directory
if(!dir.exists(deseqdir)) stop("Cannot load DESeq-files because the directory does not exist.")
setwd(deseqdir)

# We consider all DESeq result files (as opposed to selective ChIP-seq loading)
deseqfiles <- list.files(getwd())
if(length(deseqfiles) == 0) stop(paste("No DESeq files found in: ", deseqdir, "Stopping."))

# Extracts the TF target from a DESeq file according to their format
# Assumes convention column name (e.g. Th17.batf.wt -> batf)
extract.tf <- function(deseq.cols, i) {
  # regex looks for a sring of the pattern 'pre.string.suf' where string will be the extracted tf 
  # this will extract the name from the specific DESeq files in GSE40918 (e.g. th17.batf.wt) and allow anyone
  # to add custom names
  tf.col <- grep(c("(.*\\.(?!.*rpkm).+\\..*)($|\n)"), deseq.cols, ignore.case = TRUE, perl = TRUE, value = TRUE)
  
  if(length(tf.col) > 0) {
    splitcol <- strsplit(tf.col[[1]][1], "[.]")
    tf <- splitcol[[1]][2]
  } else {
    # no name column with correct pattern 'pre.string.suf'
    tf <- NA
  }
  
  if(is.na(tf)) {
    warning(paste("Could not load DESeq file:", i, "(skipped).
                  No TF name could be found, cannot assign column in interaction matrix. 
                  Add a column (header is important) to the DESeq file with the pattern 
                  pre.TFNAME.suf (e.g. th17.batf.wt)."))
    return(NA)
  }
  
  # only use target TFs
  if(!tf %in% GLOBAL[["CORE_TFS"]]) {
    warning(paste("Could not load DESeq file for:", tf, "(skipped).
                  The transcription factor is not in the allowed list of core transcription factors.
                  Edit the global variable CORE_TFS to allow the transcription factor to be loaded."))
    return(NA)
  }
  
  return(toupper(tf))
}

# Scans DESeq files for genes and transcription factors (matching CORE_TFS) to generate a pre-allocated skeleton matrix
get.skel.mat <- function() {
  # Ensure we are in correct directory
  if(!dir.exists(deseqdir)) stop("Cannot load ChIP-files because the directory does not exist.")
  setwd(deseqdir)
  println(paste("Reading", length(deseqfiles),"DESeq files to create complete list of genes."))
  
  # Vectors for row and column names of final DESeq-matrix
  all.genes <- c()
  all.tfs <- c()
  
  for(i in deseqfiles) {
    # read in the data and extract the library name
    cst <- read.table(i, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    tf <- extract.tf(colnames(cst), i)
    if(is.na(tf)) next
    
    all.tfs <- c(all.tfs, tf)
    
    # Collect all genes from current DESeq file; regex wants to capture anything like: id, *gene.id, *gene_id
    # Tested on regex101.com =; perl = TRUE option lets grep use pcre regex library
    gene.id.cols <- grep(c("\\b(.*gene[_|\\.])?id($|\n)"), colnames(cst), ignore.case = TRUE, perl = TRUE, value = TRUE)
    gene.ids <- cst[, gene.id.cols[[1]][1]] # first element if multiple matches
    file.genes <- as.character(gene.ids)
    all.genes <- append(all.genes, file.genes)
  }
  
  println("Generating zero-filled DESeq-matrix skeleton.")
  all.genes.unique <- toupper(sort(unique(all.genes)))
  all.tfs.unique <- toupper(sort(unique(all.tfs)))
  
  # 0-initialized matrix  
  scores.empty <- matrix(0, nrow = length(all.genes.unique), ncol = length(all.tfs.unique), 
                         dimnames = list(all.genes.unique, all.tfs.unique))
  return(scores.empty)
}

# Extract p-values and log2(foldchange) values from DESeq results files and 
# fill pre-allocated confidence score matrix according to formula in Computational Methods:
# score = pval * sign(log2(foldchange))
populate.deseq.scores <- function(deseq.scores) {
  println("Extract DESeq non-adjusted p-values and log2(foldchange) from files.")
  for(i in deseqfiles) {
    # read in the data and extract the library name
    # DESeq p-values have too many numbers and are makes it a factor. stringsAsFactors option prevents that
    # but then we have a 'character' p-value column which still needs conversion to 'numeric' class
    cst <- read.table(i, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    tf <- extract.tf(colnames(cst), i)
    if(is.na(tf) || nchar(tf) == 0) next
    
    # get the DESeq p-values by iterating and accessing matrix via id and TF-name (genes are not ordered by name in DESeq files!)
    gene.id.cols <- grep(c("\\b(.*gene[_|\\.])?id($|\n)"), colnames(cst), ignore.case = TRUE, perl = TRUE, value = TRUE)
    gene.ids <- toupper(cst[, gene.id.cols[[1]][1]]) # first element if multiple matches
    idx <- 1
    for(gn in gene.ids) {
      if(nchar(gn) == 0) next
      deseq.scores[gn, tf] <- -log10(cst$pval[idx]) * sign(cst$log2FoldChange[idx])
      idx <- idx + 1
    }
  }
  
  # replace NA and Inf values in matrix with 0s (for later ranking procedure)
  deseq.scores[deseq.scores == Inf] <- 0
  deseq.scores[is.na(deseq.scores)] <- 0
  
  # drop 0-only-rows
  # println(dim(deseq.scores))
  # deseq.scores <- deseq.scores[rowSums(abs(deseq.scores[, -1]))>(1e-10),]
  # println("DROPPED---------------------------------------")
  # println(dim(deseq.scores))
  
  return(deseq.scores)
}

# ----------------
# Main function: load & process DEseq data and return the confidence score matrix S(DESeq)
# ----------------
load.deseq <- function() {
  empty.score.mat <- get.skel.mat()
  deseq.scores <- populate.deseq.scores(empty.score.mat)
  return(deseq.scores)
}