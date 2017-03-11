# TODO: could be combined with generate-chipseq-mat.r to reduce code
# Read a DeSeq files from dir
load.deseq <- function(dir, CORE_TFS) {
  setwd(dir)
  # DEseq results files
  deseqfiles <- list.files(getwd())
  
  # Vectors for row and column names of final KO-matrix
  all.genes <- c()
  all.tfs <- c()
  
  print(paste("Reading", length(deseqfiles),"DESeq files to create complete list of genes."))
        
  # Iteration 1: get list of all tested genes and TFs so a matrix can be set up
  for(i in deseqfiles) {
    # read in the data and extract the library name
    cst <- read.table(i, sep="\t", header=TRUE)
    
    # get TF name and add it to list
    tfcol <- colnames(cst)[3]
    tf <- strsplit(tfcol, "[.]")[[1]][2] # assumes convention column name (e.g. Th17.batf.wt -> batf)
    
    # only use target TFs
    if(!tf %in% CORE_TFS) {
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
  ko.scores <- matrix(0, nrow = length(all.genes.unique), ncol = length(all.tfs))
  # unique gene list makes up rows
  rownames(ko.scores) <- all.genes.unique
  # unique transcription factor list makes up columns
  colnames(ko.scores) <- toupper(all.tfs)
  
  print("Extract DESeq non-adjusted p-values and log2(foldchange) from files.")
  # Iteration 2: extract p-values and log2(foldchange) values from DESeq results files and 
  # fill confidence score matrix according to formula in Computational Methods
  for(i in deseqfiles) {
    # read in the data and extract the library name
    cst <- read.table(i, sep="\t", header=TRUE)
    
    # get TF name and add it to list
    tfcol <- colnames(cst)[3]
    tf <- toupper(strsplit(tfcol, "[.]")[[1]][2]) # assumes convention column name (e.g. Th17.batf.wt -> batf)
    
    # only use target TFs
    if(!tolower(tf) %in% CORE_TFS) {
      print(paste("Transcription factor not in target group:", tf, "(skipped)"))
      next
    }
    
    # get the DESeq p-values by iterating and accessing matrix via id and TF-name
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