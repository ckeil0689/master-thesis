#!/usr/bin/env Rscript
# TODO: could be combined with generate-rnaseq-mat.r to reduce code
# Read a ChIP-seq files from dir
setwd(paste0(getwd(), "/../suppl/data/chipseq"))

args <- commandArgs(trailingOnly=TRUE)
thx <- toupper(as.character(args[1]))
thx <-'Th17'

# fixed target transcription factors to use
target_tfs <- c("batf", "irf4", "stat3", "hif1a", "maf", "fosl2", "rorc")

# ensure we have a usable argument to work with
if(thx != 'Th17' && thx != 'Th0') {
  print("Incorrect or missing argument.")
  print("Usage: script.R <(Th17|Th0)>")
  stop()
}

# contains info about ChIP experiments
ref_filepath <- paste0(getwd(), "/../../mmc4.csv")
ref_file <- read.table(ref_filepath, sep=",", header=TRUE)
thx_rows <- c()

# ChIP-seq results files
chipfiles <- list.files(getwd())

# Vectors for row and column names of final Thx (x=0/=17) matrix
all_genes_thx_wt <- c()
tfs_thx <- c()

print("Reading ChIP files to create lists of TFs and genes.")
# Iteration 1: create unique, maximal list of tested TFs and genes with MACS peaks
for(i in chipfiles) {
  
  # skip non-ChIP-seq files
  if(!grepl('^GSM[0-9]*(_SL[0-9]{1,9}){2}_genes.txt$', i)) {
    next
  }
  
  # read in the data and extract the library name
  cst <- read.table(i, sep="\t", header=TRUE)
  libname <- gsub('^GSM[0-9]*_(SL[0-9]{1,9})_.*.txt$','\\1', basename(i))
  
  # use library name to get tf from reference for each condition
  row <- match(libname, ref_file$Library_ID)
  if(ref_file$Condition[row]==thx && ref_file$Genotype[row]=='wt') {
    tf <- ref_file$Factor[row]
    # correct for hyphens and underscores in TF names
    tf <- gsub("(-|_)", "", tolower(tf))
    tf <- gsub("cmaf", "maf", tf) # only these TF names are inconsistent all the time...
    tf <- gsub("rorg", "rorc", tf)
    
    # only use target TFs
    if(!tf %in% target_tfs) {
      next
    }
    
    thx_rows <- c(thx_rows, row)

    tf <- paste0(tf, "-", libname)
    tfs_thx <- c(tfs_thx, tf)
    
    # create full gene list by adding all genes from this file
    genes_thx_wt <- as.character(cst$Gene_ID)
    all_genes_thx_wt <- append(all_genes_thx_wt, genes_thx_wt)
    
  } else {
    next
  }
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

print("Extract Poisson model p-values from ChIP files.")
# iteration 2: extract and assign associated Poisson model p-values to the matrix
for(i in chipfiles) {
  
  # skip non-ChIP-seq files
  if(!grepl('^GSM[0-9]*(_SL[0-9]{1,9}){2}_genes.txt$', i)) {
    next
  }
  
  cst <- read.table(i, sep="\t", header=TRUE)
  libname <- gsub('^GSM[0-9]*_(SL[0-9]{1,9})_.*.txt$','\\1', basename(i))

  row <- match(libname, ref_file$Library_ID)
  if(ref_file$Condition[row]==thx && ref_file$Genotype[row]=='wt') {
    # get tf from reference for Thx/wt condition
    tf <- ref_file$Factor[row]
    
    # fix hyphens and underscores in TF names
    tf <- gsub("(-|_)", "", tolower(tf))
    tf <- gsub("cmaf", "maf", tf) # only these TF names are inconsistent all the time...
    tf <- gsub("rorg", "rorc", tf)
    
    # only use target TFs
    if(!tf %in% target_tfs) {
      next
    }
    
    tf <- paste0(tf, "-", libname)
    
    # order the genes, get index to also reorder Poisson model p-values
    genes_thx_wt <- cst$Gene_ID
    
    # get the Poisson p-values by iterating and accessing matrix via Gene_ID and TF-name
    idx <- 1
    for(j in genes_thx_wt) {
      thx_mat[j, tf] <- cst$genewide_pois_model_pval[idx]
      idx <- idx + 1
    }
    
  } else {
    next
  }
}

print("Created sorted, unique TF list.")
# remove library suffix from transcription factor names and create a sorted, unique TF list
tfs_thx_unique <- gsub("-(SL[0-9]{1,9})$", "", tfs_thx)
tfs_thx_unique <- sort(unique(tfs_thx_unique))

print("Generate sorted, unique matrix skeleton.")
thx_unique_mat <- matrix(0, nrow = length(all_genes_thx_unique), ncol = length(tfs_thx_unique))
rownames(thx_unique_mat) <- all_genes_thx_unique
colnames(thx_unique_mat) <- tfs_thx_unique

cols <- colnames(thx_mat)

print("Calculate mean p-values for unique TFs...")
for(i in tfs_thx_unique) {
  
  # pattern to match
  p <- paste0(i, "-(SL[0-9]{1,9})$")
  tf_colset <- c()
  
  # extract all columns that match the current TF
  for(j in cols) {
    if(grepl(p, j)) {
      tf_colset <- c(tf_colset, j)
    }
  }
  
  subsetCols <- subset(thx_mat, select = tf_colset)
  # now matrix has a column for each library file - take the mean for each TF and write that in ONE column (final result: one column per unique TF)
  tf_meancol <- rowMeans(subset(thx_mat, select = tf_colset), na.rm = TRUE)
  #cbind(thx_unique_mat, tf_meancol)
  thx_unique_mat[,i] <- tf_meancol
}

print("Writing matrix to file...")
# debug
#filename=paste0("C_", thx, "debug_mat.txt")
#write.table(thx_mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)

# finally let labels be capital letters
rownames(thx_unique_mat) <- toupper(all_genes_thx_unique)
colnames(thx_unique_mat) <- toupper(tfs_thx_unique)

# write matrices to a tab-delimited file
filename=paste0("C_", thx, "_mat.txt")
write.table(thx_unique_mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)