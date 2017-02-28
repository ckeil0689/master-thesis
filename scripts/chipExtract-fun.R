# !/usr/bin/env Rscript
# ----------------
# Load & process ChIP-seq data and return the confidence score matrix S(ChIP)
# ----------------
load.chip <- function(dir, reflibfile, thx, CORE_TFS) {
  print("Reading ChIP files to create lists of TFs and genes.")
  setwd(dir)
  ref_file <- read.table(reflibfile, sep=",", header=TRUE)
  
  # ChIP-seq results files from GEO (assumes local copy in chipdir)
  chipfiles <- list.files(getwd())
  thx_rows <- c()
  
  # Vectors for row and column names of final Thx (x=0/=17) matrix
  all_genes_thx_wt <- c()
  tfs_thx <- c()
  
  # Iteration 1: create unique, maximal list of tested TFs and genes with MACS peaks
  print("Finding unique list of all genes tested in all ChIP files.")
  for(i in chipfiles) {
    # skip non-ChIP-seq files
    if(!grepl('^GSM[0-9]*(_SL[0-9]{1,9}){2}_genes.txt$', i)) {
      next
    }
    
    libname <- gsub('^GSM[0-9]*_(SL[0-9]{1,9})_.*.txt$','\\1', basename(i))
    
    # use library name to get tf from reference for each condition
    row <- match(libname, ref_file$Library_ID)
    if(tolower(ref_file$Condition[row])==thx && ref_file$Genotype[row]=='wt' && ref_file$Time_Point[row]=='') {
      tf <- ref_file$Factor[row]
      # correct for hyphens and underscores in TF names
      tf <- gsub("(-|_)", "", tolower(tf))
      tf <- gsub("cmaf", "maf", tf) # only these TF names are inconsistent all the time...
      tf <- gsub("rorg", "rorc", tf)
      
      # only use target TFs
      if(!tf %in% CORE_TFS) {
        next
      }
      
      thx_rows <- c(thx_rows, row)
      
      tf <- paste0(tf, "-", libname)
      tfs_thx <- c(tfs_thx, tf)
      
      # read in the data and extract the library name
      cst <- read.table(i, sep="\t", header=TRUE)
      # create full gene list by adding all genes from this file
      genes_thx_wt <- as.character(cst$Gene_ID)
      all_genes_thx_wt <- append(all_genes_thx_wt, genes_thx_wt)
      
    } else {
      next
    }
  }
  
  # generate the empty Th17 and Th0 matrices for ChIP-seq
  all_genes_thx_unique <- sort(unique(all_genes_thx_wt))
  
  print(paste("Generate a zero-filled matrix skeleton.", "(genes =", length(all_genes_thx_unique), ", TF columns =", length(tfs_thx), ")"))
  # 0-initialized matrix  
  thx_mat <- matrix(0, nrow = length(all_genes_thx_unique), ncol = length(tfs_thx))
  # unique gene list makes up rows
  rownames(thx_mat) <- all_genes_thx_unique
  # unique transcription factor list makes up columns
  colnames(thx_mat) <- tfs_thx
  
  print("Extract Poisson model p-values from ChIP-seq files.")
  # iteration 2: extract and assign associated Poisson model p-values to the matrix
  for(i in chipfiles) {
    
    # skip non-ChIP-seq files
    if(!grepl('^GSM[0-9]*(_SL[0-9]{1,9}){2}_genes.txt$', i)) {
      next
    }
    
    libname <- gsub('^GSM[0-9]*_(SL[0-9]{1,9})_.*.txt$','\\1', basename(i))
    
    row <- match(libname, ref_file$Library_ID)
    if(tolower(ref_file$Condition[row])==thx && ref_file$Genotype[row]=='wt' && ref_file$Time_Point[row]=='') {
      # get tf from reference for Thx/wt condition
      tf <- ref_file$Factor[row]
      
      # fix hyphens and underscores in TF names
      tf <- gsub("(-|_)", "", tolower(tf))
      tf <- gsub("cmaf", "maf", tf) # only these TF names are inconsistent all the time...
      tf <- gsub("rorg", "rorc", tf)
      
      # only use target TFs
      if(!tf %in% CORE_TFS) {
        next
      }
      
      tf <- paste0(tf, "-", libname)
      
      # order the genes, get index to also reorder Poisson model p-values
      cst <- read.table(i, sep="\t", header=TRUE)
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
  
  print("Create sorted, unique TF list from loaded files.")
  # remove library suffix from transcription factor names and create a sorted, unique TF list
  tfs_thx_unique <- gsub("-(SL[0-9]{1,9})$", "", tfs_thx)
  tfs_thx_unique <- sort(unique(tfs_thx_unique))
  
  print(paste("Generate a zero-filled confidence score matrix skeleton.", "(genes =", length(all_genes_thx_unique), ", TFs (unique) =", length(tfs_thx_unique), ")"))
  thx_unique_mat <- matrix(0, nrow = length(all_genes_thx_unique), ncol = length(tfs_thx_unique))
  rownames(thx_unique_mat) <- all_genes_thx_unique
  colnames(thx_unique_mat) <- tfs_thx_unique
  
  cols <- colnames(thx_mat)
  
  print("Calculate actual confidence scores (mean p-values for unique TFs)...")
  # TODO check if 0s should be counted in mean
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
    #subsetCols[subsetCols == 0] <- NA # exclude zeroes from mean (debug)
    # now matrix has a column for each library file - take the mean for each TF and write that in ONE column (final result: one column per unique TF)
    tf_meancol <- rowMeans(subsetCols, na.rm = TRUE)
    thx_unique_mat[,i] <- tf_meancol
  }
  
  # finally let labels be capital letters
  rownames(thx_unique_mat) <- toupper(all_genes_thx_unique)
  colnames(thx_unique_mat) <- toupper(tfs_thx_unique)
  
  print("Completed generation of ChIP-seq confidence score matrix.")
  return(thx_unique_mat)
}