# Read a ChIP-seq files from dir
setwd("/home/chrisk/Documents/uni/thesis/suppl/data/chipseq")

# contains info about ChIP experiments
ref_file <- read.table("/home/chrisk/Documents/uni/thesis/suppl/mmc4.csv", sep=",", header=TRUE)
th17_rows <- c()
th0_rows <- c()

# ChIP-seq results files
chipfiles <- list.files(getwd())

# Vectors for row and column names of final matrices
all_genes_th17_wt <- c()
all_genes_th0_wt <- c()
tfs_th17 <- c()
tfs_th0 <- c()

#Iteration 1: create unique, maximal list of tested TFs and genes with MACS peaks
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
  if(ref_file$Condition[row]=='Th17' && ref_file$Genotype[row]=='wt') {
    th17_rows <- c(th17_rows, row)
    tf_th17_wt <- ref_file$Factor[row]
    
    # correct for hyphens and underscores in TF names
    tf_th17_wt <- gsub("(-|_)", "", tolower(tf_th17_wt))
    tfs_th17 <- c(tfs_th17, tf_th17_wt)
    
    # create full gene list by adding all genes from this file
    genes_th17_wt <- as.character(cst$Gene_ID)
    all_genes_th17_wt <- append(all_genes_th17_wt, genes_th17_wt)
    
  } else if(ref_file$Condition[row]=='Th0' && ref_file$Genotype[row]=='wt') {
    th0_rows <- c(th0_rows, row)
    tf_th0_wt <- ref_file$Factor[row]
    
    # correct for hyphens and underscores in TF names
    tf_th0_wt <- gsub("(-|_)", "", tolower(tf_th0_wt))
    tfs_th0 <- c(tfs_th0, tf_th0_wt)
    
    # create full gene list by adding all genes from this file
    genes_th0_wt <- as.character(cst$Gene_ID)
    all_genes_th0_wt <- append(all_genes_th0_wt, genes_th0_wt) 
    
  } else {
    next
  }
}

# Generate the empty Th17 and Th0 matrices for ChIP-seq
tfs_th17_unique <- sort(unique(tfs_th17))
tfs_th0_unique <- sort(unique(tfs_th0))
all_genes_th17_unique <- sort(unique(all_genes_th17_wt))
all_genes_th0_unique <- sort(unique(all_genes_th0_wt))

th17_mat <- matrix(0, nrow = length(all_genes_th17_unique), ncol = length(tfs_th17_unique))
rownames(th17_mat) <- all_genes_th17_unique
colnames(th17_mat) <- tfs_th17_unique

th0_mat <- matrix(0, nrow = length(all_genes_th0_unique), ncol = length(tfs_th0_unique))
rownames(th0_mat) <- all_genes_th0_unique
colnames(th0_mat) <- tfs_th0_unique

#Iteration 2: extract and assign associated Poisson model p-values to the matrices
for(i in chipfiles) {
  
  # skip non-ChIP-seq files
  if(!grepl('^GSM[0-9]*(_SL[0-9]{1,9}){2}_genes.txt$', i)) {
    next
  }
  
  cst <- read.table(i, sep="\t", header=TRUE)
  libname <- gsub('^GSM[0-9]*_(SL[0-9]{1,9})_.*.txt$','\\1', basename(i))
  
  #get tf from reference for Th17/wt condition
  row <- match(libname, ref_file$Library_ID)
  if(ref_file$Condition[row]=='Th17' && ref_file$Genotype[row]=='wt') {
    tf_th17_wt <- ref_file$Factor[row]
    
    # fix hyphens and underscores in TF names
    tf_th17_wt <- gsub("(-|_)", "", tolower(tf_th17_wt))
    
    
    # order the genes, get index to also reorder Poisson model p-values
    genes_th17_wt <- cst$Gene_ID
    
    # get the Poisson p-values by iterating and accessing matrix via Gene_ID and TF-name
    idx <- 1
    for(j in genes_th17_wt) {
      th17_mat[j, tf_th17_wt] <- cst$genewide_pois_model_pval[idx]
      idx <- idx + 1
    }
    
  } else if(ref_file$Condition[row]=='Th0' && ref_file$Genotype[row]=='wt') {
    tf_th0_wt <- ref_file$Factor[row]
    
    # fix hyphens and underscores in TF names
    tf_th0_wt <- gsub("(-|_)", "", tolower(tf_th0_wt))
    
    # order the genes, get index to also reorder Poisson model p-values
    genes_th0_wt <- cst$Gene_ID
    
    # get the Poisson p-values by iterating and accessing matrix via Gene_ID and TF-name
    # TODO --> MEAN for the different experiments
    idx <- 1
    for(j in genes_th0_wt) {
      th0_mat[j, tf_th0_wt] <- cst$genewide_pois_model_pval[idx]
      idx <- idx + 1
    }
  }
}

# Write matrices to a tab-delimited file
write.table(th17_mat, file = "C_th17_mat.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table(th0_mat, file = "C_th0_mat.txt", sep = "\t", row.names = TRUE, col.names = NA)