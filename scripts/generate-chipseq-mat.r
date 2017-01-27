# Read a ChIP-seq files from dir
setwd("/home/chrisk/Documents/uni/thesis/suppl/data/chipseq")

args <- commandArgs(trailingOnly=TRUE)
thx <- toupper(as.character(args[1]))
thx <-'Th17'

# ensure we have a usable argument to work with
if(thx != 'Th17' && thx != 'Th0') {
  print("Incorrect or missing argument.")
  print("Usage: script.R <(Th17|Th0)>")
  stop()
}

# contains info about ChIP experiments
ref_file <- read.table("/home/chrisk/Documents/uni/thesis/suppl/mmc4.csv", sep=",", header=TRUE)
thx_rows <- c()

# ChIP-seq results files
chipfiles <- list.files(getwd())

# Vectors for row and column names of final Thx (x=0/=17) matrix
all_genes_thx_wt <- c()
tfs_thx <- c()

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
    thx_rows <- c(thx_rows, row)
    tf_thx_wt <- ref_file$Factor[row]
    
    # correct for hyphens and underscores in TF names
    tf_thx_wt <- gsub("(-|_)", "", tolower(tf_thx_wt))
    tf_thx_wt <- paste0(tf_thx_wt, "-", libname)
    tfs_thx <- c(tfs_thx, tf_thx_wt)
    
    # create full gene list by adding all genes from this file
    genes_thx_wt <- as.character(cst$Gene_ID)
    all_genes_thx_wt <- append(all_genes_thx_wt, genes_thx_wt)
    
  } else {
    next
  }
}

# generate the empty Th17 and Th0 matrices for ChIP-seq
#tfs_thx_unique <- sort(unique(tfs_thx))
all_genes_thx_unique <- sort(unique(all_genes_thx_wt))

# 0-initialized matrix  
#thx_mat <- matrix(0, nrow = length(all_genes_thx_unique), ncol = length(tfs_thx_unique))
thx_mat <- matrix(0, nrow = length(all_genes_thx_unique), ncol = length(tfs_thx))
# unique gene list makes up rows
rownames(thx_mat) <- all_genes_thx_unique
# unique transcription factor list makes up columns
#colnames(thx_mat) <- tfs_thx_unique
colnames(thx_mat) <- tfs_thx

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
    tf_thx_wt <- ref_file$Factor[row]
    
    # fix hyphens and underscores in TF names
    tf_thx_wt <- gsub("(-|_)", "", tolower(tf_thx_wt))
    tf_thx_wt <- paste0(tf_thx_wt, "-", libname)
    
    # order the genes, get index to also reorder Poisson model p-values
    genes_thx_wt <- cst$Gene_ID
    
    # get the Poisson p-values by iterating and accessing matrix via Gene_ID and TF-name
    idx <- 1
    for(j in genes_thx_wt) {
      thx_mat[j, tf_thx_wt] <- cst$genewide_pois_model_pval[idx]
      idx <- idx + 1
    }
    
  } else {
    next
  }
}

# now matrix has a column for each library file - take the mean for each TF and write that in ONE column (result: one column per TF)
#z$mean <- rowMeans(subset(z, select = c(x, y)), na.rm = TRUE)

# write matrices to a tab-delimited file
filename=paste("C_", thx, "_mat.txt")
write.table(thx_mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)