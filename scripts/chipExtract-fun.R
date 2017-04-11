# !/usr/bin/env Rscript
# ----------------
# Constants
# ----------------
# ChIP-seq file directory relative to /scripts/
chipdir <- paste0(getwd(), "/../suppl/data/chipseq/")

# Experiment-library reference table; relative to /scripts/ directory
reflibfile <- paste0(getwd(), "/../suppl/mmc4.csv")
if(!file.exists(reflibfile)) stop(paste("Reference file does not exist, cannot load ChIP-files:", reflibfile))
ref.table <- read.table(reflibfile, sep=",", header=TRUE)

# Selection of ChIP-seq results files from GEO (assumes local copy in chipdir)
# The chosen files here correspond with the ones which the original authors
# have found to be qualitative superior results files. Each TF has one Th0 and one Th17 file.
th0_chipfiles <- c("GSM1004785_SL3192_SL3190_genes.txt", # BATF wt
                   "GSM1004824_SL1235_SL1234_genes.txt", # IRF4 wt
                   "GSM1004798_SL4424_SL4425_genes.txt", # MAF wt
                   "GSM1004857_SL3780_SL3778_genes.txt", # STAT3 wt
                   "GSM1004853_SL3779_SL3778_genes.txt", # RORC wt
                   "GSM1004808_SL6500_SL6499_genes.txt") # FOSL2 wt

th17_chipfiles <- c("GSM1004787_SL3037_SL3036_genes.txt", # BATF wt
                    "GSM1004833_SL2872_SL2876_genes.txt", # IRF4 rorc wt
                    "GSM1004800_SL3032_SL2871_genes.txt", # MAF wt
                    "GSM1004865_SL3315_SL3319_genes.txt", # STAT3 rorc wt
                    "GSM1004856_SL2870_SL2871_genes.txt", # RORC wt
                    "GSM1004809_SL6498_SL6497_genes.txt") # FOSL2 wt

p300_th0_chipfile <- "GSM1004842_SL1948_SL1947_genes.txt" # P300 Th0 wt
p300_th17_chipfile <- "GSM1004851_SL3594_SL3592_genes.txt" # P300 Th17 rorc wt

# ----------------
# Sub-routines
# ----------------
# Correct for hyphens and underscores in TF names
fix.tf.name <- function(tf.name) {
  tf.name <- gsub("(-|_)", "", tolower(tf.name))
  tf.name <- gsub("cmaf", "maf", tf.name) # only these TF names are inconsistent all the time...
  tf.name <- gsub("rorg", "rorc", tf.name)
  return(tf.name)
}

# Use library name to get tf from reference for each condition
# exp.name - the experiment name for lookup in ref.table
extract.tf.from.ref <- function(exp.name, boost.p300) {
  if(is.null(exp.name)) stop("No experiment passed. Stopping.")
  libname <- gsub('^GSM[0-9]*_(SL[0-9]{1,9})_.*_genes\\.txt$','\\1', basename(exp.name))
  row <- match(libname, ref.table$Library_ID)
  if(is.na(row)) stop("No library match found.")
  tf <- ref.table$Factor[row]
  tf <- fix.tf.name(tf)
  if(tf == "") stop("Could not identify TF from given experiment file.")
  
  # only use target TFs or p300 (when boosting)
  if(!tf %in% GLOBAL[["CORE_TFS"]]) {
    if(!(boost.p300 && tf == "p300")) {
      stop("TF not a core TF and not p300 while filling the ChIP-seq activator matrix.")
    }
  }
  
  thx <- tolower(ref.table$Condition[row])
  tf <- paste0(tf, "-", thx)
  return(tf)
}

# Create unique, maximal list of tested TFs and genes with MACS peaks and form an empty matrix
get.skel.matrix <- function(all_chipfiles, boost.p300) {
  # Ensure we are in correct directory
  if(!dir.exists(chipdir)) stop("Cannot load ChIP-files because the directory does not exist.")
  setwd(chipdir)
  # Vectors for row and column names of final Thx (x=0/=17) matrix
  genes.all <- c()
  tf.list <- c()
  
  println("Finding unique list of all genes tested in all ChIP files.")
  for(i in all_chipfiles) {
    if(!file.exists(i)) stop("Not all required ChIP files are present. Stopping.")
    tf <- extract.tf.from.ref(i, boost.p300)
    if(is.null(tf)) next
    
    tf.list <- c(tf.list, tf)
    
    # read in the data and extract the library name
    cst <- read.table(i, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    # create full gene list by adding all genes from this file
    genes <- as.character(cst$Gene_ID)
    genes.all <- append(genes.all, genes)
  }
  
  # generate the empty Th17 and Th0 matrices for ChIP-seq
  genes.unique <- sort(unique(genes.all))
  
  println(paste0("Generate a zero-filled matrix skeleton.", " (genes = ", length(genes.unique), ", TF columns = ", length(tf.list), ")"))
  thx_mat <- matrix(0, nrow = length(genes.unique), ncol = length(tf.list), dimnames = list(genes.unique, tf.list))
  return(thx_mat)
}

# Extract and assign associated Poisson model p-values to a matrix
get.pois.vals <- function(pois.mat, all_chipfiles, boost.p300) {
  println("Extract Poisson model p-values from ChIP-seq files.")
  for(i in all_chipfiles) {
    if(!file.exists(i)) stop("Not all required ChIP files are present. Stopping.")
    tf <- extract.tf.from.ref(i, boost.p300)
    if(is.null(tf) || !(tf %in% colnames(pois.mat))) next
    
    # Order the genes, get index to also reorder Poisson model p-values
    cst <- read.table(i, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    genes <- cst$Gene_ID
    
    # Get the Poisson p-values by iterating and accessing matrix via Gene_ID and TF-name
    idx <- 1
    for(j in genes) {
      # There will be a warning if conversion produces an NA value (e.g. non-numeric string)
      pois.val <- suppressWarnings(as.numeric(cst$genewide_pois_model_pval[idx]))
      if(is.na(pois.val)) {
        pois.val <- 0
      }
      pois.mat[j, tf] <- as.numeric(pois.val)
      idx <- idx + 1
    }
  }
  return(pois.mat)
}

# Make sure all data in the data.frame is numeric - if not then abort
convert.to.numeric <- function(df) {
  for(i in 1:ncol(df)) {
    tryCatch({
      df[,i] <- as.numeric(as.character(df[,i]))
    }, error = function(e) {
      println(e)
      stop("Error when converting data frame columns to numeric type.")
    })
  }
  
  # Convert NA values to 0 --> treated as no measured activity
  df[is.na(df)] = 0
  
  return(df)
}

# Calculate ChIP-scores: score = TF_Th17 - TF_Th0 (+ p300-score if activator)
calc.chipscores <- function(pois.mat, genes.unique, tfs.list.unique, boost.p300) {
  println(paste("Generate a zero-filled confidence score matrix skeleton.", "(genes =", length(genes.unique), ", TFs (unique) =", length(tfs.list.unique), ")"))
  chipscores <- matrix(0, nrow = length(genes.unique), ncol = length(tfs.list.unique), dimnames = list(genes.unique, tfs.list.unique))
  cols <- colnames(pois.mat)
  
  if(boost.p300) {
    p300.th17.col <- pois.mat[,"p300-th17"]
    p300.th0.col <- pois.mat[,"p300-th0"]
    p300.df <- data.frame(p300.th17.col, p300.th0.col)
    p300.df <- convert.to.numeric(p300.df)
    p300.score.col <- (p300.df$p300.th17.col - p300.df$p300.th0.col)
  }
  
  println("Calculate final ChIP-seq score matrix (cs = Th17 - Th0)...")
  # The Th17 and Th0 column for each TF in pois.mat will be replaced with a single TF column
  # Th0 value will be substracted from Th17 to create the final chip score 
  for(i in tfs.list.unique) {
    # patterns to match
    p_th0 <- paste0(i, "-th0")
    p_th17 <- paste0(i, "-th17")
    
    # Assumes one Th17 and one Th0 experiment for each TF
    for(j in cols) {
      if(grepl(p_th0, j)) {
        th0_col <- pois.mat[,j]
      } else if(grepl(p_th17, j)) {
        th17_col <- pois.mat[,j]
      } else {
        next
      }
    }
    
    df <- data.frame(th17_col, th0_col)
    df <- convert.to.numeric(df)
    
    if(boost.p300) {
      th.diff.col <- (df$th17_col - df$th0_col) + p300.score.col
    } else {
      th.diff.col <- (df$th17_col - df$th0_col)
    }
    
    chipscores[,i] <- th.diff.col
  }
  
  # finally make labels capital
  rownames(chipscores) <- toupper(genes.unique)
  colnames(chipscores) <- toupper(tfs.list.unique)
  
  return(chipscores)
}

# ----------------
# Main function: load & process ChIP-seq data and return the confidence score matrix S(ChIP)
# ----------------
load.chip <- function(boost.p300 = FALSE) {
  println("Reading ChIP files to create lists of TFs and genes.")
  
  # Define files to consider
  if(boost.p300) {
    println("ChIP-scores activator with p300 boost files.")
    all_chipfiles <- c(th0_chipfiles, th17_chipfiles, p300_th0_chipfile, p300_th17_chipfile)
  } else {
    println("ChIP-scores repressor/absolute without p300 boost files.")
    all_chipfiles <- c(th0_chipfiles, th17_chipfiles)
  }
  
  # Get an empty skeleton matrix (to avoid dynamic memory reallocation in R...)
  thx_mat <- get.skel.matrix(all_chipfiles, boost.p300)
  genes.unique <- rownames(thx_mat)
  tf.list <- colnames(thx_mat)
  
  # Fill the matrix with Poisson model p-values as found in the files
  pois.mat <- get.pois.vals(thx_mat, all_chipfiles, boost.p300)
  
  # debug only --- remove
  filename <- paste0(getwd(), "/../analysis/non-compressed-chip.txt")
  write.table(pois.mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)
  
  println("Create sorted, unique TF list from loaded files.")
  # remove library suffix from transcription factor names and create a sorted, unique TF list
  tfs.list <- gsub("-(th0|th17)$", "", tf.list)
  tfs.list.unique <- sort(unique(tfs.list))
  
  # Drop p300 - it was only needed for boosting (TODO: confirm if this is correct)
  if(boost.p300) {
    tfs.list.unique <- tfs.list.unique[tfs.list.unique !="p300"]
  }
  
  # Get the final confidence score matrix for ChIP-seq data (if activator or repressor depends on p300 boost)
  chipscores <- calc.chipscores(pois.mat, genes.unique, tfs.list.unique, boost.p300)
  
  # remove all 0-only-rows
  # chipscores <- chipscores[rowSums(abs(chipscores[, -1]))>(1e-10),]
  
  println("Completed generation of ChIP-seq confidence score matrix.")
  return(chipscores)
}