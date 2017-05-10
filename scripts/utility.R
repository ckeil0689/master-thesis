# A wrapper for normal print to reduce code clutter
println <- function(string) {
  if(!GLOBAL[["TEST"]]) print(string)
}

# Write a matrix to tab-delimited .txt-file
write.mat <- function(mat, outpath, prefix, suffix = "") {
  # some input checking
  if(!(is.matrix(mat) || is.data.frame(mat) || is.table(mat))) {
    stop("First argument is not a matrix, data frame, or table.")
  } 
  
  if(!dir.exists(outpath)) {
    stop(paste("Given output path is not a directory:", outpath))
  }
  
  if(prefix != "" && suffix != "") {suffix <- paste0("_", suffix)}
  filename = paste0(outpath, prefix, suffix, ".txt")
  print(paste("Writing matrix to file:", filename)) # do not use println() here --> test will fail, output is wanted
  write.table(mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)
}

# Load z-scores located at the given path
# Returns: a 1-column matrix of z-scores with genes as row names 
load.zscores <- function(zscores.path) {
  if(file.exists(zscores.path)) {
    zscore.table <- read.csv(zscores.path)
    if(!"Gene_id" %in% colnames(zscore.table)) stop("mmc5.csv could not be properly read. 
                                             Make sure the file is really comma-separated. Stopping.")
    zscore.col.idx <- grep("th17.th0.zscores$", colnames(zscore.table), ignore.case = TRUE, perl = TRUE)
    zscore_col <- colnames(zscore.table)[zscore.col.idx]
    # Th17_vs_Th0_Zscores matches name used in Cytoscape Style for example KC.cys
    zscores.all <- as.matrix(zscore.table[, zscore_col])
    rownames(zscores.all) <- zscore.table[, "Gene_id"]
    colnames(zscores.all) <- "Th17_vs_Th0_Zscores"
    return(zscores.all)
  }
  println(paste("Z-score table not found at: ", zscores.path))
  return(NULL)
}

# Filters row names (gene ids) in a table of genes by their z-score using the globally defined z.abs.cut
filter.genes.by.zscore <- function(zscores.all, z.abs.cut) {
  if(is.null(zscores.all) || is.null(z.abs.cut)) {
    stop("A necessary parameter is missing. Cannot filter by zscores.")
  }
  if(!is.numeric(z.abs.cut)) {
    stop("zscore cut is not numeric.")
  }
  genes.final.idx <- which(abs(zscores.all) > z.abs.cut)
  genes.final <- rownames(zscores.all)[genes.final.idx]
  return(genes.final)
}