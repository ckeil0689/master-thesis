# Write a matrix to tab-delimited .txt-file
write.mat <- function(mat, outpath, prefix, suffix) {
  filename = paste0(outpath, prefix, "_", suffix, ".txt")
  print(paste("Writing matrix to file:", filename))
  write.table(mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)
}

# Load z-scores located at the given path
# Returns: a 1-column matrix of z-scores with genes as row names 
load.zscores <- function(zscores.path) {
  if(file.exists(zscores.path)) {
    zscore.table <- read.table(zscores.path, sep="\t", header=TRUE)
    zscore_col <- "Th17.Th0.zscores"
    zscores.all <- as.matrix(zscore.table[, zscore_col])
    rownames(zscores.all) <- zscore.table[, "Gene_id"]
    colnames(zscores.all) <- "Th17_vs_Th0_Zscores" # matches name used in Cytoscape Style for example KC.cys
    return(zscores.all)
  }
  print(paste("Z-score table not found at: ", zscores.path))
  print("Returning NULL.")
  return(NULL)
}

# Filters a table of genes by their z-score using the globally defined z.abs.cut
filter.genes.by.zscore <- function(zscores.all) {
  genes.final.idx <- which(abs(zscores.all) > GLOBAL[["z.abs.cut"]])
  genes.final <- rownames(zscores.all)[genes.final.idx]
  return(genes.final)
}