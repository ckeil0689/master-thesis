# Write a matrix to tab-delimited .txt-file
write.mat <- function(mat, outpath, prefix, suffix) {
  filename = paste0(outpath, prefix, "_", suffix, ".txt")
  print(paste("Writing matrix to file:", filename))
  write.table(mat, file = filename, sep = "\t", row.names = TRUE, col.names = NA)
}