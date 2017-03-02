# Laod confidence score matrices by option
write.mat <- function(mat, outpath, prefix, thx, ext) {
  filename = paste0(outpath, prefix, thx, ext)
  print(paste("Writing matrix to file:", filename))
  write.table(mat, file = filename, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# Takes a combined matrix file and transforms it to a list of node-node interactions that can be loaded into Cytoscape
create.interactions <- function(combomat, outpath, combo, thx) {
  print("Transforming matrix to node-node-value list.")
  
  pos <- length(combomat[combomat > 1.5])
  neg <- length(combomat[combomat < -1.5])
  
  tot <- pos + neg
  
  print(paste("Total edges (pos):", pos))
  print(paste("Total edges (neg):", neg))
  
  print(paste("sif.table to write:", tot))
  
  #pre-allocate data table since dimensions are known
  sif.table <- data.table("nodeA"=as.character(rep(NA, tot)), "interaction"=rep("NA_kc", tot), "nodeB"=as.character(rep(NA, tot)))
  eda.table <- data.table("Activity"=as.character(rep(NA, tot)))
  
  # Fill table with values from the combined matrix
  listrow <- 1
  for (i in 1:nrow(combomat)) {
    if(i%%100==0) cat("\r", paste0("Progress: ", round((i*100/nrow(combomat)), digits = 0), "%"))
    for (j in 1:ncol(combomat)) {
      val <- combomat[i,j]
      if(abs(val) > 1.50) {
        
        signv <- sign(val)
        if(signv == 1) {
          ia <- "positive_kc"
        } else if(signv) {
          ia <- "negative_kc"
        } else {
          ia <- "NA_kc"
        }
        
        set(sif.table, listrow, "nodeA", rownames(combomat)[i])
        set(sif.table, listrow, "interaction", ia)
        set(sif.table, listrow, "nodeB", colnames(combomat)[j])
        # set(sif.table, listrow, "value", val)
        
        eda.entry <- paste0(rownames(combomat)[i], " (", ia, ") ", colnames(combomat)[j], "=", val)
        set(eda.table, listrow, "Activity", eda.entry)
        listrow <- listrow + 1
      }
    }
  }
  cat("\n")
  write.mat(sif.table, outpath, paste0(combo, "_"), thx, ".sif")
  write.mat(eda.table, outpath, paste0(combo, "_"), thx, ".eda.attrs")
  print("Done creating .sif and .eda for Cytoscape.")
  # return(sif.table)
}