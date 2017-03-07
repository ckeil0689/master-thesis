# Load confidence score matrices by option
write.mat <- function(mat, outpath, combo, type, append = FALSE) {
  filename = paste0(outpath, combo, "_", type,".xlsx")
  print(paste("Writing matrix to file:", filename))
  write.table(mat, file = filename, append = append, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# Takes a combined matrix file and transforms it to a list of node-node interactions that can be loaded into Cytoscape
create.interactions <- function(combomat, outpath, combo, type, append = FALSE) {
  print("Transforming matrix to node-node-value list.")
  
  # Select top 20% of edges from signed combined matrix
  m.cut <- quantile(combomat, probs=0.8)
  
  print(paste("Combomat uncut:", dim(combomat)))
  print(paste("Combomat cut:", dim(m.cut)))
  
  cut <- GLOBAL[["abs.cut"]]
  pos <- length(combomat[combomat > cut])
  neg <- length(combomat[combomat < -1*cut])
  
  tot <- pos + neg
  
  print(paste("Positive edges:", pos))
  print(paste("Negative edges:", neg))
  print(paste("Edges to write:", tot))
  
  #pre-allocate data table since dimensions are known
  cyt.table <- data.table("nodeA"=as.character(rep(NA, tot)), "interaction"=as.character(rep("neutral", tot)), 
                          "nodeB"=as.character(rep(NA, tot)), "confidence_score"=as.double(rep(0, tot)))
  
  # Fill table with values from the combined matrix
  listrow <- 1
  for (i in 1:nrow(combomat)) {
    if(i%%100==0) cat("\r", paste0("Progress: ", round((i*100/nrow(combomat)), digits = 0), "%"))
    for (j in 1:ncol(combomat)) {
      val <- combomat[i,j]
      
      if(val > cut) {
        edge.type <- "positive_KC"
      } else if(val < -1*cut) {
        edge.type <- "negative_KC"
      } else {
        next # do not set any edge (list size limited to 'tot')
      }
      
      set(cyt.table, as.integer(listrow), "nodeA", colnames(combomat)[j])
      set(cyt.table, as.integer(listrow), "interaction", edge.type)
      set(cyt.table, as.integer(listrow), "nodeB", rownames(combomat)[i])
      set(cyt.table, as.integer(listrow), "confidence_score", val)
      
      listrow <- listrow + 1
    }
  }
  cat("\n")
  write.mat(cyt.table, outpath, combo, type, append)
  print("Done creating data table for Cytoscape.")
}