# Takes a combined matrix file and transforms it to a list of node-node interactions that can be loaded into Cytoscape
create.interactions <- function(combomat) {
  print("Transforming matrix to node-node-value list.")
  
  pos <- length(combomat[combomat > 1.5])
  neg <- length(combomat[combomat < -1.5])
  
  tot <- pos + neg
  
  print(paste("Total edges (pos):", pos))
  print(paste("Total edges (neg):", neg))
  
  edge_num <- nrow(combomat) * ncol(combomat)
  print(paste("Edges to write:", edge_num))
  
  #pre-allocate data table since dimensions are known
  edges <- data.table("node1"=as.character(rep(NA, tot)), "node2"=as.character(rep(NA, tot)), "value"=rep(0, tot))
  
  # Fill table with values from the combined matrix
  listrow <- 1
  for (i in 1:nrow(combomat)) {
    if(i%%100==0) cat("\r", paste0("Progress: ", round((i*100/nrow(combomat)), digits = 0), "%"))
    for (j in 1:ncol(combomat)) {
      #edges <- rbind(edges, c(rownames(combomat)[i], colnames(combomat)[j], combomat[i,j]))
      val <- combomat[i,j]
      if(abs(val) >= 1.50) {
        #listrow <- ((i-1) * ncol(combomat) + j)
        set(edges, listrow, "node1", rownames(combomat)[i])
        set(edges, listrow, "node2", colnames(combomat)[j])
        set(edges, listrow, "value", val)
        listrow <- listrow + 1
      }
    }
  }
  cat("\n")
  print("Done creating interaction list for Cytoscape.")
  return(edges)
}