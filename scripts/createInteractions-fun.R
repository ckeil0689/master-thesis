# Load confidence score matrices by option
write.mat <- function(mat, outpath, combo, type, append = FALSE) {
  filename = paste0(outpath, combo, "_", type,".xlsx")
  print(paste("Writing matrix to file:", filename))
  
  # remove header row if appending
  if(append) {
    mat <- mat[-1,]
  }
  write.table(mat, file = filename, append = append, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# Takes a combined matrix file and transforms it to a list of node-node interactions that can be loaded into Cytoscape
create.interactions <- function(combomat, outpath, combo, type, pos.edge = "positive", neg.edge = "negative", append = FALSE) {
  print("Transforming matrix to node-node-value list.")
  
  # Select top 20% of edges from signed combined matrix
  m.cut <- quantile(combomat, probs=0.8)
  print(paste("Determined cut:", m.cut))
  
  abs.cut <- GLOBAL[["abs.cut"]]
  tot <- length(combomat[abs(combomat > abs.cut)])
  
  #pre-allocate data table since dimensions are known
  cyt.table <- data.table("nodeA"=as.character(rep(NA, tot)), "interaction"=as.character(rep("neutral", tot)), 
                          "nodeB"=as.character(rep(NA, tot)), "confidence_score"=as.double(rep(0, tot)))
  
  # Fill table with values from the combined matrix
  listrow <- 1
  # Iterate over TFs (columns)
  for (i in 1:ncol(combomat)) {
    if(i%%100==0) cat("\r", paste0("Progress: ", round((i*100/nrow(combomat)), digits = 0), "%"))
    # Per TF, only look at genes with absolute interaction value over the cutoff
    target.genes <- which(abs(combomat[,i]) > abs.cut)
    for (j in 1:length(target.genes)) {
      gene <- target.genes[j]
      val <- combomat[gene, i]
      
      if(val > abs.cut) {
        edge.type <- pos.edge
      } else if(val < abs.cut) {
        edge.type <- neg.edge
      } else {
        next # do not set any edge (list size limited to 'tot')
      }
      

      set(cyt.table, as.integer(listrow), "nodeA", colnames(combomat)[i])
      set(cyt.table, as.integer(listrow), "interaction", edge.type)
      set(cyt.table, as.integer(listrow), "nodeB", rownames(combomat)[j])
      set(cyt.table, as.integer(listrow), "confidence_score", val)
      
      listrow <- listrow + 1
    }
  }
  cat("\n")
  write.mat(cyt.table, outpath, combo, type, append)
  print("Done creating data table for Cytoscape.")
}