# Load confidence score matrices by option
write.interactions <- function(mat, outpath, combo, type, used.cut, append = FALSE) {
  filename = paste0(outpath, combo, "_", type, "_", used.cut, "_cs-cut_", Sys.Date(), ".csv")
  println(paste("Writing matrix to file:", filename))
  # no column names when appending (otherwise it will be treated as random data entry by Cytoscape)
  write.table(mat, file = filename, append = append, sep = ",", row.names = FALSE, col.names = !append, quote = FALSE)
}

# Pre-allocates data table since final dimensions are known - performance is much better than dynamic resize during edge selection
create.empty.table <- function(total.edge.num) {
  # 4-column table: nodeA | nodeB | interaction | confidence_score
  empty.table <- data.table("nodeA"=as.character(rep(NA, total.edge.num)), "interaction"=as.character(rep("neutral", total.edge.num)), 
                          "nodeB"=as.character(rep(NA, total.edge.num)), "confidence_score"=as.double(rep(0, total.edge.num)))
  return(empty.table)
}

# Fill edge table for Cytoscape with TF-gene interaction values from the combined matrix (edge selection logic adapted from Aviv Madar)
select.edges <- function(combo.mat, cyt.table, used.cut, pos.edge, neg.edge) {
  listrow <- 1
  # Iterate over TFs (columns)
  for (i in 1:ncol(combo.mat)) {
    if(i%%100==0 && !GLOBAL[["TEST"]]) cat("\r", paste0("Progress: ", round((i*100/nrow(combo.mat)), digits = 0), "%"))
    # Per TF, only look at genes with absolute interaction value over the cutoff
    target.genes <- which(abs(combo.mat[,i]) > used.cut)
    if(length(target.genes) == 0) {
      # println(paste("No targets found. Skipping", colnames(combo.mat)[i]))
      next
    }
    for (j in 1:length(target.genes)) {
      gene <- target.genes[j]
      val <- combo.mat[gene, i]
      
      if(val > used.cut) {
        edge.type <- "positive_KC" #pos.edge
      } else if(val < used.cut) {
        edge.type <- "negative_KC" #neg.edge
      } else {
        next # do not set any edge (list size limited to 'tot')
      }
      
      set(cyt.table, as.integer(listrow), "nodeA", colnames(combo.mat)[i])
      set(cyt.table, as.integer(listrow), "interaction", edge.type)
      set(cyt.table, as.integer(listrow), "nodeB", rownames(combo.mat)[gene])
      set(cyt.table, as.integer(listrow), "confidence_score", val)
      
      listrow <- listrow + 1
    }
  }
  if(!GLOBAL[["TEST"]]) cat("\n")
  return(cyt.table)
}

# Takes a combined matrix file and transforms it to a list of node-node interactions that can be loaded into Cytoscape
create.interactions <- function(combo.mat, outpath, combo, type, pos.edge = "positive", neg.edge = "negative", append = FALSE) {
  println("Transforming matrix to node-node-value list.")
  
  # Select top 20% of edges from signed combined matrix
  m.cut <- quantile(combo.mat, probs=.97)
  
  # Set which cut value is used (quantile m.cut or defined cs.abs.cut)
  used.cut = GLOBAL[["cs.abs.cut"]]
  println(paste("Using absolute confidence score cut:", used.cut))
  
  # Info about the total number of edges
  tot <- length(combo.mat[abs(combo.mat) > used.cut])
  println(paste0(tot, " [", pos.edge, "]"))
  
  cyt.table.empty <- create.empty.table(tot)
  cyt.table.full <- select.edges(combo.mat, cyt.table.empty, used.cut, pos.edge, neg.edge)
  
  write.interactions(cyt.table.full, outpath, combo, type, used.cut, append)
  println("Done creating data table for Cytoscape.")
}