# Combine previously generated Q-matrices into a single matrix
combine.qmats <- function(qmatfiles) {
  
  # Vectors for row and column names of final Thx (x=0/=17) matrix
  all_genes_thx_wt <- c()
  tfs_thx <- c()
  
  print("Finding all genes and TFs to build a matrix.")
  # Iteration 1: create unique, maximal list of TFs and genes
  for(qmat in qmatfiles) {
    if(is.null(qmat)) next
    
    tfs_thx <- c(tfs_thx, colnames(qmat))
    all_genes_thx_wt <- append(all_genes_thx_wt, rownames(qmat))
  }
  
  print("Generate zero-filled matrix skeleton.")
  # generate the empty Th17 and Th0 matrices for ChIP-seq
  all_genes_thx_unique <- toupper(sort(unique(all_genes_thx_wt)))
  tfs_thx_unique <- toupper(sort(unique(tfs_thx)))
  
  # 0-initialized matrix  
  combined_mat <- matrix(0, nrow = length(all_genes_thx_unique), ncol = length(tfs_thx_unique))
  # unique gene list makes up rows
  rownames(combined_mat) <- all_genes_thx_unique
  # unique transcription factor list makes up columns
  colnames(combined_mat) <- tfs_thx_unique
  
  print("Filling combination matrix with summed values from Q-matrices.")
  # Iteration 2: filling data by matrix addition
  for(qmat in qmatfiles) {
    if(is.null(qmat)) next
    
    cAB <- union(colnames(combined_mat), colnames(qmat))
    rAB <- union(rownames(combined_mat), rownames(qmat))
    
    A1 <- matrix(0, ncol=length(cAB), nrow=length(rAB), dimnames=list(rAB, cAB))
    B1 <- A1
    
    indxA <- outer(rAB, cAB, FUN=paste) %in% outer(rownames(combined_mat), colnames(combined_mat), FUN=paste) 
    indxB <- outer(rAB, cAB, FUN=paste) %in% outer(rownames(qmat), colnames(qmat), FUN=paste)
    A1[indxA] <- combined_mat
    B1[indxB] <- as.matrix(qmat)
    
    combined_mat <- A1 + B1
  }
  print("Finished integration of data.")
  return(combined_mat)
}