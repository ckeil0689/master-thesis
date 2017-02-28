# Combine previously generated Q-matrices into a single matrix
combine.qmats <- function(kqmat, cqmat, rqmat, iqmat, CORE_TFS) {
  
  print(paste("Dim K-qmat:", dim(kqmat)))
  
  qmatlist <- list(kqmat, cqmat, rqmat, iqmat)
  
  print(paste("Qmatlist:", length(qmatlist)))
  # Vectors for row and column names of final Thx (x=0/=17) matrix
  genes <- rownames(kqmat)
  tfs_thx <- CORE_TFS
  
  print("Finding all genes and TFs to build a matrix.")
  # # Iteration 1: create unique, maximal list of TFs and genes
  # for(i in qmatlist) {
  #   if(is.null(i)) next
  #   tfs_thx <- c(tfs_thx, colnames(i))
  #   all_genes_thx_wt <- append(all_genes_thx_wt, rownames(i))
  # }
  # 
  # print("Generate zero-filled matrix skeleton.")
  # # generate the empty Th17 and Th0 matrices for ChIP-seq
  # all_genes_thx_unique <- toupper(sort(unique(all_genes_thx_wt)))
  # tfs_thx_unique <- toupper(sort(unique(tfs_thx)))
  
  # print(paste("Unique genes (total):", length(all_genes_thx_unique)))
  print(paste("Unique genes (total):", length(genes)))
  
  # # 0-initialized matrix  
  # combined_mat <- matrix(0, nrow = length(all_genes_thx_unique), ncol = length(tfs_thx_unique))
  # # unique gene list makes up rows
  # rownames(combined_mat) <- all_genes_thx_unique
  # # unique transcription factor list makes up columns
  # colnames(combined_mat) <- tfs_thx_unique
  
  # 0-initialized matrix  
  combined_mat <- matrix(0, nrow = length(genes), ncol = length(CORE_TFS))
  # unique gene list makes up rows
  rownames(combined_mat) <- toupper(sort(genes))
  # unique transcription factor list makes up columns
  colnames(combined_mat) <- toupper(sort(CORE_TFS))
  
  print("Filling combination matrix with summed values from Q-matrices.")
  # Iteration 2: performing Q-matrix addition to combine data types
  for(i in qmatlist) {
    if(is.null(i)) next

    A_df=as.data.frame(as.table(combined_mat))
    B_df=as.data.frame(as.table(i))
    
    merged_df=rbind(A_df,B_df)
    
    combined_mat=acast(merged_df, Var1 ~ Var2, sum)
  }
  
  print(paste("Combined matrix dim:", dim(combined_mat)))
  print("Finished integration of data.")
  return(combined_mat)
}