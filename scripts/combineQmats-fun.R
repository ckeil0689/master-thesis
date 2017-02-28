# Combine previously generated Q-matrices into a single matrix
combine.qmats <- function(kqmat, cqmat, rqmat, iqmat, CORE_TFS) {
  
  print(paste("Dim K-qmat:", dim(kqmat)))
  
  qmatlist <- list(kqmat, cqmat, rqmat, iqmat)
  
  print(paste("Qmatlist:", length(qmatlist)))
  # Vectors for row and column names of final Thx (x=0/=17) matrix
  genes <- toupper(sort(rownames(kqmat)))
  
  # print(paste("Unique genes (total):", length(all_genes_thx_unique)))
  print(paste("Unique genes (total):", length(genes)))

  # 0-initialized matrix  
  combined_mat <- matrix(0, nrow = length(genes), ncol = length(CORE_TFS))
  # unique gene list makes up rows
  rownames(combined_mat) <- genes
  # unique transcription factor list makes up columns
  colnames(combined_mat) <- toupper(sort(CORE_TFS))
  
  print("Filling combination matrix with summed values from Q-matrices.")
  # Iteration 2: performing Q-matrix addition to combine data types
  for(i in qmatlist) {
    if(is.null(i)) next

    # from https://stackoverflow.com/questions/26042738/r-add-matrices-based-on-row-and-column-names
    A_df=as.data.frame(as.table(combined_mat))
    B_df=as.data.frame(as.table(i))
    
    merged_df=rbind(A_df,B_df)
    
    combined_mat=acast(merged_df, Var1 ~ Var2, sum)
  }
  
  combined_mat <- combined_mat[genes,] #TODO fix
  
  print(paste("Combined matrix dim:", dim(combined_mat)))
  print("Finished integration of data.")
  return(combined_mat)
}