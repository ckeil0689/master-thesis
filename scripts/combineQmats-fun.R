# Combine previously generated Q-matrices into a single matrix
combine.qmats <- function(kqmat, cqmat, rqmat, iqmat, CORE_TFS) {
  
  qmatlist <- list(kqmat, cqmat, rqmat, iqmat)
  
  # Vectors for row and column names of final Thx (x=0/=17) matrix
  ko.genes <- toupper(sort(rownames(kqmat)))
  print(paste("Unique genes (total):", length(ko.genes)))

  # 0-initialized matrix  
  combined_mat <- matrix(0, nrow = length(ko.genes), ncol = length(CORE_TFS))
  # unique gene list makes up rows
  rownames(combined_mat) <- ko.genes
  # unique transcription factor list makes up columns
  colnames(combined_mat) <- toupper(sort(CORE_TFS))
  
  # Performing Q-matrix addition to combine data types
  for(i in qmatlist) {
    if(is.null(i)) next

    # from https://stackoverflow.com/questions/26042738/r-add-matrices-based-on-row-and-column-names
    A_df=as.data.frame(as.table(combined_mat))
    B_df=as.data.frame(as.table(i))
    
    print(B_df[1:5,])
    merged_df = rbind(A_df, B_df)
    
    combined_mat = acast(merged_df, Var1 ~ Var2, sum)
  }
  
  # Only include RNA-seq KO-ko.genes (derived from mmc4.xlsx from materials)
  combined_mat <- combined_mat[ko.genes,] 
  
  print("Finished integration of data.")
  return(combined_mat)
}