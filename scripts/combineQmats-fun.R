get.genes <- function(qmatlist) {
  # Only keep non-NULL row names (genes)
  genes <- c()
  for(i in qmatlist) {
    if(is.null(i)) next
    genes <- c(genes, rownames(i))
  }
  if(length(genes) == 0) stop("No genes could be extracted from the Q-matrices. Stopping.")
  return(toupper(sort(unique(genes))))
}

get.tfs <- function(qmatlist) {
  # Only keep non-NULL col names (tfs)
  tfs <- c()
  for(i in qmatlist) {
    if(is.null(i)) next
    tfs <- c(tfs, colnames(i))
  }
  if(length(tfs) == 0) stop("No transcription factors could be extracted from the Q-matrices. Stopping.")
  return(toupper(sort(unique(tfs))))
}

# Combine previously generated Q-matrices into a single matrix
combine.qmats <- function(kqmat, cqmat, rqmat, iqmat) {
  
  qmatlist <- list(kqmat, cqmat, rqmat, iqmat)
  
  # Vectors for row and column names of final Thx (x=0/=17) matrix
  final.genes <- get.genes(qmatlist)
  final.tfs <- get.tfs(qmatlist)

  # 0-initialized matrix  
  combo.mat <- matrix(0, nrow = length(final.genes), ncol = length(final.tfs), 
                         dimnames = list(final.genes, toupper(sort(final.tfs))))
  
  # Performing Q-matrix addition to combine data types
  for(i in qmatlist) {
    if(is.null(i)) next
    # from https://stackoverflow.com/questions/26042738/r-add-matrices-based-on-row-and-column-names
    A_df=as.data.frame(as.table(combo.mat))
    B_df=as.data.frame(as.table(i))
    merged_df = rbind(A_df, B_df)
    combo.mat = acast(merged_df, Var1 ~ Var2, sum)
  }
  
  println("Finished integration of data.")
  return(combo.mat)
}