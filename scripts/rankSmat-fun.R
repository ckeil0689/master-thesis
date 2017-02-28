# Takes a confidence score S-matrix and converts it to a rank matrix 
# (see Computational Methods paper, http://www.cell.com/action/showImagesData?pii=S0092-8674%2812%2901123-3 ) 

# applies rank function to matrix <matname> and writes it into the same directory as the source file
# wrapped in function for possible reuse or automation 
rank.smat <- function(smat) {
  print("Ranking.")
  # replace all zeroes or infinities with NA to exclude them from ranking
  smat[smat == 0] <- NA
  # apply rank to absolute value of confidence scores (descending order --> negative sign)
  ranked_mat <- matrix(rank(-abs(smat), na.last = "keep"), ncol=ncol(smat))
  # replace all NAs with zeroes again post-ranking
  # ranked_mat[is.na(ranked_mat)] <- 0
  
  rownames(ranked_mat) <- rownames(smat)
  colnames(ranked_mat) <- colnames(smat)
  
  return(ranked_mat)
}