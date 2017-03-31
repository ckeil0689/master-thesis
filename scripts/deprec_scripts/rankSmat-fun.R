# Takes a confidence score S-matrix and converts it to a rank matrix 
# (see Computational Methods paper, http://www.cell.com/action/showImagesData?pii=S0092-8674%2812%2901123-3 ) 

# applies rank function to matrix <matname> and writes it into the same directory as the source file
# wrapped in function for possible reuse or automation 
rank.smat <- function(smat) {
  # replace all zeroes or infinities with NA to exclude them from ranking
  smat[smat == 0] <- NA
  # apply rank to absolute value of confidence scores (descending rank order --> negative sign before abs())
  ranked_mat <- matrix(rank(-abs(smat), na.last = "keep"), ncol=ncol(smat), dimnames = list(rownames(smat), colnames(smat)))
  # replace all NAs with zeroes again post-ranking
  ranked_mat[is.na(ranked_mat)] <- 0
  
  print("Done ranking.")
  qmat <- calc.qmat(smat, ranked_mat)
  return(qmat)
}

# convert non-zero confidence score into quantile score that ranges 
# from zero (lowest confidence) to 1 (highest confidence) 
qscore <- function(rank, n_nonzero) {
  (1 - (rank/n_nonzero))
}

calc.qmat <- function(orig_mat, ranked_mat) {
  print("Applying quantile score function over ranked matrix.")
  n_nonzero = sum(orig_mat != 0)
  df <- as.data.frame(ranked_mat)
  qmat <- as.matrix(mapply(qscore, df, MoreArgs = list(n_nonzero)))
  qmat[is.na(qmat)] <- 0
  rownames(qmat) <- rownames(orig_mat)
  return(qmat)
}