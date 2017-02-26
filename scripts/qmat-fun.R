# Takes a ranked S-matrix and converts it to a Q-matrix 
# (see Computational Methods paper, http://www.cell.com/action/showImagesData?pii=S0092-8674%2812%2901123-3 ) 

# convert non-zero confidence score into quantile score that ranges 
# from zero (lowest confidence) to 1 (highest confidence) 
qscore <- function(rank, n_nonzero) {
  (1 - (rank/n_nonzero))
}

calc.qmat <- function(orig_mat, ranked_mat) {
  print("Calculating sign matrix.")
  sign_mat <- sign(as.data.frame(orig_mat))
  
  print("Applying quantile score function over ranked matrix.")
  n_nonzero = sum(orig_mat != 0)
  qmat <- as.matrix(mapply(qscore, as.data.frame(ranked_mat), MoreArgs = list(n_nonzero)))
  rownames(qmat) <- rownames(orig_mat)
  
  # Carry over signs from original matrix 
  # TODO move this outside later?
  qmat <- qmat %*% sign_mat
  
  return(qmat)
}