# Takes a ranked S-matrix and converts it to a Q-matrix 
# (see Computational Methods paper, http://www.cell.com/action/showImagesData?pii=S0092-8674%2812%2901123-3 ) 

calc.qmat <- function(orig_mat, ranked_mat) {
  
  print("Calculating sign matrix.")
  sign_mat <- sign(orig_mat)
  
  print("Applying quantile score function over ranked matrix.")
  qmat <- mapply(qscore, ranked_mat, sign_mat, MoreArgs = list(n_nonzero = sum(orig_mat != 0)))
  rownames(qmat) <- rownames(orig_mat)
  
  return(qmat)
}

# convert non-zero confidence score into quantile score that ranges 
# from zero (lowest confidence) to 1 (highest confidence) 
qscore <- function(rank, signval, n_nonzero) {
  (1 - (rank/n_nonzero)) * signval
}