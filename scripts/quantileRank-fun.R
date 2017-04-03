# Takes a confidence score S-matrix and converts it to a rank matrix 
# (see Computational Methods paper, http://www.cell.com/action/showImagesData?pii=S0092-8674%2812%2901123-3 ) 

# applies rank function to matrix <matname> and writes it into the same directory as the source file
# wrapped in function for possible reuse or automation 
# positiveOnly - Exclude all negative and zero values from ranking and treat them as NA (used for all but ChIP)
calc.quantile.ranks <- function(smat, positiveOnly = FALSE) {
  # replace all zeroes or infinities with NA to exclude them from ranking
  if(positiveOnly) {
    smat[smat <= 0.0] <- NA
  } else {
    smat[smat == 0.0] <- NA
  }
  
  # Squash bad input#apply
  storage.mode(smat) <- "numeric"
  
  # This should only happen when ALL values are NAs at this point
  # Return unranked matrix
  if(!is.numeric(smat)) {
    smat[is.na(smat)] <- 0.0
    # Somehow we can be character type again here (???) so ensure numeric matrix
    storage.mode(smat) <- "numeric"
    return(smat)
  }
  
  # apply rank to absolute value of confidence scores (descending rank order --> negative sign before abs())
  ranked.mat <- matrix(rank(-abs(smat), na.last = "keep"), ncol=ncol(smat), dimnames = list(rownames(smat), colnames(smat)))
  # replace all NAs with zeroes again post-ranking
  smat[is.na(smat)] <- 0.0
  ranked.mat[is.na(ranked.mat)] <- 0.0
  
  print("Done ranking.")
  qmat <- calc.qmat(smat, ranked.mat)
  return(qmat)
}

# convert non-zero confidence score into quantile score that ranges 
# from zero (lowest confidence) to 1 (highest confidence) 
qscore <- function(rank, n_nonzero) {
  (1 - (rank/n_nonzero))
}

calc.qmat <- function(orig.mat, ranked.mat) {
  print("Applying quantile score function over ranked matrix.")
  ix.zero <- which(orig.mat == 0)
  n_nonzero = sum(orig.mat != 0)
  
  if(length(ranked.mat[ranked.mat > n_nonzero]) != 0) {
    stop("No rank may be larger than the total amount of non-zero values in ranked matrix.")
  }
  
  df <- as.data.frame(ranked.mat)
  qmat <- as.matrix(mapply(qscore, df, MoreArgs = list(n_nonzero)))
  rownames(qmat) <- rownames(orig.mat)
  
  # No rank value means no confidence score
  qmat[is.infinite(qmat)] <- 0.0
  qmat[is.na(qmat)] <- 0.0
  qmat[ix.zero] <- 0.0

  return(qmat)
}