# Filter the combined matrices by zscores from mmc5
apply.zscore.filter <- function(combo.mat, genes.final, combo, type) {
  filtered.genes.idx <- which(rownames(combo.mat) %in% genes.final)
  if(length(filtered.genes.idx) == 0) stop("When filtering genes by differential expression z-score, no genes were left.")
  filtered.genes <- rownames(combo.mat)[filtered.genes.idx]
  combo.mat.filtered <- combo.mat[filtered.genes.idx,]
  rownames(combo.mat.filtered) <- filtered.genes
  if(GLOBAL[["DEBUG"]]) write.mat(combo.mat.filtered, outpath.debug, combo, paste0("mat_", type))
  return(combo.mat.filtered)
}

apply.sign.mat <- function(combo.mat.filtered, ko.scores, combo, type) {
  println("Applying signs to matrix.")
  # Empty matrix with same dimension as combined matrix
  mat.sign <- matrix(0, nc=ncol(combo.mat.filtered), nr=nrow(combo.mat.filtered), dimnames=dimnames(combo.mat.filtered))
  # Only set values which also appear in KO matrix (TF-target gene pairs)
  ko.genes <- rownames(ko.scores)[which(rownames(ko.scores) %in% rownames(combo.mat.filtered))]
  mat.sign[ko.genes, colnames(ko.scores)] <- ko.scores[ko.genes,]
  # The knockout values will give us signs, everything else treated as positive (ChIP!)
  mat.sign <- sign(mat.sign)
  mat.sign[which(mat.sign==0)] <- 1
  if(GLOBAL[["DEBUG"]]) write.mat(mat.sign, outpath.debug, combo, paste0("_", type, "_signmat"))
  
  if(!identical(dim(mat.sign), dim(combo.mat.filtered))) {
    stop("Sign matrix does not have the same dimension as combined matrix, things will break. Stopping.")
  }
  
  factor <- 1
  if(type == "repressor") {
    factor <- -1
  }
  
  # Element-wise multiplication with sign matrix
  combo.mat.signed <- combo.mat.filtered * as.vector(mat.sign * factor)
  if(GLOBAL[["DEBUG"]])  write.mat(combo.mat.signed, outpath.debug, combo, paste0("_signed_", type))
  return(combo.mat.signed)
}

# Combines the supplied Q-matrices to a single matrix and applies the determined sign matrix to it. The resulting matrix will be returned.
# Q-matrices which are passed as NULL will not be considered when combining the matrices, thus allowing for different combinations.
# Arguments
# combo - a String describing the data type combination (e.g. kc --> knockout and chip)
# type - The type (activator or repressor)
# *.qmat - The q-matrices (ranked and adjusted data to fit in range of [0-1] per data type)
# genes.final - The genes to be included in the combined matrix (rows). These could, for example, have been filtered by z-scores.
createCombinedMat <- function(combo, type, ko.qmat, chip.qmat, rna.qmat, immgen.qmat, genes.final, ko.scores) {
  println("Combining Q-matrices to a single matrix.")
  if((length(ko.scores[is.na(ko.scores)]) + length(ko.scores[is.infinite(ko.scores)])) > 0) {
    stop("NA or Inf values found in knockout score matrix. Stopping.")
  }
  combo.mat <- combine.qmats(ko.qmat, chip.qmat, rna.qmat, immgen.qmat)
  combo.mat.filtered <- apply.zscore.filter(combo.mat, genes.final, combo, type)
  combo.mat.signed <- apply.sign.mat(combo.mat.filtered, ko.scores, combo, type)
  return(combo.mat.signed)
}