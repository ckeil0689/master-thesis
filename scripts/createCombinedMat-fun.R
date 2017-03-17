# Combines the supplied Q-matrices to a single matrix and applies the determined sign matrix to it. The resulting matrix will be returned.
# Q-matrices which are passed as NULL will not be considered when combining the matrices, thus allowing for different combinations.
# Arguments
# combo - a String describing the data type combination (e.g. kc --> knockout and chip)
# type - The type (activator or repressor)
# *_qmat - The q-matrices (ranked and adjusted data to fit in range of [0-1] per data type)
# genes.final - The genes to be included in the combined matrix (rows). These could, for example, have been filtered by z-scores.
# CORE_TFS - The transcription factors to be included. RNA-seq compendium and Immgen could be up to 28 TFs.
createCombinedMat <- function(combo, type, ko_qmat, chip_qmat, rna_qmat, immgen_qmat, genes.final, CORE_TFS) {
  source(paste0(getwd(), "/" , "combineQmats-fun.R"))
  print("Combining Q-matrices to a single matrix.")
  combined_mat <- combine.qmats(ko_qmat, chip_qmat, rna_qmat, immgen_qmat, CORE_TFS)
  
  # Filter the combined matrices by zscores from mmc5
  filtered.genes.idx <- which(rownames(combined_mat) %in% genes.final)
  filtered.genes <- rownames(combined_mat)[filtered.genes.idx]
  combined_mat <- combined_mat[filtered.genes.idx,]
  rownames(combined_mat) <- filtered.genes
  if(GLOBAL[["DEBUG"]]) write.mat(combined_mat, outpath.debug, combo, "_mat_", type)
  
  # --------------
  # Apply sign matrix
  # --------------
  print("Applying signs to matrix.")
  # Empty matrix with same dimension as combined matrix
  mat.sign <- matrix(0, nc=ncol(combined_mat), nr=nrow(combined_mat), dimnames=dimnames(combined_mat))
  # Only set values which also appear in KO matrix (TF-target gene pairs)
  ko.genes <- rownames(ko.scores)[which(filtered.genes %in% rownames(ko.scores))]
  mat.sign[rownames(ko.genes), colnames(ko.scores)] <- ko.scores[ko.genes,]
  # The knockout values will give us signs, everything else treated as positive (ChIP!)
  mat.sign <- sign(mat.sign)
  mat.sign[which(mat.sign==0)] <- 1
  if(GLOBAL[["DEBUG"]]) write.mat(mat.sign, outpath.debug, combo, "_", type, "_signmat")
  
  print("Checking dimensions...")
  if(!identical(dim(mat.sign), dim(combined_mat))) {
    print(paste("Dimension sign_mat:", dim(mat.sign)))
    print(paste("Dimension combined_mat:", dim(combined_mat)))
    print(paste("Dimension ko_qmat.activator:", dim(ko_qmat)))
    print(paste("Dimension chip_qmat:", dim(chip_qmat)))
    stop("Sign matrix does not have the same dimension as combined matrix, things will break. Stopping.")
  }
  
  factor <- 1
  if(type == "repressor") {
    factor <- -1
  }
  
  print("Applying sign matrix to combined activator matrix...")
  # Element-wise multiplication with sign matrix
  combined_mat <- combined_mat * as.vector(mat.sign * factor)
  
  if(GLOBAL[["DEBUG"]])  write.mat(combined_mat, outpath.debug, combo, "_signed_", type)
}