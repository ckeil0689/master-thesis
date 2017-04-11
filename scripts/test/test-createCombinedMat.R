# !/usr/bin/env Rscript

# Testing the application of the z-score filter and the sign-matrix from DESeq differential expression data
source(paste0(getwd(), "/../" , "createCombinedMat-fun.R"), chdir = TRUE)

context("Generation of combined matrix")

test_that("Rownames/genes are filtered by z-scores as expected", {
  
  combo.genes <- c("A", "B", "C", "D", "E", "J", "K", "O", "X", "S", "L", "W", "Q", "H")
  combo.tfs <- c("BATF", "MAF", "STAT3")
  combo.mat <- matrix(0, nrow=length(combo.genes), ncol=length(combo.tfs), dimnames=list(combo.genes, combo.tfs))
  
  # Matrix dimensions
  genes.final <- c("A", "B", "X", "S", "L", "W", "Q", "H")
  expected.result <- matrix(0, nrow=length(genes.final), ncol=length(combo.tfs), dimnames=list(genes.final, combo.tfs))
  
  filtered.mat <- apply.zscore.filter(combo.mat, genes.final, "testcombo", "testtype")
  expect_that(filtered.mat, is_a("matrix"))
  expect_that(filtered.mat, is_identical_to(expected.result))
  
  combo.genes <- c("A", "B", "C")
  combo.tfs <- c("BATF", "MAF")
  combo.vals <- c(0.000, 0.949, 0.984, 
                  0.566, 0.254, 0.000)
  combo.mat <- matrix(combo.vals, nrow=length(combo.genes), ncol=length(combo.tfs), dimnames=list(combo.genes, combo.tfs))
  
  # Selection of data values
  genes.final <- c("A", "C")
  expected.vals <- c(0.000, 0.984,
                     0.566, 0.000)
  expected.result <- matrix(expected.vals, nrow=length(genes.final), ncol=length(combo.tfs), dimnames=list(genes.final, combo.tfs))
  
  filtered.mat <- apply.zscore.filter(combo.mat, genes.final, "testcombo", "testtype")
  expect_that(filtered.mat, is_a("matrix"))
  expect_that(filtered.mat, is_identical_to(expected.result))
  
  # Irregularities
  # Filtered genes not at all in combo matrix
  genes.final <- c("D", "F")
  expect_that(apply.zscore.filter(combo.mat, genes.final, "testcombo", "testtype"), throws_error())
  
  # Filtered genes are equal to non-filtered
  genes.final <- combo.genes
  expected.result <- matrix(combo.vals, nrow=length(combo.genes), ncol=length(combo.tfs), dimnames=list(combo.genes, combo.tfs))
  
  filtered.mat <- apply.zscore.filter(combo.mat, genes.final, "testcombo", "testtype")
  expect_that(filtered.mat, is_a("matrix"))
  expect_that(filtered.mat, is_identical_to(expected.result))
})

test_that("Sign matrix is applied as expected", {
  
  combo.genes <- c("A", "B", "C")
  combo.tfs <- c("BATF", "MAF")
  combo.vals <- c(0.000, 0.949, 0.984, 
                  0.566, 0.254, 0.000)
  combo.mat <- matrix(combo.vals, nrow=length(combo.genes), ncol=length(combo.tfs), dimnames=list(combo.genes, combo.tfs))
  
  ko.genes <- c("B", "C")
  ko.tfs <- c("BATF", "MAF")
  ko.vals <- c(0.000, -0.949, 
               -0.566, 0.040)
  ko.scores <- matrix(ko.vals, nrow=length(ko.genes), ncol=length(ko.tfs), dimnames=list(ko.genes, ko.tfs))
  
  # Normal data
  # activator
  expected.vals <- c(0.000, 0.949, -0.984, 
                     0.566, -0.254, 0.000)
  expected.result <- matrix(expected.vals, nrow=length(combo.genes), ncol=length(combo.tfs), dimnames=list(combo.genes, combo.tfs))
  
  sign.mat <- apply.sign.mat(combo.mat, ko.scores, "test", "activator")
  
  expect_that(sign.mat, is_a("matrix"))
  expect_that(sign.mat, is_identical_to(expected.result))
  
  # repressor
  expected.vals.r <- c(0.000, -0.949, 0.984, 
                       -0.566, 0.254, 0.000)
  expected.result.r <- matrix(expected.vals.r, nrow=length(combo.genes), ncol=length(combo.tfs), dimnames=list(combo.genes, combo.tfs))
  
  sign.mat.r <- apply.sign.mat(combo.mat, ko.scores, "test", "repressor")
  
  expect_that(sign.mat.r, is_a("matrix"))
  expect_that(sign.mat.r, is_identical_to(expected.result.r))
  
  # ko.genes not found in combo matrix --> no signs applied
  rownames(ko.scores) <- c("D", "E")
  
  expected.vals <- c(0.000, 0.949, 0.984, 
                     0.566, 0.254, 0.000)
  expected.result <- matrix(expected.vals, nrow=length(combo.genes), ncol=length(combo.tfs), dimnames=list(combo.genes, combo.tfs))
  
  sign.mat <- apply.sign.mat(combo.mat, ko.scores, "test", "activator")
  
  expect_that(sign.mat, is_a("matrix"))
  expect_that(sign.mat, is_identical_to(expected.result))
  
})

# Integration
test_that("Complete system of combining Q-matrices works as expected", {
  
  # KC ------------------------
  chip.genes <- c("A", "B", "C")
  chip.tfs <- c("BATF", "MAF")
  chip.vals <- c(0.000, 0.949, 0.984, 
                  0.566, 0.254, 0.000)
  chip.qmat <- matrix(chip.vals, nrow=length(chip.genes), ncol=length(chip.tfs), dimnames=list(chip.genes, chip.tfs))
  
  ko.genes <- c("B", "C")
  ko.tfs <- c("BATF", "MAF")
  ko.vals <- c(0.980, 0.563, 
               0.298, 0.098)
  ko.qmat <- matrix(ko.vals, nrow=length(ko.genes), ncol=length(ko.tfs), dimnames=list(ko.genes, ko.tfs))
  
  ko.score.vals <- c(-14.897, -8.987,
                     0.765, 2.123)
  ko.scores <- matrix(ko.score.vals, nrow=length(ko.genes), ncol=length(ko.tfs), dimnames=list(ko.genes, ko.tfs))
  
  genes.final <- c("A", "C")
  combo.tfs <- c("BATF", "MAF")
  
  # activator
  expected.vals <- c(0.000, -1.547,
                     0.566, 0.098)
  expected.result <- matrix(expected.vals, nrow=length(genes.final), ncol=length(combo.tfs), dimnames=list(genes.final, combo.tfs))
  
  combo.mat <- createCombinedMat("kc", "activator", ko.qmat, chip.qmat, NULL, NULL, genes.final, ko.scores)
  
  expect_that(combo.mat, is_a("matrix"))
  expect_that(combo.mat, is_identical_to(expected.result))
  
  # repressor
  expected.vals.r <- -1*expected.vals
  expected.result.r <- matrix(expected.vals.r, nrow=length(genes.final), ncol=length(combo.tfs), dimnames=list(genes.final, combo.tfs))
  combo.mat.r <- createCombinedMat("kc", "repressor", ko.qmat, chip.qmat, NULL, NULL, genes.final, ko.scores)
  
  expect_that(combo.mat.r, is_a("matrix"))
  expect_that(combo.mat.r, is_identical_to(expected.result.r))
  
  # KCRI -----------------------
  rna.genes <- c("A", "C", "D")
  rna.tfs <- c("BATF", "MAF", "RORC", "STAT3")
  rna.vals <- c(0.893, 0.000, 0.415,
                0.935, 0.000, 0.999,
                0.734, 0.239, 0.256,
                0.000, 0.000, 0.300)
  rna.qmat <- matrix(rna.vals, nrow=length(rna.genes), ncol=length(rna.tfs), dimnames=list(rna.genes, rna.tfs))
  
  immgen.genes <- c("A", "B", "D")
  immgen.tfs <- c("BATF", "MAF", "RORC", "STAT3")
  immgen.vals <- c(0.429, 0.448, 0.000,
                   0.735, 0.087, 0.567,
                   0.846, 0.000, 0.855,
                   0.253, 0.000, 0.628)
  immgen.qmat <- matrix(immgen.vals, nrow=length(immgen.genes), ncol=length(immgen.tfs), dimnames=list(immgen.genes, immgen.tfs))
  
  genes.final <- c("A", "B", "C", "D")
  combo.tfs <- c("BATF", "MAF", "RORC", "STAT3")
  
  # activator
  expected.vals <- c(1.322, -2.377, -1.547, 0.415, #BATF
                     2.236, 0.639, 0.098, 1.566, #MAF
                     1.580, 0.000, 0.239, 1.111, #RORC
                     0.253, 0.000, 0.000, 0.928) #STAT3
  expected.result <- matrix(expected.vals, nrow=length(genes.final), ncol=length(combo.tfs), dimnames=list(genes.final, combo.tfs))
  
  combo.mat <- createCombinedMat("kcri", "activator", ko.qmat, chip.qmat, rna.qmat, immgen.qmat, genes.final, ko.scores)
  
  expect_that(combo.mat, is_a("matrix"))
  expect_that(dimnames(combo.mat),is_identical_to(dimnames(expected.result)))
  expect_that(combo.mat, equals(expected.result, tolerance = 0.00001)) # floating point calc --> not using is_identical_to()
  
  # repressor
  expected.vals.r <- -1*expected.vals
  expected.result.r <- matrix(expected.vals.r, nrow=length(genes.final), ncol=length(combo.tfs), dimnames=list(genes.final, combo.tfs))
  
  combo.mat <- createCombinedMat("kcri", "repressor", ko.qmat, chip.qmat, rna.qmat, immgen.qmat, genes.final, ko.scores)
  
  expect_that(combo.mat, is_a("matrix"))
  expect_that(dimnames(combo.mat),is_identical_to(dimnames(expected.result.r)))
  expect_that(combo.mat, equals(expected.result.r, tolerance = 0.00001)) # floating point calc --> not using is_identical_to()
})