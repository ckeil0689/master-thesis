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
  
})