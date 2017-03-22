#!/usr/bin/env Rscript

# Testing utility functions in util.R
# util.R should be located in dir above /test/
source(paste0(getwd(), "/../" , "util.R"))
context("Utility functions")

# write.mat function
test_that("writ.mat writes tab-delimited matrix file", {
  # A matrix to write
  r <- 3
  c <- 3
  mat <- matrix(rbinom(r*c,1,0.5),r,c)
  # give it row and column names
  rnms <- c("a", "b", "c")
  cnms <- c("X", "Y", "Z")
  rownames(mat) <- rnms
  colnames(mat) <- cnms
  outpath <- paste0(getwd(), "/")
  prefix <- "tmp"
  suffix <- "matrix"
  
  # make sure we are not verifying the creation of a file that already exists
  expect_that(file.exists("tmp_matrix.txt"), is_false())
  # now write that file
  write.mat(mat, outpath, prefix, suffix)
  # was it written?
  expect_that(file.exists("tmp_matrix.txt"), is_true())
  
  # load it again
  loaded_mat <- as.matrix(read.table("tmp_matrix.txt", header = TRUE, sep = "\t", row.names = 1))
  
  # compare
  expect_that(colnames(loaded_mat), equals(cnms))
  expect_that(rownames(loaded_mat), equals(rnms))
  expect_that(loaded_mat, equals(mat))
  
  # delete tmp file
  file.remove("tmp_matrix.txt")
  expect_that(file.exists("tmp_matrix.txt"), is_false())
})