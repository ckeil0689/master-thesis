#!/usr/bin/env Rscript

# Testing utility functions in util.R
# util.R should be located in dir above /test/
source(paste0(getwd(), "/../" , "util.R"))

context("Write.mat")
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
  filename <- "tmp_matrix.txt"
  
  # make sure we are not verifying the creation of a file that already exists
  expect_that(file.exists(filename), is_false())
  # now write that file
  expect_that(write.mat(mat, outpath, prefix, suffix), prints_text(paste("Writing matrix to file:", paste0(outpath, filename))))
  # was it written?
  expect_that(file.exists(filename), is_true())
  
  # load it again
  loaded_mat <- as.matrix(read.table(filename, header = TRUE, sep = "\t", row.names = 1))
  
  # compare
  expect_that(colnames(loaded_mat), equals(cnms))
  expect_that(rownames(loaded_mat), equals(rnms))
  expect_that(loaded_mat, equals(mat))
  
  # test some problematic input
  expect_that(write.mat(outpath, prefix, suffix), throws_error())
  expect_that(write.mat(mat, prefix, suffix), throws_error())
  expect_that(write.mat(mat, outpath, prefix), prints_text(paste("Writing matrix to file:", paste0(outpath, prefix, ".txt"))))
  expect_that(write.mat(), throws_error())
  expect_that(write.mat(outpath, prefix, suffix, mat), throws_error())
  
  # delete tmp file
  file.remove(filename)
  expect_that(file.exists(filename), is_false())
})

context("Loading zscores")
# Loading z-score table mmc5
test_that("zscores are successfully loaded", {
  zscores.path <- paste0(getwd(), "/../../suppl/mmc5.xls" )
  expect_that(file.exists(zscores.path), is_true())
  
  # Load the scores using the function
  zscores.all <- load.zscores(zscores.path)
  
  # Check if the scores are as expected
  score.col <- "Th17_vs_Th0_Zscores"
  expect_that(zscores.all, is_a("matrix"))
  expect_that(colnames(zscores.all), matches(c(score.col)))
  expect_that(is.numeric(zscores.all[score.col]), is_true())
})