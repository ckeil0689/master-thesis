#!/usr/bin/env Rscript

# Testing utility functions in util.R
# util.R should be located in dir above /test/
source(paste0(getwd(), "/../main/" , "utility.R"), chdir = TRUE)

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
  tmp_file <- paste0(prefix, ".txt")
  expect_that(write.mat(mat, outpath, prefix), prints_text(paste("Writing matrix to file:", paste0(outpath, tmp_file))))
  expect_that(write.mat(), throws_error())
  expect_that(write.mat(outpath, prefix, suffix, mat), throws_error())
  
  # delete tmp file
  file.remove(filename)
  file.remove(tmp_file)
})

# Loading z-score table mmc5
test_that("zscores are successfully loaded", {
  zscores.path <- paste0(getwd(), "/../../suppl/mmc5.csv" )
  zscores.exist <- file.exists(zscores.path)
  expect_that(zscores.exist, is_true())
  
  if(zscores.exist) {
    # Load the scores using the function
    zscores.all <- load.zscores(zscores.path)
    
    # Check if the scores are as expected
    score.col <- "Th17_vs_Th0_Zscores"
    expect_that(zscores.all, is_a("matrix"))
    expect_that(colnames(zscores.all), matches(c(score.col)))
    expect_that(is.numeric(zscores.all[score.col]), is_true())
  }
})

test_that("Genes are correctly filtered by absolute zscore", {
  vals <- c(-9.8, -2, 5.78, 0.0, 3.46, -12.7, -1.3, 3, 2.50, -2.50, 2.49, -2.51, -2.49, 2.51)
  expected_out <- c("1", "3", "5", "6", "8", "12", "14") # indexes will be row names, so these are indexes bigger than z.abs.val
  zscores <- matrix(vals, ncol = 1, nrow = length(vals))
  rownames(zscores) <- as.character(c(1:length(vals)))
  filtered.genes <- filter.genes.by.zscore(zscores, 2.5) 
  expect_that(filtered.genes, is_identical_to(expected_out))
})

test_that("Incorrect input is handled correctly", {
  vals <- c(-9.8, -2, 5.78, 0.0, 3.46, -12.7, -1.3, 3, 2.50, -2.50, 2.49, -2.51, -2.49, 2.51)
  zscores <- matrix(vals, ncol = 1, nrow = length(vals))
  
  # character as zscore
  expect_that(filter.genes.by.zscore(zscores, "2.5"), throws_error())
  
  # undefined zscores
  expect_that(filter.genes.by.zscore(NULL, 2.5), throws_error())
  
  # non-numeric score values
  vals <- c("A", "B", "C", "D")
  zscores <- matrix(vals, ncol = 1, nrow = length(vals))
  expect_that(filter.genes.by.zscore(zscores, 2.5), throws_error())
  
  # No score cut defined
  expect_that(filter.genes.by.zscore(zscores), throws_error())
})