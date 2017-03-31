# !/usr/bin/env Rscript

# Testing creation of rank matrix
source(paste0(getwd(), "/../" , "quantileRank-fun.R"), chdir = TRUE)

context("Rank matrix creation")
test_that("Confidence score matrix is transformed to rank matrix as expected", {
  mat <- matrix(0, nrow = 3, ncol = 3, dimnames = list(c("g1", "g2", "g3"), c("BATF", "MAF", "RORC")))
  mat["g1",] <- c(1.44583, NA, -13.98272)
  mat["g2",] <- c(8.000, 21.500, -5.22111)
  mat["g3",] <- c(0, -4.56320, 3.56708)
  
  # Ranking is done on absolute values and in descending order
  expected.result <- matrix(0, nrow = 3, ncol = 3, dimnames = list(c("g1", "g2", "g3"), c("BATF", "MAF", "RORC")))
  expected.result["g1",] <- c(7, 0, 2)
  expected.result["g2",] <- c(3, 1, 4)
  expected.result["g3",] <- c(0, 5, 6)
  
  rank.result <- rank.smat(mat)
  
  expect_that(rank.result, is_a("matrix"))
  # no NA values
  expect_that(length(rank.result[is.na(rank.result)]) == 0, is_true())
  # no values outside of range
  # expect_that(length(rank.result[rank.result > 1.0 || rank.result < 0.0]) == 0, is_true())
  expect_that(rank.result, is_identical_to(expected.result))
})