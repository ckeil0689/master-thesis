# !/usr/bin/env Rscript

# Testing creation of rank matrix
source(paste0(getwd(), "/../" , "quantileRank-fun.R"), chdir = TRUE)

context("Rank matrix creation")
test_that("Confidence score matrix is transformed to rank matrix as expected", {
  # Typical matrix input
  mat <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("g1", "g2"), c("BATF", "MAF")))
  mat["g1",] <- c(1.44583, NA)
  mat["g2",] <- c(8.000, 21.500)
  
  # Ranking is done on absolute values and in descending order
  expected.result <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("g1", "g2"), c("BATF", "MAF")))
  expected.result["g1",] <- c(0.000, 0)
  expected.result["g2",] <- c(0.333, 0.666)
  
  func.result <- calc.quantile.ranks(mat, positiveOnly = FALSE)
  
  expect_that(func.result, is_a("matrix"))
  # no NA values
  expect_that(length(func.result[is.na(func.result)]) == 0, is_true())
  # no values outside of range
  expect_that(length(func.result[func.result > 1.0 || func.result < 0.0]) == 0, is_true())
  expect_that(func.result, equals(expected.result, tolerance = 0.001))
  
  # Bad input
  mat["g1",] <- c("hello", NA)
  mat["g2",] <- c(NA, NA)
  
  expected.result["g1",] <- c(0, 0)
  expected.result["g2",] <- c(0, 0)
  
  func.result <- calc.quantile.ranks(mat, positiveOnly = FALSE)
  
  expect_that(func.result, is_a("matrix"))
  # no NA values
  expect_that(length(func.result[is.na(func.result)]) == 0, is_true())
  # no values outside of range
  expect_that(length(func.result[func.result > 1.0 || func.result < 0.0]) == 0, is_true())
  expect_that(func.result, equals(expected.result, tolerance = 0.001))
})