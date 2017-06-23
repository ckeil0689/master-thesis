# !/usr/bin/env Rscript

# Testing creation of rank matrix
source(paste0(getwd(), "/../main/" , "quantileRank-fun.R"), chdir = TRUE)

context("Calculating quantile rank matrix")

test_that("Confidence score matrix is transformed to quantile rank matrix as expected", {
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
  
  expect_that(func.result <- calc.quantile.ranks(mat, positiveOnly = FALSE), gives_warning())
  
  expect_that(func.result, is_a("matrix"))
  # no NA values
  expect_that(length(func.result[is.na(func.result)]) == 0, is_true())
  # no values outside of range
  expect_that(length(func.result[func.result > 1.0 || func.result < 0.0]) == 0, is_true())
  expect_that(func.result, equals(expected.result, tolerance = 0.001))
})

# According to mmc1 paper formula
test_that("Qscore is calculated as expected", {
  
  # Easy case 1
  n_nonzero <- 10
  rank <- 1
  qs <- qscore(rank, n_nonzero)
  expected_qs <- 0.9
  expect_that(qs, equals(expected_qs, tolerance = 0.0000001))
  
  # More complex case
  n_nonzero <- 236
  rank <- 78
  qs <- qscore(rank, n_nonzero)
  expected_qs <- 0.66949152542
  expect_that(qs, equals(expected_qs, tolerance = 0.0000001))
  
  # Division by zero
  n_nonzero <- 0
  rank <- 78
  qs <- qscore(rank, n_nonzero)
  expected_qs <- -Inf
  expect_that(qs, is_identical_to(expected_qs))
  
  # No rank
  n_nonzero <- 130
  rank <- 0
  qs <- qscore(rank, n_nonzero)
  expected_qs <- 1.0
  expect_that(qs, is_identical_to(expected_qs))
})