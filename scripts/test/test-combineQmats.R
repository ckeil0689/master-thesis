# !/usr/bin/env Rscript

# Testing combination of Q-matrices
source(paste0(getwd(), "/../" , "combineQmats-fun.R"), chdir = TRUE)

context("Combination of Q-matrices")

test_that("Rownames/genes are combined as expected", {
  
  k_genes <- c("A", "B", "C", "D")
  c_genes <- c("B", "D", "E", "J", "K", "O")
  r_genes <- c("A", "B", "C",  "E", "J", "X", "S", "L")
  i_genes <- c("B", "C",  "E", "J", "X", "W", "Q", "H")
  
  # Test matrices
  kqmat <- matrix(0, nrow = length(k_genes), ncol = length(GLOBAL[["CORE_TFS"]]), dimnames = list(k_genes, GLOBAL[["CORE_TFS"]]))
  cqmat <- matrix(0, nrow = length(c_genes), ncol = length(GLOBAL[["CORE_TFS"]]), dimnames = list(c_genes, GLOBAL[["CORE_TFS"]]))
  rqmat <- matrix(0, nrow = length(r_genes), ncol = length(GLOBAL[["CORE_TFS"]]), dimnames = list(r_genes, GLOBAL[["CORE_TFS"]]))
  iqmat <- matrix(0, nrow = length(i_genes), ncol = length(GLOBAL[["CORE_TFS"]]), dimnames = list(i_genes, GLOBAL[["CORE_TFS"]]))
  
  # KCRI - normal
  qmatlist <- list(kqmat, cqmat, rqmat, iqmat)
  genes <- get.genes(qmatlist)
  expected.result <- sort(c("A", "B", "C", "D", "E", "J", "K", "O", "X", "S", "L", "W", "Q", "H"))
  expect_that(genes, is_a("character"))
  expect_that(genes, is_identical_to(expected.result))
  
  # KC -normal
  qmatlist <- list(kqmat, cqmat, NULL, NULL)
  genes <- get.genes(qmatlist)
  expected.result <- sort(c("A", "B", "C", "D", "E", "J", "K", "O"))
  expect_that(genes, is_a("character"))
  expect_that(genes, is_identical_to(expected.result))
  
  # Only one non-NULL
  qmatlist <- list(kqmat, NULL, NULL, NULL)
  genes <- get.genes(qmatlist)
  expected.result <- sort(c("A", "B", "C", "D"))
  expect_that(genes, is_a("character"))
  expect_that(genes, is_identical_to(expected.result))
  
  # NULL input
  qmatlist <- list(NULL, NULL, NULL, NULL)
  expect_that(get.genes(qmatlist), throws_error())
})

test_that("Q-score matrices are combined as expected", {
  
})