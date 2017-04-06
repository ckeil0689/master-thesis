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
  
  k_genes <- c("A", "B")
  c_genes <- c("A", "B", "C")
  r_genes <- c("C",  "E")
  i_genes <- c("B", "J")
  
  test_tfs <- toupper(sort(c("maf", "stat3", "batf")))
  
  # Test matrices
  k_values <- c(0.992, 0.346, # A, B
                0.283, 0.625, # ...
                0.526, 0.987)
  c_values <- c(0.815, 0.245, 0.555, # A B C
                0.871, 0.658, 0.514, # ...
                0.374, 0.595, 0.119)
  r_values <- c(0.000, 0.949, # C E
                0.984, 0.566, # ...
                0.254, 0.000)
  i_values <- c(0.563, 0.000, # B J
                0.722, 0.312, # ...
                0.645, 0.000)
  
  kqmat <- matrix(k_values, nrow = length(k_genes), ncol = length(test_tfs), dimnames = list(k_genes, test_tfs))
  cqmat <- matrix(c_values, nrow = length(c_genes), ncol = length(test_tfs), dimnames = list(c_genes, test_tfs))
  rqmat <- matrix(r_values, nrow = length(r_genes), ncol = length(test_tfs), dimnames = list(r_genes, test_tfs))
  iqmat <- matrix(i_values, nrow = length(i_genes), ncol = length(test_tfs), dimnames = list(i_genes, test_tfs))
  
  # KCRI
  print("KCRI combo test")
  genes <- toupper(sort(unique(c(k_genes, c_genes, r_genes, i_genes))))
  expected.result <- matrix(0, nrow = length(genes), ncol = length(test_tfs))
  
  combo.mat <- combine.qmats(kqmat, cqmat, rqmat, iqmat)
  expect_that(combo.mat, is_a("matrix"))
  expect_that(rownames(combo.mat), is_identical_to(genes))
  print(paste("COMBO COLNAMES:", colnames(combo.mat)))
  expect_that(colnames(combo.mat), is_identical_to(test_tfs))
  expect_that(dim(combo.mat), is_identical_to(dim(expected.result)))
  
  # KC
  print("KC combo test")
  genes <- toupper(sort(unique(c(k_genes, c_genes))))
  expected.result <- matrix(0, nrow = length(genes), ncol = length(test_tfs))
  
  combo.mat <- combine.qmats(kqmat, cqmat, NULL, NULL)
  expect_that(combo.mat, is_a("matrix"))
  expect_that(rownames(combo.mat), is_identical_to(genes))
  print(paste("COMBO COLNAMES:", colnames(combo.mat)))
  expect_that(colnames(combo.mat), is_identical_to(test_tfs))
  expect_that(dim(combo.mat), is_identical_to(dim(expected.result)))
  
})