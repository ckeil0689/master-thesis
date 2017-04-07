# !/usr/bin/env Rscript

# Testing combination of Q-matrices
source(paste0(getwd(), "/../" , "combineQmats-fun.R"), chdir = TRUE)

context("Combination of Q-matrices")

test_that("Rownames/genes are combined as expected", {
  
  k.genes <- c("A", "B", "C", "D")
  c.genes <- c("B", "D", "E", "J", "K", "O")
  r.genes <- c("A", "B", "C",  "E", "J", "X", "S", "L")
  i.genes <- c("B", "C",  "E", "J", "X", "W", "Q", "H")
  
  # Test matrices
  kqmat <- matrix(0, nrow = length(k.genes), ncol = length(GLOBAL[["CORE_TFS"]]), dimnames = list(k.genes, GLOBAL[["CORE_TFS"]]))
  cqmat <- matrix(0, nrow = length(c.genes), ncol = length(GLOBAL[["CORE_TFS"]]), dimnames = list(c.genes, GLOBAL[["CORE_TFS"]]))
  rqmat <- matrix(0, nrow = length(r.genes), ncol = length(GLOBAL[["CORE_TFS"]]), dimnames = list(r.genes, GLOBAL[["CORE_TFS"]]))
  iqmat <- matrix(0, nrow = length(i.genes), ncol = length(GLOBAL[["CORE_TFS"]]), dimnames = list(i.genes, GLOBAL[["CORE_TFS"]]))
  
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

test_that("Colnames/TFs are combined as expected", {
  
  k.tfs <- c("A", "B", "C", "D")
  c.tfs <- c("B", "D", "E", "J", "K", "O")
  r.tfs <- c("A", "B", "C",  "E", "J", "X", "S", "L")
  i.tfs <- c("B", "C",  "E", "J", "X", "W", "Q", "H")
  
  # Test matrices - CORE TFs used as row names
  kqmat <- matrix(0, nrow = length(GLOBAL[["CORE_TFS"]]), ncol = length(k.tfs), dimnames = list(GLOBAL[["CORE_TFS"]], k.tfs))
  cqmat <- matrix(0, nrow = length(GLOBAL[["CORE_TFS"]]), ncol = length(c.tfs), dimnames = list(GLOBAL[["CORE_TFS"]], c.tfs))
  rqmat <- matrix(0, nrow = length(GLOBAL[["CORE_TFS"]]), ncol = length(r.tfs), dimnames = list(GLOBAL[["CORE_TFS"]], r.tfs))
  iqmat <- matrix(0, nrow = length(GLOBAL[["CORE_TFS"]]), ncol = length(i.tfs), dimnames = list(GLOBAL[["CORE_TFS"]], i.tfs))
  
  # KCRI - normal
  qmatlist <- list(kqmat, cqmat, rqmat, iqmat)
  tfs <- get.tfs(qmatlist)
  expected.result <- sort(c("A", "B", "C", "D", "E", "J", "K", "O", "X", "S", "L", "W", "Q", "H"))
  expect_that(tfs, is_a("character"))
  expect_that(tfs, is_identical_to(expected.result))
  
  # KC -normal
  qmatlist <- list(kqmat, cqmat, NULL, NULL)
  tfs <- get.tfs(qmatlist)
  expected.result <- sort(c("A", "B", "C", "D", "E", "J", "K", "O"))
  expect_that(tfs, is_a("character"))
  expect_that(tfs, is_identical_to(expected.result))
  
  # Only one non-NULL
  qmatlist <- list(kqmat, NULL, NULL, NULL)
  tfs <- get.tfs(qmatlist)
  expected.result <- sort(c("A", "B", "C", "D"))
  expect_that(tfs, is_a("character"))
  expect_that(tfs, is_identical_to(expected.result))
  
  # All NULL matrices input
  qmatlist <- list(NULL, NULL, NULL, NULL)
  expect_that(get.tfs(qmatlist), throws_error())
  
  # List NULL input
  qmatlist <- NULL
  expect_that(get.tfs(qmatlist), throws_error())
})

test_that("Q-score matrices are combined as expected", {
  
  k.genes <- c("A", "B")
  c.genes <- c("A", "B", "C")
  r.genes <- c("C",  "E")
  i.genes <- c("B", "J")
  
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
  
  kqmat <- matrix(k_values, nrow = length(k.genes), ncol = length(test_tfs), dimnames = list(k.genes, test_tfs))
  cqmat <- matrix(c_values, nrow = length(c.genes), ncol = length(test_tfs), dimnames = list(c.genes, test_tfs))
  rqmat <- matrix(r_values, nrow = length(r.genes), ncol = length(test_tfs), dimnames = list(r.genes, test_tfs))
  iqmat <- matrix(i_values, nrow = length(i.genes), ncol = length(test_tfs), dimnames = list(i.genes, test_tfs))
  
  # KCRI
  genes <- toupper(sort(unique(c(k.genes, c.genes, r.genes, i.genes))))
  combo.vals <- c(1.807, 1.154, 0.555, 0.949, 0.000, # BATF | A B C E J
                  1.154, 1.595, 1.498, 0.566, 0.312, # MAF | A B C E J
                  0.900, 1.582, 0.373, 0.000, 0.000) # STAT3 | A B C E J
    
  expected.result <- matrix(combo.vals, nrow = length(genes), ncol = length(test_tfs))
  
  combo.mat <- combine.qmats(kqmat, cqmat, rqmat, iqmat)
  expect_that(combo.mat, is_a("matrix"))
  expect_that(rownames(combo.mat), is_identical_to(genes))
  expect_that(colnames(combo.mat), is_identical_to(test_tfs))
  expect_that(dim(combo.mat), is_identical_to(dim(expected.result)))
  
  # KC
  genes <- toupper(sort(unique(c(k.genes, c.genes))))
  combo.vals <- c(1.807, 0.591, 0.555, # BATF | A B C
                  1.154, 1.283, 0.514, # MAF | A B C
                  0.900, 1.582, 0.119 )# STAT3 | A B C
  expected.result <- matrix(combo.vals, nrow = length(genes), ncol = length(test_tfs))
  
  combo.mat <- combine.qmats(kqmat, cqmat, NULL, NULL)
  expect_that(combo.mat, is_a("matrix"))
  expect_that(rownames(combo.mat), is_identical_to(genes))
  expect_that(colnames(combo.mat), is_identical_to(test_tfs))
  expect_that(dim(combo.mat), is_identical_to(dim(expected.result)))
  
  # K
  genes <- toupper(sort(unique(k.genes)))
  combo.vals <- c(0.992, 0.346, # BATF | A B
                  0.283, 0.625, # MAF | A B
                  0.526, 0.987) # STAT3 | A B
  expected.result <- matrix(combo.vals, nrow = length(genes), ncol = length(test_tfs))
  
  combo.mat <- combine.qmats(kqmat, NULL, NULL, NULL)
  expect_that(combo.mat, is_a("matrix"))
  expect_that(rownames(combo.mat), is_identical_to(genes))
  expect_that(colnames(combo.mat), is_identical_to(test_tfs))
  expect_that(dim(combo.mat), is_identical_to(dim(expected.result)))
  
  # All NULL
  expect_that(combine.qmats(NULL, NULL, NULL, NULL), throws_error())
})