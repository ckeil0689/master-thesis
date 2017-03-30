#!/usr/bin/env Rscript

# Testing knockout score matrix generation
source(paste0(getwd(), "/../" , "deseqExtract-fun.R"), chdir = TRUE)

context("Testing TF extraction")

test_that("TF name is correctly extracted from DESeq file", {
  
  # As expected
  batf <- "Th0.batf.wt"
  maf <- "Th17.maf.wt"
  
  # Incomplete
  stat <- "Th17.stat3."
  rorc <- ".rorc.wt"
  hifla <- ".hif1a."
  
  # Faulty 
  fosl2 <- "Th17-fosl2-wt"
  
  expect_that(extract.tf(batf), is_identical_to("BATF"))
  expect_that(extract.tf(maf), is_identical_to("MAF"))
  expect_that(extract.tf(stat), is_identical_to("STAT3"))
  expect_that(extract.tf(rorc), is_identical_to("RORC"))
  expect_that(extract.tf(hifla), is_identical_to("HIF1A"))
  expect_that(tf <- extract.tf(fosl2), gives_warning())
  expect_that(is.na(tf), is_true())
})

context("Testing setup of skeleton matrix")

test_that("DESeq files are correctly parsed and KO skeleton matrix is setup as expected", {
  
  expect_that(skel.mat <- get.skel.mat(), gives_warning()) # skips TF ikzf3 (not in CORE)
  expect_that(skel.mat, is_a("matrix"))
  expect_that(length(rownames(skel.mat)) > 0, is_true())
  expect_that(length(colnames(skel.mat)) > 0, is_true())
  
  # All values should be zero
  expect_that(length(skel.mat[skel.mat != 0]) == 0, is_true())
  
  # Column names should all be part of CORE_TFS
  expect_that(all(tolower(colnames(skel.mat)) %in% GLOBAL[["CORE_TFS"]]), is_true())
})