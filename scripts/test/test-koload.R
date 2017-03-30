#!/usr/bin/env Rscript

# Testing knockout score matrix generation
source(paste0(getwd(), "/../" , "deseqExtract-fun.R"), chdir = TRUE)

context("Testing TF extraction")

test_that("TF name is correctly extracted from DESeq file", {
  
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
  expect_that(extract.tf(fosl2), is_identical_to(NA))
})