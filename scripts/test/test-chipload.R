#!/usr/bin/env Rscript

# Testing ChIP score matrix generation
source(paste0(getwd(), "/../" , "chipExtract-fun.R"), chdir = TRUE)

CORE_TFS <- c("batf", "irf4", "stat3", "maf", "rorc", "fosl2", "hif1a")

context("Fixing TF names")

test_that("TF names are fixed to match conventional output (matching Cell authors example matrices)", {
  expect_that(fix.tf.name("I_am_TF-GAMMA"), equals("iamtfgamma"))
  expect_that(fix.tf.name("RoRG"), equals("rorc"))
  expect_that(fix.tf.name("c-mAF"), equals("maf"))
  expect_that(fix.tf.name(""), equals(""))
  expect_that(fix.tf.name("----__-"), equals(""))
})

context("Extracting information from reference table")

test_that("TF names are extracted correctly from reference table", {
  
  # Files known as TF-Thx experiment
  rorc_th0 <- "GSM1004853_SL3779_SL3778_genes.txt"
  maf_th0 <- "GSM1004798_SL4424_SL4425_genes.txt"
  batf_th17 <- "GSM1004787_SL3037_SL3036_genes.txt"
  stat3_th17 <- "GSM1004865_SL3315_SL3319_genes.txt"
  p300_th0 <- "GSM1004842_SL1948_SL1947_genes.txt"
  p300_th17 <- "GSM1004851_SL3594_SL3592_genes.txt"
  
  non_existent <- "GSM0000000_SL00000_SL00000_genes.txt"
  gibberish <- "ceverwrgwegrfr.txt"
  empty <- ""
  
  expect_that(extract.tf.from.ref(rorc_th0, boost.p300 = TRUE, CORE_TFS), is_identical_to("rorc-th0"))
  expect_that(extract.tf.from.ref(maf_th0, boost.p300 = TRUE, CORE_TFS), is_identical_to("maf-th0"))
  expect_that(extract.tf.from.ref(batf_th17, boost.p300 = TRUE, CORE_TFS), is_identical_to("batf-th17"))
  expect_that(extract.tf.from.ref(stat3_th17, boost.p300 = TRUE, CORE_TFS), is_identical_to("stat3-th17"))
  expect_that(extract.tf.from.ref(p300_th0, boost.p300 = TRUE, CORE_TFS), is_identical_to("p300-th0"))
  expect_that(extract.tf.from.ref(p300_th17, boost.p300 = TRUE, CORE_TFS), is_identical_to("p300-th17"))
  expect_that(extract.tf.from.ref(non_existent, boost.p300 = TRUE, CORE_TFS), throws_error("No library match found."))
  expect_that(extract.tf.from.ref(gibberish, boost.p300 = TRUE, CORE_TFS), throws_error("No library match found."))
  expect_that(extract.tf.from.ref(empty, boost.p300 = TRUE, CORE_TFS), throws_error("No library match found."))
  expect_that(extract.tf.from.ref(NULL, boost.p300 = TRUE, CORE_TFS), throws_error("No experiment passed. Stopping."))
})

context("Generating the ChIP-seq confidence score matrix")

test_that("Skeleton matrix is created as expected", {
  
})