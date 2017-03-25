#!/usr/bin/env Rscript

# Testing ChIP score matrix generation
source(paste0(getwd(), "/../" , "chipExtract-fun.R"))

context("Fixing TF names")
test_that("TF names are fixed to match conventional output (matching Cell authors example matrices)", {
  
  tf1 <- "I_am_TF-GAMMA"
  tf2 <- "RoRG"
  tf3 <- "c-mAF"
  
  expect_that(fix.tf.name(tf1), equals("iamtfgamma"))
  expect_that(fix.tf.name(tf2), equals("rorc"))
  expect_that(fix.tf.name(tf3), equals("maf"))
})