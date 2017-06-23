# !/usr/bin/env Rscript

# Integration test for master script all tests focus on K + C data types as R + I have been abandoned for this project!
context("Running master script")
test_that("Master script completes a full run without errors", {
  # For master.R the source() call is wrapped in a test function because it will immediately run when sourcing
  # A main function would change that but it would have to be called every time in the command line
  source(paste0(getwd(), "/../main/" , "master.R"), chdir = TRUE)
})

context("Results from master script")
test_that("Master script finishes with directory structure as expected", {
  # relative to /test/ from here
  outpath.debug <- paste0(getwd(), "/../../suppl/data/analysis/debug/")
  outpath.cyt <- paste0(getwd(), "/../../suppl/data/analysis/cyt/")
  expect_that(dir.exists(outpath.debug), is_true())
  expect_that(dir.exists(outpath.cyt), is_true())
})

test_that("Master script finishes with all files generated as expected", {
  
  # files in debug
  k.smat <- paste0(outpath.debug, "K_smat.txt")
  c.smat.activator <- paste0(outpath.debug, "C_activator_smat.txt")
  c.smat.repressor <- paste0(outpath.debug, "C_repressor_smat.txt")
  k.qmat.activator <- paste0(outpath.debug, "K_activator_qmat.txt")
  k.qmat.repressor <- paste0(outpath.debug, "K_repressor_qmat.txt")
  c.qmat.activator <- paste0(outpath.debug, "C_activator_qmat.txt")
  c.qmat.repressor <- paste0(outpath.debug, "C_repressor_qmat.txt")
  kc.activator <- paste0(outpath.debug, "kc_mat_activator.txt")
  kc.repressor <- paste0(outpath.debug, "kc_mat_repressor.txt")
  kc.activator.signmat <- paste0(outpath.debug, "kc_activator_signmat.txt")
  kc.repressor.signmat <- paste0(outpath.debug, "kc_repressor_signmat.txt")
  kc.activator.signed <- paste0(outpath.debug, "kc_signed_activator.txt")
  kc.repressor.signed <- paste0(outpath.debug, "kc_signed_repressor.txt")
  
  # files in cyt
  zscores <- paste0(outpath.cyt, "zscores.txt")
  kc.single <- paste0(outpath.cyt, "kc_single_", GLOBAL[["cs.abs.cut"]], "_cs-cut_", Sys.Date(), ".csv")
  kc.activator <- paste0(outpath.cyt, "kc_activator_", GLOBAL[["cs.abs.cut"]], "_cs-cut_", Sys.Date(), ".csv")
  kc.repressor <- paste0(outpath.cyt, "kc_repressor_", GLOBAL[["cs.abs.cut"]], "_cs-cut_", Sys.Date(), ".csv")
  
  if(GLOBAL[["DEBUG"]]) {
    # S-matrices
    expect_that(file.exists(k.smat), is_true())
    expect_that(file.exists(c.smat.activator), is_true())
    expect_that(file.exists(c.smat.repressor), is_true())
    
    # Q-matrices
    expect_that(file.exists(k.qmat.activator), is_true())
    expect_that(file.exists(k.qmat.repressor), is_true())
    expect_that(file.exists(c.qmat.activator), is_true())
    expect_that(file.exists(c.qmat.repressor), is_true())
    
    # Combined Q-mats
    expect_that(file.exists(kc.activator), is_true())
    expect_that(file.exists(kc.repressor), is_true())
    
    # Sign matrices & signed matrices
    expect_that(file.exists(kc.activator.signmat), is_true())
    expect_that(file.exists(kc.repressor.signmat), is_true())
    expect_that(file.exists(kc.activator.signed), is_true())
    expect_that(file.exists(kc.repressor.signed), is_true())
  }
  
  # Cytoscape required files
  expect_that(file.exists(zscores), is_true())
  expect_that(file.exists(kc.single), is_true())
  expect_that(file.exists(kc.activator), is_true())
  expect_that(file.exists(kc.repressor), is_true())
})