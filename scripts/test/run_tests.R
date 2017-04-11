#!/usr/bin/env Rscript

# Test script for unit testing the single functions & components used in the project
source(paste0(getwd(), "/../addLibraries.R"))

quiet <- TRUE
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1 && args[1] == "verbose") {
  quiet <- FALSE
}

# Test directory paths
test.dir <- getwd()
results.out <- paste0(test.dir, "/results/")
if(!dir.exists(results.out)) {
  dir.create(results.out)
}

# Run the tests
source(paste0(getwd(), "/../setGlobalVars.R"))
GLOBAL[["DEBUG"]] <- FALSE # turn off for testing if it is enabled in setGlobalVars.R
GLOBAL[["TEST"]] <- quiet # if TRUE, no output from functions will be printed in addition to testthat output
source(paste0(getwd(), "/../util.R")) # for write.mat() & println()
test_dir(test.dir, reporter="summary")

print("Done testing.")