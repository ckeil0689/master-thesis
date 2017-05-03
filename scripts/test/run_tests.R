#!/usr/bin/env Rscript

# Test script for unit testing the single functions & components used in the project
source(paste0(getwd(), "/../addLibraries.R"))

# Turn verbose testing on or leave it quiet
quiet <- TRUE
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 1 && args[1] == "verbose") {
  quiet <- FALSE
}

# Run the tests
print("Beginning tests. Some tests may take a while because they load and parse large amounts of data!")
source(paste0(getwd(), "/../setGlobalVars.R"))
GLOBAL[["DEBUG"]] <- FALSE # turn off for testing if it is enabled in setGlobalVars.R
GLOBAL[["TEST"]] <- quiet # if TRUE, no output from functions will be printed in addition to testthat output
source(paste0(getwd(), "/../utility.R")) # if needed before source'd in master.R (write.mat etc.)
test_dir(getwd(), reporter="summary")

print("Done testing.")