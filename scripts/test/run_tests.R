#!/usr/bin/env Rscript

# Test script for unit testing the single functions & components used in the project
source(paste0(getwd(), "/../addLibraries.R"))

# Test directory paths
test.dir <- getwd()
results.out <- paste0(test.dir, "/results/")
if(!dir.exists(results.out)) {
  dir.create(results.out)
}

# Run the tests
source(paste0(getwd(), "/../setGlobalVars.R"))
test_dir(test.dir, reporter="summary")

print("Done testing.")