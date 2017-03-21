#!/usr/bin/env Rscript

# Test script for unit testing the single functions & components used in the project
library(testthat) 

# Test directory paths
test.dir <- getwd()
results.out <- paste0(test.dir, "/results/")

# 1) Testing 
# source("path/to/fibo.R")
test_results <- test_dir(results.out, reporter="summary")