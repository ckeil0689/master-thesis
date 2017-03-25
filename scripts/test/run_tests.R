#!/usr/bin/env Rscript

# Test script for unit testing the single functions & components used in the project
# Perform some intial setups
list.of.packages <- c("testthat")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# install if missing
if(length(new.packages)) {
  print("Installing data.table package...")
  install.packages(new.packages, repos="http://cran.rstudio.com/")
}

library(testthat) 

# Test directory paths
test.dir <- getwd()
results.out <- paste0(test.dir, "/results/")
if(!dir.exists(results.out)) {
  dir.create(results.out)
}

# Run the tests
test_dir(test.dir, reporter="summary")

print("Done testing.")