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

# 1) Testing 
# source("path/to/fibo.R")
test_results <- test_dir(results.out, reporter="summary")