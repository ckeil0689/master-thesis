#!/usr/bin/env Rscript

# Checks for and if needed installs required R-packages
# Perform some intial setups
list.of.packages <- c("data.table", "reshape2", "testthat", "ggplot2", "plyr", "igraph", "gplots")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# install if missing
if(length(new.packages)) {
  print("Installing data.table package...")
  tryCatch({
    install.packages(new.packages, repos="http://cran.rstudio.com/")
  }, error = function(err) {
    print(paste("ERROR:", err))
    stop("Try running the script as administrator (sudo).")
  })
}

library(data.table)
library(reshape2)
library(testthat)