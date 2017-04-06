#!/usr/bin/env Rscript

# Checks for and if needed installs required R-packages
# Perform some intial setups
list.of.packages <- c("data.table", "reshape2", "testthat")#, "xlsx")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# install if missing
if(length(new.packages)) {
  print("Installing data.table package...")
  install.packages(new.packages, repos="http://cran.rstudio.com/")
}

library(data.table)
library(reshape2)
library(testthat) 
# library(xlsx) // causes issues with rJava package installs which rely on specific Java dev kit installs --> unreliable mess, very probable to break on different machines