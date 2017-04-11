#!/usr/bin/env Rscript

# Tiny script to set global values. This is useful because it can be used outside of master.R, for example for testing.

GLOBAL <- list()
GLOBAL[["DEBUG"]] <- TRUE
GLOBAL[["TEST"]] <- FALSE # is set to TRUE in run-tests.R so no print statements are made when testing other than from testthat
GLOBAL[["use.nyu.rank"]] <- TRUE
# Carried over from original Aviv Madar code
GLOBAL[["z.abs.cut"]] <- 0.00
# Confidence score cut, first inferred from KC.cys example file (1.50), then altered --> optimize?
GLOBAL[["cs.abs.cut"]] <- 1.65
# Core target transcription factors
# GLOBAL[["CORE_TFS"]] <- c("batf", "irf4", "stat3", "maf", "rorc")
GLOBAL[["CORE_TFS"]] <- c("batf", "irf4", "stat3", "maf", "rorc", "fosl2", "hif1a")
