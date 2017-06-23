#!/usr/bin/env Rscript

# Tiny script to set global values. This is useful because it can be used outside of master.R, for example for testing.

GLOBAL <- list()
GLOBAL[["DEBUG"]] <- TRUE
GLOBAL[["TEST"]] <- FALSE # don't change; it is set to TRUE in run-tests.R so no print statements are made when testing other than from testthat
# at the moment, Aviv Madar's ranking method is included as an option (originally for comparison reasons when implementing the ranking from 'Computational Methods')
GLOBAL[["use.nyu.rank"]] <- FALSE
# Carried over from original Aviv Madar code / z-score cut which selects genes from DESeq data based on Th0/Th17 differential expression
GLOBAL[["z.abs.cut"]] <- 2.50
# Confidence score cut, first inferred from KC.cys example file (1.50), then altered --> optimize?
GLOBAL[["cs.abs.cut"]] <- 1.50
# Core target transcription factors
# GLOBAL[["CORE_TFS"]] <- c("batf", "irf4", "stat3", "maf", "rorc")
GLOBAL[["CORE_TFS"]] <- c("batf", "irf4", "stat3", "maf", "rorc", "fosl2", "hif1a")
