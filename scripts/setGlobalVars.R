#!/usr/bin/env Rscript

# Tiny script to set global values. This is useful because it can be used outside of master.R, for example for testing.

GLOBAL <- list()
GLOBAL[["DEBUG"]] <- TRUE
# Carried over from original Aviv Madar code
GLOBAL[["z.abs.cut"]] <- 0.00
# Core target transcription factors
# GLOBAL[["CORE_TFS"]] <- c("batf", "irf4", "stat3", "maf", "rorc")
GLOBAL[["CORE_TFS"]] <- c("batf", "irf4", "stat3", "maf", "rorc", "fosl2", "hif1a")