#!/usr/bin/env Rscript

# Analysis of Network Generation Output
scriptdir <- getwd()
analysis.dir <- paste0(getwd(), "/../../suppl/data/analysis/")
outpath.debug <- paste0(analysis.dir, "debug")
outpath.cyt <- paste0(analysis.dir, "cyt")

if(!dir.exists(outpath.cyt) || !dir.exists(outpath.debug)) {
  print("No analysis or debug output folders found.")
  print(paste("Folders expected in:"), analysis.dir)
  stop("Stopping because no files available for analysis.")
}

print("Beginning analysis of latest files.")
# Order all output by date
details = file.info(list.files(path = outpath.cyt, pattern = "*.csv", full.names=TRUE))
details = details[with(details, order(as.POSIXct(mtime))), ]
files = rownames(details)
latest.files <- tail(files, 3) # 3 because activator, repressor, single

print("Analyzed files:")
print(latest.files)

# Analysis loop (one run per file)
