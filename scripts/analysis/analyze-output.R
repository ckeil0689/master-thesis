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

# Results output directory
outpath.results <- paste0(scriptdir, "/results/")
if(!dir.exists(outpath.results)) {
  dir.create(outpath.results)
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
for(i in latest.files) {
  data <- read.csv(i)
  type <- gsub(paste0(outpath.cyt, "/kc_(activator|repressor|single)_.*csv"),'\\1', i)
  scores <- data[,"confidence_score"]
  # Histogram
  scores.hist <- hist(abs(scores), main = paste("KC score density distribution for", type), 
                      freq = FALSE, col = "blue", xlab = "abs(confidence_score)", ylim = range(c(0:3)))
  png(paste0(outpath.results, type, "-scores-hist.png"))
  plot(scores.hist)
  dev.off()
}

# Ciofani Network Example (KC)
cio.filename <- paste0(scriptdir, "/kc-edges-formatted.csv")
if(file.exists(cio.filename)) {
  print("Analyzing Ciofani et al. example KC-network for comparison.")
  cio.kc.net <- read.csv(cio.filename)
  
  type <- "single (Ciofani et al.)"
  scores <- cio.kc.net[,"confidence_score"]
  # Histogram
  scores.hist <- hist(abs(scores), main = paste("KC score density distribution for", type), 
                      freq = FALSE, col = "blue", xlab = "abs(confidence_score)", ylim = range(c(0:3)))
  png(paste0(outpath.results, "ciofani-scores-hist.png"))
  plot(scores.hist)
  dev.off()
} else {
  print("Could not detect example edge table file for Ciofani et al. KC network. Skipping analysis.")
}

# Matrix analysis
print("Analyzing intermediate matrix files stored in debug folder.")
debug.files <- list.files(path = outpath.debug, pattern = "*.csv", full.names=TRUE)
print("Analyzed files:")
print(debug.files)
