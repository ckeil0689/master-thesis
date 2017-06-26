#!/usr/bin/env Rscript
library(ggplot2)
library(gplots)
library(plyr)
library(igraph)

# Analysis of Network Generation Output
scriptdir <- getwd()
analysis.dir <- paste0(getwd(), "/../../suppl/data/analysis/")
outpath.debug <- paste0(analysis.dir, "debug")
outpath.cyt <- paste0(analysis.dir, "cyt")
outpath.results <- paste0(scriptdir, "/", "results/")

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

# Retrieve latest output files (by POSIX modification time)
details = file.info(list.files(path = outpath.cyt, pattern = "*.csv", full.names=TRUE))
details = details[with(details, order(as.POSIXct(mtime))), ]
files = rownames(details)
latest.files <- tail(files, 3) # 3 because activator, repressor, single

print("Analyzed files:")
print(basename(latest.files))

gen.kc.net <- NULL
# Confidence score density histograms for each file
for(i in latest.files) {
  data <- read.csv(i, stringsAsFactors = FALSE)
  type <- gsub(paste0(outpath.cyt, "/kc_(activator|repressor|single)_.*csv"),'\\1', i)
  if(type == "single") {
    gen.kc.net <- data
  }
  scores <- data[,"confidence_score"]
  scores.hist <- hist(abs(scores), main = paste("KC score density distribution for", type), 
                      freq = FALSE, col = "blue", xlab = "abs(confidence_score)", ylim = range(c(0:3)))
  png(paste0(outpath.results, type, "-scores-hist.png"))
  plot(scores.hist)
  dev.off()
}

# Barplots for TF - interactions
count.ia <- count(gen.kc.net, vars=c("nodeA","interaction"))
ggplot(count.ia, aes(x = nodeA, y = freq, fill = interaction)) +   
  geom_bar(position = "dodge", stat="identity") + 
  scale_fill_manual(breaks=levels("interaction"), values=c('green', 'red')) + 
  coord_cartesian(ylim=c(0,1400)) + geom_text(aes(label=freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(text = element_text(size=30))
ggsave("all-interaction-barplot-GN.pdf", width = 16, height = 9)

# Adjacency heatmaps (latest final matrices from debug directory)
print("Making interaction heatmaps. This uses the last debug output (enable DEBUG when running master.R).")
debug.matrices <- list.files(path = outpath.debug, pattern = "kc_signed_(activator|repressor).txt", full.names=TRUE)
col_breaks = c(seq(-2,1.50,length=400),  # for red
               seq(1.51,2,length=100))
for(i in debug.matrices) {
  data <- read.csv(i, sep = "\t", row.names = 1, header = TRUE, stringsAsFactors = FALSE)
  data.m <- data.matrix(data, rownames.force = TRUE)
  type <- gsub(paste0(outpath.debug, "/kc_signed_(activator|repressor).txt"),'\\1', i)
  print(paste("Drawing heatmap for", type, "..."))
  pdf(paste0(outpath.results, "/", type, "-heatmap-GN.pdf"))
  heatmap.2(data.m, Rowv = TRUE, Colv = FALSE, dendrogram = "none", 
            main = paste("Score Heatmap for", type),
            labRow = FALSE, margins = c(8, 2), trace = "none", breaks = col_breaks)
  dev.off()
  
  print(paste("Calculating sum scores for", type, "..."))
  # Ranked sum score list
  sumscores <- rowSums(abs(data.m))
  sums.genes <- data.frame(sumscores, row.names = rownames(data.m))
  sums.genes.rankidx <- order(sums.genes[,"sumscores"], decreasing = TRUE)
  ranked.genes <- sums.genes[sums.genes.rankidx,, drop = FALSE]
  colnames(ranked.genes) <- NULL
  rownames(ranked.genes) <- toupper(rownames(ranked.genes))
  write.table(ranked.genes, paste0(outpath.results, type, "-sum-genes-GN.rnk"), quote = FALSE, sep="\t", row.names = TRUE)
}

print("Doing analysis for Ciofani net.")
# --------------------------------------
# Ciofani Network Example (KC) ---------
# --------------------------------------
cio.filename <- paste0(scriptdir, "/kc-edges-formatted.csv")
cio.kc.net <- NULL
if(file.exists(cio.filename)) {
  print("Analyzing Ciofani et al. example KC-network for comparison.")
  cio.kc.net <- read.csv(cio.filename, stringsAsFactors = FALSE)
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

# Single barplot for TF interactions
count.ia <- count(cio.kc.net, vars=c("tf","interaction"))
ggplot(count.ia, aes(x = tf, y = freq, fill = interaction)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(breaks=levels("interaction"), values=c('green', 'red')) + 
  coord_cartesian(ylim=c(0,1400)) + geom_text(aes(label=freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(text = element_text(size=30))
ggsave("all-interaction-barplot-ciofani.pdf", width = 16, height = 9)

# Matrix analysis
print("Analyzing intermediate matrix files stored in debug folder.")
debug.files <- list.files(path = outpath.debug, pattern = "*.csv", full.names=TRUE)
print("Analyzed files:")
print(debug.files)
