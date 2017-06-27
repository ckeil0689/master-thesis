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

find.top5.targets.list <- function(ia.list, net.name, type) {
  ia.list <- ia.list[,c("nodeA", "nodeB", "confidence_score")]
  print(paste("Getting top scores for each TF for", type, net.name, "."))
  tfs <- unique(ia.list[,"nodeA"])
  top.num <- 5
  top5.targets.by.tf <- data.frame(matrix(ncol=3, nrow = top.num * length(tfs), 
                                          dimnames = list(NULL, c("tf", "target", "confidence_score"))))
  idx <- 1
  for(tf in tfs) {
    tf.rows <- ia.list[ia.list[,"nodeA"]==tf,] # select rows where nodeA == tf
    top5.idx <- head(order(abs(tf.rows[,"confidence_score"]), decreasing = T), n = 5)
    top5.targets <- tf.rows[top5.idx, "nodeB"]
    top5.scores <- tf.rows[top5.idx, "confidence_score"]
    print(paste("Top 5 Targets for", tf, ":", top5.targets))
    top5.targets.by.tf[idx:(idx+top.num-1),"tf"] <- rep(tf, each = top.num)
    top5.targets.by.tf[idx:(idx+top.num-1),"target"] <- top5.targets
    top5.targets.by.tf[idx:(idx+top.num-1),"confidence_score"] <- top5.scores
    idx <- idx + top.num 
  }
  print(top5.targets.by.tf)
  write.table(top5.targets.by.tf, paste0(outpath.results, type, "LIST-top5-targets-by-tf-GN.txt"), 
              quote = FALSE, sep="\t", row.names = FALSE)
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
print("Generating confidence score density histograms...")
for(i in latest.files) {
  # --------------------------------------
  # Confidence score density histograms --
  # --------------------------------------
  kc.el <- read.csv(i, stringsAsFactors = FALSE)
  type <- gsub(paste0(outpath.cyt, "/kc_(activator|repressor|single)_.*csv"),'\\1', i)
  if(type == "single") {
    gen.kc.net <- kc.el
  }
  scores <- kc.el[,"confidence_score"]
  scores.hist <- hist(abs(scores), main = paste("KC score density distribution for", type), 
                      freq = FALSE, col = "blue", xlab = "abs(confidence_score)", ylim = range(c(0:3)))
  png(paste0(outpath.results, type, "-scores-hist-GN.png"))
  plot(scores.hist)
  dev.off()
  
  # --------------------------------------
  # Top scoring targets per TF ---------
  # --------------------------------------
  if(type == "activator" || type == "repressor") {
    find.top5.targets.list(kc.el, "GN", type)
  }
}

# --------------------------------------
# Barplots for TF interactions ---------
# --------------------------------------
print("Creating barplots for TF interaction types...")
count.ia <- count(gen.kc.net, vars=c("nodeA","interaction"))
ggplot(count.ia, aes(x = nodeA, y = freq, fill = interaction)) +   
  geom_bar(position = "dodge", stat="identity") + 
  scale_fill_manual(breaks=levels("interaction"), values=c('green', 'red')) + 
  coord_cartesian(ylim=c(0,1400)) + geom_text(aes(label=freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(text = element_text(size=30))
ggsave("tf-interaction-barplots-GN.pdf", width = 16, height = 9)

file.pattern = "kc_signed_(activator|repressor).txt"
debug.matrices <- list.files(path = outpath.debug, pattern = file.pattern, full.names=TRUE)
col_breaks = c(seq(-2,1.50,length=400),  # for red
               seq(1.51,2,length=100))
for(i in debug.matrices) {
  # --------------------------------------
  # Confidence score heatmaps ------------
  # --------------------------------------
  print(paste("Making interaction heatmap for", type,
              ". This uses the last debug output (enable DEBUG when running master.R)."))
  kc.mat <- read.csv(i, sep = "\t", row.names = 1, header = TRUE, stringsAsFactors = FALSE)
  kc.mat.num <- data.matrix(kc.mat, rownames.force = TRUE)
  type <- gsub(paste0(outpath.debug, "/", file.pattern),'\\1', i)
  pdf(paste0(outpath.results, "/", type, "-heatmap-GN.pdf"))
  heatmap.2(kc.mat.num, Rowv = TRUE, Colv = FALSE, dendrogram = "none", 
            main = paste("Score Heatmap for", type),
            labRow = FALSE, margins = c(8, 2), trace = "none", breaks = col_breaks)
  dev.off()
  
  # --------------------------------------
  # TF sum scores for GSEA ---------------
  # --------------------------------------
  print(paste("Calculating sum scores for", type, "to be used as GSEA ranks."))
  # Ranked sum score list
  sumscores <- rowSums(abs(kc.mat.num))
  sums.genes <- data.frame(sumscores, row.names = rownames(kc.mat.num))
  sums.genes.rankidx <- order(sums.genes[,"sumscores"], decreasing = TRUE)
  ranked.genes <- sums.genes[sums.genes.rankidx,, drop = FALSE]
  colnames(ranked.genes) <- NULL
  rownames(ranked.genes) <- toupper(rownames(ranked.genes))
  write.table(ranked.genes, paste0(outpath.results, type, "-sum-genes-GN.rnk"), quote = FALSE, sep="\t", row.names = TRUE)
  
  # --------------------------------------
  # Top scoring targets per TF ---------
  # --------------------------------------
  print(paste("Getting top scores for each TF for", type, "."))
  tfs <- colnames(kc.mat)
  top.num <- 5
  top5.targets.by.tf <- data.frame(matrix(ncol=3, nrow = top.num * length(tfs), 
                                          dimnames = list(NULL, c("tf", "target", "confidence_score"))))
  idx <- 1
  for(tf in tfs) {
    tf.col <- kc.mat[,tf]
    top5.idx <- head(order(abs(tf.col), decreasing = T), n = 5)
    top5.targets <- rownames(kc.mat)[top5.idx]
    top5.scores <- tf.col[top5.idx]
    print(paste("Top 5 Targets for", tf, ":", top5.targets))
    top5.targets.by.tf[idx:(idx+top.num-1),"tf"] <- rep(tf, each = top.num)
    top5.targets.by.tf[idx:(idx+top.num-1),"target"] <- top5.targets
    top5.targets.by.tf[idx:(idx+top.num-1),"confidence_score"] <- top5.scores
    idx <- idx + top.num 
  }
  print(top5.targets.by.tf)
  write.table(top5.targets.by.tf, paste0(outpath.results, type, "-top5-targets-by-tf-GN.txt"), 
              quote = FALSE, sep="\t", row.names = FALSE)
}

print("Doing analysis for Ciofani net.")
# --------------------------------------
# Ciofani Network Example (KC) ---------
# --------------------------------------
cio.filename <- paste0(scriptdir, "/kc-edges-formatted.csv")
cio.kc.el <- NULL
if(file.exists(cio.filename)) {
  print("Analyzing Ciofani et al. example KC-network for comparison.")
  cio.kc.el <- read.csv(cio.filename, stringsAsFactors = FALSE)
  type <- "single (Ciofani et al.)"
  scores <- cio.kc.el[,"confidence_score"]
  # Histogram
  scores.hist <- hist(abs(scores), main = paste("KC score density distribution for", type), 
                      freq = FALSE, col = "blue", xlab = "abs(confidence_score)", ylim = range(c(0:3)))
  png(paste0(outpath.results, "-scores-hist-Ciofani.png"))
  plot(scores.hist)
  dev.off()
} else {
  print("Could not detect example edge table file for Ciofani et al. KC network. Skipping analysis.")
}

# Single barplot for TF interactions
count.ia <- count(cio.kc.el, vars=c("tf","interaction"))
ggplot(count.ia, aes(x = tf, y = freq, fill = interaction)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(breaks=levels("interaction"), values=c('green', 'red')) + 
  coord_cartesian(ylim=c(0,1400)) + geom_text(aes(label=freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(text = element_text(size=30))
ggsave("all-interaction-barplot-Ciofani.pdf", width = 16, height = 9)

# Top TF targets
cio.act <- cio.kc.el[cio.kc.el[,"interaction"] == "positive_KC",]
cio.rep <- cio.kc.el[cio.kc.el[,"interaction"] == "negative_KC",]
