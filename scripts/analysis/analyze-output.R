#!/usr/bin/env Rscript
library(ggplot2)
library(gplots)
library(plyr)
library(igraph)
library(RColorBrewer)

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

# --------------------------------------
# Confidence score density histograms --
# --------------------------------------
draw.score.distr.hist <- function(el, net.name, type) {
  scores <- el[,"confidence_score"]
  scores.hist <- hist(abs(scores), main = paste("KC score density distribution for", type), 
                      freq = FALSE, col = "blue", xlab = "abs(confidence_score)", ylim = range(c(0:3)))
  png(paste0(outpath.results, type, "-scores-hist-", net.name, ".png"))
  plot(scores.hist)
  dev.off()
}

# --------------------------------------
# Top scoring targets per TF -----------
# --------------------------------------
find.top5.targets.list <- function(ia.list, net.name, type) {
  ia.list <- ia.list[, c("nodeA", "nodeB", "confidence_score")]
  print(paste("Getting top scores for each TF for", type, net.name, "."))
  tfs <- unique(ia.list[,"nodeA"])
  tfs <- tfs[tfs != "HIF1"]
  top.num <- 5
  top5.targets.by.tf <- data.frame(matrix(ncol=3, nrow = top.num * length(tfs), 
                                          dimnames = list(NULL, c("tf", "target", "confidence_score"))))
  idx <- 1
  all.targets <- c()
  for(tf in tfs) {
    tf.rows <- ia.list[ia.list[,"nodeA"]==tf,] # select rows where nodeA == tf
    top5.idx <- head(order(abs(tf.rows[,"confidence_score"]), decreasing = T), n = 5)
    top5.targets <- tf.rows[top5.idx, "nodeB"]
    top5.scores <- tf.rows[top5.idx, "confidence_score"]
    top5.targets.by.tf[idx:(idx+top.num-1),"tf"] <- rep(tf, each = top.num)
    top5.targets.by.tf[idx:(idx+top.num-1),"target"] <- top5.targets
    top5.targets.by.tf[idx:(idx+top.num-1),"confidence_score"] <- top5.scores
    idx <- idx + top.num
    all.targets <- c(all.targets, top5.targets)
  }
  
  unique.targets <- unique(all.targets)
  print(paste(net.name, type, "has", length(unique.targets), "unique targets."))
  write.table(top5.targets.by.tf, paste0(outpath.results, type, "-top5-targets-by-tf-", net.name,".csv"), 
              quote = FALSE, sep=",", row.names = FALSE)
}

# --------------------------------------
# Barplots for TF interactions ---------
# --------------------------------------
draw.all.tf.ia.barplot <- function(el, net.name) {
  print(paste("Creating barplots for TF interaction types for", net.name))
  count.ia <- count(el, vars=c("nodeA","interaction"))
  ggplot(count.ia, aes(x = nodeA, y = freq, fill = interaction)) +   
    geom_bar(position = "dodge", stat="identity") + 
    scale_fill_manual(breaks=levels("interaction"), values=c('green', 'red')) + 
    coord_cartesian(ylim=c(0,1400)) + geom_text(aes(label=freq), size=15, 
                                                position=position_dodge(width=0.9), vjust=-0.25) +
    theme_bw() + theme(text = element_text(size=40), axis.title.x = element_blank())
  ggsave(paste0(outpath.results, "tf-ia-barplots-", net.name, ".pdf"), width = 16, height = 9)
}

# --------------------------------------
# Score heatmaps for VIP genes ---------
# --------------------------------------
th17.vip.targets <- c("IL21", "IL10", "IL17A", "IL23R", "IL17F", "FOXP3", 
                      "GATA3", "IFNG", "TGFB1", "IL2", "TNFA")
draw.heatmap.from.el <- function(el, net.name, type) {
  
  print(paste0("Making selected TF -> target gene interaction heatmap for ", type, " ", net.name,
               ". Using a supplied edge list..."))
  el.tfs <- sort(unique(el[,"nodeA"]))
  kc.mat.vip <- matrix(0L, nrow = length(th17.vip.targets), ncol = length(el.tfs), 
                       dimnames = list(th17.vip.targets, el.tfs))
  el.filtered.idx <- which((el[,"nodeA"] %in% el.tfs & el[,"nodeB"] %in% th17.vip.targets), arr.ind = T)
  el.filtered <- el[el.filtered.idx,]
  kc.mat.vip[el.filtered[,"nodeB"], el.filtered[,"nodeA"]] <- el.filtered[,"confidence_score"]

  filepath <- paste0(outpath.results, "/", type, "-el-heatmap-VIP-", net.name, ".pdf")
  draw.heatmap(kc.mat.vip, filepath, type)
}
draw.heatmap.from.mat <- function(kc.mat.num, net.name, type) {

  print(paste0("Making selected TF -> target gene interaction heatmap for ", type, " ", net.name,
               ". Using a supplied matrix (likely debug matrix, enable DEBUG when running master.R)..."))
  th17.vip.targets.idx <- which(rownames(kc.mat.num) %in% th17.vip.targets, arr.ind = T)
  kc.mat.num.vip <- kc.mat.num[th17.vip.targets.idx,]

  filepath <- paste0(outpath.results, "/", type, "-heatmap-VIP-", net.name, ".pdf")
  draw.heatmap(kc.mat.num.vip, filepath, type)
}
draw.heatmap <- function(vip.mat, filepath, type) {
  
  colors.heat <- rev(brewer.pal(9,"YlOrRd"))
  colors.breaks <- c(seq(0, 1, length=100),seq(1, 2, length=100))
  
  pdf(filepath)
  heatmap.2(abs(vip.mat), Rowv = TRUE, Colv = FALSE, dendrogram = "row",
            trace = "none", symkey=F, symbreaks=F, keysize=1,
            col = colors.heat,
            #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
            key.par=list(mar=c(3.5,0,3,0)), key.title = paste("Score", toupper(type)),
            cexRow = 1.7, cexCol = 1.7
            # lmat -- added 2 lattice sections (5 and 6) for padding
            # 1. Heatmap, 2. Row Dend, 3. Col Dend, 4. Key
            , lmat=rbind(c(5, 4, 3), c(2, 1, 6)), lwid = c(1, 5, 1)
  )
  dev.off()
}

# --------------------------------------
# TF sum scores for GSEA ---------------
# --------------------------------------
rank.genes.by.tf.sumscores <- function(kc.mat.num, type) {
  
  print(paste("Calculating sum scores for", type, "to be used as GSEA ranks."))
  # Ranked sum score list
  sumscores <- rowSums(abs(kc.mat.num))
  sums.genes <- data.frame(sumscores, row.names = rownames(kc.mat.num))
  sums.genes.rankidx <- order(sums.genes[,"sumscores"], decreasing = TRUE)
  ranked.genes <- sums.genes[sums.genes.rankidx,, drop = FALSE]
  colnames(ranked.genes) <- NULL
  rownames(ranked.genes) <- toupper(rownames(ranked.genes))
  write.table(ranked.genes, paste0(outpath.results, type, "-sum-genes-GN.rnk"), 
              quote = FALSE, sep="\t", row.names = TRUE)
}

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Run analysis <<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Beginning analysis of latest generated files.")
# Retrieve latest output files (by POSIX modification time)
details = file.info(list.files(path = outpath.cyt, pattern = "*.csv", full.names=TRUE))
details = details[with(details, order(as.POSIXct(mtime))), ]
files = rownames(details)
latest.files <- tail(files, 3) # 3 because activator, repressor, single
print(paste("Analyzed:", basename(latest.files)))

print("Generating confidence score density histograms...")
for(i in latest.files) {
  kc.el <- read.csv(i, stringsAsFactors = FALSE)
  type <- gsub(paste0(outpath.cyt, "/kc_(activator|repressor|single)_.*csv"),'\\1', i)

  draw.score.distr.hist(kc.el, "GN", type)
  
  if(type == "activator" || type == "repressor") {
    find.top5.targets.list(kc.el, "GN", type)
  } else {
    draw.all.tf.ia.barplot(kc.el, "GN")
  }
}

file.pattern <- "kc_signed_(activator|repressor).txt"
debug.matrices <- list.files(path = outpath.debug, pattern = file.pattern, full.names=TRUE)
for(i in debug.matrices) {
  type <- gsub(paste0(outpath.debug, "/", file.pattern),'\\1', i)
  kc.mat <- read.csv(i, sep = "\t", row.names = 1, header = TRUE, stringsAsFactors = FALSE)
  kc.mat.num <- data.matrix(kc.mat, rownames.force = TRUE)
  
  draw.heatmap.from.mat(kc.mat.num, "GN", type)
  rank.genes.by.tf.sumscores(kc.mat.num, type)
}

# --------------------------------------
# Ciofani Network Example (KC) ---------
# --------------------------------------
print("Analyzing example Ciofani net.")
# Load edge list / interaction list
cio.filename <- paste0(scriptdir, "/kc-edges-formatted.csv")
cio.kc.el <- NULL
if(file.exists(cio.filename)) {
  print("Analyzing Ciofani et al. example KC-network for comparison.")
  cio.kc.el <- read.csv(cio.filename, stringsAsFactors = FALSE)
  colnames(cio.kc.el) <- c("suid", "nodeA", "interaction", "nodeB", "confidence_score") # to make code reuse easier
} else {
  print("Could not detect example edge table file for Ciofani et al. KC network. Skipping analysis.")
}

# Single score distribution
draw.score.distr.hist(cio.kc.el, "Ciofani", "single (Ciofani et al.)")

# Split into activator and repressor
cio.act <- cio.kc.el[cio.kc.el[,"interaction"] == "positive_KC",]
cio.rep <- cio.kc.el[cio.kc.el[,"interaction"] == "negative_KC",]

# Cio VIP Heatmaps
draw.heatmap.from.el(cio.act, "Ciofani", "activator")
draw.heatmap.from.el(cio.rep, "Ciofani", "repressor")

# Top TF targets
find.top5.targets.list(cio.act, "Ciofani", "activator")
find.top5.targets.list(cio.rep, "Ciofani", "repressor")

# Interaction type barplot
draw.all.tf.ia.barplot(cio.kc.el, "Ciofani")