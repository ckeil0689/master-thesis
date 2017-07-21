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
  title <- expression(paste("KC absolute score density distribution ", "C" ["net"]))
  if(type == "single") {
    title <- expression(paste("KC absolute score density distribution ", "P" ["net"]))
  }
  scores.hist <- hist(abs(scores), main = title, 
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
            key.par=list(mar=c(3.5,0,3,0)), key.title = paste("Score", tools::toTitleCase(type)),
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
print("Creating an example Poisson distribution image.")
x <- runif(5, min=1, max=10)
lambdas <- runif(5, min=100, max=500)
# Vector of densities for ggplot
dens <- c()
for(l in lambdas) {
   pois.dist <- ppois(x, l, lower.tail = FALSE)
   dens <- c(dens, pois.dist)
}
lines <- rep(LETTERS[seq(1, 5)], each = 5)
df <- data.frame(dens, lines)
ggplot(df, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5) + theme_bw() + 
  labs(x = "x", y = "Density") + ggtitle("Poisson Cumulative Density Functions")
ggsave(paste0(outpath.results, "pois-dist-example.pdf"))
       
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

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Recreate Cytoscape Charts (HARDCODED data... I know...)<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Recreating Cytoscape charts for improved visuals.")
# 1) Cnet vs Pnet Average Clustering Coefficient
acc.vals <- c(0.521, 0.528, 0.527, 0.528, 0.526, 0.524, 0.005, 0.004, 0.001, 0.002)
neighbors.num <- c(2, 3, 4, 5, 6, 7, 360, 500, 683, 739)
net <- rep("cnet", length(acc.vals))
acc.df.cnet <- data.frame(c(acc.vals, neighbors.num, net))

acc.vals <- c(0.6880, 0.6769, 0.6664, 0.6767, 0.6666, 0.0026, 0.0015, 0.0014, 0.0017, 0.0014, 9.5846e-4)
neighbors.num <- c(2, 3, 4, 5, 6, 899, 1103, 1149, 1254, 1421, 1682)
net <- rep("pnet", length(acc.vals))
acc.df.pnet <- data.frame(c(acc.vals, neighbors.num, net))

p <- ggplot(NULL, aes(y = acc.vals, x = neighbors.num, colour = net), log10 = "x")
p + geom_point(data = acc.df.cnet) + geom_step(data = acc.df.pnet) + 
  theme_bw() + theme(text = element_text(size=20), axis.title.x = "Number of neighbors", 
                     axis.title.y = "Avg. clustering coefficient")
ggsave(paste0(outpath.results, "avg-clus-coeff-pnet-vs-cnet.pdf"))

# 2) Shortest path length distribution histogram for Cnet and Pnet
obs = c(1, 2, 3, 4)
freq = c(4237, 7354, 3371, 557)
spld.cnet.df <- data.frame(fill="blue", obs, freq)

obs = c(1, 2, 3, 4)
freq = c(7500, 10747, 2153, 720)
spld.pnet.df <- data.frame(fill="green", obs, freq)

spld.df <- rbind(spld.cnet.df, spld.pnet.df)
ggplot(spld.df, aes(x=obs, y=freq, fill=fill)) +
  geom_histogram(binwidth=1, colour="black", position="dodge") +
  scale_fill_identity()
ggsave(paste0(outpath.results, "shortest-path-dist-pnet-vs-cnet.pdf"))

# 3) Cnet vs Pnet BCV
bcv.vals <- c(1.4415e-4, 4.866e-4, 3.5041e-4, 8.0744e-4, 6.7927e-4, 3.0552e-4, 4.3588e-4)
neighbors.num <- c(360, 500, 500, 683, 732, 735, 739)
net <- rep("cnet", length(bcv.vals))
bcv.df.cnet <- data.frame(c(bcv.vals, neighbors.num, net))

bcv.vals <- c(1.3175e-4, 5.5606e-4, 1.8849e-4, 1.0849e-4, 1.8057e-4, 2.2423e-4)
neighbors.num <- c(899, 1682, 1103, 1149, 1254, 1421)
net <- rep("pnet", length(bcv.vals))
bcv.df.pnet <- data.frame(c(bcv.vals, neighbors.num, net))

p <- ggplot(NULL, aes(y = bcv.vals, x = neighbors.num, colour = net), log10 = "x")
p + geom_point(data = bcv.df.cnet) + geom_step(data = bcv.df.pnet) + 
  theme_bw() + theme(text = element_text(size=20), axis.title.x = "Number of neighbors", 
                     axis.title.y = "Betweenness centrality")
ggsave(paste0(outpath.results, "bcv-pnet-vs-cnet.pdf"))

# 4) In-degree distribution histogram for Cnet and Pnet
obs = c(1, 2, 3, 4, 5, 6, 7)
freq = c(1202, 451, 281, 166, 86, 28, 4)
indeg.cnet.df <- data.frame(fill="blue", obs, freq)

obs = c(1, 2, 3, 4, 5, 6)
freq = c(1661, 719, 504, 362, 204, 71)
indeg.pnet.df <- data.frame(fill="green", obs, freq)

indeg.df <- rbind(indeg.cnet.df, indeg.pnet.df)
ggplot(indeg.df, aes(x=obs, y=freq, fill=fill)) +
  geom_histogram(binwidth=1, colour="black", position="dodge") +
  scale_fill_identity()
ggsave(paste0(outpath.results, "indeg-dist-pnet-vs-cnet.pdf"))