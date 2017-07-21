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
cnet.color <- "#F8766D"
pnet.color <- "#00BFC4"
znet.color <- "#0e85f2"
mafnet.color <- "#0ef2e1"
chart.textsize <- 16
# 1) Cnet vs Pnet Average Clustering Coefficient
print("ACC distribution...")
acc.vals <- c(0.521, 0.528, 0.527, 0.528, 0.526, 0.524, 0.005, 0.004, 0.001, 0.002)
neighbors.num <- c(2, 3, 4, 5, 6, 7, 360, 500, 683, 739)
acc.df.cnet <- data.frame(acc.vals, neighbors.num)
name.cnet <- "C [net]"

acc.vals <- c(0.6880, 0.6769, 0.6664, 0.6767, 0.6666, 0.0026, 0.0015, 0.0014, 0.0017, 0.0014, 9.5846e-4)
neighbors.num <- c(2, 3, 4, 5, 6, 899, 1103, 1149, 1254, 1421, 1682)
net.pnet <- rep("pnet", length(acc.vals))
acc.df.pnet <- data.frame(acc.vals, neighbors.num)
name.pnet <- "P [net]"

ggplot() + 
  geom_point(data = acc.df.cnet, aes(x = neighbors.num, y = acc.vals, colour = name.cnet)) + 
  geom_point(data = acc.df.pnet, aes(x = neighbors.num, y = acc.vals, colour = name.pnet)) + 
  theme_bw() + coord_trans(x="log10") + labs(x="Number of Neighbors", y="Avg. Clustering Coefficient", 
                                             title="Avg. Clustering Coefficient Distribution", colour="Network") +
  scale_colour_discrete(labels = function(x) parse(text=x)) +
  theme(text = element_text(size=chart.textsize))
ggsave(paste0(outpath.results, "acc-pnet-vs-cnet.pdf"))

# 2) Shortest path length distribution histogram for Cnet and Pnet
print("SPL distribution...")
obs = c(1, 2, 3, 4)
freq = c(4237, 7354, 3371, 557)
spld.cnet.df <- data.frame(fill=cnet.color, obs, freq)

obs = c(1, 2, 3, 4)
freq = c(7500, 10747, 2153, 720)
spld.pnet.df <- data.frame(fill=pnet.color, obs, freq)

spld.df <- rbind(spld.cnet.df, spld.pnet.df)
ggplot(spld.df, aes(x=obs, y=freq, fill=fill)) +
  geom_histogram(position="dodge", stat = "identity") +
  scale_fill_identity(guide = "legend", name = "Network", 
                      labels = c(expression("C"["net"]), expression("P"["net"]))) + theme_bw() + 
  labs(x="Path Length", y="Frequency", title="Shortest Path Length Distribution", colour="Network") +
  theme(text = element_text(size=chart.textsize))
ggsave(paste0(outpath.results, "spl-pnet-vs-cnet.pdf"))

# 3) Cnet vs Pnet BCV
print("BCV distribution...")
bcv.vals <- c(1.4415e-4, 4.866e-4, 3.5041e-4, 8.0744e-4, 6.7927e-4, 3.0552e-4, 4.3588e-4)
neighbors.num <- c(360, 500, 500, 683, 732, 735, 739)
bcv.df.cnet <- data.frame(bcv.vals, neighbors.num)

bcv.vals <- c(1.3175e-4, 5.5606e-4, 1.8849e-4, 1.0849e-4, 1.8057e-4, 2.2423e-4)
neighbors.num <- c(899, 1682, 1103, 1149, 1254, 1421)
bcv.df.pnet <- data.frame(bcv.vals, neighbors.num)

ggplot() + 
  geom_point(data = bcv.df.cnet, aes(y = bcv.vals, x = neighbors.num, color = name.cnet)) + 
  geom_point(data = bcv.df.pnet, aes(y = bcv.vals, x = neighbors.num, color = name.pnet)) + 
  theme_bw() + labs(x="Number of Neighbors", y="Betweenness Centrality", 
                    title="Betweenness Centrality", colour="Network") +
  scale_colour_discrete(labels = function(x) parse(text=x)) +
  theme(text = element_text(size=chart.textsize))
ggsave(paste0(outpath.results, "bcv-pnet-vs-cnet.pdf"))

# 4) In-degree distribution histogram for Cnet and Pnet
print("In-degree distribution...")
obs = c(1, 2, 3, 4, 5, 6, 7)
freq = c(1202, 451, 281, 166, 86, 28, 4)
indeg.cnet.df <- data.frame(fill=cnet.color, obs, freq)

obs = c(1, 2, 3, 4, 5, 6)
freq = c(1661, 719, 504, 362, 204, 71)
indeg.pnet.df <- data.frame(fill=pnet.color, obs, freq)

indeg.df <- rbind(indeg.cnet.df, indeg.pnet.df)
ggplot(indeg.df, aes(x=obs, y=freq, fill=fill)) +
  geom_histogram(position="dodge", stat = 'identity') +
  scale_fill_identity(guide = "legend", name = "Network", 
                      labels = c(expression("C"["net"]), expression("P"["net"]))) + theme_bw() +
  labs(x="In-degree", y="Number of Nodes", title="In-degree Distribution", colour="Network") +
  theme(text = element_text(size=chart.textsize))
ggsave(paste0(outpath.results, "indeg-dist-pnet-vs-cnet.pdf"))

# 5) In-degree distribution histogram for Pnet @ zscore = 2.50
print("In-degree distribution (zscore = 2.50)...")
obs = c(1, 2, 3, 4, 5, 6)
freq = c(306, 212, 197, 191, 127, 57)
indeg.znet.df <- data.frame(fill=znet.color, obs, freq)

obs = c(1, 2, 3, 4, 5, 6)
freq = c(1661, 719, 504, 362, 204, 71)
indeg.pnet.df <- data.frame(fill=pnet.color, obs, freq)

indeg.df <- rbind(indeg.pnet.df, indeg.znet.df)
ggplot(indeg.df, aes(x=obs, y=freq, fill=fill)) +
  geom_histogram(position="dodge", stat = 'identity') +
  scale_fill_identity(guide = "legend", name = "Network", 
                      labels = c(expression("P"["net"]), expression("P"["zcut"]))) + theme_bw() +
  labs(x="In-degree", y="Number of Nodes", title="In-degree Distribution", colour="Network") +
  theme(text = element_text(size=chart.textsize))
ggsave(paste0(outpath.results, "indeg-dist-pnet-vs-znet.pdf"))

# 6) Pnet vs PMAF Average Clustering Coefficient
print("ACC distribution (PMAF)...")
acc.vals <- c(0.0027, 0.6626, 0.5977, 0.5870, 0.5870, 0.6000, 0.0018, 0.0015, 4.9187e-4, 0.0018, 0.0016)
neighbors.num <- c(802, 2, 3, 4, 5, 6, 903, 1257, 2042, 1115, 974)
acc.df.mafnet <- data.frame(acc.vals, neighbors.num)
name.mafnet <- "P [MAF]"

acc.vals <- c(0.6880, 0.6769, 0.6664, 0.6767, 0.6666, 0.0026, 0.0015, 0.0014, 0.0017, 0.0014, 9.5846e-4)
neighbors.num <- c(2, 3, 4, 5, 6, 899, 1103, 1149, 1254, 1421, 1682)
net.pnet <- rep("pnet", length(acc.vals))
acc.df.pnet <- data.frame(acc.vals, neighbors.num)
name.pnet <- "P [net]"

ggplot() + 
  geom_point(data = acc.df.mafnet, aes(x = neighbors.num, y = acc.vals, colour = name.mafnet)) + 
  geom_point(data = acc.df.pnet, aes(x = neighbors.num, y = acc.vals, colour = name.pnet)) + 
  theme_bw() + coord_trans(x="log10") + labs(x="Number of Neighbors", y="Avg. Clustering Coefficient", 
                                             title="Avg. Clustering Coefficient Distribution", colour="Network") +
  scale_colour_discrete(labels = function(x) parse(text=x)) +
  theme(text = element_text(size=chart.textsize))
ggsave(paste0(outpath.results, "acc-pnet-vs-mafnet.pdf"))