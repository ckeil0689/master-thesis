##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jun 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

cat("reading data\n")
source("r_scripts/th17/used_for_paper/fig2_util.R")

## set paths
path.input <- "input/th17/used_for_paper/"
path.input.save.results <- "input/th17/used_for_paper/chipClusterAnalysis/"
path.sam.small <- "input/th17/used_for_paper/rawData/sam_files_small_Nov_7_2011/"

# path.output <- "results/"
# path.input.counts <- "input/th17/used_for_paper/rawData/histone_marks_counts/"

f.nm.th17.counts <- paste(sep="",path.output,"tf_clusters_th17_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.is,".xls")
f.nm.th0.counts <- paste(sep="",path.output,"tf_clusters_th0_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.is,".xls")
th17.counts <- read.delim(f.nm.th17.counts,header=T)
th0.counts <- read.delim(f.nm.th0.counts,header=T)
macs.id.to.expt.name <- read.delim(paste(sep="",path.input.refseq,"macs_id_to_expt_map.txt"))
## get th17 histone counts in each cluster and total counts in experiment
# th17.counts <- read.delim(paste(sep="",path.input.counts,"counts_th17.txt"),header=T)
# th0.counts <- read.delim(paste(sep="",path.input.counts,"counts_th0.txt"),header=T)
# total.counts <- read.delim(paste(sep="",path.input.counts,"aligned_read_counts.txt"),header=F)
# colnames(total.counts) <- c("expt","counts")
# x <- total.counts[,"counts"]
# names(x) <- total.counts[,"expt"]
# total.counts <- x
# rm(x)

## I only keep one replicate for each experiment, the one with more peaks
## keepers <- c("SL3191_CTCF_Th0","SL3034_CTCF_Th17","SL1948_P300_Th0","SL1041_P300_Th17","SL3594_P300_Th17rorcwt","SL3595_P300_Th17rorcko",
##              "SL1952_POL2_Th0","SL1946_POL2_Th17","SL4030_H3K4me1_Th17rorcwt","SL4031_H3K4me1_Th17rorcko",
##              "SL4034_H3K4me3_Th17rorcwt","SL4035_H3K4me3_Th17rorcko","SL2993_H3K4me2_Th17rorcwt","SL2994_H3K4me2_Th17rorcko",
##              "SL4038_H3AcK914_Th17rorcwt","SL4039_H3AcK914_Th17rorcko","SL2997_H3K27me3_Th17rorcwt","SL2998_H3K27me3_Th17rorcko")
keepers.th17 <- c("CTCF_Th17_SL3034_SL3035_1","P300_Th17_SL1041_SL972_1","POL2_Th17_SL1946_SL1945_1","P300_Th17rorcwt_SL3594_SL3592_1",
             "P300_Th17rorcko_SL3595_SL3593_1","H3K4ME1_Th17rorcwt_SL4030_SL4028_1","H3K4ME3_Th17rorcwt_SL4034_SL4028_2","H3K4ME2_Th17rorcwt_SL2991_SL2999_1",
             "H3ACK9_Th17rorcwt_SL4038_SL4028_2","H3K27ME3_Th17rorcwt_SL2997_SL2999_1")
keepers.th0 <- c()
keepers.th0 <- c("IRF4_Th0_SL1235_SL1234",
             "BATF_Th0_SL3192_SL3190",
             "MAF_Th0_SL4424_SL4425",
             "STAT3_Th0_SL3780_SL3778",
             "RORC_Th0_SL3779_SL3778")
             # "P300_Th0_SL1948_SL1947",
             # "FOSL2_Th0_SL6500_SL6499")
# keepers <- c("SL3034_CTCF_Th17","SL1041_P300_Th17","SL1946_POL2_Th17","SL3594_P300_Th17rorcwt",
#              "SL3595_P300_Th17rorcko","SL4030_H3K4me1_Th17rorcwt","SL4034_H3K4me3_Th17rorcwt","SL2993_H3K4me2_Th17rorcwt",
#              "SL4038_H3AcK914_Th17rorcwt","SL2997_H3K27me3_Th17rorcwt")
# keepers.ctrl <- c("SL3035_CTCF_Th17","SL972_P300_Th17","SL1945_POL2_Th17","SL3592_P300_Th17rorcwt",
#                   "SL3593_P300_Th17rorcko","SL4028_H3K4me1_Th17rorcwt","SL4028_H3K4me3_Th17rorcwt","SL2999_H3K4me2_Th17rorcwt",
#                   "SL4028_H3AcK914_Th17rorcwt","SL2999_H3K27me3_Th17rorcwt")



keepers.sl.num <- sapply(strsplit(keepers,"_"),function(i) i[3])
keepers.ctrl.sl.num <- sapply(strsplit(keepers,"_"),function(i) i[4])


expts <- c(keepers.sl.num,keepers.ctrl.sl.num)
total.mapped.reads.th17 <- numeric(length=length(expts))
names(total.mapped.reads.th17) <- expts

for(i in 1:length(expts)){
	e <- paste(sep="",expts[i],".sam")
	f.nm <- paste(sep="",path.sam.small,e)
	# reading takes time so I added some optimization like use scan vs. read.table and give number of lines to save mem
	n.lines <- as.numeric(strsplit(system(paste(sep="","wc -l ",f.nm),intern=T),split=" ")[[1]][2])
	cat("adding total reads counts for experiment ", e,"\n")
	cat("reading data\n")
	m.s <- scan(file=f.nm,what <- list(chr=character(),read_5_prime_bp_location=numeric()),sep="\t",nlines=n.lines)
	# l is a list with each element being the indices in m.s where a read on that chr is found (sorted from begining to end of chr)
	cat("spliting reads per chr\n")
	l <- split.matrix.to.chrmosomes(x=m.s,chr=chr,n.p=n.p)
	# o is the ruler for reads per chr
	o <- list()
	for(j in 1:length(l)){
		o[[chr[j]]]	 <- m.s[["read_5_prime_bp_location"]][l[[j]]]
	}
	rm(l,m.s)
	gc()
	cat("adding counts to ll.th17\n")
	read.counts <- count.reads.to.ll.list(ll.th17,o)
	# add counts to each MTBS
	for(j in 1:length(ll.th17)){
		ll.th17[[j]][[expts[i]]] <- read.counts[j]
	}
	total.mapped.reads.th17[expts[i]] <- n.lines
}







lambda.keepers.per.bp <- total.mapped.reads.th17/mm9.genome.size
## the righ terms normalize the ctrl library to have the same number of reads as experiments
norm.ctrl.to.expt <- total.mapped.reads.th17[keepers.sl.num]/total.mapped.reads.th17[keepers.ctrl.sl.num]
names(norm.ctrl.to.expt) <- keepers.ctrl.sl.num

## load cluster stack for th17 --- named ll.th17
## load(paste(sep="",path.input,"chipClusterAnalysis/tf_clusters_th17.RData"))

## put results in matrices and then into ll.th17
m.counts <- as.matrix(th17.counts[,keepers.sl.num ])
m.ctrl.counts <- as.matrix(th17.counts[,keepers.ctrl.sl.num ])
m.pvals <- matrix(0,nr=nrow(m.counts),nc=ncol(m.counts),dimnames=dimnames(m.counts))
m.fc <- matrix(0,nr=nrow(m.counts),nc=ncol(m.counts),dimnames=dimnames(m.counts)) ## not implemented atm

## get each cluster length into a vector
clusts.length.bp <- sapply(ll.th17,function(i) i$end-i$start)
for(i in 1:length(keepers.sl.num)){
  e <- keepers.sl.num[i] ## get expt name
  e.c <- keepers.ctrl.sl.num[i] ## get ctrl name
  ## get global lambda by rescaling the global labmda by cluster length (using th17 observation)
  lambda.global <- lambda.keepers.per.bp[e]*clusts.length.bp ## get lambda for experiment (expected reads per length of cluster)
  ## get lambda local by (using chip ctrl library)
  lambda.local <- norm.ctrl.to.expt[e.c]*m.ctrl.counts[,e.c] ## get lambda for experiment (expected reads per length of cluster)
  ## keep maximum lambda (as in MACS) to capture clusters which just sit on a region with high #tags in ctrl
  lambda <- pmax(lambda.global,lambda.local) ## lambda local seems to be bigger for all clusters...
  m.pvals[,e] <- -10*ppois(m.counts[,e],lambda=lambda,lower.tail=FALSE,log.p = TRUE)/log(10) ## covert from ln to log10
  m.fc[,e] <- log2((m.counts[,e]+1)/(m.ctrl.counts[,e.c]+1))
}

## cap -10log10(pval) at 3100
m.pvals[which(m.pvals>3100)] <- 3100
colnames(m.pvals) <- colnames(m.fc) <- keepers

## get values into ll.th17[["pvals"]]
for(i in 1:length(ll.th17)){
  ll.th17[[i]][["h.pvals"]] <- m.pvals[i,]
  ll.th17[[i]][["h.fc"]] <- m.fc[i,]
}

########################
## repeat for th0
## I only keep one replicate for each experiment, the one with more peaks
keepers <- c("SL3191_CTCF_Th0","SL1948_P300_Th0","SL1952_POL2_Th0","SL3594_P300_Th17rorcwt",
             "SL3595_P300_Th17rorcko",
             "SL4031_H3K4me1_Th17rorcko","SL4035_H3K4me3_Th17rorcko",
             "SL2994_H3K4me2_Th17rorcko","SL4039_H3AcK914_Th17rorcko","SL2998_H3K27me3_Th17rorcko")
keepers.ctrl <- c("SL3190_CTCF_Th0","SL1947_P300_Th0","SL1951_POL2_Th0","SL3592_P300_Th17rorcwt",
                  "SL3593_P300_Th17rorcko",
                  "SL4029_H3K4me1_Th17rorcko","SL4029_H3K4me3_Th17rorcko",
                  "SL3000_H3K4me2_Th17rorcko","SL4029_H3AcK914_Th17rorcko","SL3000_H3K27me3_Th17rorcko")

keepers.sl.num <- sapply(strsplit(keepers,"_"),function(i) i[1])
keepers.ctrl.sl.num <- sapply(strsplit(keepers.ctrl,"_"),function(i) i[1])
lambda.keepers.per.bp <- total.counts[keepers.sl.num]/mm9.genome.size
## the righ terms normalize the ctrl library to have the same number of reads as experiments
norm.ctrl.to.expt <- total.counts[keepers.sl.num]/total.counts[keepers.ctrl.sl.num]
names(norm.ctrl.to.expt) <- keepers.ctrl.sl.num

## load cluster stack for th0 --- named ll.th0
## load(paste(sep="",path.input,"chipClusterAnalysis/tf_clusters_th0.RData"))

## put results in matrices and then into ll.th0
m.counts <- as.matrix(th0.counts[,keepers.sl.num ])
m.ctrl.counts <- as.matrix(th0.counts[,keepers.ctrl.sl.num ])
m.pvals <- matrix(0,nr=nrow(m.counts),nc=ncol(m.counts),dimnames=dimnames(m.counts))
m.fc <- matrix(0,nr=nrow(m.counts),nc=ncol(m.counts),dimnames=dimnames(m.counts)) ## not implemented atm

## get each cluster length into a vector
clusts.length.bp <- sapply(ll.th0,function(i) i$end-i$start)
for(i in 1:length(keepers.sl.num)){
  e <- keepers.sl.num[i] ## get expt name
  e.c <- keepers.ctrl.sl.num[i] ## get ctrl name
  ## get global lambda by rescaling the global labmda by cluster length (using th0 observation)
  lambda.global <- lambda.keepers.per.bp[e]*clusts.length.bp ## get lambda for experiment (expected reads per length of cluster)
  ## get lambda local by (using chip ctrl library)
  lambda.local <- norm.ctrl.to.expt[e.c]*m.ctrl.counts[,e.c] ## get lambda for experiment (expected reads per length of cluster)
  ## keep maximum lambda (as in MACS) to capture clusters which just sit on a region with high #tags in ctrl
  lambda <- pmax(lambda.global,lambda.local) ## lambda local seems to be bigger for all clusters...
  m.pvals[,e] <- -10*ppois(m.counts[,e],lambda=lambda,lower.tail=FALSE,log.p = TRUE)/log(10) ## covert from ln to log10
  m.fc[,e] <- log2((m.counts[,e]+1)/(m.ctrl.counts[,e.c]+1))
}

## cap -10log10(pval) at 3100 like MACS ChIP peaks calling algorithm
m.pvals[which(m.pvals>3100)] <- 3100
colnames(m.pvals) <- colnames(m.fc) <- keepers

## get values into ll.th0[["pvals"]]
for(i in 1:length(ll.th0)){
  ll.th0[[i]][["h.pvals"]] <- m.pvals[i,]
  ll.th0[[i]][["h.fc"]] <- m.fc[i,]
}
