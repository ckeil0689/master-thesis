##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2012 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
source("r_scripts/th17/used_for_paper/fig2_util.R")

cat("adding total reads to ll.th17 and ll.th0 for")
# this script is supposed to be run right after cluster_peaks_for_fig2_v2.R
## set paths
path.sam.small <- "/data/th17/data/shrink_sam/"

path.input <- "input/th17/used_for_paper/diff_chip_volcanoplots/"
path.results <- "results/diff_chip_volcanoplots/"

# load mtls without counts
f.nm <- paste(sep="",path.input,"extended_tf_clusters_",paste(core.tfs,collapse="_"),"_th17_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.last.count,".RData")
load(f.nm)
f.nm <- paste(sep="",path.input,"extended_tf_clusters_",paste(core.tfs,collapse="_"),"_th0_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.last.count,".RData")
load(f.nm)

# load file annotations
x <- read.table(expts.file.nm,header=T,sep="\t",colClasses = "character")
# we will add total reads to ll.th17 and ll.th0 for different experiments
# get all genetic backgrounds
expts.wt <- sapply(strsplit(x[,"expt"],"::"),function(i) i[1])
expts.ko <- sapply(strsplit(x[,"expt"],"::"),function(i) i[2])
expts <- c(expts.wt,expts.ko)
expts <- sort(unlist(sapply(strsplit(expts,"_"),function(i) i)))
chr <- c(paste(sep="","chr",1:19),paste(sep="","chr",c("X","Y")))
n.reads.per.expt <- numeric(length(expts))
names(n.reads.per.expt) <- expts

# already calculated experiments perviously
prev.expts <- intersect(expts,names(ll.th17[[1]]))
# only work on new experiments (don't calculate reads on each mtls again...)
expts <- setdiff(expts,prev.expts)


for(i in 1:length(expts)){
	e <- expts[i]
	# get .sam file from th17 repository
  	cat("get read counts for expt SL numbers:",e,"\n")
	cmd.line <- paste(sep="","scp am2654@bot.bio.nyu.edu:",path.sam.small,e,"/",e,".txt", " tmp/")
	system(cmd.line)
	f.nm <- paste(sep="","tmp/",e,".txt")	
	# reading takes time so I added some optimization like use scan vs. read.table and give number of lines to save mem
	n.lines <- as.numeric(strsplit(system(paste(sep="","wc -l ",f.nm),intern=T),split=" ")[[1]][2])
	# keep record of number of mappable reads
	n.reads.per.expt[e] <- n.lines
	cat("adding total reads counts to mtl list for experiment ", e,"\n")
	cat("reading data\n")
	m.s <- scan(file=f.nm,what <- list(chr=character(),read_5_prime_bp_location=numeric()),sep="\t",nlines=n.lines)
	# l is a list with each element being the indices in m.s where a read on that chr is found (sorted from begining to end of chr)
	cat("spliting reads per chr\n")
	o <- split.matrix.to.chrmosomes(x=m.s,chr=chr,n.p=n.p)
	rm(m.s)
	gc()
	cat("adding counts to ll.th17\n")
	read.counts <- count.reads.to.ll.list(ll.th17,o)
	# add counts to each MTLS th17
	for(j in 1:length(ll.th17)){
		ll.th17[[j]][[expts[i]]] <- read.counts[j]
	}
	cat("adding counts to ll.th0\n")
	read.counts <- count.reads.to.ll.list(ll.th0,o)
	# add counts to each MTLS th0
	for(j in 1:length(ll.th0)){
		ll.th0[[j]][[expts[i]]] <- read.counts[j]
	}		
	# remove .sam file from th17 repository
  	cat("get read counts for expt SL numbers:",e,"\n")
	cmd.line <- paste(sep="","rm ","tmp/",e,".txt")
	system(cmd.line)
	f.nm <- paste(sep="",path.input,"extended_tf_clusters_",paste(core.tfs,collapse="_"),"_th17_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.is,".RData")
	save(ll.th17,file=f.nm)
	f.nm <- paste(sep="",path.input,"extended_tf_clusters_",paste(core.tfs,collapse="_"),"_th0_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.is,".RData")
	save(ll.th0,file=f.nm)
}

f.nm <- paste(sep="",path.input,"mapped_reads_per_expt_",date.is,".RData")
save(n.reads.per.expt,file=f.nm)
