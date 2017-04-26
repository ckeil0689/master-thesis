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
# path.sam.small <- "input/th17/used_for_paper/rawData/sam_files_small_Nov_7_2011/"
# path.sam.small <- "/Volumes/Drobo/aviv/th17/shrinked_sam_Apr_3_2012/"
path.sam.small <- "/Volumes/Drobo/aviv/th17/shrinked_sam_Apr_8_2012/"
path.input <- "input/th17/used_for_paper/fig2/"

path.results <- "results/fig2/"

x <- read.table(paste(sep="",path.input,"map_sam_files_to_expt_May_27_2012.txt"),header=T,sep="\t",colClasses = "character")
# we will add total reads to ll.th17 and ll.th0 for different experiments
# get all genetic backgrounds
bgs <- sapply(sapply(x[,"name"],strsplit, "_"),function(i) i[2])
ix <- numeric()
# if(length(tfs.to.cluster.t)==2){
# 	genetic.bg <- tolower(strsplit(tfs.to.cluster.t,"_")[[1]][1])
# } else {
	# genetic.bg <- tolower(tfs.to.cluster.t)
# }
# for(i in 1:length(genetic.bg)){
#  	ix <- c(ix,which(bgs == genetic.bg[i]))
# }
ix <- 1:dim(x)[1]

prev.expts <- unique(grep("SL",names(ll.th17[[1]]),value=T),grep("SL",names(ll.th0[[1]]),value=T))

expts <- sort(unique(sapply(strsplit(x[ix,"expt"],"::"),function(i) i[1])))
expts <- sort(unlist(sapply(strsplit(expts,"_"),function(i) i)))
stop("AM")
expts <- setdiff(expts,prev.expts)
for(i in 1:length(expts)){
	e <- paste(sep="",expts[i],".sam")
	ix <- grep(expts[i],x[,"expt"])
	lineage <- strsplit(x[ix,"name"],"_")[[1]][1]
	f.nm <- paste(sep="",path.sam.small,e)
	# reading takes time so I added some optimization like use scan vs. read.table and give number of lines to save mem
	n.lines <- as.numeric(strsplit(system(paste(sep="","wc -l ",f.nm),intern=T),split=" ")[[1]][2])
	cat("adding total reads counts for experiment ", e,"\n")
	cat("reading data\n")
	m.s <- scan(file=f.nm,what <- list(chr=character(),read_5_prime_bp_location=numeric()),sep="\t",nlines=n.lines)
	# l is a list with each element being the indices in m.s where a read on that chr is found (sorted from begining to end of chr)
	cat("spliting reads per chr\n")
	o <- split.matrix.to.chrmosomes(x=m.s,chr=chr,n.p=n.p)
	# o is the ruler for reads per chr
	# o <- list()
	# for(j in 1:length(l)){
	# 	o[[chr[j]]]	 <- m.s[["read_5_prime_bp_location"]][l[[j]]]
	# }
	rm(m.s)
	gc()
	if(lineage=="th17"){
		cat("adding counts to ll.th17\n")
		read.counts <- count.reads.to.ll.list(ll.th17,o)
		# add counts to each MTBS
		for(j in 1:length(ll.th17)){
			ll.th17[[j]][[expts[i]]] <- read.counts[j]
		}
	} else if(lineage=="th0"){
		cat("adding counts to ll.th0\n")
		read.counts <- count.reads.to.ll.list(ll.th0,o)
		# add counts to each MTBS
		for(j in 1:length(ll.th0)){
			ll.th0[[j]][[expts[i]]] <- read.counts[j]
		}
	} else {
		stop("my logic was flawed. lineage does not equal th17 or th0, lineage = ",lineage,".\n")
	}
}
f.nm <- paste(sep="",path.input,"extended_tf_clusters_",paste(tfs,collapse="_"),"_th17_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.is,".RData")
save(ll.th17,file=f.nm)
f.nm <- paste(sep="",path.input,"extended_tf_clusters_",paste(tfs,collapse="_"),"_th0_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.is,".RData")
save(ll.th0,file=f.nm)