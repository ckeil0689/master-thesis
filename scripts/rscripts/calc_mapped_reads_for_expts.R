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
	# remove .sam file from th17 repository
	cmd.line <- paste(sep="","rm ","tmp/",e,".txt")
	system(cmd.line)
}

f.nm <- paste(sep="",path.input,"mapped_reads_per_expt_",date.is,".RData")
save(n.reads.per.expt,file=f.nm)


