##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jun 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

#####################################################################
# 0. set up variables
#####################################################################
cat("producing volcano plots\n")
source("r_scripts/th17/used_for_paper/cluster_peaks_util.R")

# genetic.bg.trgt.g <- genetic.bg.trgt[g]
# genetic.bg.g <- strsplit(genetic.bg.trgt.g,"_")[[1]][1]
# set paths
path.input <- "input/th17/used_for_paper/diff_chip_volcanoplots/"
path.results <- paste(sep="","results/diff_chip_volcanoplots/",date.is,"/")
system(paste(sep="","mkdir ", path.results))			

# read experimental map (sl# to experiment)
x <- read.table(paste(sep="",path.input,"map_sam_files_to_expt_Jul_22_2012.txt"),header=T,sep="\t",colClasses = "character")

# load mtls with counts
f.nm <- paste(sep="",path.input,"extended_tf_clusters_",paste(core.tfs,collapse="_"),"_th17_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.last.count,".RData")
load(f.nm)
f.nm <- paste(sep="",path.input,"extended_tf_clusters_",paste(core.tfs,collapse="_"),"_th0_d.tfs_",d.cut,"_d.tss_",mtbs.tss.dist,"_",date.last.count,".RData")
load(f.nm)
ll.th17.all.mtls <- ll.th17
ll.th0.all.mtls <- ll.th0

# get the number of reads per expt (for normalizing to rpms)
f.nm <- paste(sep="",path.input,"mapped_reads_per_expt_",date.last.count,".RData")
load(f.nm)


# filter low pval mtls for ll.th17, ll.th0, ll.th17.all, and ll.th0.all
max.pval.mtls <- sapply(ll.th17.all.mtls,function(i) max(i$pval))
ix.sig.mtls <- which(max.pval.mtls > mtl.min.pval)
ll.th17.all.mtls <- ll.th17.all.mtls[ix.sig.mtls]
max.pval.mtls <- sapply(ll.th0.all.mtls,function(i) max(i$pval))
ix.sig.mtls <- which(max.pval.mtls > mtl.min.pval)
ll.th0.all.mtls <- ll.th0.all.mtls[ix.sig.mtls]

# get list of diff chip experiments to perform
expts <- x[,"name"] # exp names
lin <- sapply(sapply(expts,strsplit, "\\."),function(i) i[1]) # linage tested
ko.tfs <- toupper(sapply(sapply(expts,strsplit, "\\."),function(i) i[2])) # genetic bg (tf deficiency or tf wt)
chiped.tfs <- toupper(sapply(sapply(expts,strsplit, "\\."),function(i) i[4])) # which tf was chiped
sl.wt.all <- sapply(strsplit(x[,"expt"],"::"),function(i) i[1]) # unique identifier for wt chip samples
sl.ko.all <- sapply(strsplit(x[,"expt"],"::"),function(i) i[2]) # unique identifier for ko chip samples
reads.per.expt <- n.reads.per.expt # tatal mappable reads per experiment

#####################################################################
# 1. get fold change and pval for each mtl for each diff chip we do #
#####################################################################
# for(i in 1:length(expts)){ # for each diff chip expt
	for(i in 1:3){ # for each diff chip expt
	# get lineage of comparison (so we know what mtl list to use)
	cond <- lin[i]
	if(tolower(cond)==tolower("Th17")){
		ll <- ll.th17.all.mtls
	} else {
		ll <- ll.th0.all.mtls
	}
	n <- length(ll)
	tf.chip <- chiped.tfs[i]
	sl.wt <- sl.wt.all[i]
	sl.ko <- sl.ko.all[i]

	rep.num <- length(strsplit(sl.wt,"_")[[1]]) 
	for(j in 1:length(rep.num)){
		# get sl numbers for j'th repeat
		sl.wt <- strsplit(sl.wt,"_")[[1]][j]
		sl.ko <- strsplit(sl.ko,"_")[[1]][j]
		# get read total numbers
		# read.tot.wt <- reads.per.expt[sl.wt]
		# read.tot.ko <- reads.per.expt[sl.ko]
		read.tot.wt <- 1.2*10^6
		read.tot.ko <- 1*10^6
		read.total <- c(read.tot.wt,read.tot.ko)
		names(read.total) <- c(sl.ko,sl.wt)
		# get wt and ko read counts per mtl in a vector form
		wt <- (sapply(ll,function(i) i[[sl.wt]]) + 1)
		ko <- (sapply(ll,function(i) i[[sl.ko]]) + 1)
		# calc fold change for each mtl
		fc <- log2((read.total[sl.ko]/read.total[sl.wt])*wt/ko)
		# get the span (in bp) of each mtl
		span <- sapply(ll,function(i) i$span.l)
		# calc lambda: how many RPMs expected by chance for each MTL in wt/ko expts (based on whole genome read distribution)
		lambda.global.wt <- read.total[sl.wt]/mm9.effective.genome.size*span
		lambda.global.ko <- read.total[sl.ko]/mm9.effective.genome.size*span
		# to be more stringent calculate an alternative lambda (local) 
		# (based on the alternate experiment: for ko use lambda based on wt, for wt use lambda based on ko)
		lambda.local.ko <- wt*(read.total[sl.ko]/read.total[sl.wt])
		lambda.local.wt <- ko*(read.total[sl.wt]/read.total[sl.ko])
		# chose the largest lambda based on local and global (more stringent)
		lambda.wt <- pmax(lambda.global.wt,lambda.local.wt)
		lambda.ko <- pmax(lambda.global.ko,lambda.local.ko)
		# report fold change and pval for chip diff in each mtl
		e.fc <- paste(tf.chip,"_FC_",sl.wt,"_",sl.ko,"_rep_",j,sep="")
		e.pval <- paste(tf.chip,"_pval_",sl.wt,"_",sl.ko,"_rep_",j,sep="")
		for(l in 1:length(ll)){
			ll[[l]][[e.fc]] <- fc[l]
			if(wt[l]>=ko[l]){
				ll[[l]][[e.pval]] <- -1*ppois(wt[l],lambda=lambda.wt[l],lower.tail=FALSE,log.p = TRUE)/log(10)
			} else {
				ll[[l]][[e.pval]] <- -1*ppois(ko[l],lambda=lambda.ko[l],lower.tail=FALSE,log.p = TRUE)/log(10)
			}
		}	
	}
	# overwrite the mtl list with mtl list that has pval and fc
	if(tolower(cond)==tolower("Th17")){
		ll.th17.all.mtls <- ll
	} else {
		ll.th0.all.mtls <- ll		
	}	
}
stop("AM")
#####################################################################
# 2. create figures for diff chip (volcanoplots, boxplots)          #
#####################################################################
for(i in 1:length(expts)){ # for each diff chip expt
	# get lineage of comparison (so we know what mtl list to use)
	cond <- lin[i]
	if(tolower(cond)==tolower("Th17")){
		ll <- ll.th17.all.mtls
	} else {
		ll <- ll.th0.all.mtls
	}
	n <- length(ll)
	tf.chip <- chiped.tfs[i]
	tf.ko <- ko.tfs[i]
	sl.wt <- sl.wt.all[i]
	sl.ko <- sl.ko.all[i]

	####################################	
	## now we seperate into two cases ##
	####################################
	# 1, where both chiped tf and deficient tf are tfs that were clustered into mtls (e.g. irf4 and batf, or rorc and stat3)
	# 2, where chiped tf is not part of the mtls (e.g. p300, histone modifications) but the deficient tf is
	## for case 1 we look at three types of mtls
	# a. mtls where chiped tf was singleton
	# b. mtls where only the chiped tf and ko tf where part (two ways mtls)
	# c. mtls where the deficient tf and ko tf were part together with at least one other t
	## for case 2 we look at three types of mtls
	# a. mtls where the deficient tf was not part, 
	# b. mtls where the deficient tf was singleton
	# c. mtls where the deficient tf was with 1 or more other tfs

	# get all mtl tf combinations we have
	mtl.type <- sapply(ll,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_")) 
	# handle case 1
	if( (tf.chip %in% core.tfs) & (tf.ko %in% core.tfs) ){
		ix.mtls <- list()
		# get 1.a
		ix.mtls[[paste(sep="","Singleton (",tf.chip,")")]] <- which(mtl.type %in% tf.chip)
		# get 1.b
		search.term <- paste(sep="",sort(c(tf.ko,tf.chip)),collapse="_")
		ix.mtls[[paste(sep="",search.term)]] <- which(mtl.type %in% search.term)
		# get 1.c
		ix1 <- grep(tf.ko,mtl.type)
		ix2 <- grep(tf.chip,mtl.type)
		ix.mtls[[paste(sep="",search.term,"_+")]] <- setdiff(intersect(ix1,ix2),ix.mtls[[paste(sep="",search.term)]])
		# get all combined
		ix.all <- c(ix.mtls[[paste(sep="","Singleton (",tf.chip,")")]],ix.mtls[[paste(sep="",search.term)]],ix.mtls[[paste(sep="",search.term,"_+")]])
	} else if( !(tf.chip %in% core.tfs) & (tf.ko %in% core.tfs) ) {
		ix.mtls <- list()
		# get 2.a
		ix.mtls[[paste(sep="","Singleton (other than ",tf.ko,")")]] <-  which(mtl.type %in% core.tfs[which(core.tfs!=tf.ko)])
		# get 2.b
		ix.mtls[[paste(sep="","Singleton (",tf.ko,")")]] <- which(mtl.type %in% tf.ko)
		# get 2.c
		ix.mtls[[paste(sep="",tf.ko,"_+")]] <- setdiff(grep(tf.ko,mtl.type),ix.mtls[[paste(sep="",tf.ko)]])
		# get all combined
		ix.all <- c(ix.mtls[[paste(sep="","Singleton (other than ",tf.ko,")")]],
					ix.mtls[[paste(sep="","Singleton (",tf.ko,")")]],
					ix.mtls[[paste(sep="",tf.ko,"_+")]])
	} else {
		stop("conceptual error in diff chip!!! bailing out...")
	}
	
	####################################	
	## pretty plot                    ##
	####################################
	
	# get trgt genes for each mtl (says if they are proximal or distal)
	trgts <- sapply(ll,function(i) paste(sep="",sort(unique(i$trgts)),collapse="_"))
	# some cutoffs for plotting volcano plot (where to draw quadrant lines)
	pval.cut <- 20
	fc.cut <- 1

	# there may be repeats plot for each a volcanoplot and boxplot
	sl.wt <- sl.wt.all[i]
	sl.ko <- sl.ko.all[i]
	pdf(paste(sep="",path.results,cond,"_",tf.ko,"_wt_vs_ko_",tf.chip,"_ChIP_volcanos_min_mtl_pval_",mtl.min.pval,"_",date.is,".pdf"))
	for(k in 1:length(sl.wt)){
		sl.wt.k <- sl.wt[k]
		sl.ko.k <- sl.ko[k]
		# get the index in mtl list where fc and pval results are stored for each mtl
		ix.fc <- grep( paste(sep="","FC","_",sl.wt.k,"_",sl.ko.k),names(ll[[1]]) )
		ix.pval <- grep( paste(sep="","pval","_",sl.wt.k,"_",sl.ko.k),names(ll[[1]]) )
		# get the fc and pval for diff chip i
		fc <- sapply(ll,function(i) i[[ ix.fc ]])
		pval <- sapply(ll,function(i) i[[ ix.pval ]])
		xlab <- paste(sep="","Fold Change [log2(RPKM_in_",tf.ko,"_wt/RPKM_in_",tf.ko,"_ko)]")
		ylab <- paste(sep="","pvalue [-log10]")
		main <- paste(sep="",tf.chip," ChIP ","(in ",tf.ko," WT vs. KO (",cond,") ",sl.wt.k,"/",sl.ko.k,")")
		# mark some genes to highlight
		ix.highlight <- numeric()
		for(j in 1:length(highlight.genes)){
		 	ix.highlight <- unique(c(ix.highlight,grep(highlight.genes[j],trgts,ignore.case=T)))
		}
		# keep highlight genes that were significant
		ix.sig.volcano.points <- which(pval>pval.cut & fc>fc.cut)
		ix.highlight <- intersect(intersect(ix.highlight,ix.sig.volcano.points),ix.all)
		# define proximal vs. distal genes
		ix.proximal <- which(trgts!="")
		ix.distal <- which(trgts=="")
		# cap the pvals and fc (to make display with less white space)
		ylim <- c(0,200)
		xlim <- c(-3,6)
		pval[which(pval>ylim[2])] <- ylim[2]
		fc[which(fc>xlim[2])] <- xlim[2]
		pval[which(pval<ylim[1])] <- ylim[1]
		fc[which(fc<xlim[1])] <- xlim[1]
		# define color scheme
		col.vec.rgb <- col2rgb(c("darkgreen","darkorange","purple"))
		col.vec.rgb.trnsprnt <- c(rgb(red=col.vec.rgb["red",1], green=col.vec.rgb["green",1], blue=col.vec.rgb["blue",1], alpha=255,maxColorValue=255),
			rgb(red=col.vec.rgb["red",2], green=col.vec.rgb["green",2], blue=col.vec.rgb["blue",2], alpha=255,maxColorValue=255),
			rgb(red=col.vec.rgb["red",3], green=col.vec.rgb["green",3], blue=col.vec.rgb["blue",3], alpha=255,maxColorValue=255)
		)
		col.vec <- col.vec.rgb.trnsprnt
		# plot volcano plot
		for(l in 1:2){
			if(l == 1){ix.focus <- ix.distal;pch=20;cex=.3} else {ix.focus <- ix.proximal;pch=20;cex=.3}
			for(j in length(ix.mtls):1){
				ix <- intersect(ix.mtls[[j]],ix.focus)
				if(l==1&j==length(ix.mtls)){
					plot(fc[ix],pval[ix],cex=cex,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,main=main,pch=pch,col=col.vec[j],xaxt="n")
					axis(1,at=seq(xlim[1],xlim[2],by=1))	
					lines(x=c(min(fc)-1,max(fc)+1),y=c(pval.cut,pval.cut))
					lines(x=c(-1*fc.cut,-1*fc.cut),y=c(min(pval)-100,max(pval)+100))
					lines(x=c(fc.cut,fc.cut),y=c(min(pval)-100,max(pval)+100))
				} else {
					points(fc[ix],pval[ix],cex=cex,pch=pch,col=col.vec[j])
				}
			}
		}
		if(length(ix.highlight)>0 & show.highlights){
			arrows(x0=fc[ix.highlight],y0=pval[ix.highlight]-3,x1=fc[ix.highlight],y1=pval[ix.highlight],col="black",code=2,length=0.05)
			text(x=fc[ix.highlight],y=pval[ix.highlight]-4,trgts[ix.highlight],col="black",cex=.7)
		}
		legend(x=xlim[1]-0.3,y=ylim[2]+6,paste(sep="",names(ix.mtls)," (n=",sapply(ix.mtls,length),")"),fill=col.vec,cex=.6,title="MTL type")
		# plot boxplots
		pval.list <- list()
		fc.list <- list()
		for(l in 1:2){
			if(l == 1){ix.focus <- ix.proximal;type="Prox"} else {ix.focus <- ix.distal;type="Dist"}
			for(j in 1:length(ix.mtls)){
				ix <- intersect(ix.mtls[[j]],ix.focus)	
				e <- paste(sep="",type,"_",names(ix.mtls)[j],"\n(n=",length(ix),")")
				pval.list[[e]] <- pval[ix]
				fc.list[[e]] <- fc[ix]		
			}
		}			
		ylim <- c(0,100)
		boxplot(pval.list,outline=FALSE,col=col.vec,xaxt="n",yaxt="n",notch=T,ylim=ylim)
		axis(2,at=seq(ylim[1],ylim[2],by=10))	
		ylim <- c(-3,5)
		boxplot(fc.list,outline=FALSE,col=col.vec ,xaxt="n",yaxt="n",notch=T,ylim=ylim)	
		axis(2,at=seq(ylim[1],ylim[2],by=1))
	}
	dev.off()
}
