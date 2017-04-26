##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jun 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

cat("producing volcano plots for figure 2\n")
source("r_scripts/th17/used_for_paper/cluster_peaks_util.R")

genetic.bg.trgt.g <- genetic.bg.trgt[g]
genetic.bg.g <- strsplit(genetic.bg.trgt.g,"_")[[1]][1]
# set paths
path.input <- "input/th17/used_for_paper/fig2/"
path.results <- "results/fig2/"

# read experimental map (sl# to experiment)
# x <- read.table(paste(sep="",path.input,"map_sam_files_to_expt.txt"),header=T,sep="\t",colClasses = "character")
x <- read.table(paste(sep="",path.input,"map_sam_files_to_expt_Apr_8_2012.txt"),header=T,sep="\t",colClasses = "character")


# get all genetic backgrounds

bgs <- sapply(sapply(x[,"name"],strsplit, "_"),function(i) i[2])
ix <- numeric()
ix <- which(bgs == genetic.bg.g)


expts <- sort(unique(sapply(strsplit(x[ix,"expt"],"::"),function(i) i[1])))
expts <- sort(unlist(sapply(strsplit(expts,"_"),function(i) i)))

ix.list.th17 <- list()
ix.list.th0 <- list()
load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_8_2012.RData")
ll.th17.all.mtls <- ll.th17
load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th0_d.tfs_100_d.tss_5000_Apr_8_2012.RData")
ll.th0.all.mtls <- ll.th0

if(genetic.bg.trgt.g=="batf_irf4"){
	lineages <- c("Th17","Th0")
	ix.list.th17[["IRF4"]][["i.ko.th17"]] <- grep("th17_batf_ko_irf_chip",x[,"name"])
	ix.list.th17[["IRF4"]][["i.wt.th17"]] <- grep("th17_batf_wt_irf_chip",x[,"name"])
	ix.list.th0[["IRF4"]][["i.ko.th0"]] <- grep("th0_batf_ko_irf_chip",x[,"name"])
	ix.list.th0[["IRF4"]][["i.wt.th0"]] <- grep("th0_batf_wt_irf_chip",x[,"name"])
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_BATF_IRF4_th17_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_BATF_IRF4_th0_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
} else if (genetic.bg.trgt.g=="irf4_batf"){
	lineages <- c("Th17","Th0")	
	ix.list.th17[["BATF"]][["i.ko.th17"]] <- grep("th17_irf4_ko_batf_chip",x[,"name"])
	ix.list.th17[["BATF"]][["i.wt.th17"]] <- grep("th17_irf4_wt_batf_chip",x[,"name"])
	ix.list.th0[["BATF"]][["i.ko.th0"]] <- grep("th0_irf4_ko_batf_chip",x[,"name"])
	ix.list.th0[["BATF"]][["i.wt.th0"]] <- grep("th0_irf4_wt_batf_chip",x[,"name"])
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_BATF_IRF4_th17_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_BATF_IRF4_th0_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
} else if (genetic.bg.trgt.g=="rorc_irf4"){
	lineages <- c("Th17")
	ix.list.th17 <- list()
	ix.list.th17[["IRF4"]][["i.ko.th17"]] <- grep("th17_rorc_ko_irf4_chip",x[,"name"])
	ix.list.th17[["IRF4"]][["i.wt.th17"]]  <- grep("th17_rorc_wt_irf4_chip",x[,"name"])
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_RORC_IRF4_th17_d.tfs_100_d.tss_5000_Jan_13_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_RORC_IRF4_th0_d.tfs_100_d.tss_5000_Jan_13_2012.RData")	
} 	else if (genetic.bg.trgt.g=="rorc_stat3"){
	lineages <- c("Th17")
	ix.list.th17 <- list()
	ix.list.th17[["STAT3"]][["i.ko.th17"]] <- grep("th17_rorc_ko_stat3_chip",x[,"name"])
	ix.list.th17[["STAT3"]][["i.wt.th17"]] <- grep("th17_rorc_wt_stat3_chip",x[,"name"])
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_RORC_STAT3_th17_d.tfs_100_d.tss_5000_Jan_13_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_RORC_STAT3_th0_d.tfs_100_d.tss_5000_Jan_13_2012.RData")			
} else if (genetic.bg.trgt.g=="rorc_all"){
	lineages <- c("Th17")
	ix.list.th17 <- list()
	ix.list.th17[["P300"]][["i.ko.th17"]] <- grep("th17_rorc_ko_p300_chip",x[,"name"])
	ix.list.th17[["P300"]][["i.wt.th17"]] <- grep("th17_rorc_wt_p300_chip",x[,"name"])
	ix.list.th17[["h3k4me1"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k4me1_chip",x[,"name"])
	ix.list.th17[["h3k4me1"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k4me1_chip",x[,"name"])	
	ix.list.th17[["h3k4me2"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k4me2_chip",x[,"name"])
	ix.list.th17[["h3k4me2"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k4me2_chip",x[,"name"])	
	ix.list.th17[["h3k4me3"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k4me3_chip",x[,"name"])
	ix.list.th17[["h3k4me3"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k4me3_chip",x[,"name"])	
	ix.list.th17[["h3ack914"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3ack914_chip",x[,"name"])
	ix.list.th17[["h3ack914"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3ack914_chip",x[,"name"])
	ix.list.th17[["h3k27me3"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k27me3_chip",x[,"name"])
	ix.list.th17[["h3k27me3"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k27me3_chip",x[,"name"])
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th0_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
} else if (genetic.bg.trgt.g=="irf4_stat3"){
	lineages <- c("Th17")
	ix.list.th17 <- list()
	ix.list.th17[["STAT3"]][["i.ko.th17"]] <- grep("th17_irf4_ko_stat3_chip",x[,"name"])
	ix.list.th17[["STAT3"]][["i.wt.th17"]]  <- grep("th17_irf4_wt_stat3_chip",x[,"name"])
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_3_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_3_2012.RData")
} else if (genetic.bg.trgt.g=="batf_stat3"){
	lineages <- c("Th17")
	ix.list.th17 <- list()
	ix.list.th17[["STAT3"]][["i.ko.th17"]] <- grep("th17_batf_ko_stat3_chip",x[,"name"])
	ix.list.th17[["STAT3"]][["i.wt.th17"]]  <- grep("th17_batf_wt_stat3_chip",x[,"name"])
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_3_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_3_2012.RData")
} else if (genetic.bg.trgt.g=="stat3_irf4"){
	lineages <- c("Th17")
	ix.list.th17 <- list()
	ix.list.th17[["IRF4"]][["i.ko.th17"]] <- grep("th17_stat3_ko_irf4_chip",x[,"name"])
	ix.list.th17[["IRF4"]][["i.wt.th17"]]  <- grep("th17_stat3_wt_irf4_chip",x[,"name"])
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_3_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_3_2012.RData")
} else if (genetic.bg.trgt.g=="stat3_batf"){
	lineages <- c("Th17")
	ix.list.th17 <- list()
	ix.list.th17[["BATF"]][["i.ko.th17"]] <- grep("th17_stat3_ko_batf_chip",x[,"name"])
	ix.list.th17[["BATF"]][["i.wt.th17"]]  <- grep("th17_stat3_wt_batf_chip",x[,"name"])
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_3_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_3_2012.RData")
} else if (genetic.bg.trgt.g=="irf4_p300"){
	lineages <- c("Th17")
	ix.list.th17 <- list()
	ix.list.th17[["P300"]][["i.ko.th17"]] <- grep("th17_irf4_ko_p300_chip",x[,"name"])
	ix.list.th17[["P300"]][["i.wt.th17"]]  <- grep("th17_irf4_wt_p300_chip",x[,"name"])
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_8_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_8_2012.RData")
} else if (genetic.bg.trgt.g=="stat3_p300"){
	lineages <- c("Th17")
	ix.list.th17 <- list()
	ix.list.th17[["P300"]][["i.ko.th17"]] <- grep("th17_stat3_ko_p300_chip",x[,"name"])
	ix.list.th17[["P300"]][["i.wt.th17"]]  <- grep("th17_stat3_wt_p300_chip",x[,"name"])
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_8_2012.RData")
	load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Apr_8_2012.RData")
} else {
	stop("unspecified genetic bg (",genetic.bg.g,"), bailing out...")
}

# filter low pval mtls for ll.th17,ll.th0,ll.th17.all,and ll.th0.all
max.pval.mtls <- sapply(ll.th17,function(i) max(i$pval))
ix.sig.mtls <- which(max.pval.mtls > mtl.min.pval)
ll.th17 <- ll.th17[ix.sig.mtls]
max.pval.mtls <- sapply(ll.th0,function(i) max(i$pval))
ix.sig.mtls <- which(max.pval.mtls > mtl.min.pval)
ll.th0 <- ll.th0[ix.sig.mtls]
max.pval.mtls <- sapply(ll.th17.all.mtls,function(i) max(i$pval))
ix.sig.mtls <- which(max.pval.mtls > mtl.min.pval)
ll.th17.all.mtls <- ll.th17.all.mtls[ix.sig.mtls]
max.pval.mtls <- sapply(ll.th0.all.mtls,function(i) max(i$pval))
ix.sig.mtls <- which(max.pval.mtls > mtl.min.pval)
ll.th0.all.mtls <- ll.th0.all.mtls[ix.sig.mtls]

# calc fold changes and pvals for singleton mtls from ll.x.all.mtls
for(lin in 1:length(lineages)){
	cond <- lineages[lin]
	if(cond=="Th17"){
		ix.list <- ix.list.th17
		ll <- ll.th17.all.mtls
	} else {
		ix.list <- ix.list.th0
		ll <- ll.th0.all.mtls
	}
	for(i in 1:length(ix.list)){
		tf <- names(ix.list[i])
		ix <- ix.list[[i]]
		# switch name to SL number
		if(cond=="Th17"){
			names(ix) <- c(strsplit(x[ix["i.ko.th17"],"expt"],"::")[[1]][1],strsplit(x[ix["i.wt.th17"],"expt"],"::")[[1]][1])
		} else {
			names(ix) <- c(strsplit(x[ix["i.ko.th0"],"expt"],"::")[[1]][1],strsplit(x[ix["i.wt.th0"],"expt"],"::")[[1]][1])
		}
		rep.num <- length(strsplit(names(ix)[1],"_")[[1]])  # AM not working for repeats yet, need to implement
		for(j in 1:length(rep.num)){
			sl.ko <- strsplit(names(ix)[1],"_")[[1]][j]
			sl.wt <- strsplit(names(ix)[2],"_")[[1]][j]
			read.total <- as.numeric(sapply(strsplit(sapply(strsplit(x[,"read_counts"],"::")[ix],function(i) i[1]),"_"),function(i) i[j]))
			names(read.total) <- c(sl.ko,sl.wt)
			wt <- (sapply(ll,function(i) i[[sl.wt]]) + 1)
			ko <- (sapply(ll,function(i) i[[sl.ko]]) + 1)
			fc <- log2((read.total[sl.ko]/read.total[sl.wt])*wt/ko)
			span <- sapply(ll,function(i) i$span.l)
			lambda.global.ko <- read.total[sl.ko]/mm9.effective.genome.size*span # how many RPMs expected by chance for each MTL in wt expt
			lambda.global.wt <- read.total[sl.wt]/mm9.effective.genome.size*span # how many RPMs expected by chance for each MTL in ko expt
			lambda.local.ko <- wt*(read.total[sl.ko]/read.total[sl.wt])
			lambda.local.wt <- ko*(read.total[sl.wt]/read.total[sl.ko])
			lambda.wt <- pmax(lambda.global.wt,lambda.local.wt)
			lambda.ko <- pmax(lambda.global.ko,lambda.local.ko)
			e.fc <- paste(tf,"_FC_",sl.wt,"_",sl.ko,sep="")
			e.pval <- paste(tf,"_pval_",sl.wt,"_",sl.ko,sep="")
			for(l in 1:length(ll)){
				ll[[l]][[e.fc]] <- fc[l]
				if(wt[l]>=ko[l]){
					ll[[l]][[e.pval]] <- -1*ppois(wt[l],lambda=lambda.wt[l],lower.tail=FALSE,log.p = TRUE)/log(10)
				} else {
					ll[[l]][[e.pval]] <- -1*ppois(ko[l],lambda=lambda.ko[l],lower.tail=FALSE,log.p = TRUE)/log(10)
				}
			}
		}
	}
	if(cond=="Th17"){
		ll.th17.all.mtls <- ll
	} else {
		ll.th0.all.mtls <- ll		
	}
}

# calc fold changes and pvals
for(lin in 1:length(lineages)){
	cond <- lineages[lin]
	if(cond=="Th17"){
		ix.list <- ix.list.th17
		ll <- ll.th17.all.mtls
	} else {
		ix.list <- ix.list.th0
		ll <- ll.th0.all.mtls
	}
		mtl.type <- sapply(ll,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_")) 
		if(genetic.bg.trgt.g=="rorc_all"){
			ix.mtls <- list()
			ix.mtls[["Singleton (-Rorc)"]] <- which(mtl.type %in% core.tfs[which(core.tfs!="RORC")])
			search.term <- "RORC"
			ix.mtls[[paste(sep="",search.term)]] <- which(mtl.type %in% search.term)
			ix.mtls[[paste(sep="",search.term,"_+")]] <- setdiff(grep(search.term,mtl.type),ix.mtls[[paste(sep="",search.term)]])
			ix.all <- c(ix.mtls[["Singleton (-Rorc)"]],ix.mtls[[paste(sep="",search.term)]],ix.mtls[[paste(sep="",search.term,"_+")]])
		} else if( length(grep("p300",genetic.bg.trgt.g)) > 0 ){
			# genetic.bg.trgt.g=="stat3_p300" | genetic.bg.trgt.g=="irf4_p300"
			search.term <- toupper(strsplit(genetic.bg.trgt.g,"_")[[1]][1])
			ix.mtls <- list()
			ix.mtls[[paste(sep="","Singleton (-",search.term,")")]] <- which(mtl.type %in% core.tfs[which(core.tfs!=search.term)])
			ix.mtls[[paste(sep="",search.term)]] <- which(mtl.type %in% search.term)
			ix.mtls[[paste(sep="",search.term,"_+")]] <- setdiff(grep(search.term,mtl.type),ix.mtls[[paste(sep="",search.term)]])
			ix.all <- c(ix.mtls[[paste(sep="","Singleton (-",search.term,")")]],ix.mtls[[paste(sep="",search.term)]],ix.mtls[[paste(sep="",search.term,"_+")]])
		} else {
			ix.mtls <- list()
			# ix.mtls[["Singleton"]] <- which(mtl.type %in% core.tfs)
			tfs <- strsplit(toupper(genetic.bg.trgt.g),"_")[[1]]
			tf.chipped <- strsplit(toupper(genetic.bg.trgt.g),"_")[[1]][2]
			# ix.mtls[["Singleton"]] <- which(mtl.type %in% tf.chipped)
			ix.mtls[[paste(sep="","Singleton (",tf.chipped,")")]] <- which(mtl.type %in% tf.chipped)
			search.term <- paste(sep="",sort(tfs),collapse="_")
			ix.mtls[[paste(sep="",search.term)]] <- which(mtl.type %in% search.term)
			ix1 <- grep(tfs[1],mtl.type)
			ix2 <- grep(tfs[2],mtl.type)
			ix.mtls[[paste(sep="",search.term,"_+")]] <- setdiff(intersect(ix1,ix2),ix.mtls[[paste(sep="",search.term)]])
			# ix.all <- c(ix.mtls[["Singleton"]],ix.mtls[[paste(sep="",search.term)]],ix.mtls[[paste(sep="",search.term,"_+")]])
			ix.all <- c(ix.mtls[[paste(sep="","Singleton (",tf.chipped,")")]],ix.mtls[[paste(sep="",search.term)]],ix.mtls[[paste(sep="",search.term,"_+")]])
		}
		
		trgts <- sapply(ll,function(i) paste(sep="",sort(unique(i$trgts)),collapse="_"))
		pval.cut <- 20
		fc.cut <- 1
		for(i in 1:length(ix.list)){
			tf <- names(ix.list[i])
			ix <- grep(tf,names(ll[[1]]))
			fc <- sapply(ll,function(i) i[[ix[1]]])
			pval <- sapply(ll,function(i) i[[ix[2]]])
			xlab <- paste(sep="","Fold Change [log2(RPKM_in_",toupper(genetic.bg.g),"_wt/RPKM_in_",toupper(genetic.bg.g),"_ko)]")
			ylab <- paste(sep="","pvalue [-log10]")
			main <- paste(sep="",tf," ChIP ","(genetic bg: ",cond," ",toupper(genetic.bg.g)," wt/ko ",sl.wt,"/",sl.ko,")")
			if(genetic.bg.g == "rorc"){
				highlight.genes <- c("IL17A","IL17F","IL23R","FURIN","CXCL3","CCL20","IL1R1","IL12RB2","IL22")			
				ix.highlight <- numeric()
				for(j in 1:length(highlight.genes)){
				 	ix.highlight <- unique(c(ix.highlight,grep(highlight.genes[j],trgts,ignore.case=T)))
				}
			} else {
				# ix.highlight <- numeric()			
				highlight.genes <- c("IL17A","IL17F","IL23R","FURIN","CXCL3","CCL20","IL1R1","IL12RB2","IL22")			
				ix.highlight <- numeric()
				for(j in 1:length(highlight.genes)){
				 	ix.highlight <- unique(c(ix.highlight,grep(highlight.genes[j],trgts,ignore.case=T)))
				}
			}
			ix.sig.volcano.points <- which(pval>pval.cut & fc>fc.cut)
			ix.highlight <- intersect(intersect(ix.highlight,ix.sig.volcano.points),ix.all)
			ix.proximal <- which(trgts!="")
			ix.distal <- which(trgts=="")
			# ylim <- c(0,min(max(pval),200))
			# xlim <- c(min(min(fc),-3),min(max(fc),6))
			ylim <- c(0,200)
			xlim <- c(-3,6)
			pval[which(pval>ylim[2])] <- ylim[2]
			fc[which(fc>xlim[2])] <- xlim[2]
			pval[which(pval<ylim[1])] <- ylim[1]
			fc[which(fc<xlim[1])] <- xlim[1]
			col.vec.rgb <- col2rgb(c("darkgreen","darkorange","purple"))
			col.vec.rgb.trnsprnt <- c(rgb(red=col.vec.rgb["red",1], green=col.vec.rgb["green",1], blue=col.vec.rgb["blue",1], alpha=255,maxColorValue=255),
				rgb(red=col.vec.rgb["red",2], green=col.vec.rgb["green",2], blue=col.vec.rgb["blue",2], alpha=255,maxColorValue=255),
				rgb(red=col.vec.rgb["red",3], green=col.vec.rgb["green",3], blue=col.vec.rgb["blue",3], alpha=255,maxColorValue=255)
			)
			col.vec <- col.vec.rgb.trnsprnt
			pdf(paste(sep="",path.results,cond,"_",genetic.bg.g,"_wt_vs_ko_",tf,"_ChIP_volcanos_min_mtl_pval_",mtl.min.pval,"_",date.is,".pdf"))
			# plot volcano plot
			for(l in 1:2){
				if(l == 1){ix.focus <- ix.distal;pch=20;cex=.3} else {ix.focus <- ix.proximal;pch=3;cex=.5}
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
			if(length(ix.highlight)>0){
				arrows(x0=fc[ix.highlight],y0=pval[ix.highlight]-3,x1=fc[ix.highlight],y1=pval[ix.highlight],col="black",code=2,length=0.05)
				text(x=fc[ix.highlight],y=pval[ix.highlight]-4,trgts[ix.highlight],col="black",cex=.5)
			}		
			# legend(x="topleft",paste(sep="",names(ix.mtls)," (n=",sapply(ix.mtls,length),")"),fill=col.vec,cex=.5,title="MTL type")
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
			# boxplot(pval.list,outline=FALSE,ylab=ylab,col=col.vec,xaxt="n",yaxt="n",main=main,notch=T,ylim=ylim)
			boxplot(pval.list,outline=FALSE,col=col.vec,xaxt="n",yaxt="n",notch=T,ylim=ylim)
			axis(2,at=seq(ylim[1],ylim[2],by=10))	
			# text(1:length(pval.list), par("usr")[3], labels = names(pval.list), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
			# legend(x="topleft",paste(sep="",names(ix.mtls)," (n=",sapply(ix.mtls,length),")"),fill=col.vec,cex=.5,title="MTL type")
			ylim <- c(-3,5)
			# boxplot(fc.list,outline=FALSE,ylab=ylab,col=col.vec ,xaxt="n",yaxt="n",main=main,notch=T,ylim=ylim)	
			boxplot(fc.list,outline=FALSE,col=col.vec ,xaxt="n",yaxt="n",notch=T,ylim=ylim)	
			axis(2,at=seq(ylim[1],ylim[2],by=1))
			# text(1:length(fc.list), par("usr")[3], labels = names(fc.list), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
			# legend(x="topleft",paste(sep="",names(ix.mtls)," (n=",sapply(ix.mtls,length),")"),fill=col.vec,cex=.5,title="MTL type")
			dev.off()
		}
}
# f.nm <- paste(sep="",path.results,cond,"_",genetic.bg.g,"_wt_vs_ko_min_mtl_pval_",mtl.min.pval,"_",date.is,".xls")
# print.ll.verbose.all(ll,f.nm=f.nm)


