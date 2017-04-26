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

genetic.bg.trgt.g <- genetic.bg.trgt[g]
genetic.bg.g <- strsplit(genetic.bg.trgt.g,"_")[[1]][1]
# set paths
path.input <- "input/th17/used_for_paper/fig2/"
path.results <- "results/fig2/"

# read experimental map (sl# to experiment)
x <- read.table(paste(sep="",path.input,"map_sam_files_to_expt.txt"),header=T,sep="\t",colClasses = "character")
# get all genetic backgrounds

bgs <- sapply(sapply(x[,"name"],strsplit, "_"),function(i) i[2])
ix <- numeric()
ix <- which(bgs == genetic.bg.g)


expts <- sort(unique(sapply(strsplit(x[ix,"expt"],"::"),function(i) i[1])))
expts <- sort(unlist(sapply(strsplit(expts,"_"),function(i) i)))

ix.list.th17 <- list()
ix.list.th0 <- list()
load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th17_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
ll.th17.all.mtls <- ll.th17
load("input/th17/used_for_paper/fig2/extended_tf_clusters_IRF4_STAT3_RORC_BATF_MAF_th0_d.tfs_100_d.tss_5000_Jan_16_2012.RData")
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
	ix.list.th17[["IRF4"]][["i.ko.th17"]] <- grep("th17_rorc_ko_irf4_chip",x[,"name"])
	ix.list.th17[["IRF4"]][["i.wt.th17"]]  <- grep("th17_rorc_wt_irf4_chip",x[,"name"])
	ix.list.th17[["STAT3"]][["i.ko.th17"]] <- grep("th17_rorc_ko_stat3_chip",x[,"name"])
	ix.list.th17[["STAT3"]][["i.wt.th17"]] <- grep("th17_rorc_wt_stat3_chip",x[,"name"])
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
# keep only singleton mtls from ll.x.
# mtl.type <- sapply(ll.th17,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
# ix <- which(sapply(mtl.type,function(i) length(strsplit(i,"_")[[1]]))==1)
# ll.th17 <- ll.th17[ix]
# mtl.type <- sapply(ll.th0,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
# ix <- which(sapply(mtl.type,function(i) length(strsplit(i,"_")[[1]]))==1)
# ll.th0 <- ll.th0[ix]
mtl.type <- sapply(ll.th17.all.mtls,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
ix <- which(sapply(mtl.type,function(i) length(strsplit(i,"_")[[1]]))==1)
ll.th17.all.mtls <- ll.th17.all.mtls[ix]
mtl.type <- sapply(ll.th0.all.mtls,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
ix <- which(sapply(mtl.type,function(i) length(strsplit(i,"_")[[1]]))==1)
ll.th0.all.mtls <- ll.th0.all.mtls[ix]

# stop("AM")
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
		ll <- ll.th17
		mtl.type <- sapply(ll,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
		if(genetic.bg.trgt.g=="rorc_all"){
			mtl.type.uniq <- table(mtl.type)
		} else {
			mtl.type.uniq <- table(mtl.type)[c(1,3,2)]			
		}
		ll.all <- ll.th17.all.mtls
		mtl.type.all <- sapply(ll.all,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
		mtl.type.uniq.all <- table(mtl.type.all)
		mtl.type.uniq.all <- mtl.type.uniq.all[which(!names(mtl.type.uniq.all) %in% names(mtl.type.uniq))]
		ix <- which(mtl.type.all %in% names(mtl.type.uniq.all))
		ll.all <- ll.th17.all.mtls[ix]
	} else {
		ix.list <- ix.list.th0
		ll <- ll.th0
		mtl.type <- sapply(ll,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
		if(genetic.bg.trgt.g=="rorc_all"){
			mtl.type.uniq <- table(mtl.type)
		} else {
			mtl.type.uniq <- table(mtl.type)[c(1,3,2)]			
		}
		ll.all <- ll.th0.all.mtls
		mtl.type.all <- sapply(ll.all,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
		mtl.type.uniq.all <- table(mtl.type.all)
		mtl.type.uniq.all <- mtl.type.uniq.all[which(!names(mtl.type.uniq.all) %in% names(mtl.type.uniq))]
		ix <- which(mtl.type.all %in% names(mtl.type.uniq.all))
		ll.all <- ll.th0.all.mtls[ix]
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
# stop("AM")
	trgts <- sapply(ll,function(i) paste(sep="",sort(unique(i$trgts)),collapse="_"))
	# trgts[which(trgts=="")] <- ""
	pval.cut <- 20
	fc.cut <- 1
	# n <- 2000		
	for(i in 1:length(ix.list)){
		tf <- names(ix.list[i])
		ix <- grep(tf,names(ll[[1]]))
		fc <- sapply(ll,function(i) i[[ix[1]]])
		pval <- sapply(ll,function(i) i[[ix[2]]])
		xlab <- paste(sep="","Fold Change [log2(reads_in_",toupper(genetic.bg.g),"_wt/reads_in_",toupper(genetic.bg.g),"_ko)]")
		ylab <- paste(sep="","pvalue [-log10]")
		main <- paste(sep="",tf," ChIP ","(genetic bg: ",cond," ",toupper(genetic.bg.g)," wt/ko ",sl.wt,"/",sl.ko,")")
		if(genetic.bg.g == "rorc"){
			highlight.genes <- c("IL17A","IL17F","IL23R","FURIN","CXCL3","CCL20","IL1R1","IL12RB2","IL22")			
			ix.highlight <- numeric()
			for(j in 1:length(highlight.genes)){
			 	ix.highlight <- unique(c(ix.highlight,grep(highlight.genes[j],trgts,ignore.case=T)))
			}
		} else {
			ix.highlight <- numeric()			
		}
		ix.sig.volcano.points <- which(pval>pval.cut & fc>fc.cut)
		ix.highlight <- intersect(ix.highlight,ix.sig.volcano.points)
		ix.proximal <- which(trgts!="")
		ix.distal <- which(trgts=="")
		ylim <- c(min(pval),min(max(pval),200))
		xlim <- c(min(fc),min(max(fc),6))		
		pval[which(pval>ylim[2])] <- ylim[2]
		fc[which(fc>xlim[2])] <- xlim[2]
		col.vec <- rainbow(length(mtl.type.uniq))
		col.vec.rgb <- col2rgb(c("darkgreen","darkorange","purple"))
		# col.vec.rgb <- c(col2rgb("darkgreen"),col2rgb("darkorange"),col2rgb("purple"))
		col.vec.rgb.trnsprnt <- c(rgb(red=col.vec.rgb["red",1], green=col.vec.rgb["green",1], blue=col.vec.rgb["blue",1], alpha=255,maxColorValue=255),
			rgb(red=col.vec.rgb["red",2], green=col.vec.rgb["green",2], blue=col.vec.rgb["blue",2], alpha=255,maxColorValue=255),
			rgb(red=col.vec.rgb["red",3], green=col.vec.rgb["green",3], blue=col.vec.rgb["blue",3], alpha=255,maxColorValue=255)
		)
		col.vec <- col.vec.rgb.trnsprnt
		col.vec <- rainbow(length(mtl.type.uniq))
		pdf(paste(sep="",path.results,cond,"_",genetic.bg.g,"_wt_vs_ko_",tf,"_ChIP_volcanos_min_mtl_pval_",mtl.min.pval,"_",date.is,".pdf"))
		# plot volcano plot
		for(l in 1:2){
			if(l == 1){ix.focus <- ix.distal;pch=20;cex=.3} else {ix.focus <- ix.proximal;pch=3;cex=.5}
			for(j in length(mtl.type.uniq):1){
				ix <- which(mtl.type==names(mtl.type.uniq)[j])
				ix <- intersect(ix,ix.focus)
				if(l == 1){mtl.type.uniq[j] <- length(ix)} else {mtl.type.uniq[j] <- mtl.type.uniq[j]+length(ix)}
				if(l==1&j==length(mtl.type.uniq)){
					plot(fc[ix],pval[ix],cex=cex,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,main=main,pch=pch,col=col.vec[j])
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
		legend(x="topleft",paste(sep="",names(mtl.type.uniq)," (n=",mtl.type.uniq,")"),fill=col.vec,cex=.5,title="MTL type")
		# plot boxplots
		trgts <- sapply(ll,function(i) paste(sep="",sort(unique(i$trgts)),collapse="_"))
		ix.proximal.all <- which(trgts!="")
		ix.distal.all <- which(trgts=="")
		# mtl.type.all <- sapply(ll.all,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
		# mtl.type.uniq.all <- table(mtl.type.all)
		# mtl.type.uniq.all <- mtl.type.uniq.all[which(!names(mtl.type.uniq.all) %in% names(mtl.type.uniq))]
		ix <- grep(tf,names(ll.all[[1]]))
		fc.all <- sapply(ll.all,function(i) i[[ix[1]]])
		pval.all <- sapply(ll.all,function(i) i[[ix[2]]])
		pval.list.all <- list()
		fc.list.all <- list()
		for(l in 1:2){
			if(l == 1) {
				ix.focus <- ix.proximal.all;pch=3;type="p"
			} else {
				ix.focus <- ix.distal.all;pch=20;type="d"
			}		
			for(j in 1:length(mtl.type.uniq.all)){
				mtl <- names(mtl.type.uniq.all)[j]
				ix <- intersect(which(mtl.type.all==mtl),ix.focus)
				e <- paste(sep="",type,"_",mtl," (n=",length(ix),")")
				pval.list.all[[e]] <- pval.all[ix]
				fc.list.all[[e]] <- fc.all[ix]
			}
		}		
		pval.list <- list()
		fc.list <- list()
		for(l in 1:2){
			if(l == 1)
			{ix.focus <- ix.proximal;pch=3;type="p"} 
			else 
			{ix.focus <- ix.distal;pch=20;type="d"}		
			for(j in 1:length(mtl.type.uniq)){
				mtl <- names(mtl.type.uniq)[j]
				ix <- intersect(which(mtl.type==mtl),ix.focus)
				e <- paste(sep="",type,"_",mtl," (n=",length(ix),")")
				pval.list[[e]] <- pval[ix]
				fc.list[[e]] <- fc[ix]
			}
		}
		p.l <- c(pval.list.all[1:3],pval.list.all[4:6],pval.list[1:3],pval.list[4:6])
		f.l <- c(fc.list.all[1:3],fc.list.all[4:6],fc.list[1:3],fc.list[4:6])
		# len <- length(p.l)*2
		# boxplot(p.l,outline=FALSE,ylab=ylab,col="gray",xaxt="n",main=main,notch=T,at=seq(from=1,to=len,by=2))	
		# text(seq(from=1,to=len,by=2),xlim=c(1,len), par("usr")[3], labels = names(p.l), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
		# par(new=TRUE)
		# boxplot(f.l,outline=FALSE,ylab=ylab,col="gray",xaxt="n",main=main,notch=T,at=seq(from=2,to=len,by=2),xlim=c(1,len))	
		# text(seq(from=2,to=len,by=2), par("usr")[3], labels = names(f.l), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
		# pf.l <- list()
		# for(i in 1:length(p.l)){
		# 	pf.l[[paste(sep="",names(f.l)[i],"_fc")]] <- f.l[[i]]
		# 	pf.l[[paste(sep="",names(f.l)[i],"_pval")]] <- p.l[[i]]
		# }		
		ix <- 1:length(p.l)
		mid.point <- length(p.l)/2+0.5
		col.vec.boxes <- c(rep("white",6),rep(col.vec,2))
		boxplot(p.l,outline=FALSE,ylab=ylab,col=col.vec.boxes,xaxt="n",main=main,notch=T)	
		text(ix, par("usr")[3], labels = names(p.l), srt = 45, adj = c(1,2), xpd = TRUE,cex=.7)
		lines(x=c(mid.point,mid.point),y=c(-10^6,10^6))
		boxplot(f.l,outline=FALSE,ylab=xlab,col=col.vec.boxes,xaxt="n",main=main,notch=T)	
		text(ix, par("usr")[3], labels = names(f.l), srt = 45, adj = c(1,2), xpd = TRUE,cex=.7)		
		lines(x=c(mid.point,mid.point),y=c(-10^6,10^6))


		
		boxplot(pval.list,outline=FALSE,ylab=ylab,col="gray",xaxt="n",main=main,notch=T)	
		text(1:length(pval.list), par("usr")[3], labels = names(pval.list), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
		boxplot(fc.list,outline=FALSE,ylab=xlab,col="gray",xaxt="n",main=main,notch=T)	
		text(1:length(fc.list), par("usr")[3], labels = names(fc.list), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
		dev.off()

		# f.nm <- paste(sep="",path.results,cond,"_",genetic.bg.g,"_wt_vs_ko_",tf,"_ChIP_volcanos_min_mtl_pval_",mtl.min.pval,"_",date.is,".pdf")
		# print.ll.verbose.fc(ll,f.nm)
	}
	
} 
# else if (genetic.bg.g=="rorc"){
# 	g  <-  1 # genetic bg from loop above (for consistancy of code)
# 	cond <- "Th17"
# 	ix.list <- list()
# 	ix.list[["P300"]][["i.ko.th17"]] <- grep("th17_rorc_ko_p300_chip",x[,"name"])
# 	ix.list[["P300"]][["i.wt.th17"]] <- grep("th17_rorc_wt_p300_chip",x[,"name"])
# 	ix.list[["IRF4"]][["i.ko.th17"]] <- grep("th17_rorc_ko_irf4_chip",x[,"name"])
# 	ix.list[["IRF4"]][["i.wt.th17"]]  <- grep("th17_rorc_wt_irf4_chip",x[,"name"])
# 	ix.list[["STAT3"]][["i.ko.th17"]] <- grep("th17_rorc_ko_stat3_chip",x[,"name"])
# 	ix.list[["STAT3"]][["i.wt.th17"]] <- grep("th17_rorc_wt_stat3_chip",x[,"name"])
# 	ix.list[["h3k4me1"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k4me1_chip",x[,"name"])
# 	ix.list[["h3k4me1"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k4me1_chip",x[,"name"])	
# 	ix.list[["h3k4me2"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k4me2_chip",x[,"name"])
# 	ix.list[["h3k4me2"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k4me2_chip",x[,"name"])	
# 	ix.list[["h3k4me3"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k4me3_chip",x[,"name"])
# 	ix.list[["h3k4me3"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k4me3_chip",x[,"name"])	
# 	ix.list[["h3ack914"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3ack914_chip",x[,"name"])
# 	ix.list[["h3ack914"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3ack914_chip",x[,"name"])
# 	ix.list[["h3k27me3"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k27me3_chip",x[,"name"])
# 	ix.list[["h3k27me3"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k27me3_chip",x[,"name"])
# 
# 
# 	# calc fold changes and pvals
# 	for(i in 1:length(ix.list)){
# 		tf <- names(ix.list[i])
# 		ix <- ix.list[[i]]
# 		# switch name to SL number
# 		names(ix) <- c(strsplit(x[ix["i.ko.th17"],"expt"],"::")[[1]][1],strsplit(x[ix["i.wt.th17"],"expt"],"::")[[1]][1])
# 		rep.num <- length(strsplit(names(ix)[1],"_")[[1]])
# 		for(j in 1:length(rep.num)){
# 			sl.ko <- strsplit(names(ix)[1],"_")[[1]][j]
# 			sl.wt <- strsplit(names(ix)[2],"_")[[1]][j]
# 			read.total <- as.numeric(sapply(strsplit(sapply(strsplit(x[,"read_counts"],"::")[ix],function(i) i[1]),"_"),function(i) i[j]))
# 			names(read.total) <- c(sl.ko,sl.wt)
# 			wt <- (sapply(ll.th17,function(i) i[[sl.wt]]) + 1)
# 			ko <- (sapply(ll.th17,function(i) i[[sl.ko]]) + 1)
# 			fc <- log2((read.total[sl.ko]/read.total[sl.wt])*wt/ko)
# 			span <- sapply(ll.th17,function(i) i$span.l)
# 			lambda.global.ko <- read.total[sl.ko]/mm9.effective.genome.size*span # how many RPMs expected by chance for each MTL in wt expt
# 			lambda.global.wt <- read.total[sl.wt]/mm9.effective.genome.size*span # how many RPMs expected by chance for each MTL in ko expt
# 			lambda.local.ko <- wt*(read.total[sl.ko]/read.total[sl.wt])
# 			lambda.local.wt <- ko*(read.total[sl.wt]/read.total[sl.ko])
# 			lambda.wt <- pmax(lambda.global.wt,lambda.local.wt)
# 			lambda.ko <- pmax(lambda.global.ko,lambda.local.ko)
# 			e.fc <- paste(tf,"_FC_",sl.wt,"_",sl.ko,sep="")
# 			e.pval <- paste(tf,"_pval_",sl.wt,"_",sl.ko,sep="")
# 			for(l in 1:length(ll.th17)){
# 				ll.th17[[l]][[e.fc]] <- fc[l]
# 				if(wt[l]>=ko[l]){
# 					ll.th17[[l]][[e.pval]] <- -1*ppois(wt[l],lambda=lambda.wt[l],lower.tail=FALSE,log.p = TRUE)/log(10)
# 				} else {
# 					ll.th17[[l]][[e.pval]] <- -1*ppois(ko[l],lambda=lambda.ko[l],lower.tail=FALSE,log.p = TRUE)/log(10)
# 				}
# 			}
# 		}
# 	}
# 
# 	# print ll.th17
# 	f.nm <- paste(path.results,"mtls.th17.rorc.wt.vs.ko.fc.xls",sep="")
# 	print.ll.verbose.fc(ll.th17,f.nm)	
# 
# 	trgts <- sapply(ll.th17,function(i) paste(sep="",sort(unique(i$trgts)),collapse="_"))
# 	trgts[which(trgts=="")] <- ""
# 	pval.cut <- 20
# 	fc.cut <- 1
# 	n <- # 2000
# 	# plot volcano plots
# 	for(i in 1:length(ix.list)){
# 		tf <- names(ix.list[i])
# 		ix <- grep(tf,names(ll.th17[[1]]))
# 		fc <- sapply(ll.th17,function(i) i[[ix[1]]])
# 		pval <- sapply(ll.th17,function(i) i[[ix[2]]])
# 		mtl.type <- sapply(ll.th17,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
# 		mtl.type.uniq <- table(mtl.type)[c(1,7,4,5,2,6,3)]
# 
# 		xlab <- paste(sep="","Fold Change [log2(reads_in_",toupper(genetic.bg.g[g]),"_wt/reads_in_",toupper(genetic.bg.g[g]),"_ko)]")
# 		ylab <- paste(sep="","pvalue [-log10]")
# 		main <- paste(sep="",tf," ChIP ","(genetic bg: ",cond," ",toupper(genetic.bg.g[g])," wt/ko ",sl.wt,"/",sl.ko,")")
# 	
# 		# xlab <- paste(sep="",names(ll.th17[[1]][ix[1]])," (rorc wt/ko)")
# 		# ylab <- paste(sep="",names(ll.th17[[1]][ix[2]])," (rorc wt/ko)")
# 		# main <- paste(sep="",tf," ChIP ","(genetic bg: Th17 ",toupper("rorc")," wt/ko)")
# 	
# 		max.pval.mtls <- sapply(ll.th17,function(i) max(i$pval))
# 		ix.sig.mtls <- which(max.pval.mtls > mtl.min.pval)
# 		# highlight.genes <- c("IL17A","IL17F","IL23R","FURIN","SMAD3","CYSLTR1","CCR6","RORC","IL21R")
# 		highlight.genes <- c("IL17A","IL17F","IL23R","FURIN","CXCL3","CCL20","IL1R1","IL12RB2","IL22")
# 		ix.highlight <- numeric()
# 		for(j in 1:length(highlight.genes)){
# 			ix.highlight <- unique(c(ix.highlight,grep(highlight.genes[j],trgts,ignore.case=T)))
# 		}
# 		ix.sig.volcano.points <- which(pval>pval.cut & fc>fc.cut)
# 		ix.highlight <- intersect(ix.highlight,intersect(ix.sig.volcano.points,ix.sig.mtls))
# 		ix.proximal <- which(trgts!="")
# 		ix.distal <- which(trgts=="")
# 		ylim <- c(min(pval),min(max(pval),200))
# 		xlim <- c(min(fc),min(max(fc),4))		
# 		pval[which(pval>ylim[2])] <- ylim[2]
# 		fc[which(fc>xlim[2])] <- xlim[2]
# 		col.vec <- rainbow(length(mtl.type.uniq))
# 		pdf(paste(sep="",path.output,genetic.bg.g[g],"_wt_vs_ko_",tf,"_ChIP_volcanos_",date.is,".pdf"))
# 		for(l in 1:2){
# 			if(l == 1){ix.focus <- intersect(ix.distal,ix.sig.mtls);pch=20;cex=.3} else {ix.focus <- intersect(ix.proximal,ix.sig.mtls);pch=3;cex=.5}
# 			for(j in 1:length(mtl.type.uniq)){
# 				ix <- which(mtl.type==names(mtl.type.uniq)[j])
# 				ix <- intersect(ix,ix.focus)
# 				if(l == 1){mtl.type.uniq[j] <- length(ix)} else {mtl.type.uniq[j] <- mtl.type.uniq[j]+length(ix)}
# 				if(l==1&j==1){
# 					plot(fc[ix],pval[ix],cex=cex,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,main=main,pch=pch,col=col.vec[j])
# 					lines(x=c(min(fc)-1,max(fc)+1),y=c(pval.cut,pval.cut))
# 					lines(x=c(-1*fc.cut,-1*fc.cut),y=c(min(pval)-100,max(pval)+100))
# 					lines(x=c(fc.cut,fc.cut),y=c(min(pval)-100,max(pval)+100))
# 				} else {
# 					points(fc[ix],pval[ix],cex=cex,pch=pch,col=col.vec[j])
# 				}
# 			}
# 		}
# 		if(length(ix.highlight)>0){
# 			arrows(x0=fc[ix.highlight],y0=pval[ix.highlight]-3,x1=fc[ix.highlight],y1=pval[ix.highlight],col="black",code=2,length=0.05)
# 			text(x=fc[ix.highlight],y=pval[ix.highlight]-4,trgts[ix.highlight],col="black",cex=.5)
# 		}
# 		legend(x="topleft",paste(sep="",names(mtl.type.uniq)," (n=",mtl.type.uniq,")"),fill=col.vec,cex=.5,title="MTL type")
# 		# get pvals and fc per mtl
# 		pval.list <- list()
# 		fc.list <- list()
# 		for(l in 1:2){
# 			if(l == 1)
# 			{ix.focus <- intersect(ix.proximal,ix.sig.mtls);pch=3;type="p"} 
# 			else 
# 			{ix.focus <- intersect(ix.distal,ix.sig.mtls);pch=20;type="d"}		
# 			for(j in 1:length(mtl.type.uniq)){
# 				mtl <- names(mtl.type.uniq)[j]
# 				ix <- intersect(which(mtl.type==mtl),ix.focus)
# 				e <- paste(sep="",type,"_",mtl,"\n(n=",length(ix),")")
# 				pval.list[[e]] <- pval[ix]
# 				fc.list[[e]] <- fc[ix]
# 			}
# 		}
# 		boxplot(pval.list,outline=FALSE,ylab=ylab,col="gray",xaxt="n",main=main,notch=T)	
# 		text(1:length(pval.list), par("usr")[3], labels = names(pval.list), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
# 		boxplot(fc.list,outline=FALSE,ylab=xlab,col="gray",xaxt="n",main=main,notch=T)	
# 		text(1:length(fc.list), par("usr")[3], labels = names(fc.list), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
# 		dev.off()
# 	}			
# 		
# }		
