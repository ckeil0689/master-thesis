##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jun 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '


cat("adding total reads to ll.th17 and ll.th0 for")
# this script is supposed to be run right after cluster_peaks_for_fig2.R
## set paths
path.sam.small <- "input/th17/used_for_paper/rawData/sam_files_small_Nov_7_2011/"
path.results <- "results/fig2/"

x <- read.table(paste(sep="",path.input.save.results,"map_sam_files_to_expt.txt"),header=T,sep="\t",colClasses = "character")
# we will add total reads to ll.th17 and ll.th0 for different experiments
# get all genetic backgrounds
bgs <- sapply(sapply(x[,"name"],strsplit, "_"),function(i) i[2])
ix <- numeric()
for(i in 1:length(genetic.bg)){
 	ix <- c(ix,which(bgs == genetic.bg[i]))
}

expts <- sort(unique(sapply(strsplit(x[ix,"expt"],"::"),function(i) i[1])))
expts <- sort(unlist(sapply(strsplit(expts,"_"),function(i) i)))

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
		stop("my logic is flawed. lineage does not equal th17 or th0, lineage = ",lineage,".\n")
	}
}
save(ll.th0,file=paste(sep="",path.input.save.results,"ll_th0_bg_", paste(sep="",genetic.bg,collapse="_"),"_tf_cluster_",
	paste(sep="",tfs.to.cluster,collapse="_"),"_",date.is,".RData"))
save(ll.th17,file=paste(sep="",path.input.save.results,"ll_th17_bg_", paste(sep="",genetic.bg,collapse="_"),"_tf_cluster_",
		paste(sep="",tfs.to.cluster,collapse="_"),"_",date.is,".RData"))
stop("AM")

if (genetic.bg[1]=="batf" | genetic.bg[1]=="irf4"){
	# for(i in 1:length(genetic.bg)){
	# 	if(genetic.bg[i]=="batf"){
	# 		i.ko.th17 <- grep("th17_batf_ko_irf_chip",x[,"name"])
	# 		i.wt.th17 <- grep("th17_batf_wt_irf_chip",x[,"name"])
	# 		i.ko.th0 <- grep("th0_batf_ko_irf_chip",x[,"name"])
	# 		i.wt.th0 <- grep("th0_batf_wt_irf_chip",x[,"name"])
	# 		main1="ChIP peaks pvalues (Th17 wt)"
	# 		main2="IRF4 ChIP - BATF (wt/ko) genetic background";
	# 		# main3="IRF4 ChIP Pval>500 (BATF wt vs. ko)";main4="ChIP peaks pvalues Pval<200 (BATF wt)"
	# 		ylab1="log2( reads_batf_wt / reads_batf_ko )";ylab2="pval"
	# 	} else if (genetic.bg[i]=="irf4"){
	# 		i.ko.th17 <- grep("th17_irf4_ko_batf_chip",x[,"name"])
	# 		i.wt.th17 <- grep("th17_irf4_wt_batf_chip",x[,"name"])
	# 		i.ko.th0 <- grep("th0_irf4_ko_batf_chip",x[,"name"])
	# 		i.wt.th0 <- grep("th0_irf4_wt_batf_chip",x[,"name"])
	# 		main1="ChIP peaks pvalues (Th17 wt)"				
	# 		main2="BATF ChIP - IRF4 (wt/ko) genetic background"
	# 		# main3="BATF ChIP Pval>500 (IRF4 wt vs. ko)";main4="ChIP peaks pvalues Pval<200 (IRF4 wt)"	
	# 		ylab1="log2( reads_irf4_wt / reads_irf4_ko )";ylab2="pval"
	# 	}	
	# 	ix <- c(i.ko.th17,i.wt.th17,i.ko.th0,i.wt.th0)
	# 	names(ix) <- c(strsplit(x[i.ko.th17,"expt"],"::")[[1]][1],strsplit(x[i.wt.th17,"expt"],"::")[[1]][1],
	# 	strsplit(x[i.ko.th0,"expt"],"::")[[1]][1],strsplit(x[i.wt.th0,"expt"],"::")[[1]][1])
	# 	sl.ko.th17 <- names(ix)[1]
	# 	sl.wt.th17 <- names(ix)[2]
	# 	sl.ko.th0 <- names(ix)[3]
	# 	sl.wt.th0 <- names(ix)[4]
	# 
	# 	read.total <- as.numeric(sapply(strsplit(x[,"read_counts"],"::")[ix],function(i) i[1]))
	# 	names(read.total) <- names(ix)
	# 
	# 	pdf(paste(sep="",path.output,genetic.bg[i],"_wt_vs_ko_",date.is,".pdf"))	
	# 	# par(mfrow=c(3,2))
	# 	## plot for th17
	# 	boxplots.wt.vs.ko(ll=ll.th17,sl.wt=sl.wt.th17,sl.ko=sl.ko.th17,tf.to.s=tfs.to.cluster,
	# 	main1=paste("Th17",main1),main2=paste("Th17",main2),
	# 	ylab1=ylab1,ylab2=ylab2,read.total=read.total)
	# 	## plot for th0
	# 	boxplots.wt.vs.ko(ll=ll.th0,sl.wt=sl.wt.th0,sl.ko=sl.ko.th0,tf.to.s=tfs.to.cluster,
	# 	main1=paste("Th0",main1),main2=paste("Th0",main2),
	# 	# main3=paste("Th0",main3),main4=paste("Th0",main4),
	# 	ylab1=ylab1,ylab2=ylab2,read.total=read.total)	
	# 	dev.off()
	# }
	for(g in 1:length(genetic.bg)){
		ix.list.th17 <- list()
		ix.list.th0 <- list()
		if(genetic.bg[g]=="batf"){
			ix.list.th17[["IRF4"]][["i.ko.th17"]] <- grep("th17_batf_ko_irf_chip",x[,"name"])
			ix.list.th17[["IRF4"]][["i.wt.th17"]] <- grep("th17_batf_wt_irf_chip",x[,"name"])
			ix.list.th0[["IRF4"]][["i.ko.th0"]] <- grep("th0_batf_ko_irf_chip",x[,"name"])
			ix.list.th0[["IRF4"]][["i.wt.th0"]] <- grep("th0_batf_wt_irf_chip",x[,"name"])
		} else if (genetic.bg[g]=="irf4"){
			ix.list.th17[["BATF"]][["i.ko.th17"]] <- grep("th17_irf4_ko_batf_chip",x[,"name"])
			ix.list.th17[["BATF"]][["i.wt.th17"]] <- grep("th17_irf4_wt_batf_chip",x[,"name"])
			ix.list.th0[["BATF"]][["i.ko.th0"]] <- grep("th0_irf4_ko_batf_chip",x[,"name"])
			ix.list.th0[["BATF"]][["i.wt.th0"]] <- grep("th0_irf4_wt_batf_chip",x[,"name"])
		}
		# calc fold changes and pvals
		for(lineage in 1:2){
			if(lineage==1){
				cond <- "Th17"
				ix.list <- ix.list.th17
				ll <- ll.th17
			} else {
				cond <- "Th0"					
				ix.list <- ix.list.th0
				ll <- ll.th0
			}
			for(i in 1:length(ix.list)){
				tf <- names(ix.list[i])
				ix <- ix.list[[i]]
				# switch name to SL number
				if(lineage==1){
					names(ix) <- c(strsplit(x[ix["i.ko.th17"],"expt"],"::")[[1]][1],strsplit(x[ix["i.wt.th17"],"expt"],"::")[[1]][1])
				} else {
					names(ix) <- c(strsplit(x[ix["i.ko.th0"],"expt"],"::")[[1]][1],strsplit(x[ix["i.wt.th0"],"expt"],"::")[[1]][1])
				}
				rep.num <- length(strsplit(names(ix)[1],"_")[[1]])
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

			trgts <- sapply(ll,function(i) paste(sep="",sort(unique(i$trgts)),collapse="_"))
			trgts[which(trgts=="")] <- ""
			pval.cut <- 20
			fc.cut <- 1
			# n <- 2000		
			for(i in 1:length(ix.list)){
				tf <- names(ix.list[i])
				ix <- grep(tf,names(ll[[1]]))
				fc <- sapply(ll,function(i) i[[ix[1]]])
				pval <- sapply(ll,function(i) i[[ix[2]]])
				mtl.type <- sapply(ll,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
				mtl.type.uniq <- table(mtl.type)[c(1,3,2)]
				# xlab <- paste(sep="",names(ll[[1]][ix[1]])," (",toupper(genetic.bg[g])," wt/ko)")
				xlab <- paste(sep="","Fold Change [log2(reads_in_",toupper(genetic.bg[g]),"_wt/reads_in_",toupper(genetic.bg[g]),"_ko)]")
				ylab <- paste(sep="","pvalue [-log10]")
				main <- paste(sep="",tf," ChIP ","(genetic bg: ",cond," ",toupper(genetic.bg[g])," wt/ko ",sl.wt,"/",sl.ko,")")
				max.pval.mtls <- sapply(ll,function(i) max(i$pval))
				ix.sig.mtls <- which(max.pval.mtls > 500)
				highlight.genes <- c("IL17A","IL17F","IL23R","FURIN","CXCL3","CCL20","IL1R1","IL12RB2","IL22")
				ix.highlight <- numeric()
				# for(j in 1:length(highlight.genes)){
				# 	ix.highlight <- unique(c(ix.highlight,grep(highlight.genes[j],trgts,ignore.case=T)))
				# }
				ix.sig.volcano.points <- which(pval>pval.cut & fc>fc.cut)
				ix.highlight <- intersect(ix.highlight,intersect(ix.sig.volcano.points,ix.sig.mtls))
				ix.proximal <- which(trgts!="")
				ix.distal <- which(trgts=="")
				ylim <- c(min(pval),min(max(pval),200))
				xlim <- c(min(fc),min(max(fc),6))		
				pval[which(pval>ylim[2])] <- ylim[2]
				fc[which(fc>xlim[2])] <- xlim[2]
				col.vec <- rainbow(length(mtl.type.uniq))
				pdf(paste(sep="",path.output,cond,"_",genetic.bg[g],"_wt_vs_ko_",tf,"_ChIP_volcanos_",date.is,".pdf"))
				for(l in 1:2){
					if(l == 1){ix.focus <- intersect(ix.distal,ix.sig.mtls);pch=20;cex=.3} else {ix.focus <- intersect(ix.proximal,ix.sig.mtls);pch=3;cex=.5}
					for(j in length(mtl.type.uniq):1){
						ix <- which(mtl.type==names(mtl.type.uniq)[j])
						ix <- intersect(ix,ix.focus)
						if(l == 1){mtl.type.uniq[j] <- length(ix)} else {mtl.type.uniq[j] <- mtl.type.uniq[j]+length(ix)}
						if(l==1&j==3){
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
				pval.list <- list()
				fc.list <- list()
				for(l in 1:2){
					if(l == 1)
					{ix.focus <- intersect(ix.proximal,ix.sig.mtls);pch=3;type="p"} 
					else 
					{ix.focus <- intersect(ix.distal,ix.sig.mtls);pch=20;type="d"}		
					for(j in 1:length(mtl.type.uniq)){
						mtl <- names(mtl.type.uniq)[j]
						ix <- intersect(which(mtl.type==mtl),ix.focus)
						e <- paste(sep="",type,"_",mtl,"\n(n=",length(ix),")")
						pval.list[[e]] <- pval[ix]
						fc.list[[e]] <- fc[ix]
					}
				}
				boxplot(pval.list,outline=FALSE,ylab=ylab,col="gray",xaxt="n",main=main,notch=T)	
				text(1:length(pval.list), par("usr")[3], labels = names(pval.list), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
				boxplot(fc.list,outline=FALSE,ylab=xlab,col="gray",xaxt="n",main=main,notch=T)	
				text(1:length(fc.list), par("usr")[3], labels = names(fc.list), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
				dev.off()
				if(cond=="Th17"){
					ll.th17 <- ll
					f.nm <- paste(path.results,"mtls.th17.",genetic.bg[g],".wt.vs.ko.fc.xls",sep="")
					print.ll.verbose.fc(ll.th17,f.nm)
				} else {
					ll.th0 <- ll
					f.nm <- paste(path.results,"mtls.th0.",genetic.bg[g],".wt.vs.ko.fc.xls",sep="")
					print.ll.verbose.fc(ll.th0,f.nm)
				}
			}
		}
	}
	
} else if (genetic.bg[1]=="rorc"){
	g  <-  1 # genetic bg from loop above (for consistancy of code)
	cond <- "Th17"
	ix.list <- list()
	ix.list[["P300"]][["i.ko.th17"]] <- grep("th17_rorc_ko_p300_chip",x[,"name"])
	ix.list[["P300"]][["i.wt.th17"]] <- grep("th17_rorc_wt_p300_chip",x[,"name"])
	ix.list[["IRF4"]][["i.ko.th17"]] <- grep("th17_rorc_ko_irf4_chip",x[,"name"])
	ix.list[["IRF4"]][["i.wt.th17"]]  <- grep("th17_rorc_wt_irf4_chip",x[,"name"])
	ix.list[["STAT3"]][["i.ko.th17"]] <- grep("th17_rorc_ko_stat3_chip",x[,"name"])
	ix.list[["STAT3"]][["i.wt.th17"]] <- grep("th17_rorc_wt_stat3_chip",x[,"name"])
	ix.list[["h3k4me1"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k4me1_chip",x[,"name"])
	ix.list[["h3k4me1"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k4me1_chip",x[,"name"])	
	ix.list[["h3k4me2"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k4me2_chip",x[,"name"])
	ix.list[["h3k4me2"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k4me2_chip",x[,"name"])	
	ix.list[["h3k4me3"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k4me3_chip",x[,"name"])
	ix.list[["h3k4me3"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k4me3_chip",x[,"name"])	
	ix.list[["h3ack914"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3ack914_chip",x[,"name"])
	ix.list[["h3ack914"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3ack914_chip",x[,"name"])
	ix.list[["h3k27me3"]][["i.ko.th17"]] <- grep("th17_rorc_ko_h3k27me3_chip",x[,"name"])
	ix.list[["h3k27me3"]][["i.wt.th17"]] <- grep("th17_rorc_wt_h3k27me3_chip",x[,"name"])
	
	
	# calc fold changes and pvals
	for(i in 1:length(ix.list)){
		tf <- names(ix.list[i])
		ix <- ix.list[[i]]
		# switch name to SL number
		names(ix) <- c(strsplit(x[ix["i.ko.th17"],"expt"],"::")[[1]][1],strsplit(x[ix["i.wt.th17"],"expt"],"::")[[1]][1])
		rep.num <- length(strsplit(names(ix)[1],"_")[[1]])
		for(j in 1:length(rep.num)){
			sl.ko <- strsplit(names(ix)[1],"_")[[1]][j]
			sl.wt <- strsplit(names(ix)[2],"_")[[1]][j]
			read.total <- as.numeric(sapply(strsplit(sapply(strsplit(x[,"read_counts"],"::")[ix],function(i) i[1]),"_"),function(i) i[j]))
			names(read.total) <- c(sl.ko,sl.wt)
			wt <- (sapply(ll.th17,function(i) i[[sl.wt]]) + 1)
			ko <- (sapply(ll.th17,function(i) i[[sl.ko]]) + 1)
			fc <- log2((read.total[sl.ko]/read.total[sl.wt])*wt/ko)
			span <- sapply(ll.th17,function(i) i$span.l)
			lambda.global.ko <- read.total[sl.ko]/mm9.effective.genome.size*span # how many RPMs expected by chance for each MTL in wt expt
			lambda.global.wt <- read.total[sl.wt]/mm9.effective.genome.size*span # how many RPMs expected by chance for each MTL in ko expt
			lambda.local.ko <- wt*(read.total[sl.ko]/read.total[sl.wt])
			lambda.local.wt <- ko*(read.total[sl.wt]/read.total[sl.ko])
			lambda.wt <- pmax(lambda.global.wt,lambda.local.wt)
			lambda.ko <- pmax(lambda.global.ko,lambda.local.ko)
			e.fc <- paste(tf,"_FC_",sl.wt,"_",sl.ko,sep="")
			e.pval <- paste(tf,"_pval_",sl.wt,"_",sl.ko,sep="")
			for(l in 1:length(ll.th17)){
				ll.th17[[l]][[e.fc]] <- fc[l]
				if(wt[l]>=ko[l]){
					ll.th17[[l]][[e.pval]] <- -1*ppois(wt[l],lambda=lambda.wt[l],lower.tail=FALSE,log.p = TRUE)/log(10)
				} else {
					ll.th17[[l]][[e.pval]] <- -1*ppois(ko[l],lambda=lambda.ko[l],lower.tail=FALSE,log.p = TRUE)/log(10)
				}
			}
		}
	}
	
	# print ll.th17
	f.nm <- paste(path.results,"mtls.th17.rorc.wt.vs.ko.fc.xls",sep="")
	print.ll.verbose.fc(ll.th17,f.nm)	

	trgts <- sapply(ll.th17,function(i) paste(sep="",sort(unique(i$trgts)),collapse="_"))
	trgts[which(trgts=="")] <- ""
	pval.cut <- 20
	fc.cut <- 1
	n <- # 2000
	# plot volcano plots
	for(i in 1:length(ix.list)){
		tf <- names(ix.list[i])
		ix <- grep(tf,names(ll.th17[[1]]))
		fc <- sapply(ll.th17,function(i) i[[ix[1]]])
		pval <- sapply(ll.th17,function(i) i[[ix[2]]])
		mtl.type <- sapply(ll.th17,function(i) paste(sep="",sort(unique(i$tf.to.s)),collapse="_"))
		mtl.type.uniq <- table(mtl.type)[c(1,7,4,5,2,6,3)]

		xlab <- paste(sep="","Fold Change [log2(reads_in_",toupper(genetic.bg[g]),"_wt/reads_in_",toupper(genetic.bg[g]),"_ko)]")
		ylab <- paste(sep="","pvalue [-log10]")
		main <- paste(sep="",tf," ChIP ","(genetic bg: ",cond," ",toupper(genetic.bg[g])," wt/ko ",sl.wt,"/",sl.ko,")")
		
		# xlab <- paste(sep="",names(ll.th17[[1]][ix[1]])," (rorc wt/ko)")
		# ylab <- paste(sep="",names(ll.th17[[1]][ix[2]])," (rorc wt/ko)")
		# main <- paste(sep="",tf," ChIP ","(genetic bg: Th17 ",toupper("rorc")," wt/ko)")
		
		max.pval.mtls <- sapply(ll.th17,function(i) max(i$pval))
		ix.sig.mtls <- which(max.pval.mtls > 500)
		# highlight.genes <- c("IL17A","IL17F","IL23R","FURIN","SMAD3","CYSLTR1","CCR6","RORC","IL21R")
		highlight.genes <- c("IL17A","IL17F","IL23R","FURIN","CXCL3","CCL20","IL1R1","IL12RB2","IL22")
		ix.highlight <- numeric()
		for(j in 1:length(highlight.genes)){
			ix.highlight <- unique(c(ix.highlight,grep(highlight.genes[j],trgts,ignore.case=T)))
		}
		ix.sig.volcano.points <- which(pval>pval.cut & fc>fc.cut)
		ix.highlight <- intersect(ix.highlight,intersect(ix.sig.volcano.points,ix.sig.mtls))
		ix.proximal <- which(trgts!="")
		ix.distal <- which(trgts=="")
		ylim <- c(min(pval),min(max(pval),200))
		xlim <- c(min(fc),min(max(fc),4))		
		pval[which(pval>ylim[2])] <- ylim[2]
		fc[which(fc>xlim[2])] <- xlim[2]
		col.vec <- rainbow(length(mtl.type.uniq))
		pdf(paste(sep="",path.output,genetic.bg[g],"_wt_vs_ko_",tf,"_ChIP_volcanos_",date.is,".pdf"))
		for(l in 1:2){
			if(l == 1){ix.focus <- intersect(ix.distal,ix.sig.mtls);pch=20;cex=.3} else {ix.focus <- intersect(ix.proximal,ix.sig.mtls);pch=3;cex=.5}
			for(j in 1:length(mtl.type.uniq)){
				ix <- which(mtl.type==names(mtl.type.uniq)[j])
				ix <- intersect(ix,ix.focus)
				if(l == 1){mtl.type.uniq[j] <- length(ix)} else {mtl.type.uniq[j] <- mtl.type.uniq[j]+length(ix)}
				if(l==1&j==1){
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
		# get pvals and fc per mtl
		pval.list <- list()
		fc.list <- list()
		for(l in 1:2){
			if(l == 1)
			{ix.focus <- intersect(ix.proximal,ix.sig.mtls);pch=3;type="p"} 
			else 
			{ix.focus <- intersect(ix.distal,ix.sig.mtls);pch=20;type="d"}		
			for(j in 1:length(mtl.type.uniq)){
				mtl <- names(mtl.type.uniq)[j]
				ix <- intersect(which(mtl.type==mtl),ix.focus)
				e <- paste(sep="",type,"_",mtl,"\n(n=",length(ix),")")
				pval.list[[e]] <- pval[ix]
				fc.list[[e]] <- fc[ix]
			}
		}
		boxplot(pval.list,outline=FALSE,ylab=ylab,col="gray",xaxt="n",main=main,notch=T)	
		text(1:length(pval.list), par("usr")[3], labels = names(pval.list), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
		boxplot(fc.list,outline=FALSE,ylab=xlab,col="gray",xaxt="n",main=main,notch=T)	
		text(1:length(fc.list), par("usr")[3], labels = names(fc.list), srt = 45, adj = c(1,1), xpd = TRUE,cex=.7)
		dev.off()
	}			
			
}		
