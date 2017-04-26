##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# R code to create cyto-nets (after running combineData.R)
# output is:
#   core.net syf file
#   core.net eda file
#   dataset tab delim file with accompanying extra information helpful for visulization
source("r_scripts/th17/used_for_paper/util.R")

## paths
path.input <- "input/th17/used_for_paper/"
path.input.combined <- paste(sep="","results/combinedAnalysis/",date.combine.data.run,"/")

if(filter.by.sam==TRUE){
#	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",z.abs.cut)
   add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",z.abs.cut,"_deseq_cut_",deseq.pval.cut)
	add.str.2 <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_cut_deseq_pval_",deseq.pval.cut,"_sam_",z.abs.cut)
} else {
#	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",0)
   add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",0,"_deseq_cut_",deseq.pval.cut)
	add.str.2 <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_cut_deseq_pval_",deseq.pval.cut,"_sam_",0)
}
add.str.3 <- combine.type.2
## path for output
path.output <- paste(sep="","results/usedForPaper/cytoscape_",date.is,add.str.2,"/")
# path.output <- paste(sep="","results/usedForPaper/cytoscape_",date.is,"/")
system(paste(sep="","mkdir ",path.output,add.str.3,"_heatmpas/"))


## define input file names
infl.nm <- paste(path.input.combined,"/results_combine_data_exprsn",add.str,".Rdata",sep="")
# infl.nm <- paste(path.input.combined,"/results.combine.data.Rdata",sep="")
fl.combos <- paste(path.input,"heatmap_combos_Jun_21_2012.txt",sep="")

## define global output file names
# fl.annot <- paste(sep="",path.output,cut.abs,"_","node_annot_",date.is,".txt")

# load input files
load(infl.nm)
# get deseq values for one TF (it has the same rownames for all tfs)
gns.all <- names(res[[1]][[1]][[1]][,"DEseq"])

# get list of genes to highlight in network
f.nm.rpkm <- paste(sep="",path.input,"gn_nms_to_show_in_cytoscape.txt")
highlight.gns <- as.character(read.delim(f.nm.rpkm,header=F,sep="\n")[[1]])

# get combinations to create heatmap for
combos <- read.delim(fl.combos,header=F,sep="\n",stringsAsFactors=F)

##################################################
# get rpkm values
f.nm.rpkm <- paste(sep="",path.input,"ranseqDatasetNoQuartNorm.RData") # which experiment to compare

## read rpkm data
load(f.nm.rpkm)
d <- rnaseq.complete.matrix
expt.names <-colnames(d)
th17.wt.ix <- grep("th17_ss_48hr_.{1,6}_wt",colnames(d),value=T,perl=T)
th17.wt.ix <- c(th17.wt.ix,"SL1858_th17_ts_48hr#1")
th17.mean <- apply(d[,th17.wt.ix],1,mean)
th0.wt.ix <- grep("th0_ss_48hr_.{1,6}_wt",colnames(d),value=T,perl=T)
th0.wt.ix <- c(th0.wt.ix,"SL2673_th0_ts_48hr#1")
th0.mean <- apply(d[,th0.wt.ix],1,mean)
rpkm.mean <- pmax(th17.mean,th0.mean)


# th17
m <- numeric(length=length(gns.all))
names(m) <- gns.all
ix <- which(gns.all %in% names(th17.mean))
m[gns.all[ix]] <- th17.mean[gns.all[ix]]
th17.rpkm <- m
# th0
m <- numeric(length=length(gns.all))
names(m) <- gns.all
ix <- which(gns.all %in% names(th0.mean))
m[gns.all[ix]] <- th0.mean[gns.all[ix]]
th0.rpkm <- m
##################################################
##################################################
# get SAM differential expression z scores
x <- as.matrix(read.table(paste(path.input,"samTh17VsTh0Zscores.xls",sep="")))
x.sam <- x[,"Score_d"]
##################################################
# keep genes with z score above cutoff and genes that are expressed in either th17 or th0
ix.sam <- which(abs(x.sam) > GLOBAL[["z.abs.cut"]] )
# AND expressed in either th17 or th0
max.rpkm <- pmax(th0.rpkm,th17.rpkm)
ix.rpkm <- which(max.rpkm > GLOBAL[["min.th17.or.th0.rpkm"]] )
# gns <- intersect(names(ix.sam),names(ix.rpkm))
if(filter.by.sam==T){
	ix.rpkm <- which(max.rpkm > GLOBAL[["min.th17.or.th0.rpkm"]] )
	gns <- intersect(names(ix.sam),names(ix.rpkm))
} else {
	gns <- gns.all
}
##################################################
#33333333333333333333333333333333
# Now write additional info files
#33333333333333333333333333333333
# get all tf names
load(paste(sep="",path.input,"infDataStructures/immgen/z_abs_cut_2.5/tfNames.RData"))
tfNames <- tfNames[-grep("hist",tfNames,ignore.case=T)]
## adding cell surface markers taken as a union of gene names in our dataset from three lists taken from human biomart, mouse biomart and GO mouse for
## external side of plasma membrane GO:0009897
srfc.gns.go <- toupper(unique(as.character(read.table(paste(sep="",path.input,"cytoscape/plsma_mmbrn_GO_0009897_Feb_21_2011.txt"),header=F,sep="\t")$V1)))
srfc.gns.biomart.mm <- toupper(unique(as.character(read.table(paste(sep="",path.input,"cytoscape/mart_export_plsma_mmbrn_mm_Feb_21_2011.txt"),header=F,sep="\t")$V1)))
srfc.gns.biomart.hg <- toupper(unique(as.character(read.table(paste(sep="",path.input,"cytoscape/mart_export_plsma_mmbrn_hg_Feb_21_2011.txt"),header=F,sep="\t")$V1)))
srfc.gns <- union(c(srfc.gns.biomart.mm,srfc.gns.biomart.hg),srfc.gns.go)
## catalytic activity GO:0003824
enzyme.gns <- unique(toupper(unique(as.character(read.table(paste(sep="",path.input,"cytoscape/catalytic_activity_GO_0003824_Mar_12_2011.txt"),header=F,sep="\t")$V1))))
## protein kinase activity GO:0004672
kinase.gns <- unique(toupper(unique(as.character(read.table(paste(sep="",path.input,"cytoscape/protein_kinase_activity_GO_0004672_Mar_12_2011.txt"),header=F,sep="\t")$V1))))

# tf=1, surface marker=2, kinase=3, enzyme=4, other=0 (if tf has priority than enzyme than srfc mrkr)
gn.function <- character(length=length(gns))
for (i in 1:length(gns)){
  if (gns[i] %in% tfNames){
    gn.function[i] <- "tf"
  } else if (gns[i] %in% srfc.gns){
    gn.function[i] <- "srfc.mrkr"
  } else if (gns[i] %in% kinase.gns){
    gn.function[i] <- "kinase"
  } else if (gns[i] %in% enzyme.gns){
    gn.function[i] <- "enzyme"
  } else {
    gn.function[i] <- "other"
  }
}





# cases <- c("activation","repression","absolute",whole)
cases <- c("whole")
for(i in 1:length(cases)){
	cs <- cases[i]
	cat(sep="","working on heatmaps for ",cs,"\n")
	if (i==1) {
		# get node annotation data (stays the same for each case)
		m <- matrix(0,nr=nrow(res[[cs]][[tfs[1]]][[combine.type.1]]),nc=(length(tfs)+4))
		colnames(m) <- c(tfs,"SAM","th0_rpkm","th17_rpkm","highlight.gns")
		rownames(m) <- rownames(res[[cs]][[tfs[1]]][[combine.type.1]])
		m[,c("SAM","th0_rpkm","th17_rpkm")] <- res[[cs]][[tfs[1]]][[combine.type.1]][,c("SAM","th0_rpkm","th17_rpkm")]
		m[which(rownames(m) %in% highlight.gns),"highlight.gns"] <- 1
		# m[gns,"gn.function"] <- gn.function
		for(j in 1:length(tfs)){
			m[,tfs[j]] <- -1*res[[cs]][[tfs[j]]][[combine.type.1]][,"DEseq"]
		}
		gns <- gns[which(gns %in% rownames(m))]
		m.annot <- m[gns,]
	}

	combine.type.2 <- combine.type.2
	m <- matrix(0,nr=length(gns),nc=length(tfs))
	colnames(m) <- core.2.tfs
	rownames(m) <- gns
	for(j in 1:length(tfs)){
		if(cs=="absolute"){
			if(tfs[j]!="P300"){
				m[,j] <- res[[cs]][[tfs[j]]][[combine.type.1]][gns,combine.type.2]*sign(res[[cs]][[tfs[j]]][[combine.type.1]][gns,"DEseq"])
			} else {
				m[,j] <- 2*res[[cs]][[tfs[j]]][[combine.type.1]][gns,combine.type.2]*sign(res[[cs]][[tfs[j]]][[combine.type.1]][gns,"SAM"])
			}
		} else {
			m[,j] <- res[[cs]][[tfs[j]]][[combine.type.1]][gns,combine.type.2]
		}
	}
	# remove self regulation (we can't make it out with KO it will look like they all autoregulate)
	for(j in 1:length(tfs)){
		if(tfs[j]!="P300"){
			m[colnames(m)[j],j] <- 0
		}
	}
	m.kc <- m
	# get matrix for DEseq pvals	
	combine.type.2 <- "DEseq"
	m <- matrix(0,nr=length(gns),nc=length(tfs))
	colnames(m) <- core.2.tfs
	rownames(m) <- gns
	for(j in 1:length(tfs)){
		m[,j] <- res[[cs]][[tfs[j]]][[combine.type.1]][gns,combine.type.2]
	}
	# remove self regulation (we can't make it out with KO it will look like they all autoregulate)
	for(j in 1:length(tfs)){
		if(tfs[j]!="P300"){
			m[colnames(m)[j],j] <- 0
		}
	}
	m.deseq <- m
	# get matrix for FC
	# get matrix for DEseq pvals	
	combine.type.2 <- "FC"
	m <- matrix(0,nr=length(gns),nc=length(tfs))
	colnames(m) <- core.2.tfs
	rownames(m) <- gns
	for(j in 1:length(tfs)){
		m[,j] <- res[[cs]][[tfs[j]]][[combine.type.1]][gns,combine.type.2]
	}
	# remove self regulation (we can't make it out with KO it will look like they all autoregulate)
	for(j in 1:length(tfs)){
		if(tfs[j]!="P300"){
			m[colnames(m)[j],j] <- 0
		}
	}
	m.fc <- m

													##### plot heatmap for combinatorics argument #####
													
													
	#stop("AM")
	library(pheatmap)
	library(gplots)
	col.to.white.gradient <- function(n,r.st,g.st,b.st){
    if (n <   0){ n = 0}
    if (n > 255){ n = 255}
		rgb.grad.vec <- numeric(length(n))
		red.step <- floor((255-r.st)/n)
		green.step <- floor((255-g.st)/n)
		blue.step <- floor((255-b.st)/n)				
		for(i in 1:n){
	    r = r.st+red.step*i
	    g = g.st+green.step*i
	    b = b.st+blue.step*i
			rgb.grad.vec[i] <- rgb (r,g,b,max = 255)
		}
    return (rgb.grad.vec)
	}
	
	# filter based on deseq pvals
	ix <- which(abs(m.deseq)<abs(log10(deseq.pval.cut)))

	m <- m.kc
	m.cut <- m.kc
	m.cut[ix] <- 0
	ix <- which(abs(m.cut)<cut.abs)
	m.cut[ix] <- 0
	m.num.reg <- apply(m.cut,1,function(i) length(which(abs(i[1:5])>cut.abs)))
	x.way <- 4
	ix.show <- which(m.num.reg>=x.way)
	n <- 10
	or.col <- c(255,102,0)
	orange.to.white <- col.to.white.gradient(n=n,r.st=or.col[1],g.st=or.col[2],b.st=or.col[3])
	bl.col <- c(0,51,255)
	blue.to.white <- col.to.white.gradient(n=n,r.st=bl.col[1],g.st=bl.col[2],b.st=bl.col[3])
	col <- c(orange.to.white,rev(blue.to.white))
	fl.nm <- paste(sep="",path.output,add.str.3,"_heatmpas/",x.way,"way_or_more_",cs,"_heatmap",add.str.2,"_",date.is,".pdf")
	pdf(fl.nm)	

		tot <- n*4
		blue.col <- rgb(bl.col[1],bl.col[2],bl.col[3],max=255) 
		orange.col <- rgb(or.col[1],or.col[2],or.col[3],max=255) 
		
		# cluster on genes
		d <- dist(as.matrix(m.cut[ix.show,]))
		hc <- hclust(d)
		gn.order <- hc$labels[hc$order]
		# define the layout for the heatmap
		# layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(4,1),heights=c(1,1));layout.show(2)
		
		# plot regulation heatmap
		pheatmap(m.cut[gn.order,],scale="none",color=colorRampPalette(c(orange.col,rep("white",4),
		blue.col))(tot),breaks=seq(-2,2,by=0.1),cellwidth = 12, cellheight = 4,treeheight_col=0,
		cex=0.3,cluster_cols=T,main=paste(sep="","kc=",cut.abs," pval=",deseq.pval.cut),cluster_rows=F)

		
		# make sam vals symmetric around zero
		w <- x.sam[gn.order]
		eff.max.val <- min( abs( c( max(w), min(w) ) ) )
		ix <- which(w < -eff.max.val)
		if(length(ix)){w[ix] <- -eff.max.val}
		ix <- which(w > eff.max.val)
		if(length(ix.sam)){w[ix] <- eff.max.val}		
		# plot sam heatmap		
		pheatmap(w[gn.order],scale="none",cellheight = 4,main="SAM",
			colorRampPalette(c(orange.col,"white",blue.col))(tot),	
		cex=0.3,main=paste(sep="","kc=",cut.abs," pval=",deseq.pval.cut),cluster_rows=F,cluster_cols=F)

		# plot fc heatmap
		x.fc <- log2((th17.rpkm+1)/(th0.rpkm+1))
		w <- x.fc[gn.order]
		eff.max.val <- min( abs( c( max(w), min(w) ) ) )
		ix <- which(w < -eff.max.val)
		if(length(ix)){w[ix] <- -eff.max.val}
		ix <- which(w > eff.max.val)
		if(length(ix)){w[ix] <- eff.max.val}
		pheatmap(w[gn.order],scale="none",cellheight = 4,main="FC",
			colorRampPalette(c(orange.col,"white",blue.col))(tot),	
		cex=0.3,main=paste(sep="","kc=",cut.abs," pval=",deseq.pval.cut),cluster_rows=F,cluster_cols=F)

		# cluster on genes
		d <- dist(as.matrix(m[ix.show,]))
		hc <- hclust(d)
		gn.order <- hc$labels[hc$order]
		
		tot <- n*4
		blue.col <- rgb(bl.col[1],bl.col[2],bl.col[3],max=255) 
		orange.col <- rgb(or.col[1],or.col[2],or.col[3],max=255) 
		pheatmap(m[gn.order,],scale="none",color=colorRampPalette(c(orange.col,rep("white",4),
		blue.col))(tot),breaks=seq(-2,2,by=0.1),cellwidth = 12, cellheight = 4,treeheight_col=0,
		cex=0.4,cluster_cols=T,main="without_pval_cut",cluster_rows=F)

	dev.off()

													##### Done! plot heatmap for combinatorics argument #####	


													##### plot heatmap for select gene groups #####														
	f.nm <- paste(sep="",path.input,"Signature_list_for_heatmap.txt")
	x <- read.table(f.nm,sep="\t",header=T,as.is=T)

	x.fc <- log2((th17.rpkm+1)/(th0.rpkm+1))
	eff.max.val <- min( abs( c( max(x.fc), min(x.fc) ) ) )
	ix <- which(x.fc < -eff.max.val)
	if(length(ix)){x.fc[ix] <- -eff.max.val}
	ix <- which(x.fc > eff.max.val)
	if(length(ix)){x.fc[ix] <- eff.max.val}

	fl.nm <- paste(sep="",path.output,add.str.3,"_heatmpas/selected_gns_",cs,"_heatmap_",add.str.2,"_",date.is,".pdf")
	pdf(fl.nm)
	m <- m.cut
	for(j in 2:ncol(x)){
		nm <- colnames(x)[j]
		gns <- toupper(x$gene[which(x[,j]==1)])
		d <- dist(as.matrix(m[gns,]))
		hc <- hclust(d)
		gn.order <- hc$labels[hc$order]
		
		if(length(gns)>1){
			pheatmap(m[gn.order,],scale="none",color=colorRampPalette(c(orange.col,rep("white",4),
			blue.col))(tot),breaks=seq(-2,2,by=0.1),cellwidth = 12, cellheight = 4,treeheight_col=0,
			cex=0.35,cluster_cols=F,cluster_rows=F)			
		}
		
		pheatmap(x.fc[gn.order],scale="none",cellheight = 4,main="FC",breaks=seq(-eff.max.val,eff.max.val,by=(2*eff.max.val)/tot ),
			colorRampPalette(c(orange.col,"white",blue.col))(tot),	
		cex=0.3,main=paste(sep="","kc=",cut.abs," pval=",deseq.pval.cut),cluster_rows=F,cluster_cols=F)
		
	}
	dev.off()							
													
	if((show.perm==T) & (cs=="absolute" | cs=="whole")){
		# stop("AM")
		m.ref <- m.kc
		for(j in 1:nrow(combos)){
			combo <- toupper(strsplit(combos[j,],"::")[[1]])
			if( any(!(combo %in% colnames(m.ref))) ){
				next
			}
			m <- m.ref[,combo]
			m.1.2 <- m
			m.1.4 <- m
			m.1.5 <- m
			m.1.6 <- m
			m.1.7 <- m
			m.1.8 <- m

			ix.1.2 <- which(abs(m)>1.2)
			ix.1.4 <- which(abs(m)>1.4)
			ix.1.5 <- which(abs(m)>1.5)
			ix.1.6 <- which(abs(m)>1.6)
			ix.1.7 <- which(abs(m)>1.7)
			ix.1.8 <- which(abs(m)>1.8)

			ix.1.2.pos <- which(m>1.2)
			ix.1.4.pos <- which(m>1.4)
			ix.1.5.pos <- which(m>1.5)			
			ix.1.6.pos <- which(m>1.6)
			ix.1.7.pos <- which(m>1.7)
			ix.1.8.pos <- which(m>1.8)

			ix.1.2.neg <- which(m< -1.2)
			ix.1.4.neg <- which(m< -1.4)
			ix.1.5.neg <- which(m< -1.5)			
			ix.1.6.neg <- which(m< -1.6)
			ix.1.7.neg <- which(m< -1.7)
			ix.1.8.neg <- which(m< -1.8)

			if(make.binary.heatmap == T){
				m.1.2[-ix.1.2] <- 0
				m.1.2[ix.1.2.pos] <- 1
				m.1.2[ix.1.2.neg] <- -1

				m.1.4[-ix.1.4] <- 0
				m.1.4[ix.1.4.pos] <- 1
				m.1.4[ix.1.4.neg] <- -1

				m.1.5[-ix.1.5] <- 0
				m.1.5[ix.1.5.pos] <- 1
				m.1.5[ix.1.5.neg] <- -1

				m.1.6[-ix.1.6] <- 0
				m.1.6[ix.1.6.pos] <- 1
				m.1.6[ix.1.6.neg] <- -1

				m.1.7[-ix.1.7] <- 0
				m.1.7[ix.1.7.pos] <- 1
				m.1.7[ix.1.7.neg] <- -1

				m.1.8[-ix.1.8] <- 0
				m.1.8[ix.1.8.pos] <- 1
				m.1.8[ix.1.8.neg] <- -1
		    } else {
				m.1.2[-ix.1.2] <- 0
				m.1.4[-ix.1.4] <- 0
				m.1.5[-ix.1.5] <- 0				
				m.1.6[-ix.1.6] <- 0
				m.1.7[-ix.1.7] <- 0
				m.1.8[-ix.1.8] <- 0
			}
			ix.0.zero <- which(rowSums(abs(m))==0)
			ix.1.2.zero <- which(rowSums(abs(m.1.2))==0)
			ix.1.4.zero <- which(rowSums(abs(m.1.4))==0)
			ix.1.5.zero <- which(rowSums(abs(m.1.5))==0)
			ix.1.6.zero <- which(rowSums(abs(m.1.6))==0)
			ix.1.7.zero <- which(rowSums(abs(m.1.7))==0)
			ix.1.8.zero <- which(rowSums(abs(m.1.8))==0)


			# set colors
			# create a color wheel to decide on params for hsv
			# hue <- seq(0.0, 1.0, by=1/40) 
			# x <- seq(0.0, 1.0, by=1/40) 
			# 
			# pie(rep(1,length(x)), 
			# labels=formatC(hue, digits=3, 
			# format="f"), 
			# cex=0.75, 
			# col=hsv(0.075, 1.0, x), 
			# radius=1.0, 
			# main="HSV (S=1, V=1)" )
			# x <- seq(0.0, 1.0, by=1/40) 
			# 		pie(rep(1,length(x)), 
			# 		labels=formatC(hue, digits=3, 
			# 		format="f"), 
			# 		cex=0.75, 
			# 		col=hsv(0.625, 1.0, x), 
			# 		radius=1.0, 
			# 		main="HSV (S=1, V=1)" )
			# test color wheel
			x <- seq(0.0, 1.0, by=1/1000) 
			orange.hsv  <-  0.075
			orange.2.black  <-  hsv(orange.hsv,1.0,rev(x))
			blue.hsv  <-  0.65
			black.2.blue  <-  hsv(blue.hsv,1.0,x)
			col <- c(orange.2.black,black.2.blue)
			orange.2.black.staurated.top = c(orange.2.black,rep(orange.2.black[length(x)],length(x)/3))
			black.2.blue.staurated.top = c(black.2.blue,rep(black.2.blue[length(x)],length(x)/3))
			col <- c(orange.2.black.staurated.top,black.2.blue.staurated.top)
			
			# system(paste(sep="","mkdir ",path.output,"heatmpas/"))
			if(filter.by.sam==T){
				fl.nm <- paste(sep="",path.output,add.str.3,"_heatmpas/",cs,"_sam_filtered_heatmap_",paste(sep="",sort(combo),collapse="_"),"_",date.is,".pdf")
			} else {
				fl.nm <- paste(sep="",path.output,add.str.3,"_heatmpas/",cs,"_heatmap_",paste(sep="",sort(combo),collapse="_"),"_",date.is,".pdf")
			}
			# fl.nm <- paste(sep="",path.output,add.str.3,"_heatmpas/",cs,"_sam_filtered_heatmap_",paste(sep="",sort(combo),collapse="_"),"_",date.is,".pdf")
			
			if(length(combo)>3){
				cexRow=0.051
			} else {
				cexRow=0.2
			}
			pdf(fl.nm)
			if(filter.by.sam==T){	# when filtering by sam heatmaps are much smaller so plot many options	
				heatmap(m[-ix.0.zero,],scale="none",col = col,cexRow=cexRow,main="no filter (all data with any entry not equal zero)")

				heatmap(m.1.2[-ix.1.2.zero,],scale="none",col = col,cexRow=cexRow,main="1.2")

				# heatmap(m.1.4[-ix.1.4.zero,],scale="none",col = col,cexRow=cexRow,main="1.4")

				heatmap(m.1.5[-ix.1.5.zero,],scale="none",col = col,cexRow=cexRow,main="1.5")

				# heatmap(m.1.6[-ix.1.6.zero,],scale="none",col = col,cexRow=cexRow,main="1.6")

				heatmap(m.1.7[-ix.1.7.zero,],scale="none",col = col,cexRow=cexRow,main="1.7")
			
				# heatmap(m.1.8[-ix.1.8.zero,],scale="none",col = col,cexRow=cexRow,main="1.8")
			} else {
				heatmap(m.1.5[-ix.1.5.zero,],scale="none",col = col,cexRow=cexRow,main="1.5")				
			}
			dev.off()

		}
	}
}
