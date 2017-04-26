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

# make output directory
if(filter.by.sam==TRUE){
	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",z.abs.cut,"_deseq_cut_",deseq.pval.cut)
	add.str.2 <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_cut_deseq_pval_",deseq.pval.cut,"_sam_",z.abs.cut)
} else {
	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",0,"_deseq_cut_",deseq.pval.cut)
	add.str.2 <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_cut_deseq_pval_",deseq.pval.cut,"_sam_",0)
}

## path for output
# path.output <- paste(sep="","results/usedForPaper/cytoscape_",date.is,"_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_cut_sam_",z.abs.cut,"/")
path.output <- paste(sep="","results/usedForPaper/cytoscape_",date.is,add.str.2,"/")
system(paste(sep="","mkdir ",path.output))

## define input file names
# infl.nm <- paste(path.input.combined,"/results_combine_data_exprsn_prcnt_diff_",prcnt.chng.cut,"_num_tfs_",num.tfs,".Rdata",sep="")
infl.nm <- paste(path.input.combined,"/results_combine_data_exprsn",add.str,".Rdata",sep="")

## define global output file names
fl.annot <- paste(sep="",path.output,cut.abs,"_","node_annot_",date.is,".txt")

# load input files
load(infl.nm)

# get deseq values for one TF (it has the same rownames for all tfs)
gns.all <- names(res[[1]][[1]][[1]][,"DEseq"])

# get list of genes to highlight in network
f.nm.rpkm <- paste(sep="",path.input,"gn_nms_to_show_in_cytoscape.txt")
f.nm.qc <- paste(sep="",path.input,"qc",add.str,"_",date.is,".pdf")

highlight.gns <- as.character(read.delim(f.nm.rpkm,header=F,sep="\n")[[1]])

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




# stop("AM")
# cases <- c("activation","repression","absolute","whole")
cases <- c("whole")
for(i in 1:length(cases)){
	cs <- cases[i]
	if (i==1) {
		# get node annotation data (stays the same for each case)
		# c.names <- c(tfs,"SAM","th0_rpkm","th17_rpkm","highlight.gns","gn.function",
		# 			"P300_th17","P300_th0","P300_th17_minus_th0")
		c.names <- c(tfs,"SAM","th0_rpkm","th17_rpkm","highlight.gns","gn.function")
					
		m <- matrix(0,nr=nrow(res[[cs]][[tfs[1]]][[combine.type.1]]),nc=length(c.names))
		colnames(m) <- c.names
		rownames(m) <- rownames(res[[cs]][[tfs[1]]][[combine.type.1]])
		# m[,c("SAM","th0_rpkm","th17_rpkm","P300_th17","P300_th0","P300_th17_minus_th0")] <- 
		# 				res[[cs]][[tfs[1]]][[combine.type.1]][,c("SAM","th0_rpkm","th17_rpkm","P300_th17","P300_th0","P300_th17_minus_th0")]
		m[,c("SAM","th0_rpkm","th17_rpkm")] <- 
						res[[cs]][[tfs[1]]][[combine.type.1]][,c("SAM","th0_rpkm","th17_rpkm")]
		m[which(rownames(m) %in% highlight.gns),"highlight.gns"] <- 1
		m[gns,"gn.function"] <- gn.function
		for(j in 1:length(tfs)){
			m[,tfs[j]] <- -1*res[[cs]][[tfs[j]]][[combine.type.1]][,"DEseq"]
		}
		gns <- gns[which(gns %in% rownames(m))]
		m.annot <- m[gns,]
	}

	m <- matrix(0,nr=length(gns),nc=length(tfs))
	colnames(m) <- core.2.tfs
	rownames(m) <- gns
	for(j in 1:length(tfs)){
		if(cs=="absolute"){
			m[,j] <- res[[cs]][[tfs[j]]][[combine.type.1]][gns,combine.type.2]*sign(res[[cs]][[tfs[j]]][[combine.type.1]][gns,"DEseq"])
		} else {
			m[,j] <- res[[cs]][[tfs[j]]][[combine.type.1]][gns,combine.type.2]
		}
	}
	# remove self regulation (we can't make it out with KO it will look like they all autoregulate)
	for(j in 1:length(tfs)){
		m[colnames(m)[j],j] <- 0
	}

	# filter based on deseq pvals
	for(j in 1:length(tfs)){
		deseq.pvals.tf <- res[[cs]][[tfs[j]]][[combine.type.1]][gns,"DEseq"]
		ix.low.pval <- which(abs(deseq.pvals.tf)<abs(log10(deseq.pval.cut)))
		m[ix.low.pval,j] <- 0
	}


	## define ouput file names
	f.eda <- paste(sep="",path.output,combine.type.2,"_",cut.abs,"_",cs,add.str.2,"_",date.is,".eda")
	f.sif <- paste(sep="",path.output,combine.type.2,"_",cut.abs,"_",cs,add.str.2,"_",date.is,".sif")
	f.core.1.txt <- paste(sep="",path.output,cut.abs,"_core_5way_",cs,"_",date.is,".txt")
	f.core.2.txt <- paste(sep="",path.output,cut.abs,"_core_11way_",cs,"_",date.is,".txt")
	f.extend.txt <- paste(sep="",path.output,cut.abs,"_extend_alltfs_",cs,"_",date.is,".txt")
	m.core.1.tfs <- m[core.1.tfs,core.1.tfs]
	m.core.2.tfs <- m[core.2.tfs,core.2.tfs]
	m.all.tfs <- m[which(rownames(m) %in% tfNames),core.2.tfs]
	# m.cut <- quantile(m[which(m!=0)],cut.qunt)
	# m.cut <- quantile(m,cut.qunt)
	m.cut <- cut.abs
	cat("writing network files: ",paste(sep="",cs,"_",date.is,".eda"), "and",paste(sep="",cs,"_",date.is,".eda"),"\n")
	cat(sep="","m.cut = ",m.cut,"\n")
	if(cs=="repression"){
		pos.edge="negative"
	}else {
		pos.edge="positive"
	}
	write.cyto.net(m=m,m.cut=m.cut,file.eda=f.eda,file.sif=f.sif,append=FALSE,pos.edge=pos.edge)
}
# write node annotation file
cat("gene_id",colnames(m.annot),sep="\t",append=FALSE,file=fl.annot)
cat("\n",append=TRUE,file=fl.annot)
write.table(sep="\t",m.annot,quote=FALSE,append=TRUE,col.names=FALSE,file=fl.annot)


# write KC scores annotation file
cases <- c("activation","repression","absolute","whole")
for(i in 1:length(cases)){
	cs <- cases[i]
	m <- matrix(0,nr=nrow(res[[cs]][[1]][[1]]),nc=length(res[[cs]]))
	colnames(m) <- names(res[[cs]])
	rownames(m) <- rownames(res[[cs]][[1]][[1]])
	for(j in 1:length(res[[cs]])){
		tf <- names(res[[cs]])[j]
		m[,j] <- res[[cs]][[tf]][[combine.type.1]][,combine.type.2]
		# filter low pval inters if lower than pval cut
		deseq.pvals.tf <- res[[cs]][[tf]][[combine.type.1]][,"DEseq"]
		ix.low.pval <- which(abs(deseq.pvals.tf)<abs(log10(deseq.pval.cut)))
		m[ix.low.pval,j] <- 0
	}
	f.nm <- paste(sep="",path.output,combine.type.2,"_",cs,"_",add.str.2,"_",date.is,".xls")	
	cat("gene_id",colnames(m),sep="\t",append=FALSE,file=f.nm)
	cat("\n",append=TRUE,file=f.nm)
	write.table(sep="\t",m,quote=FALSE,append=TRUE,col.names=FALSE,file=f.nm)	
}


if(plot.qc==T){
	f.nm.qc <- paste(sep="",path.output,"qc",add.str.2,"_",date.is,".pdf")
	pdf(f.nm.qc)
	case <- "whole"
	chip <- combine.type.1
	comb <- combine.type.2
	x.fc <- numeric()
	x.deseq <- numeric()
	x.comb.score <- numeric()
	x.fc.all <- numeric()
	x.deseq.all <- numeric()
	x.comb.score.all <- numeric()
	for(i in 1:length(tfs)){
		tf <- tfs[i]
		m <- res[[case]][[tf]][[chip]][gns,]
		# filter based on deseq pvals
		deseq.pvals.tf <- m[,"DEseq"]
		ix.low.pval <- which(abs(deseq.pvals.tf)<abs(log10(deseq.pval.cut)))
		m[ix.low.pval,comb] <- 0
		ix <- which(abs(m[,comb])>cut.abs)
		x.fc <- c(x.fc,m[ix,"FC"])
		x.deseq <- c(x.deseq,m[ix,"DEseq"])
		x.comb.score <- c(x.comb.score,m[ix,comb])
		x.fc.all <- c(x.fc.all,m[,"FC"])
		x.deseq.all <- c(x.deseq.all,m[,"DEseq"])
		x.comb.score.all <- c(x.comb.score.all,m[,comb])
	}

	plot(x.comb.score.all,x.deseq.all,ylim=c(-3,3),pch=20,cex=0.3,col="gray")
	points(x.comb.score,x.deseq,col="blue",pch=20,cex=0.3)
	plot(x.comb.score.all,x.fc.all,ylim=c(-3,3),pch=20,cex=0.3,col="gray")
	points(x.comb.score,x.fc,col="blue",pch=20,cex=0.3)


	plot(x.deseq.all,x.fc.all,ylim=c(-3,3),xlim=c(-3,3),pch=20,cex=0.3,col="gray")
	points(x.deseq,x.fc,col="blue",pch=20,cex=0.3)

	hist(10^-abs(x.deseq.all[which(x.deseq.all!=0)]),nclass=1000)
	dev.off()
}
















