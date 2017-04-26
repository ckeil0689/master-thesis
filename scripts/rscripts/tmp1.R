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


## path for output
path.output <- paste(sep="","results/usedForPaper/cytoscape_",date.is,"/")
system(paste(sep="","mkdir ",path.output))

## define input file names
infl.nm <- paste(path.input.combined,"/results.combine.data.Rdata",sep="")

## define global output file names
fl.annot <- paste(sep="",path.output,cut.qunt,"_","node_annot_",date.is,".txt")

# load input files
load(infl.nm)

# get deseq values for one TF (it has the same rownames for all tfs)
gns.all <- names(res[[1]][[1]][[1]][,"DEseq"])

# get list of genes to highlight in network
f.nm.rpkm <- paste(sep="",path.input,"gn_nms_to_show_in_cytoscape.txt")
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
ix.rpkm <- which(max.rpkm > GLOBAL[["min.th17.or.th0.rpkm"]] )
gns <- intersect(names(ix.sam),names(ix.rpkm))
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





# cases <- c("activation","repression","absolute")
cases <- c("activation")
for(i in 1:length(cases)){
	cs <- cases[i]
	if (i==1) {
		# get node annotation data (stays the same for each case)
		# m <- matrix(0,nr=nrow(res[[cs]][[tfs[1]]][[combine.type.1]]),nc=(length(tfs)+5))
		m <- matrix(0,nr=nrow(res[[cs]][[tfs[1]]][[combine.type.1]]),nc=(length(tfs)+4))
		# colnames(m) <- c(tfs,"SAM","th0_rpkm","th17_rpkm","highlight.gns","gn.function")
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

	combine.type.2 <- "KC"
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
	m.kc <- m
	combine.type.2 <- "DEseq"
	m <- matrix(0,nr=length(gns),nc=length(tfs))
	colnames(m) <- core.2.tfs
	rownames(m) <- gns
	for(j in 1:length(tfs)){
		if(cs=="absolute"){
			m[,j] <- res[[cs]][[tfs[j]]][[combine.type.1]][gns,combine.type.2]
		} else {
			m[,j] <- res[[cs]][[tfs[j]]][[combine.type.1]][gns,combine.type.2]
		}
	}
	# remove self regulation (we can't make it out with KO it will look like they all autoregulate)
	for(j in 1:length(tfs)){
		m[colnames(m)[j],j] <- 0
	}
	m.deseq <- m
}


## plot KC scores vs. fold change or sam scores
# sum.kc <- apply(m[,1:5],1,sum)
# ix <- sort(sum.kc,index.return=T,decreasing=T)$ix

# plot(sum.kc[ix],log2((m.annot[ix,"th17_rpkm"]+1)/(m.annot[ix,"th0_rpkm"]+1)))
pdf(paste(sep="",path.output,"kc.vs.deseq.pdf"))
for(j in 1:length(tfs)){
	tf <- core.2.tfs[j]
	ix.kc <- sort(m.kc[,tf],index.return=T,decreasing=T)$ix
	ix.deseq <- sort(m.deseq[,tf],index.return=T,decreasing=T)$ix
	ix <- ix.deseq
	plot(m.kc[ix,tf],m.deseq[ix,tf], main=tf,xlab="KC score",ylab="KO significance [log10(DEseq_pval)*sign(FC)]",ylim=c(-5,5))
	text(m.kc[ix,tf],m.deseq[ix,tf], labels=rownames(m.deseq)[ix],cex=.5)
}
dev.off()

m.deseq[ix,tf][1:10]
m.kc[ix,tf][1:10]

m.ref <- m

m.1.2 <- m
m.1.4 <- m
m.1.6 <- m
m.1.7 <- m
m.1.8 <- m

ix.1.2 <- which(abs(m)>1.2)
ix.1.4 <- which(abs(m)>1.4)
ix.1.6 <- which(abs(m)>1.6)
ix.1.7 <- which(abs(m)>1.7)
ix.1.8 <- which(abs(m)>1.8)

ix.1.2.pos <- which(m>1.2)
ix.1.4.pos <- which(m>1.4)
ix.1.6.pos <- which(m>1.6)
ix.1.7.pos <- which(m>1.7)
ix.1.8.pos <- which(m>1.8)

ix.1.2.neg <- which(m< -1.2)
ix.1.4.neg <- which(m< -1.4)
ix.1.6.neg <- which(m< -1.6)
ix.1.7.neg <- which(m< -1.7)
ix.1.8.neg <- which(m< -1.8)

m.1.2[-ix.1.2] <- 0
m.1.2[ix.1.2.pos] <- 1
m.1.2[ix.1.2.neg] <- -1

m.1.4[-ix.1.4] <- 0
m.1.4[ix.1.4.pos] <- 1
m.1.4[ix.1.4.neg] <- -1

m.1.6[-ix.1.6] <- 0
m.1.6[ix.1.6.pos] <- 1
m.1.6[ix.1.6.neg] <- -1

m.1.7[-ix.1.7] <- 0
m.1.7[ix.1.7.pos] <- 1
m.1.7[ix.1.7.neg] <- -1

m.1.8[-ix.1.8] <- 0
m.1.8[ix.1.8.pos] <- 1
m.1.8[ix.1.8.neg] <- -1

ix.1.2.zero <- which(rowSums(abs(m.1.4))==0)
ix.1.4.zero <- which(rowSums(abs(m.1.4))==0)
ix.1.6.zero <- which(rowSums(abs(m.1.6))==0)
ix.1.7.zero <- which(rowSums(abs(m.1.7))==0)
ix.1.8.zero <- which(rowSums(abs(m.1.8))==0)


# stop("AM")
orange <- rgb(255,102,0,max=255)
black <- rgb(0,0,0,max=255)
blue <- rgb(0,51,255,max=255)
pdf("~/Desktop/heatmaps_reg_trends.pdf")
heatmap(m.1.2[-ix.1.2.zero,],scale="none",col = c(orange,"black",blue),cexRow=0.051,main="1.2")

heatmap(m.1.4[-ix.1.4.zero,],scale="none",col = c(orange,"black",blue),cexRow=0.051,main="1.4")

heatmap(m.1.6[-ix.1.6.zero,],scale="none",col = c(orange,"black",blue),cexRow=0.051,main="1.6")

heatmap(m.1.7[-ix.1.7.zero,],scale="none",col = c(orange,"black",blue),cexRow=0.051,main="1.7")

heatmap(m.1.8[-ix.1.8.zero,],scale="none",col = c(orange,"black",blue),cexRow=0.051,main="1.8")
dev.off()


core.tfs <- c("BATF","IRF4","MAF","STAT3","RORC")

pdf("~/Desktop/heatmaps_reg_trends_core.pdf")
heatmap(m.1.2[-ix.1.2.zero,core.tfs],scale="none",col = c(orange,"black",blue),cexRow=0.051,main="1.2")

heatmap(m.1.4[-ix.1.4.zero,core.tfs],scale="none",col = c(orange,"black",blue),cexRow=0.051,main="1.4")

heatmap(m.1.6[-ix.1.6.zero,core.tfs],scale="none",col = c(orange,"black",blue),cexRow=0.051,main="1.6")

heatmap(m.1.7[-ix.1.7.zero,core.tfs],scale="none",col = c(orange,"black",blue),cexRow=0.051,main="1.7")

heatmap(m.1.8[-ix.1.8.zero,core.tfs],scale="none",col = c(orange,"black",blue),cexRow=0.051,main="1.8")
dev.off()






