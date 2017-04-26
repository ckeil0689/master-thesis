##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
cat("\n- setting global variables\n")
rm(list=ls())
GLOBAL <- list()
#GLOBAL[["run.these.steps"]] <- c("12.1.1","14.1.1","14.2.1") # for a complete run c(1:last_step)
GLOBAL[["run.these.steps"]] <- c(12.2)
# stamp the date on this run
x <- unlist(strsplit(date()," +",perl=TRUE))
GLOBAL[["date.is"]] <- paste(x[2],x[3],x[5],sep="_")
# data set size: take genes with abs(zscore) based on diff expression btwn th17_48hr and th0_48hr
GLOBAL[["z.abs.cut"]] <- 2.5 #zscore
# Only used to determine rnaseq data set size
# GLOBAL[["median.abs.cut"]] <- 3 #rpkm
GLOBAL[["min.th17.or.th0.rpkm"]] <- 3 #rpkm
# The number of bootstraps for inferelator rnaseq and immgen
GLOBAL[["num.boots.rnaseq"]] <- 200
GLOBAL[["num.boots.immgen"]] <- 200
GLOBAL[["num.perms.enrichment.analysis"]] <- 10
# do we want to output the sequences around the peaks?
GLOBAL[["get.sequence"]] <- TRUE
# what is the minimum distance between peaks to cluster them together
## GLOBAL[["tfs.min.dist"]] <- c(50,100,150,200,250,300,500,1000)
GLOBAL[["tfs.min.dist"]] <- c(100)
# multi TF binding sites distance (up/downstrem) from TSS for gene MTBS mapping
GLOBAL[["mtbs.tss.dist"]] <- 5000
GLOBAL[["mm9.effective.genome.size"]] <- 1.87e9

# genes that we want in as a rule even if they have low median expression or z.score that kicks them out
GLOBAL[["known.tfs"]] <- c("BATF","MAF","IRF4","STAT3","RORC")
GLOBAL[["known.gns"]] <- toupper(c("Il17a","Il17f","Il23r","Il22","Il21","Foxp3","Rora"))
GLOBAL[["use.multicore"]] <- TRUE
if(GLOBAL[["use.multicore"]]==TRUE){
  library(multicore)
}

	rm(list=ls()[-which(ls() %in% c("GLOBAL","iter"))])
	date.is <- GLOBAL[["date.is"]]
	num.boots.immgen <- GLOBAL[["num.boots.immgen"]]
	num.boots.rnaseq <- GLOBAL[["num.boots.rnaseq"]]
	# date.inf.run <- "Sep_16_2011"
	date.inf.run <- "Aug_6_2012"		
	# date.deseq.run <- "Jun_18_2012"
	date.deseq.run <- "Aug_2_2012"
	chip.integration <- "genewide_pois_model_pval"
	# chip.integration <- "prox_pois_model_pval"
	z.abs.cut.old <- GLOBAL[["z.abs.cut"]]
	z.abs.cut <- 0
	prcnt.chng.cut <- 0;filter.by.sam <- F;num.tfs <- 5;deseq.pval.cut <- 1;
	combine.file.map <- paste(sep="","combine_kcri_core_",num.tfs,"_tfs.xls") # or combine_kcri_core_5tfs.xls



source("r_scripts/th17/used_for_paper/util.R")
# setup paths
path.input <- "input/th17/used_for_paper/"
path.input.immgen <- "input/th17/used_for_paper/infResults/immgen/"
path.input.rnaseq <- "input/th17/used_for_paper/infResults/rnaseq/"
path.input.chip <- "/data/th17/data/map_ChIP_peaks_to_genes/"
path.input.deseq <- paste(sep="","input/th17/used_for_paper/DEseq/",date.deseq.run,"/")
path.output <- paste(sep="","results/combinedAnalysis/",date.is,"/")
path.input.immgen.inf.data <-  paste(sep="","input/th17/used_for_paper/infDataStructures/",date.inf.run,"/immgen/z_abs_cut_",z.abs.cut.old,"/")
path.input.rnaseq.inf.data <-  paste(sep="","input/th17/used_for_paper/infDataStructures/",date.inf.run,"/rnaseq/z_abs_cut_",z.abs.cut.old,"/")

# make output directory
# if(filter.by.sam==TRUE){
add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",z.abs.cut,"_deseq_cut_",deseq.pval.cut)
# } else {
# 	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",0,"_deseq_cut_",deseq.pval.cut)
# }
system(paste(sep="","mkdir ", path.output))
system(paste(sep="","mkdir ", path.output, "activation",add.str,"/"))
system(paste(sep="","mkdir ", path.output, "repression",add.str,"/"))
# system(paste(sep="","mkdir ", path.output, "absolute",add.str,"/"))
system(paste(sep="","mkdir ", path.output, "whole",add.str,"/"))


# setup file names
file.nm.map <- paste(sep="",path.input,"combinedAnalysis/",combine.file.map)
file.nm.immgen <- paste(sep="","clrinf_cnfdnc_zcut_",z.abs.cut.old,"_Nb_",num.boots.immgen,"_",  date.inf.run,".xls")
file.nm.rnaseq <- paste(sep="","clrinf_cnfdnc_zcut_",z.abs.cut.old,"_Nb_",num.boots.rnaseq,"_",date.inf.run,".xls")
file.nm.ko <- paste(sep="","DEseq_pval_signed_",date.deseq.run,".xls")
file.nm.fc <- paste(sep="","DEseq_log2fc_",date.deseq.run,".xls")
file.nm.prcnt.chng <- paste(sep="","DEseq_prcnt_chng_",date.deseq.run,".xls")
file.nm.sam <- "samTh17VsTh0Zscores.xls"

# load mapping file
f.nm <- file.nm.map 
cat("loading SAM scores file: ",f.nm,"\n")
# scores.sam <- as.matrix(read.table(f.nm,sep="\t"))
map <- read.table(f.nm,sep="\t",colClasses = "character",header=T)
core.tfs <- toupper(map$tf)
other.tfs <- character() # other tfs will be determined based on tfs in immgen and rnaseq

# load scores
#load sam diff expression
f.nm <- paste(sep="",path.input,file.nm.sam)
cat("loading SAM scores file: ",f.nm,"\n")
# scores.sam <- as.matrix(read.table(f.nm,sep="\t"))
scores.sam <- read.table(f.nm,sep="\t",colClasses = "character")
c.nms.sam <- as.character(scores.sam[1,1:(dim(scores.sam)[2]-1)])
scores.sam <- scores.sam[-1,]
r.nms.sam <- scores.sam[,1]
scores.sam <- scores.sam[,-1]
m <- matrix(0,nr=nrow(scores.sam),nc=ncol(scores.sam))
rownames(m) <- r.nms.sam
colnames(m) <- c.nms.sam
for(j in 1:ncol(m)){
	m[,j] <- as.numeric(scores.sam[,j])
}
scores.sam <- m

# get experiments mapping
f.nm <- file.nm.map
x <- as.matrix(read.delim(f.nm,sep="\t"),header=T)
#get immgen
f.nm <- paste(sep="",path.input.immgen,file.nm.immgen)
cat("loading immgen scores file: ",f.nm,"\n")
scores.immgen <- as.matrix(read.table(f.nm,sep="\t"))
#get rnaseq
f.nm <- paste(sep="",path.input.rnaseq,file.nm.rnaseq)
cat("loading rnaseq scores file: ",f.nm,"\n")
scores.rnaseq <- as.matrix(read.table(f.nm,sep="\t"))
#######################################################
## put sign on confidence scores from immgen and ranseq
load(paste(sep="",path.input.immgen.inf.data,"ratios.RData"))
tfs.immgen <- colnames(scores.immgen)
sign.cor.mat.immgen <- sign(cor(t(ratios))[,tfs.immgen])
load(paste(sep="",path.input.rnaseq.inf.data,"ratios.RData"))
tfs.rnaseq <- colnames(scores.rnaseq)
sign.cor.mat.rnaseq <- sign(cor(t(ratios))[,tfs.rnaseq])
scores.immgen <- scores.immgen*sign(sign.cor.mat.immgen)
scores.rnaseq <- scores.rnaseq*sign(sign.cor.mat.rnaseq)
# other.tfs <- unique(tfs.immgen,tfs.rnaseq) # other tfs will be determined based on tfs in immgen and rnaseq
# other.tfs <- other.tfs[-which(other.tfs %in% core.tfs)]
rm(tfs.immgen,tfs.rnaseq,ratios)

#######################################################
#get KO
f.nm <- paste(sep="",path.input.deseq,file.nm.ko)
cat("loading KO scores file: ",f.nm,"\n")
scores.ko <- as.matrix(read.table(f.nm,sep="\t"))
#get FC
f.nm <- paste(sep="",path.input.deseq,file.nm.fc)
cat("loading FC scores file: ",f.nm,"\n")
scores.fc <- as.matrix(read.table(f.nm,sep="\t"))
#get prcnt change from wt
f.nm <- paste(sep="",path.input.deseq,file.nm.prcnt.chng)
cat("loading FC scores file: ",f.nm,"\n")
scores.prcnt.chng <- as.matrix(read.table(f.nm,sep="\t"))
# get relevant chip_experiments
scores.chip.th17 <- matrix(0,nr=nrow(scores.ko),nc=nrow(x))
rownames(scores.chip.th17) <- rownames(scores.ko)
colnames(scores.chip.th17) <- toupper(x[,"tf"])

for(i in 1:nrow(x)){
	tf <- toupper(x[i,"tf"])
	chip.th17.nm <- x[i,"chip_th17"]
	# read chip th17
	# get th17 chip scores
	if(!is.na(chip.th17.nm)){	
		cat("loading chip data for",tf,"\n")
		cmd.line <- paste(sep="","scp am2654@bot.bio.nyu.edu:",path.input.chip,chip.th17.nm,"/",chip.th17.nm,"_genes.xls" , " tmp/")
		system(cmd.line)
		tmp <- as.matrix(read.delim(paste(sep="","tmp/",chip.th17.nm,"_genes.xls"),sep="\t"))
		ix <- which(toupper(tmp[,"Gene_ID"]) %in% rownames(scores.ko))
		scores.chip.th17[toupper(tmp[,"Gene_ID"])[ix],tf] <- as.numeric(tmp[ix,chip.integration])
	}
}

# put all datasets on same gene names (adding/substracting genes based on genes discovered by KO data, refseq genes)
gns <- rownames(scores.ko)
# deal with KO data
r.scores.ko <- scores.ko
# deal with KO data
r.scores.fc <- scores.fc
# deal with KO data
r.scores.prcnt.chng <- scores.prcnt.chng
# deal with SAM data
m <- matrix(0,nr=nrow(scores.ko),nc=ncol(scores.sam))
rownames(m) <- gns
colnames(m) <- colnames(scores.sam)
ix <- which(gns %in% rownames(scores.sam))
m[gns[ix],] <- scores.sam[gns[ix],]
r.scores.sam <- m
# deal with RNAseq data
m <- matrix(0,nr=nrow(scores.ko),nc=ncol(scores.rnaseq))
rownames(m) <- gns
colnames(m) <- colnames(scores.rnaseq)
ix <- which(gns %in% rownames(scores.rnaseq))
m[gns[ix],] <- scores.rnaseq[gns[ix],]
r.scores.rnaseq <- m
# deal with immgen data
m <- matrix(0,nr=nrow(scores.ko),nc=ncol(scores.immgen))
rownames(m) <- gns
colnames(m) <- colnames(scores.immgen)
ix <- which(gns %in% rownames(scores.immgen))
m[gns[ix],] <- scores.immgen[gns[ix],]
r.scores.immgen <- m

# scores.chip.th17
m <- matrix(0,nr=nrow(scores.ko),nc=ncol(scores.chip.th17))
rownames(m) <- gns
colnames(m) <- colnames(scores.chip.th17)
ix <- which(gns %in% rownames(scores.chip.th17))
m[gns[ix],] <- scores.chip.th17[gns[ix],]
r.scores.chip.th17 <- m
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
m <- numeric(length=nrow(scores.ko))
names(m) <- gns
ix <- which(gns %in% names(th17.mean))
m[gns[ix]] <- th17.mean[gns[ix]]
th17.rpkm <- m
# th0
m <- numeric(length=nrow(scores.ko))
names(m) <- gns
ix <- which(gns %in% names(th0.mean))
m[gns[ix]] <- th0.mean[gns[ix]]
th0.rpkm <- m

#########################################################
# plot distributions for kc pvalues
f.nm <- paste(sep="","results/validation/",date.is,"/","Score_distributions_per_data.pdf")
pdf(f.nm)
order.ix <- c(3,5,4,6,7)
scores.ko.core.tfs <- abs(scores.ko[,order.ix])
scores.ko.core.tfs[which(scores.ko.core.tfs == Inf)] <- 300
colnames(scores.ko.core.tfs) <- colnames(scores.chip.th17)
boxplot(scores.ko.core.tfs,ylim=c(0,20),ylab="-log10(pvalue)",cex=0.3,pch=20,main="knockouts [RNA-seq]")
boxplot(scores.chip.th17,ylim=c(0,20),ylab="-log10(pvalue)",cex=0.3,pch=20,main="TF binding [ChIP-seq]")
boxplot(abs(scores.immgen[,core.tfs]),ylim=c(0,10),ylab="absolute(pseudo-zscore)",cex=0.3,pch=20,main="Inferelator [RNA-seq]")
boxplot(abs(scores.rnaseq[,core.tfs]),ylim=c(0,10),ylab="absolote(pseudo-zscore)",cex=0.3,pch=20,main="Inferelator [Immgen]")
dev.off()
# scores.ko.core.tfs <- scores.ko.core.tfs*log(10)
# scores.chip.th17 <- scores.chip.th17*log(10)
# quartz()
# boxplot(cbind(scores.ko.core.tfs,scores.chip.th17),ylim=c(0,20),ylab="-log2(pvalue)",axes = FALSE,cex=0.3,pch=20)
# axis(2)
# text(1:10, par("usr")[3], labels =c(colnames(scores.ko.core.tfs),colnames(scores.chip.th17)) , srt = 45, adj = 1, xpd = TRUE,cex=1)

# w <- cbind(scores.ko.core.tfs,scores.chip.th17)
# kc.comb <- scores.ko.core.tfs
# kc.comb[which(res!=0)] <- 0
# for(j in 1:length(core.tfs)){
# 	tf <- core.tfs[j]
# 	tmp <- w[,which(colnames(w)==tf)]
# 	w.chi <- apply(tmp,1,function(i) 2*sum(i) )
# 	kc.comb[,tf]  <-  pchisq(w.chi, df=2*ncol(tmp), lower.tail = F)
# }
##############################################################

##################################################

#########################################
# deal with activation
#########################################
cat("Dealing with activatory regulatory interactions:\n")
conds <- c("K","C","R","I","KC","KR","KI","CR","CI","RI","KCR","KCI","KRI","CRI","KCRI","SAM","DEseq","FC","th0_rpkm","th17_rpkm","prcnt.chng")
			# "P300_th17","P300_th0","P300_th17_minus_th0")
ix.filter.conds <- 1:15
cases <- c("th17")
m.tmp <- matrix(0,nr=nrow(scores.ko),nc=length(conds))
rownames(m.tmp) <- rownames(scores.ko)
colnames(m.tmp) <- conds
rel.rank.immgen <- convert.scores.to.relative.ranks.pos(r.scores.immgen)
rel.rank.rnaseq <- convert.scores.to.relative.ranks.pos(r.scores.rnaseq)
rel.rank.ko <- convert.scores.to.relative.ranks.pos(r.scores.ko) # convert to scores based on all tfs
# ix.ko.available <- which(!is.na(x[,"KO_deseq_file"]))
# rel.rank.ko <- rel.rank.ko[,ix.ko.available] # remove unecessary columns

rel.rank.chip.th17 <- convert.scores.to.relative.ranks(r.scores.chip.th17)

# here i keep a version that is editable (to filter by cutoffs)
m.immgen <- rel.rank.immgen
m.rnaseq <- rel.rank.rnaseq
m.ko <- rel.rank.ko
m.chip <- rel.rank.chip.th17

fisher.method <- function(w){
	tmp <- abs(w)
	w.chi <- apply(tmp,1,function(i) 2*sum(i) )
	return(pchisq(w.chi, df=2*ncol(tmp), lower.tail = F))
}


res <- list()
comb.type <- "activation"
res[[comb.type]] <- list()
for(i in 1:nrow(x)){
	m <- m.tmp
	tf <- toupper(x[i,"tf"])
	nm <- toupper(x[i,"name"])
	res[[comb.type]][[nm]] <- list()
	ko.expt.nm <- x[i,"KO_deseq_file"]
	cat("combining results for ",nm,"\n")
	m[,"SAM"] <- r.scores.sam[,"Score_d"]
	ix.bad.2 <- which(abs(m[,"SAM"]) < z.abs.cut)
	if(!is.na(ko.expt.nm)){
		m[,"K"] <- scores.ko[,ko.expt.nm]
		m[,"DEseq"] <- r.scores.ko[,ko.expt.nm]
		m[,"FC"] <- r.scores.fc[,ko.expt.nm]
		m[,"prcnt.chng"] <- r.scores.prcnt.chng[,ko.expt.nm]
		ix.bad.1 <- which(abs(m[,"prcnt.chng"]) < prcnt.chng.cut)	
		ix.bad.3 <- which(abs(m[,"DEseq"]) < -log10(deseq.pval.cut))
		ix.all <- c(ix.bad.1,ix.bad.2,ix.bad.3)
	} else {
		ix.all <- ix.bad.2
	}
	if(tf %in% colnames(rel.rank.rnaseq)){
		# m[,"R"] <- rel.rank.rnaseq[,tf]
		m[,"R"] <- -pnorm(r.scores.rnaseq[,tf],lower.tail = F,log.p = T)*log10(exp(1))
		ix <- which(gns %in% rownames(scores.rnaseq))
		m[gns[ix],"R"] <- m[gns[ix],"R"]*sign.cor.mat.rnaseq[gns[ix],tf]
	}
	if(tf %in% colnames(rel.rank.immgen)){
		# m[,"I"] <- rel.rank.immgen[,tf]
		m[,"I"] <- -pnorm(r.scores.immgen[,tf],lower.tail = F,log.p = T)*log10(exp(1))
		ix <- which(gns %in% rownames(scores.immgen))
		m[gns[ix],"I"] <- m[gns[ix],"I"]*sign.cor.mat.immgen[gns[ix],tf]
	}
	m[,"th0_rpkm"] <- th0.rpkm
	m[,"th17_rpkm"] <- th17.rpkm
	m[,"KR"] <- -log10(fisher.method(m[,c("K","R")]*log(10)))
	m[,"KI"] <- -log10(fisher.method(m[,c("K","I")]*log(10)))
	m[,"RI"] <- -log10(fisher.method(m[,c("R","I")]*log(10)))
	m[,"KRI"] <- -log10(fisher.method(m[,c("K","R","I")]*log(10)))
	# m[,"KR"] <- m[,"K"]+m[,"R"]
	# m[,"KI"] <- m[,"K"]+m[,"I"]
	# m[,"RI"] <- m[,"R"]+m[,"I"]
	# m[,"KRI"] <- m[,"K"]+m[,"R"]+m[,"I"]
	if(length(ix.all)>0){
		ix.bad <- unique(ix.all)
	}
	m[,"C"] <- scores.chip.th17[,tf]
	m[,"KC"]  <-  -log10(fisher.method(m[,c("K","C")]*log(10)))
	# m[,"KC"]  <- pnorm(sum(qnorm(x)) / sqrt(length(x)) stoufer's method, something like that, no need for now
	m[,"CR"] <- -log10(fisher.method(m[,c("C","R")]*log(10)))
	m[,"CI"] <- -log10(fisher.method(m[,c("C","I")]*log(10)))
	m[,"KCR"] <- -log10(fisher.method(m[,c("K","C","R")]*log(10)))
	m[,"KCI"] <- -log10(fisher.method(m[,c("K","C","I")]*log(10)))
	m[,"CRI"] <- -log10(fisher.method(m[,c("C","R","I")]*log(10)))
	m[,"KCRI"] <- -log10(fisher.method(m[,c("K","C","R","I")]*log(10)))
	# apply prcnt change in expression filter
	if(length(ix.all)>0){
		m[ix.bad,ix.filter.conds] <- 0
	}
	res[[comb.type]][[nm]][[cases[1]]] <- m
	path.nm <- paste(sep="",path.output,comb.type,add.str,"/")
	if(! (nm %in% list.files(path.nm))){
		system(paste(sep="","mkdir ", path.nm,nm))			
	}
	fl.nm <- paste(sep="",path.nm,nm,"/",nm,"_",cases[1],".xls")
	cat(sep="\t","Gene_ID",colnames(m),"\n",file=fl.nm)
	write.table(m,sep="\t",file=fl.nm,col.names=F,append=T)
}

fl.nm <- paste(sep="",path.output,"results_combine_data_exprsn",add.str,".Rdata")
save(res,file=fl.nm)















