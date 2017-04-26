##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

source("r_scripts/th17/used_for_paper/util.R")
# setup paths
path.input <- "input/th17/used_for_paper/"
path.input.immgen <- "input/th17/used_for_paper/infResults/immgen/"
path.input.rnaseq <- "input/th17/used_for_paper/infResults/rnaseq/"
path.input.chip <- "/data/th17/data/map_ChIP_peaks_to_genes/"
path.input.deseq <- paste(sep="","input/th17/used_for_paper/DEseq/",date.deseq.run,"/")
path.output <- paste(sep="","results/combinedAnalysis/",date.is,"/")
path.input.immgen.inf.data <-  paste(sep="","input/th17/used_for_paper/infDataStructures/immgen/z_abs_cut_",z.abs.cut.old,"/")
path.input.rnaseq.inf.data <-  paste(sep="","input/th17/used_for_paper/infDataStructures/rnaseq/z_abs_cut_",z.abs.cut.old,"/")

# make output directory
# if(filter.by.sam==TRUE){
add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",z.abs.cut,"_deseq_cut_",deseq.pval.cut)
# } else {
# 	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",0,"_deseq_cut_",deseq.pval.cut)
# }
system(paste(sep="","mkdir ", path.output))
system(paste(sep="","mkdir ", path.output, "activation",add.str,"/"))
system(paste(sep="","mkdir ", path.output, "repression",add.str,"/"))
system(paste(sep="","mkdir ", path.output, "absolute",add.str,"/"))
system(paste(sep="","mkdir ", path.output, "whole",add.str,"/"))


# setup file names
file.nm.map <- paste(sep="",path.input,"combinedAnalysis/",combine.file.map)
file.nm.immgen <- paste(sep="","clrinf_cnfdnc_zcut_",z.abs.cut.old,"_Nb_",num.boots.immgen,"_",  date.inf.run,".xls")
file.nm.rnaseq <- paste(sep="","clrinf_cnfdnc_zcut_",z.abs.cut.old,"_Nb_",num.boots.rnaseq,"_",date.inf.run,".xls")
file.nm.ko <- paste(sep="","DEseq_pval_signed_",date.deseq.run,".xls")
file.nm.fc <- paste(sep="","DEseq_log2fc_",date.deseq.run,".xls")
file.nm.prcnt.chng <- paste(sep="","DEseq_prcnt_chng_",date.deseq.run,".xls")
file.nm.sam <- "samTh17VsTh0Zscores.xls"


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

# here we keep focusing on genes with abs diff expression zscore greater than z.abs.cut.old
# gn.nms.diff <- rownames(scores.sam)[ which( abs(scores.sam[,"Score_d"]) > z.abs.cut ) ]
# scores.sam <- scores.sam[gn.nms.diff,"Score_d"]

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
tfs <- colnames(scores.immgen)
sign.cor.mat.immgen <- sign(cor(t(ratios))[,tfs])
load(paste(sep="",path.input.rnaseq.inf.data,"ratios.RData"))
tfs <- colnames(scores.rnaseq)
sign.cor.mat.rnaseq <- sign(cor(t(ratios))[,tfs])
rm(tfs,ratios)
scores.immgen <- scores.immgen*sign(sign.cor.mat.immgen)
scores.rnaseq <- scores.rnaseq*sign(sign.cor.mat.rnaseq)
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
scores.chip.th0 <- matrix(0,nr=nrow(scores.ko),nc=nrow(x))
scores.chip.th17 <- matrix(0,nr=nrow(scores.ko),nc=nrow(x))
scores.p300.chip.th0 <- matrix(0,nr=nrow(scores.ko),nc=nrow(x))
scores.p300.chip.th17 <- matrix(0,nr=nrow(scores.ko),nc=nrow(x))
rownames(scores.chip.th0) <- rownames(scores.chip.th17) <- rownames(scores.p300.chip.th0) <- rownames(scores.p300.chip.th17) <- rownames(scores.ko)
colnames(scores.chip.th0) <- colnames(scores.chip.th17) <- colnames(scores.p300.chip.th0) <- colnames(scores.p300.chip.th17) <- toupper(x[,"tf"])

for(i in 1:nrow(x)){
	tf <- toupper(x[i,"tf"])
	cat("loading chip data for",tf,"\n")
	chip.th0.nm <- x[i,"chip_th0"]
	chip.th17.nm <- x[i,"chip_th17"]
	p300.chip.th0.nm <- x[i,"p300_chip_th0"]
	p300.chip.th17.nm <- x[i,"p300_chip_th17"]
	# read chip th0
	cat("th0:\n")
	cmd.line <- paste(sep="","scp am2654@bot.bio.nyu.edu:",path.input.chip,chip.th0.nm,"/",chip.th0.nm,"_genes.xls" , " tmp/")
	system(cmd.line)
	cat("th17:\n")
	cmd.line <- paste(sep="","scp am2654@bot.bio.nyu.edu:",path.input.chip,chip.th17.nm,"/",chip.th17.nm,"_genes.xls" , " tmp/")
	system(cmd.line)
	cat("p300 th0:\n")	
	cmd.line <- paste(sep="","scp am2654@bot.bio.nyu.edu:",path.input.chip,p300.chip.th0.nm,"/",p300.chip.th0.nm,"_genes.xls" , " tmp/")
	system(cmd.line)
	cat("p300 th17:\n")
	cmd.line <- paste(sep="","scp am2654@bot.bio.nyu.edu:",path.input.chip,p300.chip.th17.nm,"/",p300.chip.th17.nm,"_genes.xls" , " tmp/")
	system(cmd.line)
	# get th0 scores
	if(!is.na(chip.th0.nm)){
		tmp <- as.matrix(read.delim(paste(sep="","tmp/",chip.th0.nm,"_genes.xls"),sep="\t"))
		ix <- which(toupper(tmp[,"Gene_ID"]) %in% rownames(scores.ko))
		scores.chip.th0[toupper(tmp[,"Gene_ID"])[ix],tf] <- as.numeric(tmp[ix,chip.integration])
	}
	# get th17 scores
	if(!is.na(chip.th17.nm)){	
		tmp <- as.matrix(read.delim(paste(sep="","tmp/",chip.th17.nm,"_genes.xls"),sep="\t"))
		ix <- which(toupper(tmp[,"Gene_ID"]) %in% rownames(scores.ko))
		scores.chip.th17[toupper(tmp[,"Gene_ID"])[ix],tf] <- as.numeric(tmp[ix,chip.integration])
	}
	# get p300 th0 scores		
	if(!is.na(p300.chip.th0.nm)){
		tmp <- as.matrix(read.delim(paste(sep="","tmp/",p300.chip.th0.nm,"_genes.xls"),sep="\t"))
		ix <- which(toupper(tmp[,"Gene_ID"]) %in% rownames(scores.ko))
		scores.p300.chip.th0[toupper(tmp[,"Gene_ID"])[ix],tf] <- as.numeric(tmp[ix,chip.integration])
	}
	# get p300 th17 scores			
	if(!is.na(p300.chip.th17.nm)){
		tmp <- as.matrix(read.delim(paste(sep="","tmp/",p300.chip.th17.nm,"_genes.xls"),sep="\t"))
		ix <- which(toupper(tmp[,"Gene_ID"]) %in% rownames(scores.ko))
		scores.p300.chip.th17[toupper(tmp[,"Gene_ID"])[ix],tf] <- as.numeric(tmp[ix,chip.integration])
	}
}
scores.chip.th17.minus.th0 <- scores.chip.th17-scores.chip.th0
scores.chip.th0.w.p300 <- scores.chip.th0 + scores.p300.chip.th0
scores.chip.th17.w.p300 <- scores.chip.th17 + scores.p300.chip.th17
scores.chip.th17.w.p300.minus.th0.w.p300 <- scores.chip.th17.w.p300-scores.chip.th0.w.p300

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

# deal with chip data
rownames(scores.chip.th0) <- rownames(scores.chip.th17) <- rownames(scores.p300.chip.th0) <- rownames(scores.p300.chip.th17)
# scores.chip.th0
m <- matrix(0,nr=nrow(scores.ko),nc=ncol(scores.chip.th0))
rownames(m) <- gns
colnames(m) <- colnames(scores.chip.th0)
ix <- which(gns %in% rownames(scores.chip.th0))
m[gns[ix],] <- scores.chip.th0[gns[ix],]
r.scores.chip.th0 <- m
# scores.chip.th17
m <- matrix(0,nr=nrow(scores.ko),nc=ncol(scores.chip.th17))
rownames(m) <- gns
colnames(m) <- colnames(scores.chip.th17)
ix <- which(gns %in% rownames(scores.chip.th17))
m[gns[ix],] <- scores.chip.th17[gns[ix],]
r.scores.chip.th17 <- m
# scores.chip.th17.minus.th0
m <- matrix(0,nr=nrow(scores.ko),nc=ncol(scores.chip.th17.minus.th0))
rownames(m) <- gns
colnames(m) <- colnames(scores.chip.th17.minus.th0)
ix <- which(gns %in% rownames(scores.chip.th17.minus.th0))
m[gns[ix],] <- scores.chip.th17.minus.th0[gns[ix],]
r.scores.chip.th17.minus.th0 <- m
# scores.chip.th0.w.p300
m <- matrix(0,nr=nrow(scores.ko),nc=ncol(scores.chip.th0.w.p300))
rownames(m) <- gns
colnames(m) <- colnames(scores.chip.th0.w.p300)
ix <- which(gns %in% rownames(scores.chip.th0.w.p300))
m[gns[ix],] <- scores.chip.th0.w.p300[gns[ix],]
r.scores.chip.th0.w.p300 <- m
# scores.chip.th17.w.p300
m <- matrix(0,nr=nrow(scores.ko),nc=ncol(scores.chip.th17.w.p300))
rownames(m) <- gns
colnames(m) <- colnames(scores.chip.th17.w.p300)
ix <- which(gns %in% rownames(scores.chip.th17.w.p300))
m[gns[ix],] <- scores.chip.th17.w.p300[gns[ix],]
r.scores.chip.th17.w.p300 <- m
# scores.chip.th17.w.p300.minus.th0.w.p300
m <- matrix(0,nr=nrow(scores.ko),nc=ncol(scores.chip.th17.w.p300.minus.th0.w.p300))
rownames(m) <- gns
colnames(m) <- colnames(scores.chip.th17.w.p300.minus.th0.w.p300)
ix <- which(gns %in% rownames(scores.chip.th17.w.p300.minus.th0.w.p300))
m[gns[ix],] <- scores.chip.th17.w.p300.minus.th0.w.p300[gns[ix],]
r.scores.chip.th17.w.p300.minus.th0.w.p300 <- m

# combine data for activation case helper function
chip.score <- function(tf,type) {
  switch(type,
         th0 = rel.rank.chip.th0[,tf],
         th17 = rel.rank.chip.th17[,tf],
         th17_minus_th0 = rel.rank.chip.th17.minus.th0[,tf],
         th0_w_p300 = rel.rank.chip.th0.w.p300[,tf],
         th17_w_p300 = rel.rank.chip.th17.w.p300[,tf],
         th17_w_p300_minus_th0_w_p300 = rel.rank.chip.th17.w.p300.minus.th0.w.p300[,tf])
}

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

##################################################

#remove from scores ko comparisons that are not used (e.g. th17 vs th0 etc.)
ix.cols <- numeric()
file.nm.map <- paste(sep="",path.input,"combinedAnalysis/",combine.file.map.all.tfs)
w <- as.matrix(read.delim(file.nm.map ,sep="\t"),header=T)
for(j in 1:nrow(w)){
	ix.cols <- c(ix.cols,which(colnames(r.scores.ko) == w[j,"KO_deseq_file"]))
}
r.scores.ko <- r.scores.ko[,ix.cols]



#########################################
# deal with activation
#########################################
cat("Dealing with activatory regulatory interactions:\n")
conds <- c("K","C","R","I","KC","KR","KI","CR","CI","RI","KCR","KCI","KRI","CRI","KCRI","SAM","DEseq","FC","th0_rpkm","th17_rpkm","prcnt.chng")
			# "P300_th17","P300_th0","P300_th17_minus_th0")
ix.filter.conds <- 1:15
cases <- c("th0","th17","th17_minus_th0","th0_w_p300","th17_w_p300","th17_w_p300_minus_th0_w_p300")
m.tmp <- matrix(0,nr=nrow(scores.ko),nc=length(conds))
rownames(m.tmp) <- rownames(scores.ko)
colnames(m.tmp) <- conds
rel.rank.immgen <- convert.scores.to.relative.ranks.pos(r.scores.immgen)
rel.rank.rnaseq <- convert.scores.to.relative.ranks.pos(r.scores.rnaseq)
rel.rank.ko <- convert.scores.to.relative.ranks.pos(r.scores.ko) # convert to scores based on all tfs
rel.rank.ko <- rel.rank.ko[,x[,"KO_deseq_file"]] # remove unecessary columns

rel.rank.chip.th0 <- convert.scores.to.relative.ranks(r.scores.chip.th0)
rel.rank.chip.th17 <- convert.scores.to.relative.ranks(r.scores.chip.th17)
rel.rank.chip.th17.minus.th0 <- convert.scores.to.relative.ranks(r.scores.chip.th17.minus.th0)
rel.rank.chip.th0.w.p300 <- convert.scores.to.relative.ranks(r.scores.chip.th0.w.p300)
rel.rank.chip.th17.w.p300 <- convert.scores.to.relative.ranks(r.scores.chip.th17.w.p300)
rel.rank.chip.th17.w.p300.minus.th0.w.p300 <- convert.scores.to.relative.ranks(r.scores.chip.th17.w.p300.minus.th0.w.p300)

# here i keep a version that is editable (to filter by cutoffs)
m.immgen <- rel.rank.immgen
m.rnaseq <- rel.rank.rnaseq
m.ko <- rel.rank.ko
m.chip <- rel.rank.chip.th17

# stop("AM")
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
	if(!is.na(ko.expt.nm)){
		m[,"K"] <- rel.rank.ko[,ko.expt.nm]
		m[,"DEseq"] <- r.scores.ko[,ko.expt.nm]
		m[,"FC"] <- r.scores.fc[,ko.expt.nm]
		m[,"prcnt.chng"] <- r.scores.prcnt.chng[,ko.expt.nm]
	}
	if(tf %in% colnames(rel.rank.rnaseq)){
		m[,"R"] <- rel.rank.rnaseq[,tf]
	}
	if(tf %in% colnames(rel.rank.immgen)){
		m[,"I"] <- rel.rank.immgen[,tf]
	}
	
	m[,"SAM"] <- r.scores.sam[,"Score_d"]
	m[,"th0_rpkm"] <- th0.rpkm
	m[,"th17_rpkm"] <- th17.rpkm
	m[,"KR"] <- m[,"K"]+m[,"R"]
	m[,"KI"] <- m[,"K"]+m[,"I"]
	m[,"RI"] <- m[,"R"]+m[,"I"]
	m[,"KRI"] <- m[,"K"]+m[,"R"]+m[,"I"]
	ix.bad.1 <- which(abs(m[,"prcnt.chng"]) < prcnt.chng.cut)	
	ix.bad.2 <- which(abs(m[,"SAM"]) < z.abs.cut)
	ix.bad.3 <- which(abs(m[,"DEseq"]) < -log10(deseq.pval.cut))
	ix.all <- c(ix.bad.1,ix.bad.2,ix.bad.3)
	if(length(ix.all)>0){
		ix.bad <- unique(c(ix.bad.1,ix.bad.2,ix.bad.3))
	}
	for(j in 1:length(cases)){
		m[,"C"] <- chip.score(tf,cases[j])
		m[,"KC"] <- m[,"K"]+m[,"C"]			
		m[,"CR"] <- m[,"C"]+m[,"R"]	
		m[,"CI"] <- m[,"C"]+m[,"I"]	
		m[,"KCR"] <- m[,"K"]+m[,"C"]+m[,"R"]	
		m[,"KCI"] <- m[,"K"]+m[,"C"]+m[,"I"]		
		m[,"CRI"] <- m[,"C"]+m[,"R"]+m[,"I"]
		m[,"KCRI"] <- m[,"K"]+m[,"C"]+m[,"R"]+m[,"I"]
		# apply prcnt change in expression filter
		if(length(ix.all)>0){
			m[ix.bad,ix.filter.conds] <- 0
		}
		res[[comb.type]][[nm]][[cases[j]]] <- m
		path.nm <- paste(sep="",path.output,comb.type,add.str,"/")
		# path.nm <- paste(sep="",path.output, "activation_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"/")
		# path.nm <- paste(sep="",path.output,"activation_cut_prcnt_",prcnt.chng.cut,"/")
		if(! (nm %in% list.files(path.nm))){
			system(paste(sep="","mkdir ", path.nm,nm))			
		}
		fl.nm <- paste(sep="",path.nm,nm,"/",nm,"_",cases[j],".xls")
		cat(sep="\t","Gene_ID",colnames(m),"\n",file=fl.nm)
		write.table(m,sep="\t",file=fl.nm,col.names=F,append=T)
	}
	if(length(ix.all)>0){
		m.immgen[ix.bad,tf] <- 0
		m.rnaseq[ix.bad,tf] <- 0
		m.ko[ix.bad,ko.expt.nm] <- 0
	  	m.chip[ix.bad,tf] <- 0
	}
}

#################################################

# output I, R, K, C, 
cs <- "th17"
fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_I",add.str,"_",date.is,".xls")
write.table(m.immgen,sep="\t",file=fl.nm)
rel.rank.immgen.activation <- m.immgen # for later to put as whole immgen (minus repression)

fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_R",add.str,"_",date.is,".xls")
write.table(m.rnaseq,sep="\t",file=fl.nm)
rel.rank.rnaseq.activation <- m.rnaseq # for later to put as whole rnaseq (minus repression)


fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_K",add.str,"_",date.is,".xls")
colnames(m.ko) <- sapply(colnames(m.ko), function(i) toupper(strsplit(i,"\\.")[[1]][6]))
ix <- grep("SIETV6",colnames(m.ko))
if(length(ix)>0){
	colnames(m.ko)[ix] <- "ETV6"
}
write.table(m.ko,sep="\t",file=fl.nm)
rel.rank.ko.activation <- m.ko

fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_C",add.str,"_",date.is,".xls")
write.table(m.chip,sep="\t",file=fl.nm)

m <- (m.chip+m.ko)
fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_KC",add.str,"_",date.is,".xls")
write.table(m,sep="\t",file=fl.nm)

#################################################

#########################################
# deal with repression
#########################################
cat("Dealing with repressive regulatory interactions:\n")
rel.rank.immgen <- convert.scores.to.relative.ranks.pos(-1*r.scores.immgen)
rel.rank.rnaseq <- convert.scores.to.relative.ranks.pos(-1*r.scores.rnaseq)
rel.rank.ko <- convert.scores.to.relative.ranks.pos(-1*r.scores.ko) # convert to scores based on all tfs
rel.rank.ko <- rel.rank.ko[,x[,"KO_deseq_file"]] # remove unecessary columns

# here i keep a version that is editable (to filter by cutoffs)
m.immgen <- rel.rank.immgen
m.rnaseq <- rel.rank.rnaseq
m.ko <- rel.rank.ko
m.chip <- rel.rank.chip.th17

comb.type <- "repression"
res[[comb.type]] <- list()
for(i in 1:nrow(x)){
	m <- m.tmp
	tf <- toupper(x[i,"tf"])
	nm <- toupper(x[i,"name"])
	res[[comb.type]][[nm]] <- list()
	ko.expt.nm <- x[i,"KO_deseq_file"]
	cat("combining results for ",nm,"\n")
	if(!is.na(ko.expt.nm)){
		m[,"K"] <- rel.rank.ko[,ko.expt.nm]
		m[,"DEseq"] <- r.scores.ko[,ko.expt.nm]
		m[,"FC"] <- r.scores.fc[,ko.expt.nm]
		m[,"prcnt.chng"] <- r.scores.prcnt.chng[,ko.expt.nm]
	}
	if(tf %in% colnames(rel.rank.rnaseq)){
		m[,"R"] <- rel.rank.rnaseq[,tf]
	}
	if(tf %in% colnames(rel.rank.immgen)){
		m[,"I"] <- rel.rank.immgen[,tf]
	}
	
	m[,"SAM"] <- r.scores.sam[,"Score_d"]
	m[,"th0_rpkm"] <- th0.rpkm
	m[,"th17_rpkm"] <- th17.rpkm
	m[,"KR"] <- m[,"K"]+m[,"R"]
	m[,"KI"] <- m[,"K"]+m[,"I"]
	m[,"RI"] <- m[,"R"]+m[,"I"]
	m[,"KRI"] <- m[,"K"]+m[,"R"]+m[,"I"]
	ix.bad.1 <- which(abs(m[,"prcnt.chng"]) < prcnt.chng.cut)	
	ix.bad.2 <- which(abs(m[,"SAM"]) < z.abs.cut)
	ix.bad.3 <- which(abs(m[,"DEseq"]) < -log10(deseq.pval.cut))
	ix.all <- c(ix.bad.1,ix.bad.2,ix.bad.3)
	if(length(ix.all)>0){
		ix.bad <- unique(c(ix.bad.1,ix.bad.2,ix.bad.3))
	}
	for(j in 1:length(cases)){
		m[,"C"] <- chip.score(tf,cases[j])
		m[,"KC"] <- m[,"K"]+m[,"C"]			
		m[,"CR"] <- m[,"C"]+m[,"R"]	
		m[,"CI"] <- m[,"C"]+m[,"I"]	
		m[,"KCR"] <- m[,"K"]+m[,"C"]+m[,"R"]	
		m[,"KCI"] <- m[,"K"]+m[,"C"]+m[,"I"]		
		m[,"CRI"] <- m[,"C"]+m[,"R"]+m[,"I"]
		m[,"KCRI"] <- m[,"K"]+m[,"C"]+m[,"R"]+m[,"I"]
		# apply prcnt change in expression filter
		if(length(ix.all)>0){
			m[ix.bad,ix.filter.conds] <- 0
		}
		res[[comb.type]][[nm]][[cases[j]]] <- m
		path.nm <- paste(sep="",path.output,comb.type,add.str,"/")
		# path.nm <- paste(sep="",path.output, "activation_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"/")
		# path.nm <- paste(sep="",path.output,"activation_cut_prcnt_",prcnt.chng.cut,"/")
		if(! (nm %in% list.files(path.nm))){
			system(paste(sep="","mkdir ", path.nm,nm))			
		}
		fl.nm <- paste(sep="",path.nm,nm,"/",nm,"_",cases[j],".xls")
		cat(sep="\t","Gene_ID",colnames(m),"\n",file=fl.nm)
		write.table(m,sep="\t",file=fl.nm,col.names=F,append=T)
	}
	if(length(ix.all)>0){
		m.immgen[ix.bad,tf] <- 0
		m.rnaseq[ix.bad,tf] <- 0
		m.ko[ix.bad,ko.expt.nm] <- 0
	  m.chip[ix.bad,tf] <- 0
	}
}

#################################################

# output I, R, K, C, 
cs <- "th17"
fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_I",add.str,"_",date.is,".xls")
write.table(m.immgen,sep="\t",file=fl.nm)
rel.rank.immgen.repression <- m.immgen # for later to put as whole immgen (minus repression)

fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_R",add.str,"_",date.is,".xls")
write.table(m.rnaseq,sep="\t",file=fl.nm)
rel.rank.rnaseq.repression <- m.rnaseq # for later to put as whole rnaseq (minus repression)


fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_K",add.str,"_",date.is,".xls")
colnames(m.ko) <- sapply(colnames(m.ko), function(i) toupper(strsplit(i,"\\.")[[1]][6]))
ix <- grep("SIETV6",colnames(m.ko))
if(length(ix)>0){
	colnames(m.ko)[ix] <- "ETV6"
}
write.table(m.ko,sep="\t",file=fl.nm)
rel.rank.ko.repression <- m.ko

fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_C",add.str,"_",date.is,".xls")
write.table(m.chip,sep="\t",file=fl.nm)

m <- (m.chip+m.ko)
fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_KC",add.str,"_",date.is,".xls")
write.table(m,sep="\t",file=fl.nm)
#################################################

#########################################
# deal with absolute regulation (repression and activation)
#########################################
cat("Dealing with absolute regulatory interactions:\n")
rel.rank.immgen <- convert.scores.to.relative.ranks.pos(abs(r.scores.immgen))
rel.rank.rnaseq <- convert.scores.to.relative.ranks.pos(abs(r.scores.rnaseq))
rel.rank.ko <- convert.scores.to.relative.ranks.pos(abs(r.scores.ko)) # convert to scores based on all tfs
rel.rank.ko <- rel.rank.ko[,x[,"KO_deseq_file"]] # remove unecessary columns

# here i keep a version that is editable (to filter by cutoffs)
m.immgen <- rel.rank.immgen
m.rnaseq <- rel.rank.rnaseq
m.ko <- rel.rank.ko
m.chip <- rel.rank.chip.th17

comb.type <- "absolute"
res[[comb.type]] <- list()
for(i in 1:nrow(x)){
	m <- m.tmp
	tf <- toupper(x[i,"tf"])
	nm <- toupper(x[i,"name"])
	res[[comb.type]][[nm]] <- list()
	ko.expt.nm <- x[i,"KO_deseq_file"]
	cat("combining results for ",nm,"\n")
	if(!is.na(ko.expt.nm)){
		m[,"K"] <- rel.rank.ko[,ko.expt.nm]
		m[,"DEseq"] <- r.scores.ko[,ko.expt.nm]
		m[,"FC"] <- r.scores.fc[,ko.expt.nm]
		m[,"prcnt.chng"] <- r.scores.prcnt.chng[,ko.expt.nm]
	}
	if(tf %in% colnames(rel.rank.rnaseq)){
		m[,"R"] <- rel.rank.rnaseq[,tf]
	}
	if(tf %in% colnames(rel.rank.immgen)){
		m[,"I"] <- rel.rank.immgen[,tf]
	}
	
	m[,"SAM"] <- r.scores.sam[,"Score_d"]
	m[,"th0_rpkm"] <- th0.rpkm
	m[,"th17_rpkm"] <- th17.rpkm
	m[,"KR"] <- m[,"K"]+m[,"R"]
	m[,"KI"] <- m[,"K"]+m[,"I"]
	m[,"RI"] <- m[,"R"]+m[,"I"]
	m[,"KRI"] <- m[,"K"]+m[,"R"]+m[,"I"]
	ix.bad.1 <- which(abs(m[,"prcnt.chng"]) < prcnt.chng.cut)	
	ix.bad.2 <- which(abs(m[,"SAM"]) < z.abs.cut)
	ix.bad.3 <- which(abs(m[,"DEseq"]) < -log10(deseq.pval.cut))
	ix.all <- c(ix.bad.1,ix.bad.2,ix.bad.3)
	if(length(ix.all)>0){
		ix.bad <- unique(c(ix.bad.1,ix.bad.2,ix.bad.3))
	}
	for(j in 1:length(cases)){
		m[,"C"] <- chip.score(tf,cases[j])
		m[,"KC"] <- m[,"K"]+m[,"C"]			
		m[,"CR"] <- m[,"C"]+m[,"R"]	
		m[,"CI"] <- m[,"C"]+m[,"I"]	
		m[,"KCR"] <- m[,"K"]+m[,"C"]+m[,"R"]	
		m[,"KCI"] <- m[,"K"]+m[,"C"]+m[,"I"]		
		m[,"CRI"] <- m[,"C"]+m[,"R"]+m[,"I"]
		m[,"KCRI"] <- m[,"K"]+m[,"C"]+m[,"R"]+m[,"I"]
		# apply prcnt change in expression filter
		if(length(ix.all)>0){
			m[ix.bad,ix.filter.conds] <- 0
		}
		res[[comb.type]][[nm]][[cases[j]]] <- m
		path.nm <- paste(sep="",path.output,comb.type,add.str,"/")
		# path.nm <- paste(sep="",path.output, "activation_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"/")
		# path.nm <- paste(sep="",path.output,"activation_cut_prcnt_",prcnt.chng.cut,"/")
		if(! (nm %in% list.files(path.nm))){
			system(paste(sep="","mkdir ", path.nm,nm))			
		}
		fl.nm <- paste(sep="",path.nm,nm,"/",nm,"_",cases[j],".xls")
		cat(sep="\t","Gene_ID",colnames(m),"\n",file=fl.nm)
		write.table(m,sep="\t",file=fl.nm,col.names=F,append=T)
	}
	if(length(ix.all)>0){
		m.immgen[ix.bad,tf] <- 0
		m.rnaseq[ix.bad,tf] <- 0
		m.ko[ix.bad,ko.expt.nm] <- 0
	  m.chip[ix.bad,tf] <- 0
	}
}

#################################################

# output I, R, K, C, 
cs <- "th17"
fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_I",add.str,"_",date.is,".xls")
write.table(m.immgen,sep="\t",file=fl.nm)

fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_R",add.str,"_",date.is,".xls")
write.table(m.rnaseq,sep="\t",file=fl.nm)


fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_K",add.str,"_",date.is,".xls")
colnames(m.ko) <- sapply(colnames(m.ko), function(i) toupper(strsplit(i,"\\.")[[1]][6]))
ix <- grep("SIETV6",colnames(m.ko))
if(length(ix)>0){
	colnames(m.ko)[ix] <- "ETV6"
}
write.table(m.ko,sep="\t",file=fl.nm)

fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_C",add.str,"_",date.is,".xls")
write.table(m.chip,sep="\t",file=fl.nm)

m <- (m.chip+m.ko)
fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_KC",add.str,"_",date.is,".xls")
write.table(m,sep="\t",file=fl.nm)
######################################################
cat("Dealing with 'whole' regulatory interactions:\n")
comb.type <- "whole"
# add another combination based on both repression and activation but keeping sign (unlike absolute which is all positive)
#############
# get fold change sign so we can multiply chip with (give sign to chip)
# c.names <- x[,"KO_deseq_file"]
# fc <- r.scores.fc[,c.names]
# colnames(fc) <- sapply(colnames(fc), function(i) toupper(strsplit(i,"\\.")[[1]][6]))
# ix <- grep("SIETV6",colnames(fc))
# if(length(ix)>0){
# 	colnames(fc)[ix] <- "ETV6"
# }
# rel.rank.fc <- fc
#############

ix.keep <- 1:15
res[[comb.type]]	 <- list()
for(j in 1:length(names(res[[1]]))){ # go over tfs
	tf <- names(res[[1]][j])
	res[[comb.type]][[tf]] <- list()
	for(k in 1:length(res[[1]][[tf]])){ # go over chip types (this is same for any tf)
		c.tp <- names(res[[1]][[tf]])[k]
		# make whole combo of act and rep
		res[[comb.type]][[tf]][[c.tp]] <- res[["activation"]][[tf]][[c.tp]]
		for(l in 1:length(ix.keep)){
			ix  <- which(res[[comb.type]][[tf]][[c.tp]][,ix.keep[l]] < abs(res[["repression"]][[tf]][[c.tp]][,ix.keep[l]]))
			res[[comb.type]][[tf]][[c.tp]][ix,ix.keep[l]] <- -1*res[["repression"]][[tf]][[c.tp]][ix,ix.keep[l]]
		}

		path.nm <- paste(sep="",path.output,comb.type,add.str,"/")
		if(! (tf %in% list.files(path.nm))){
			system(paste(sep="","mkdir ", path.nm,tf))			
		}
		m <- res[[comb.type]][[tf]][[c.tp]]
		fl.nm <- paste(sep="",path.nm,tf,"/",tf,"_",c.tp,".xls")
		cat(sep="\t","Gene_ID",colnames(m),"\n",file=fl.nm)
		write.table(m,sep="\t",file=fl.nm,col.names=F,append=T)		
	}
}
# get matrix for whole KC scores
comb.type <- "whole"
combine.type.2 <- "KC"
combine.type.1 <- "th17"
tfs <- names(res[[comb.type]])
m.kc <- matrix(0,nr=nrow(res[[comb.type]][[tf]][[c.tp]]),nc=length(tfs))

colnames(m.kc) <- tfs
rownames(m.kc) <- rownames(res[[comb.type]][[tf]][[c.tp]])
for(j in 1:length(tfs)){
	tf <- tfs[j]
	m.kc[,tf] <- res[[comb.type]][[tf]][[combine.type.1]][,combine.type.2]
}
#################################################
rel.rank.immgen <- rel.rank.immgen.activation - rel.rank.immgen.repression
rel.rank.rnaseq <- rel.rank.rnaseq.activation - rel.rank.rnaseq.repression
rel.rank.ko <- rel.rank.ko.activation - rel.rank.ko.repression
# output I, R, K, C, 

fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_I",add.str,"_",date.is,".xls")
m <- rel.rank.immgen
if(length(ix.all)>0){
	m[ix.bad,] <- 0			
}
write.table(m,sep="\t",file=fl.nm)

fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_R",add.str,"_",date.is,".xls")
m <- rel.rank.rnaseq
if(length(ix.all)>0){
	m[ix.bad,] <- 0			
}
write.table(m,sep="\t",file=fl.nm)

fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_K",add.str,"_",date.is,".xls")
m <- rel.rank.ko # already filtered b/c activation and repression ko scores are already filtered
write.table(m,sep="\t",file=fl.nm)

fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_C",add.str,"_",date.is,".xls")
m <- rel.rank.chip.th17 # already filtered 
write.table(m,sep="\t",file=fl.nm)


m <- m.kc
fl.nm <- paste(sep="",path.output,cs,"_",comb.type,"_KC",add.str,"_",date.is,".xls")
write.table(m,sep="\t",file=fl.nm)
######################################################

fl.nm <- paste(sep="",path.output,"results_combine_data_exprsn",add.str,".Rdata")
save(res,file=fl.nm)
#################################################


