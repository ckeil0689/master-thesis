##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## May 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

## add 1 read psudo count
## take pvals instead of adjusted pvals

library(gregmisc) 
## helper function
plotDE <- function( res ) {
 plot( 
 res$baseMean, 
 res$log2FoldChange, 
 log="x", pch=20, cex=.1, 
 col = ifelse( res$padj < .1, "red", "black" ) ) 
}
########## read data ##########
library(DESeq)
cat("reading data\n")
# set paths
path.input <- "input/th17/used_for_paper/"
path.input.deseq <- "input/th17/used_for_paper/DEseq/"
path.input.counts <- paste(sep="","input/th17/used_for_paper/rawData/rnaseq_ht_counts_",date.htseq.data.compiled,"/")
# path.input.counts <- "/data/th17/data/htseq-count/"
path.input.fpkm <- "/data/th17/data/cufflinks/"
path.results <- paste(sep="","results/diff_expression/", date.is,"/")
system(paste(sep="","mkdir ",path.results))

# set file names
f.nm.info <- paste(sep="",path.input.deseq,"diffExpressionExpts_",date.htseq.data.compiled,".txt") # which experiment to compare
f.nm.rpkm <- paste(sep="",path.input,"ranseqDatasetNoQuartNorm.RData") # which experiment to compare

## get count files (counting for each gene how many tags hit it)
count.files <- grep("counts",list.files(path.input.counts),value=T)

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
gns <- names(rpkm.mean)
# gns.low.rpkm <- names(rpkm.mean)[which(rpkm.mean < fpkm.cut)]


## read tab delim file with experiments we want to compare (format below):
## name	expts
## Th17.rorawt.rorcwt.vs.Th17.rorako.rorcko	SL5874_SL5878::SL5876_SL5880
## Th17.dmso.vs.Th0.dmso.human	SL5252_SL5255::SL5251_SL5254
cat("read expt mat\n")
expt.mat <- read.delim(file=f.nm.info,colClasses = "character")

for(i in 1:nrow(expt.mat)){
  e.wt <- strsplit(strsplit(expt.mat[i,2],"::")[[1]],"_")[[1]]
  e.ko <- strsplit(strsplit(expt.mat[i,2],"::")[[1]],"_")[[2]]
  e.nm <- expt.mat[i,1]
  e.nm.wt <- strsplit(e.nm,".vs.")[[1]][1]
  e.nm.ko <- strsplit(e.nm,".vs.")[[1]][2]
  cat("working on",e.nm,"\n")
  f.nm.wt<- sapply(e.wt,function(i) grep(paste(sep="",i,"\\."),count.files,value=T))
  f.nm.ko<- sapply(e.ko,function(i) grep(paste(sep="",i,"\\."),count.files,value=T))
  is.kd <- as.logical(as.numeric(expt.mat[i,3])) # val here is either 1 or zero converted to T/F
  if(i == 1){
    x <- read.delim(file=paste(sep="",path.input.counts,f.nm.wt[1]),header=T)
    n.row <- dim(x)[1]
    gn.nms <- toupper(as.character(x[,1]))
	gn.nms.2 <- gn.nms[-which(gn.nms %in% feature.rm)]
    s.mat.pval <- matrix(0,nc=nrow(expt.mat),nr=length(gn.nms.2),dimnames=list(gn.nms.2,expt.mat[,1]))
    s.mat.padj <- matrix(0,nc=nrow(expt.mat),nr=length(gn.nms.2),dimnames=list(gn.nms.2,expt.mat[,1]))
	s.mat.qval <- matrix(0,nc=nrow(expt.mat),nr=length(gn.nms.2),dimnames=list(gn.nms.2,expt.mat[,1]))
	s.mat.fdr <- matrix(0,nc=nrow(expt.mat),nr=length(gn.nms.2),dimnames=list(gn.nms.2,expt.mat[,1]))
    s.mat.log2fc <- matrix(0,nc=nrow(expt.mat),nr=length(gn.nms.2),dimnames=list(gn.nms.2,expt.mat[,1]))
	s.mat.prcnt.chng <- matrix(0,nc=nrow(expt.mat),nr=length(gn.nms.2),dimnames=list(gn.nms.2,expt.mat[,1]))
  }
  ## get ht experiments for e.wt into a matrix
  m.wt <- matrix(0,nr=n.row,nc=length(f.nm.wt))
  rownames(m.wt) <- gn.nms
  for(j in 1:length(f.nm.wt)){
    x <- read.delim(file=paste(sep="",path.input.counts,f.nm.wt[j]),header=T)
    m.wt[,j] <- as.numeric(x[,2])
  }
  m.wt <- m.wt[-which(rownames(m.wt) %in%   feature.rm),]
  ## get ht experiments for e.ko into a matrix
  m.ko <- matrix(0,nr=n.row,nc=length(f.nm.ko))
  rownames(m.ko) <- gn.nms
  for(j in 1:length(f.nm.ko)){
    x <- read.delim(file=paste(sep="",path.input.counts,f.nm.ko[j]),header=T)
    m.ko[,j] <- as.numeric(x[,2])
  }
  m.ko <- m.ko[-which(rownames(m.ko) %in%   feature.rm),]
  countData <- cbind(m.wt,m.ko) + pseudo.count
  colnames(countData) <- c(e.wt,e.ko)
  conditions <- factor(c(rep(e.nm.wt,length(e.wt)),rep(e.nm.ko,length(e.ko))),levels=c(e.nm.ko,e.nm.wt))
  cds <- newCountDataSet(countData, conditions, sizeFactors = NULL, phenoData = NULL, featureData=NULL)
  cds <- estimateSizeFactors( cds )
  ## if we have no replicates estimate gene's variance from both conditions as biological repeats
  if(dim(cds)[2]==2){
    cds <- estimateVarianceFunctions( cds,pool=T)
  } else {
    cds <- estimateVarianceFunctions( cds )
  }
  res <- nbinomTest( cds, levels(conditions)[1],levels(conditions)[2])
  colnames(res)[3:4] <- levels(conditions)
	################ get QC measures
	cor.btwn.reps <- cor(counts(cds))
	################ Estimate FDRs
	# number of reps for experiments a and b
	exp.a.n <- length(which(conditions==levels(conditions)[1]))
	exp.b.n <- length(which(conditions==levels(conditions)[2]))
	if(exp.a.n >=2 &
		 exp.b.n>=2){
		all.perms.conds <- permutations(2,length(conditions),v=c(1,2),repeats.allowed=TRUE)
		# choose only permutations that are balanced (like the real experiment)
		l <- apply(all.perms.conds,1, function(i) length(which(i==1)) )
		balanced.perms.conds <- all.perms.conds[which(l==exp.a.n),]
		# remove permutatinos that are the same as real experiment
		search <- paste(sep="",as.numeric(conditions),collapse="_")
		search.rev <- paste(sep="",rev(as.numeric(conditions)),collapse="_")
		perms <- apply(balanced.perms.conds,1,function(i) paste(sep="",i,collapse="_"))
		ix.rm <- c(which(perms==search),which(perms==search.rev))
		balanced.perms.conds <- balanced.perms.conds[-ix.rm,] # this may have repeats, ahhh, not fixed yet
		# remove perms that are the same but switch expt a with b
		balanced.perms.conds.switch <- balanced.perms.conds
		ix.1 <- which(balanced.perms.conds.switch==1)
		ix.2 <- which(balanced.perms.conds.switch==2)
		balanced.perms.conds.switch[ix.1] <- 2
		balanced.perms.conds.switch[ix.2] <- 1
		ix.rm <- numeric()
		for(p.1 in 1:floor(nrow(balanced.perms.conds)/2)){
			str.1 <- paste(sep="",balanced.perms.conds[p.1,],collapse="_")
			for(p.2 in 1:nrow(balanced.perms.conds.switch)){
				str.2 <- paste(sep="",balanced.perms.conds.switch[p.2,],collapse="_")
				if(str.2==str.1){
					ix.rm <- c(ix.rm,p.2)
				}
			}
		}
		if(length(ix.rm)>0){
			balanced.perms.conds <- balanced.perms.conds[-ix.rm,]
		}
		# put actual expt names (instead of the generic 1 or 2)
		ix.1 <- which(balanced.perms.conds==1)
		ix.2 <- which(balanced.perms.conds==2)
		balanced.perms.conds[ix.1] <- levels(conditions)[1]
		balanced.perms.conds[ix.2] <- levels(conditions)[2]

		pvals.perms <- matrix(1,nr=nrow(balanced.perms.conds),nc=nrow(res))
		for(p in 1:nrow(balanced.perms.conds)){
			# for(p in 1:3){
			conditions.perm <- as.factor(balanced.perms.conds[p,])
			cds <- newCountDataSet(countData, conditions.perm, sizeFactors = NULL, phenoData = NULL, featureData=NULL)
			cds <- estimateSizeFactors( cds )
			## if we have no replicates estimate gene's variance from both conditions as biological repeats
			if(dim(cds)[2]==2){
			  cds <- estimateVarianceFunctions( cds,pool=T)
			} else {
			  cds <- estimateVarianceFunctions( cds )
			}
			
			res.perm <- nbinomTest( cds, levels(conditions.perm)[1],levels(conditions.perm)[2])						
			pvals.perms[p,] <- res.perm$pval
		}
		# pval.p <- sort(apply(pvals.perms,2,mean))
		pval.p <- sort(apply(pvals.perms,2,function(i) prod(i)^(1/length(i))))
		pval.r <- sort(res$pval)
		pval.r.ix <- sort(res$pval,index.return=T)$ix
		fdr <- rep(NA,length(pval.p))
		qval <- rep(NA,length(pval.p))
		n <- length(pval.r)
		cnt.p <- 0 # counter over current pos in perm pvals
		max.fdr <- 0
		for(cnt.r in 1:n){ # cnt.r counter over current pos in real pvals
			while( (pval.p[cnt.p+1]<=pval.r[cnt.r]) & cnt.p<n ){
				cnt.p <- cnt.p + 1
			} 
			fdr[cnt.r] <- cnt.p/(cnt.p+cnt.r)
			if(cnt.r>20){
				if(max.fdr<fdr[cnt.r]){
					max.fdr <- fdr[cnt.r]
				}
				qval[cnt.r] <- max.fdr
			}
		}
		qval[1:20] <- qval[21] # guaranty qval is always increasing, but not noisy in first few genes
		# plot the fdr and qvals for this experiment
		f.nm <- paste(path.results,"fdr_",e.nm,"_",date.is,".pdf",sep="")
		pdf(f.nm)
			plot(fdr,main="fdr",type="l",xlab="gene.index")
			plot(qval,main="qval",type="l",xlab="gene.index")
		dev.off()
		# write correlation matrix for experiments in fastq format
		# tmp <- cor(counts(cds))
		f.nm.cor <- paste(path.results,"cor_structure_",date.is,".txt",sep="")
		if(i==1){
			cat(sep="",file=f.nm.cor,"> ",expt.mat[i,1],"_",expt.mat[i,2],"\n",append=F)
		} else {
			cat(sep="",file=f.nm.cor,"> ",expt.mat[i,1],"_",expt.mat[i,2],"\n",append=T)
		}
		cat(file=f.nm.cor,"      ",rownames(cor.btwn.reps),"\n",append=T)
		for(sl in 1:nrow(cor.btwn.reps)){
			cat(file=f.nm.cor,colnames(cor.btwn.reps)[sl],cor.btwn.reps[sl,],"\n",append=T)
		}
		res$fdr[pval.r.ix] <- fdr
		res$qval[pval.r.ix] <- qval
	} else {
		fdr <- rep(NA,length(res$pval))
		qval <- rep(NA,length(res$pval))		
		res$fdr <- fdr
		res$qval <- qval
	}
  #################### FIND genes with low expression in wt and KO ######################
  # get fpkm files from server (bot)
  f.nm.wt<- paste(sep="",e.wt,"_genes.expr")
  cat("get rpkm for expt:",e.nm.wt,"SL numbers:",e.wt,"\n")
  for(j in 1:length(f.nm.wt)){
	cmd.line <- paste(sep="","scp am2654@bot.bio.nyu.edu:",path.input.fpkm,e.wt[j],"/cufflinks_out/",f.nm.wt[j], " tmp/")
	system(cmd.line)
  }
	# read fpkm files from local (tmp/)
  for(j in 1:length(f.nm.wt)){
		# fields.type.list <- list(gene_id=character(),bundle_id=numeric(),chr=character(),left=numeric(),right=numeric(),
		# FPKM=numeric(),FPKM_conf_lo=numeric(),FPKM_conf_hi=numeric(),status=character())
		# w=as.data.frame(scan(file=paste(sep="","tmp/",f.nm.wt[j]),sep="\t",what=fields.type.list,skip = 1 ))
		tmp <- as.matrix(read.delim(paste(sep="","tmp/",f.nm.wt[j]),sep="\t"))
		if(j==1){
		  fpkm.wt <- matrix(0,nr=length(unique(tmp[,"gene_id"])),nc=length(f.nm.wt) )
		  rownames(fpkm.wt) <- unique(tmp[,"gene_id"])
		}
		tmp.non.unique <- table(tmp[,1])
		ix.non.unique <- which(tmp.non.unique>1)
		gns.unique <- names(which(tmp.non.unique==1))
		ix.tmp <- which(tmp[,"gene_id"] %in% gns.unique)
  	fpkm.wt[tmp[ix.tmp,"gene_id"],j] <- as.numeric(tmp[ix.tmp,"FPKM"])
  	for(k in 1:length(ix.non.unique)){
		gn.id <- names(ix.non.unique[k])
		ix <- which(tmp[,"gene_id"]==gn.id)
		fpkm.wt[gn.id,j] <- max(as.numeric(tmp[ix,"FPKM"]))
  	}
  }

  f.nm.ko<- paste(sep="",e.ko,"_genes.expr")
  cat("get rpkm for expt:",e.nm.ko,"SL numbers:",e.ko,"\n")
  for(j in 1:length(f.nm.ko)){
	cmd.line <- paste(sep="","scp am2654@bot.bio.nyu.edu:",path.input.fpkm,e.ko[j],"/cufflinks_out/",f.nm.ko[j], " tmp/")
	system(cmd.line)
  }

  for(j in 1:length(f.nm.ko)){
	tmp <- as.matrix(read.delim(paste(sep="","tmp/",f.nm.ko[j]),sep="\t"))
	if(j==1){
	  fpkm.ko <- matrix(0,nr=length(unique(tmp[,"gene_id"])),nc=length(f.nm.ko) )
	  rownames(fpkm.ko) <- unique(tmp[,"gene_id"])
	}		
	# tmp.non.unique <- sort(table(tmp[,1]),decreasing=T)
	tmp.non.unique <- table(tmp[,1])
	ix.non.unique <- which(tmp.non.unique>1)
	# ix.unique <- which(tmp.non.unique==1)
	gns.unique <- names(which(tmp.non.unique==1))
	ix.tmp <- which(tmp[,"gene_id"] %in% gns.unique)
  	fpkm.ko[tmp[ix.tmp,"gene_id"],j] <- as.numeric(tmp[ix.tmp,"FPKM"])
  	for(k in 1:length(ix.non.unique)){
		gn.id <- names(ix.non.unique[k])
		ix <- which(tmp[,"gene_id"]==gn.id)
		fpkm.ko[gn.id,j] <- max(as.numeric(tmp[ix,"FPKM"]))
  	}
  }

  rownames(fpkm.wt) <- toupper(rownames(fpkm.wt))
  rownames(fpkm.ko) <- toupper(rownames(fpkm.ko))
  fpkm.wt <- apply(fpkm.wt,1,mean)
  fpkm.ko <- apply(fpkm.ko,1,mean)

  # find genes that are not expressed by more than fc.cut rpkm in th0,th17, or ko condition of DEseq comparison
  gns.low.ko <- names(which((fpkm.ko < fpkm.cut)))
  gns.low.th0 <- names(which(th0.mean < fpkm.cut))
  gns.low.th17 <- names(which(th17.mean < fpkm.cut))
  gns.low.expressed <- intersect(intersect(gns.low.th17,gns.low.th0),gns.low.ko)

  ###################### END:FIND genes with low expression in wt and KO ######################
  
  ###################### START:find genes with low FC btwn wt and KO ######################
	#   if(is.kd==TRUE){
	# fc.cut <- fc.cut.kd
	#   } else {
	# fc.cut <- fc.cut.ko
	#   }
	# 
	prcnt.change <- abs(res[,4]-res[,3])/res[,4]*100 # col 4 is wt, col 3 is ko
	#   ix.low.fc <- which(prcnt.change<=fc.cut)
	#   res[ix.low.fc,c("pval","padj")] <- 1
  ###################### END:find genes with low FC btwn wt and KO ######################

  ix.low.expressed <- which(res[,"id"] %in% gns.low.expressed)
  #res[ix.low.expressed,c("pval","padj")] <- 1
  res[ix.low.expressed,c("pval","padj","fdr","qval")] <- 1

  # res$mean.rpkm.wt <- res$mean.rpkm.ko <- rep(0,length(res$id))
  res$mean.rpkm.ko <- fpkm.ko[res$id]
  res$mean.rpkm.wt <- fpkm.wt[res$id]
  # also if no fpkm value was found remove (NAs)
  ix <- which(is.na(res[,"mean.rpkm.wt"]))
  res[ix,c("pval","padj","fdr","qval")] <- 1
  # res[ix.low.expressed,c("pval","padj")] <- 1
  gns.with.rpkm <- res$id[which(res$id %in% gns)]
  res$mean.rpkm.th17 <- res$mean.rpkm.th0 <- rep(0,length(res$id))
  res$mean.rpkm.th17[which(res$id %in% gns.with.rpkm)] <- th17.mean[gns.with.rpkm]
  res$mean.rpkm.th0[which(res$id %in% gns.with.rpkm)] <- th0.mean[gns.with.rpkm]
  res$prcnt.chng <- prcnt.change
  s.mat.pval[,i] <- res$pval
  s.mat.padj[,i] <- res$padj
  s.mat.log2fc[,i] <- res$log2FoldChange
  s.mat.prcnt.chng[,i] <- res$prcnt.chng
  s.mat.qval[,i] <- res$qval
  s.mat.fdr[,i] <- res$fdr
  write.table(res,file=paste(path.results,e.nm,"_",date.is,".xls",sep=""),sep="\t",row.names=FALSE)
}

## here we 'multiply' (element wise, A_i_j times B_i_j for all i and j) the pval matrix by the sign of the log2 fold change matrix
## this step allows us to view only the pvalue and say if we have down regulation or up regulation

# this matrix will be used for ChIP scores under combine data script
path.nm <- paste(sep="",path.input,"activation/")

if(! (date.is %in% list.files(path.input.deseq))){
	system(paste(sep="","mkdir ", path.input.deseq,date.is,"/"))			
}

path.input.dated <- paste(sep="", path.input.deseq,date.is,"/")
s.mat.pval.signed <- -log10(s.mat.pval) * sign(s.mat.log2fc)
write.table(s.mat.pval.signed,file=paste(path.input.dated,"DEseq_pval_signed_",date.is,".xls",sep=""),sep="\t")
write.table(s.mat.log2fc,file=paste(path.input.dated,"DEseq_log2fc_",date.is,".xls",sep=""),sep="\t")
write.table(s.mat.prcnt.chng,file=paste(path.input.dated,"DEseq_prcnt_chng_",date.is,".xls",sep=""),sep="\t")
## other results...
write.table(s.mat.pval.signed,file=paste(path.results,"summary.pval.signed.xls",sep=""),sep="\t")
write.table(s.mat.pval,file=paste(path.results,"summary.pval.xls",sep=""),sep="\t")
write.table(s.mat.padj,file=paste(path.results,"summary.pval.adjusted.xls",sep=""),sep="\t")
write.table(s.mat.log2fc,file=paste(path.results,"summary.log2fc.xls",sep=""),sep="\t")
write.table(s.mat.prcnt.chng,file=paste(path.results,"summary.prcnt.chng.xls",sep=""),sep="\t")
write.table(s.mat.fdr,file=paste(path.results,"summary.fdr.xls",sep=""),sep="\t")
write.table(s.mat.qval,file=paste(path.results,"summary.qval.xls",sep=""),sep="\t")

n <- ncol(s.mat.pval.signed)
ix.invivo <- grep("IL17a.gfp.plusve.SI.vs.IL17a.gfp.minusve.SI",colnames(s.mat.pval.signed))
fc.invivo <-  s.mat.log2fc[,ix.invivo]
pval.invivo <-  s.mat.pval.signed[,ix.invivo]
ix <- (n-2):n
s.mat.pval.signed.small <- s.mat.pval.signed[,ix]
s.mat.log2fc.small <- s.mat.log2fc[,ix]

average.specificity <- apply(s.mat.pval.signed.small,1,mean)
sum.specificity <- apply(s.mat.pval.signed.small,1,sum)
max.specificity <- apply(abs(s.mat.pval.signed.small),1,max)*sign(apply(s.mat.pval.signed.small,1,max))
min.specificity <- apply(abs(s.mat.pval.signed.small),1,min)*sign(apply(s.mat.pval.signed.small,1,min))
s.mat.pval.signed.small <- cbind(s.mat.log2fc.small,s.mat.pval.signed.small,fc.invivo,pval.invivo,
	average.specificity,sum.specificity,max.specificity,min.specificity)
x <- colnames(s.mat.pval.signed.small)
colnames(s.mat.pval.signed.small)[1:3] <- paste(sep="","fc.",x[1:3])
colnames(s.mat.pval.signed.small)[4:6] <- paste(sep="","pval.",x[4:6])
s.mat.pval.signed.small <- s.mat.pval.signed.small[order(abs(s.mat.pval.signed.small[,"min.specificity"]),decreasing=T),]
write.table(s.mat.pval.signed.small,file=paste(path.input.dated,"specificity_mat_",date.is,".xls",sep=""),sep="\t")
write.table(s.mat.pval.signed.small,file=paste(path.results,"specificity_mat_",date.is,".xls",sep=""),sep="\t")


