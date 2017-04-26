##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Apr 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# R code to create quality control plots for entire network pipeline
library(caTools) # to use integration under curve (Trapezoid Rule Numerical Integration)

path.input.gold.standard <- "input/th17/used_for_paper/gold_standard_lists/"
## path.output <- paste(sep="","results/validation_",date.is,"/")
path.output <- paste(sep="","results/validation/",date.is,"/")
system(paste(sep="","mkdir ", path.output))			
path.input.sam <- "input/th17/used_for_paper/"
path.input.combined <- paste(sep="","results/combinedAnalysis/",date.combine.data.run,"/")



# make output directory
if(any("P300" %in% tfs)){
	num.tfs.tmp <- num.tfs-1
} else {
	num.tfs.tmp <- num.tfs
}
if(filter.by.sam==TRUE){
	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs.tmp,"_sam_",z.abs.cut,"_deseq_cut_",deseq.pval.cut)
} else {
	add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs.tmp,"_sam_",0,"_deseq_cut_",deseq.pval.cut)
}
# add.str.2 <- paste(sep="","_add_p300_",any("P300" %in% tfs))
infl.nm <- paste(path.input.combined,"results_combine_data_exprsn",add.str,".Rdata",sep="")
# load input files
load(infl.nm)
keepers <- colnames(res[[1]][[1]][[1]])[1:16]
## gold standard files
th17.gns.f.nm <- paste(sep="",path.input.gold.standard,fl.nm.gs)
gs.gns <- read.delim(sep="\t",header=T,th17.gns.f.nm, as.is=T)
ix <- which(gs.gns$distance<cut.dist)
gs.gns <- toupper(gs.gns[,"gene_id"])[ix]

# gns.dataset <- rownames(res[[1]][[1]][[1]])
# gs.gns <- gs.gns$Gene[which(gs.gns$Gene %in% gns.dataset)]

#################################
for(i in 1:length(res)){ # go over activation,repression,or absolute
	tp <- names(res)[i]
	for(j in 1:length(tfs)){ # go over tfs
		tf <- tfs[j]
		for(k in 1:length(res[[tp]][[tf]])){ # go over chip types
			c.tp <- names(res[[tp]][[tf]])[k]
			if(tp=="repression"){
				res[[tp]][[tf]][[c.tp]][,"SAM"] <- -1*res[[tp]][[tf]][[c.tp]][,"SAM"]
			} else if (tp=="absolute"){
				res[[tp]][[tf]][[c.tp]][,"SAM"] <- abs(res[[tp]][[tf]][[c.tp]][,"SAM"])
			}
		}		
	}
}
# combine scores for multiple tfs
ix.add <- 1:15 # only add for K,C,R,I or combination (not for sam, deseq, etc)
o.sum.list <- list() # here we store the summed scores over all tfs
for(i in 1:length(res)){ # go over activation,repression,or absolute, or whole
	tp <- names(res)[i]
	o.sum.list[[tp]] <- list()
	for(k in 1:length(res[[tp]][[1]])){ # go over chip types (this is same for any tf)
		c.tp <- names(res[[tp]][[1]])[k]
		for(j in 1:length(tfs)){ # go over tfs
			tf <- tfs[j]
			if(j==1){ # if first tf start the summation
				o.sum.list[[tp]][[c.tp]] <- res[[tp]][[tf]][[c.tp]]
			} else {
				o.sum.list[[tp]][[c.tp]][,ix.add] <- o.sum.list[[tp]][[c.tp]][,ix.add] + res[[tp]][[tf]][[c.tp]][,ix.add]
			}
		}
	}
}



# convert combined scores for multiple tfs into relative ranks for ploting
o.sum.rel.list <- list() # here we store the summed scores over all tfs
for(i in 1:length(res)){ # go over activation,repression,or absolute
	tp <- names(res)[i]
	o.sum.rel.list[[tp]] <- list()
	for(k in 1:length(res[[tp]][[1]])){ # go over chip types (this is same for any tf)
		c.tp <- names(res[[tp]][[1]])[k]
		# prepare list of matrices to store data for relative ranks
		gns <- rownames(o.sum.list[[tp]][[c.tp]])
		o.sum.rel.list[[tp]][[c.tp]] <- matrix(0,nr=nrow(o.sum.list[[tp]][[c.tp]]),nc=ncol(o.sum.list[[tp]][[c.tp]]))
		colnames(o.sum.rel.list[[tp]][[c.tp]]) <- colnames(o.sum.list[[tp]][[c.tp]])
		for(j in 1:ncol(o.sum.list[[tp]][[c.tp]])){
			o.sum.rel.list[[tp]][[c.tp]][,j] <- gns[order(o.sum.list[[tp]][[c.tp]][,j],decreasing=T)]
		}
	}
}

prec <- list()
rec <- list()
tpr <- list()
fpr <- list()
AUCPR <- list()
AUCROC <- list()
for(type in names(o.sum.rel.list)){
	prec[[type]] <- list()
	rec[[type]] <- list()
	tpr[[type]] <- list()
	fpr[[type]] <- list()
	AUCPR[[type]] <- list()
	AUCROC[[type]] <- list()	
	for(c.tp in names(o.sum.rel.list[[type]])){
		prec[[type]][[c.tp]] <- list()
		rec[[type]][[c.tp]] <- list()
		tpr[[type]][[c.tp]] <- list()
		fpr[[type]][[c.tp]] <- list()
		AUCPR[[type]][[c.tp]] <- list()
		AUCROC[[type]][[c.tp]] <- list()		

		x.rel.nms <- o.sum.rel.list[[type]][[c.tp]]
		x.rel <- seq(from=0,to=1,length.out=dim(x.rel.nms)[1])
		gold.gns <-  gs.gns

		positive.rel.ranks <- list() # find what relative ranks gold std gns received
		for(j in 1:length(keepers)){
			positive.rel.ranks[[ keepers[j] ]] <- x.rel[which(x.rel.nms[,keepers[j]]%in%gs.gns)]
		}
		# calc precision vs. recal and roc curve

		for(i in 1:length(keepers)){
			x <- positive.rel.ranks[[ keepers[i] ]]
			neg.tot <- length(x.rel)-length(x)
			# fpr[[ keepers[i] ]] <- numeric(length(positive.rel.ranks[[i]]))
			tp <- numeric(length(x)) # true positive so far (elem j)
			fp <- numeric(length(x)) # false positive so far, fp.j= total positive minus true positive so far 
			fn <-  numeric(length(x)) # true negative so far, tn.j= total so far minus positive so far
			tn <- numeric(length(x))	
			for(j in 1:length(positive.rel.ranks[[i]])){
				tp[j] <- j # true positive so far (elem j)
				fp[j] <- (which(x.rel==x[j])-j) # false positive so far, fp.j= total positive minus true positive so far 
				fn[j] <-  length(x) - j # true negative so far, tn.j= total so far minus positive so far
				tn[j] <- neg.tot-fp[j]
			}
			prec[[type]][[c.tp]][[ keepers[i] ]] <- tp/(tp+fp)
			rec[[type]][[c.tp]][[ keepers[i] ]] <- tp/(tp+fn)
			tpr[[type]][[c.tp]][[ keepers[i] ]] <- c(tp/(tp+fn),1)
			fpr[[type]][[c.tp]][[ keepers[i] ]] <- c(fp/(fp+tn),1)
			AUCPR[[type]][[c.tp]][[ keepers[i] ]] <- trapz(x=rec[[type]][[c.tp]][[ keepers[i] ]], 
													y=prec[[type]][[c.tp]][[ keepers[i] ]])
			AUCROC[[type]][[c.tp]][[ keepers[i] ]] <- trapz(x=fpr[[type]][[c.tp]][[ keepers[i] ]], 
													y=tpr[[type]][[c.tp]][[ keepers[i] ]])
		}
	}
}

#################################
# convert single tf scores into relative ranks for ploting
o.rel.list <- list() # here we store the relative scores for each tfs
c.tp <- c.type.all.in.one.plot
# tp <- "activation"
tp <- comb.case
type <- comb.case
for(i in 1:length(tfs)){ # go over tfs
	tf <- tfs[i]	
	# prepare list of matrices to store data for relative ranks
	gns <- rownames(res[[tp]][[tf]][[c.tp]])
	o.rel.list[[tf]] <- matrix(0,nr=nrow(res[[tp]][[tf]][[c.tp]]),nc=ncol(res[[tp]][[tf]][[c.tp]]))
	colnames(o.rel.list[[tf]]) <- colnames(res[[tp]][[tf]][[c.tp]])
	for(j in 1:ncol(o.rel.list[[tf]])){
		o.rel.list[[tf]][,j] <- gns[order(res[[tp]][[tf]][[c.tp]][,j],decreasing=T)]
	}
}

# keep AUC values for combined scores
# c.tp <- "th17_minus_th0"
# type <- "activation"
AUCROC.5.way <- AUCROC[[type]][[c.tp]]
AUCPR.5.way <- AUCPR[[type]][[c.tp]]

prec <- list()
rec <- list()
tpr <- list()
fpr <- list()
AUCPR <- list()
AUCROC <- list()
for(k in 1:length(tfs)){ # go over tfs
	tf <- tfs[k]	
	x.rel.nms <- o.rel.list[[tf]]
	x.rel <- seq(from=0,to=1,length.out=dim(x.rel.nms)[1])
	gold.gns <- gs.gns

	positive.rel.ranks <- list() # find what relative ranks gold std gns received
	for(j in 1:length(keepers)){
		positive.rel.ranks[[ keepers[j] ]] <- x.rel[which(x.rel.nms[,keepers[j]]%in%gs.gns)]
	}	
	prec[[tf]] <- list()
	rec[[tf]] <- list()
	tpr[[tf]] <- list()
	fpr[[tf]] <- list()
	AUCPR[[tf]] <- list()
	AUCROC[[tf]] <- list()
	for(i in 1:length(keepers)){
		x <- positive.rel.ranks[[ keepers[i] ]]
		neg.tot <- length(x.rel)-length(x)
		# fpr[[ keepers[i] ]] <- numeric(length(positive.rel.ranks[[i]]))
		tp <- numeric(length(x)) # true positive so far (elem j)
		fp <- numeric(length(x)) # false positive so far, fp.j= total positive minus true positive so far 
		fn <-  numeric(length(x)) # true negative so far, tn.j= total so far minus positive so far
		tn <- numeric(length(x))	
		for(j in 1:length(positive.rel.ranks[[i]])){
			tp[j] <- j # true positive so far (elem j)
			fp[j] <- (which(x.rel==x[j])-j) # false positive so far, fp.j= total positive minus true positive so far 
			fn[j] <-  length(x) - j # true negative so far, tn.j= total so far minus positive so far
			tn[j] <- neg.tot-fp[j]
		}
		prec[[tf]][[ keepers[i] ]] <- tp/(tp+fp)
		rec[[tf]][[ keepers[i] ]] <- tp/(tp+fn)
		tpr[[tf]][[ keepers[i] ]] <- c(tp/(tp+fn),1)
		fpr[[tf]][[ keepers[i] ]] <- c(fp/(fp+tn),1)
		AUCPR[[tf]][[ keepers[i] ]] <- trapz(x=rec[[tf]][[ keepers[i] ]], 
												y=prec[[tf]][[ keepers[i] ]])
		AUCROC[[tf]][[ keepers[i] ]] <- trapz(x=fpr[[tf]][[ keepers[i] ]], 
												y=tpr[[tf]][[ keepers[i] ]])
	}
}

# helper function to capitalize tfs name (only first letter is capped)
capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s,1,1)),
                  {s <- substring(s,2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


f.nm.res.pdf <- paste(sep="",path.output, "vldtn_perTF_AllInOne2_",comb.case,"_",
						c.type.all.in.one.plot,add.str,"_gs_data_",gs.list,"_",gold.stdrd.date,".pdf")
# cls <- c(rainbow(length(AUCPR)+1))[-2]

# colors for tfs: "BATF"  "MAF"   "IRF4"  "STAT3" "RORC"
cls <- c(colors()[566],#BATF="royalblue1"
		 colors()[642],#MAF=violetred1
		 colors()[496],#IRF4=olivedrab3
		 colors()[512],#STAT3=orchid4
		 colors()[566])#RORC=royalblue1
pchs <-	c(18,#BATF=diamond (full)
		 15,#MAF=square (full)
		 17,#IRF4=triangle (full)
		 16)#STAT3=circle full)
rorc.pch <- "-" #RORC=dodgerblue	
# colors for sum scores bar plots: 
cls.bars <- c(rep(colors()[616],4),#k,c,r,i=steelblue1
		 rep(colors()[652],6),#two ways=yellow
		 rep(colors()[136],5))#three-fourways=firebrick3
# ix <- c("I","R","K","C","KC","RI","CRI","KRI","KCR","KCI","KCRI") 
ix <- colnames(res[[type]][[1]][[1]])[1:15]
# reorganize ix to make trends between comb levels clear
ix=c("I","R","C","K","CI","CR","KR","RI","KI","KC","CRI","KCR","KRI","KCI","KCRI")
  
# pch=20
# tfs <- tfs
# tfs <- tfs[-which(tfs=="P300")]
pdf(f.nm.res.pdf )
# plot results AUCPR
for(i in 1:length(ix)){ # go over datatypes
	d.type <- ix[i]
	for(k in 1:length(tfs)){ # go over tfs
		tf <- tfs[k]
		x <- AUCPR
		ylim <- c(0,0.4)
		xlim <- c(1,length(ix))
		if(k==1 & i==1){
			plot(y=x[[tf]][ix][1],x=k,pch=pchs[k],main=paste(sep="---",type,c.tp),xlim=xlim,ylim=ylim,col=cls[k],
			xlab="Datatype",ylab="Area Under Curve PR",cex=1.5,type="n",axes = FALSE)
			axis(2,cex.axis=2)
			text(1:length(ix), par("usr")[3], labels = ix, srt = 90, adj = 1, xpd = TRUE,cex=1)
		}
		if(k==1){
			barplot(AUCPR.5.way[[d.type]],space=i-0.5,width=1,col=cls.bars[i],add=T,axes = FALSE)
		}
		if(tf=="RORC"){
			points(y=x[[tf]][d.type],x=i,pch=rorc.pch,cex=1.5,col=cls[k])			
		} else {
			points(y=x[[tf]][d.type],x=i,pch=pchs[k],cex=1.2,col=cls[k])			
		}
	}
}
legend("topleft",capwords(rev(c(tolower(tfs),"combined"))),col=rev(c(cls,"black")),pch=rev(c(pchs,3,0)))

# legend("topleft",capwords(rev(c(tolower(tfs),"combined"))),col=rev(c(cls,"black")),pch=rev(c(rep(pch,length(tfs)),8)))
for(i in 1:length(ix)){ # go over datatypes
	d.type <- ix[i]
	for(k in 1:length(tfs)){ # go over tfs
		tf <- tfs[k]
		x <- AUCROC
		ylim <- c(0.5,1)
		xlim <- c(1,length(ix))
		if(k==1 & i==1){
			plot(y=x[[tf]][ix][1],x=k,pch=pchs[k],main=paste(sep="---",type,c.tp),xlim=xlim,ylim=ylim,col=cls[k],
			xlab="Datatype",ylab="Area Under Curve PR",cex=1.5,type="n",axes = FALSE)
			axis(2,cex.axis=2)
			text(1:length(ix), par("usr")[3], labels = ix, srt = 90, adj = 1, xpd = TRUE,cex=1)
		}
		if(k==1){
			barplot(AUCROC.5.way[[d.type]],space=i-0.5,width=1,col=cls.bars[i],add=T,axes = FALSE)
		}
		if(tf=="RORC"){
			points(y=x[[tf]][d.type],x=i,pch=rorc.pch,cex=1.5,col=cls[k])			
		} else {
			points(y=x[[tf]][d.type],x=i,pch=pchs[k],cex=1.2,col=cls[k])			
		}
	}
}
# legend("topleft",capwords(rev(c(tolower(tfs),"combined"))),col=rev(c(cls,"black")),pch=rev(c(pchs,3,0)))
dev.off()













