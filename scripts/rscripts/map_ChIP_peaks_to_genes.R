##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Nov 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

rm(list=ls())
debug=F
script.version=0.1
print.error <- function(){
	cat("
DESCRIPTIION:	
	map_ChIP_peaks_to_genes.R takes a MACS tab delimited file as input and produces two files: 1. genes.xls-a gene centered table 
	that stores the peaks that are proximal (within x[kb] of TSS) and distal 
	(within the gene's region + y[kb] upstream or downstream of TSS and TES respectively), and 2. peaks.xls
	a peak centered table that equals the MACS input table plus colums for proximal and distal targets 
	of each peak and their distance from TSS (if exist). 
	For gene.xls script also gives a poisson model based -log10(pval) score for proximal and genewide enrichment of peaks for each gene
	with at least one prox/dist peak associated with it. 
	
INPUT:
	1.input_file=path to MACS file\n
	2.refseq_table=path to refseq table (gives TSS/TES locations for all genes)\n 
	3.path_output=path to save genes.xls and peaks.xls\names
	4.refseq_table=path to refseq table (gives TSS/TES locations for all genes)\n 
	5.tss.dist: the absolute distance from TSS where we connect a peak to gene (for proximal peaks)\n
	6.gene.wide.dist: the absolute distance from TSS or TES where we connect a peak to gene (for peaks hitting anywhere in gene)\n
	7.effective.genome.size= the effective mappable genome size (for mm9 the value is 1.87e9)
	
EXAMPLE RUN: 
	Rscript map_ChIP_peaks_to_genes_v.1.R
	input_file=input/th17/used_for_paper/rawData/MACS_tab_delim_Sep_15_2011/SL971_SL970_peaks.xls
	refseq_table=rawData/UCSC_mm9_refseq_genes_Sep_15_2011.txt
	path_output=~/Desktop/
	tss.dist=5000
	gene.wide.dist=10000
	effective.genome.size=1.87e9
	macs.skip.lines=23
	
Please cite us if you used this script: 
	The transcription factor network regulating Th17 lineage specification and function.
	Maria Ciofani, Aviv Madar, Carolina Galan, Kieran Mace, Agarwal, Kim Newberry, Richard M. Myers, 
	Richard Bonneau and Dan R. Littman et. al. (in preperation)\n\n")
}

# retrieve args
if(debug==T){
	cmd.args <- c(
		# "input_file=input/th17/used_for_paper/rawData/MACS_Sep_15_2011/SL971_SL970_peaks.xls",
		"input_file=input/th17/used_for_paper/rawData/MACS_Sep_15_2011/SL3594_SL3592_peaks.xls",		
		"refseq_table=input/th17/used_for_paper/rawData/UCSC_mm9_refseq_genes_Sep_15_2011.txt",
		"path_output=/Users/aviv/Desktop/test_script/",
		"tss.dist=5000",
		"tes.dist=10000",
		"effective.genome.size=1.87e9",
		"gene.wide.dist=10000",
		"macs.skip.lines=23"
	)
} else {
	cmd.args <- commandArgs();      
}

if(length(grep("--version",cmd.args))){
	cat("version",script.version,"\n")
	q()
}

arg.names.cmd.line <- sapply(strsplit(cmd.args,"="),function(i) i[1])
args.val.cmd.line <- sapply(strsplit(cmd.args,"="),function(i) i[2])

arg.nms <- c("input_file","refseq_table","path_output","tss.dist",
			 "gene.wide.dist","effective.genome.size","macs.skip.lines")
arg.val <- character(length=length(arg.nms))
for(i in 1:length(arg.nms)){
	ix <- which(arg.names.cmd.line==arg.nms[i])
	if(length(ix)==1){
		arg.val[i] <- args.val.cmd.line[ix]
	} else {
		stop("######could not find ",arg.nms[i]," arg######\n\n",print.error())
		
	} 
}
if(debug==T){
	print(paste(arg.nms,"=",arg.val))
}
# the files here adhere to tab delim format
path.input.macs <- arg.val[1]
path.input.refseq <- arg.val[2]
path.output <- arg.val[3]
tss.dist <- as.numeric(arg.val[4])
gene.wide.dist <- as.numeric(arg.val[5])
effective.genome.size <- as.numeric(arg.val[6])
macs.skip.lines=as.numeric(arg.val[7])

if(debug==T){
# if(1){
	load("/Users/aviv/Desktop/r.u.RData")
	gn.nms <- rownames(r.u)
} else {
	##################
	### step 1 handle refseq table file (read, make unique)
	cat("reading refseq table",path.input.refseq,"\n")
	refseq <- read.delim(file=path.input.refseq)
	# refseq can have many transcripts for each gene
	# here i make it have only one transcript for each gene (the longest one)
	cat("for transcripts matching more than one gene keeping only the longest transcript\n")
	refseq$name2 <- as.character(refseq$name2)
	refseq$chrom <- as.character(refseq$chrom)
	refseq$strand <- as.character(refseq$strand)
	gn.nms <- unique(refseq$name2)
	#x <- sort(table(refseq$name2),decreasing=T)
	#gn.nms.non.unique <- names(x)[which(x>1)]
	#gn.nms.unique <- names(x)[which(x==1)]
	# create refseq unique r.u
	n <- length(gn.nms)
	r.u <- data.frame(cbind(chrom=rep("NA",n),strand=rep("NA",n),txStart=rep(0,n),txEnd=rep(0,n)),stringsAsFactors=FALSE,row.names=gn.nms)
	for(i in 1:n){
	  ix <- which(refseq$name2==gn.nms[i])
	  if(length(ix)==0) {
	    error("could not find gene", ng.nms[i], "in refseq table.  Bailing out...\n")
	  } else if (length(ix)>1){
	    l <- apply(refseq[ix,c("txStart","txEnd")],1,function(i) abs(i[1]-i[2]) )
	    l.max <- max(l)
	    ix <- ix[which(l==l.max)[1]]
	  }
	  r.u[gn.nms[i],c("chrom","strand","txStart","txEnd")] <- refseq[ix,c("chrom","strand","txStart","txEnd")]
	}
	r.u[,"txStart"] <- as.numeric(r.u[,"txStart"])
	r.u[,"txEnd"] <- as.numeric(r.u[,"txEnd"])

	# switch TSS and TES if we have chr "-"
	cat("correcting TSS, TES assiments in refseq table based on strand\n")
	ix <- which(r.u$strand=="-")
	tmp.tes <- r.u$txStart[ix]
	tmp.tss <- r.u$txEnd[ix]
	r.u[ix,c("txStart","txEnd")] <- cbind(tmp.tss,tmp.tes)
	# ############
}


#### step 2 read macs file
cat("reading Macs file",path.input.macs,"\n")

# read macs file
M <- read.delim(file=path.input.macs,skip=macs.skip.lines)
trgt.prox <- NA
trgt.dist <- NA
dtss.prox <- NA
dtss.dist <- NA
M.annot <- cbind(M,trgt.prox,trgt.dist,dtss.prox,dtss.dist)
M.annot[is.na(M.annot)] <- ""
# get the number of genes
m <- length(gn.nms)
if(debug==T){
	# m <- 100
}
f.nm.peak <- paste(sep="",path.output,"peaks.xls")
f.nm.gene <- paste(sep="",path.output,"genes.xls")
# for each genes in the expt (j goes over genes)
cat(file=f.nm.gene,sep="\t","Gene_ID","prox_n_peak","genewide_n_peak","gn_length","strand","prox_mean_pval","genewide_mean_pval",
			"prox_max_pval","genewide_max_pval",
			"prox_pois_model_pval","genewide_pois_model_pval","(peak_id,summit,d_TSS,d_TES,class,pval,fold_enrich,FDR),(..)\n")
peaks.summit <- (M[,"start"]+M[,"summit"])
cat(sep="","mapping MACS peaks to genes\n")
for (j in 1:m){
  if(j%%20==0){cat(".")}
  tss <- r.u[j,"txStart"]
  tes <- r.u[j,"txEnd"]
  gn.lngth <- abs(tss-tes)
  strand <- r.u[j,"strand"]
  chr <- r.u[j,"chrom"]
  if(strand=="+"){
    d.tss.all <- peaks.summit-tss
   	d.tes.all <- peaks.summit-tes
  } else {
    d.tss.all <- tss-peaks.summit
   	d.tes.all <- tes-peaks.summit
  }
  ix.distal <- sort(union( c(which( abs(d.tss.all) < gene.wide.dist ),which( abs(d.tes.all) < gene.wide.dist )),
               which( sign(d.tss.all) != sign(d.tes.all) )),decreasing=TRUE)
  ix.prox <- sort(which( abs(d.tss.all) < tss.dist ),decreasing=TRUE)	
  ix.prox <- ix.prox[which(M[ix.prox,"chr"]==chr)]    
  ix.distal <- ix.distal[which(M[ix.distal,"chr"]==chr)]
  n.peaks.prox <- length(ix.prox)
  n.peaks.distal <- length(ix.distal)
	# if there is at least one peak hitting gene j
  if(n.peaks.distal > 0){
    # for each peak (l goes over peaks)
    peaks.line <- paste(sep="","(")
    for (k in 1:n.peaks.distal){
      if(k>1){
        peaks.line <- paste(sep="",peaks.line,";(")
      }
      d.tss <- d.tss.all[ ix.distal[k] ]
      d.tes <- d.tes.all[ ix.distal[k] ]
      if(sign(d.tss)!=sign(d.tes)){
        class="intra"
      } else if(d.tss<=0){
        class="upstream"
      } else if(d.tss>0){
        class="downstream"
      }
      peaks.line <- paste(sep="",peaks.line,paste(sep=",", ix.distal[k], peaks.summit[ ix.distal[k] ], 
					d.tss, d.tes, class, M[ ix.distal[k] ,7], M[ ix.distal[k] ,8], M[ ix.distal[k] ,9] ),")")
	  # handle M.annot
	  dtss <- d.tss.all[ ix.distal[k] ]
	  prev.trgt <- M.annot[ix.distal[k],"trgt.dist"]
	  if(prev.trgt!=""){ # if peak already is with trgt choose one with min dtss
		if(abs(dtss)<abs(as.numeric(M.annot[ix.distal[k],"dtss.dist"]))){
			M.annot[ix.distal[k],"dtss.dist"] <- dtss
			M.annot[ix.distal[k],"trgt.dist"] <- gn.nms[j]
		}
	  } else { # if peak does not have trgt add current trgt
		M.annot[ix.distal[k],"dtss.dist"] <- dtss
		M.annot[ix.distal[k],"trgt.dist"] <- gn.nms[j]
	  }
	# calc the pval for these many peaks in proximal region and in gene-wide region
	}
	if(n.peaks.prox > 0){
		mean.pval.prox <- mean(M[ ix.prox ,7])
		max.pval.prox <- max(M[ ix.prox ,7])
		expected.num.peaks.genome.wide.prox <- length(which(M[,7]>=mean.pval.prox))
		# lambda = num peaks / genome size * searched region
		lambda.prox <- expected.num.peaks.genome.wide.prox/effective.genome.size*(tss.dist*2)		
		pval.pois.prox <- -log10(ppois(n.peaks.prox,lambda.prox,lower.tail=FALSE))
		# handle M.annot
		for(k in 1:n.peaks.prox){
			dtss <- d.tss.all[ ix.prox[k] ]
			prev.trgt <- M.annot[ix.prox[k],"trgt.prox"]
			if(prev.trgt!=""){ # if peak already is with trgt choose one with min dtss
				if(abs(dtss)<abs(as.numeric(M.annot[ix.prox[k],"dtss.prox"]))){
					M.annot[ix.prox[k],"dtss.prox"] <- dtss
					M.annot[ix.prox[k],"trgt.prox"] <- gn.nms[j]
				}
			} else { # if peak does not have trgt add current trgt
				M.annot[ix.prox[k],"dtss.prox"] <- dtss
				M.annot[ix.prox[k],"trgt.prox"] <- gn.nms[j]
			}
		}
	} else {
		mean.pval.prox <- "NA"
		max.pval.prox <- "NA"
		pval.pois.prox <- "NA"
	}
	mean.pval.distal <- mean(M[ ix.distal ,7])
	max.pval.distal <- max(M[ ix.distal ,7])
	expected.num.peaks.genome.wide.distal <- length(which(M[,7]>=mean.pval.distal))
	# lambda = num peaks / genome size * searched region
	lambda.distal <- expected.num.peaks.genome.wide.distal/effective.genome.size*(gn.lngth+gene.wide.dist*2)
    pval.pois.distal <- -log10(ppois(n.peaks.distal,lambda.distal,lower.tail=FALSE))
    # pring peaks for gene j    
    core.line <- paste(sep="\t",gn.nms[ j ],n.peaks.prox,n.peaks.distal,gn.lngth,strand,
							mean.pval.prox,mean.pval.distal,max.pval.prox,max.pval.distal,
							pval.pois.prox,pval.pois.distal)
    cat(file=f.nm.gene,append=TRUE,sep="",core.line,"\t",peaks.line,"\n")
  }
}
cat("Done!\n")

cat(file=f.nm.peak,colnames(M.annot),"\n",sep="\t")
write.table(file=f.nm.peak,append=T,col.names=F,row.names=F,M.annot,sep="\t")







