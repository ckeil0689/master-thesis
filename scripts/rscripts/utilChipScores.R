##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Dream3/4 pipeline (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

##################helper functions for chip_seq###################

# map from peak centered chip scores (possible multiple rows for each gene) to gne centered scores (one row per gene)
# input:
#       -pml is a list of matrices where rows are peaks and columns are peak stats
#       -gns are the genes names that we want to match with peaks
#       -peaks.num.per.expt is a vector giving the number of peaks recovered by MACS for each experiment
#       -mm9.effective.genome.size is the effective mappable genome (1.87e9 for mouse)
#       -upstream.gene.region,downstream.gene.region are the regions we search for peaks upstream/downstream of a gene (in kbp)
#       -c.names.m are the cloname for the output matrix m
# output:
#       - m: a matrix where rows are gns and columns are averaged stats for peaks that hit each gene
create.gene.based.chip.scores <- function(i,pml,gns,peaks.dist.per.expt,mm9.effective.genome.size,upstream.gene.region,downstream.gene.region,c.names.m){
  # for simplicity
  w <- pml[[i]]
  # a place holder matrix in which i keep results for one experiment until moved to list of matrices
  m <- matrix(0,nr=length(gns),nc=length(c.names.m))
  colnames(m) <- c.names.m
  rownames(m) <- gns
  ## lambda.per.bp <- peaks.num.per.expt[i]/mm9.effective.genome.size
  # get gene names for experiment i
  gn.nms.i <- unique(rownames(w))

  for (j in 1:length(gn.nms.i)){
    if((j%%500) == 0){cat(sep=""," ",round(j/length(gn.nms.i)*100,1),"% ")}
    ix <- which(rownames(w) == gn.nms.i[j])
    m[gn.nms.i[j],"gn_length"] <- w[ix[1],"gn_length"]
    ## here i add 20kb to represent the effective region of a gene where peaks where searched
    m[gn.nms.i[j],"kb_searched"] <- m[gn.nms.i[j],"gn_length"] + upstream.gene.region + downstream.gene.region
    ## lambda.per.kbp.searched <- m[gn.nms.i[j],"kb_searched"]*lambda.per.bp*1000
    m[gn.nms.i[j],"num_peaks"] <- length(ix)
    m[gn.nms.i[j],"pval_macs_mean"] <- mean(w[ix,"pval"])
    ## convert pval=10*log_10(pval) into pval=ln(pval)
    ## m[gn.nms.i[j],"pval_macs_mean"] <- m[gn.nms.i[j],"pval_macs_mean"]/10*log(10) #ie. -10log10(pval)/10*ln(10)
    # get lambda for poisson distribution (lambda = #peaks with at least as high mean pval over genome)/effective genome size* searched gene region (in bp)
    num.peaks.genome.wide <- length(which(peaks.dist.per.expt[[i]]>=m[gn.nms.i[j],"pval_macs_mean"]))
    lambda <- num.peaks.genome.wide/mm9.effective.genome.size*m[gn.nms.i[j],"kb_searched"]*1000
	m[gn.nms.i[j],"pval_poisson"] <- -log10(ppois(m[gn.nms.i[j],"num_peaks"],lambda,lower.tail=FALSE))
    # m[gn.nms.i[j],"pval_poisson"] <- -log10(ppois(m[gn.nms.i[j],"num_peaks"]-1,lambda,lower.tail=FALSE))
    ## final score for gene in this chip experiment is (#peaks)*(mean peak score)/(kb_searched)
    ## m[gn.nms.i[j],"pval_gene"] <- 2*(m[gn.nms.i[j],"pval_poisson"] + m[gn.nms.i[j],"pval_macs_mean"])
  }
  cat("\n")
  return(m)
}

# map from peak centered chip scores (possible multiple rows for each gene) to gne centered scores (one row per gene)
# input:
#       -pml.mat is a matrix where rows are peaks and columns are peak stats
#       -gns are the genes names that we want to match with peaks
# output:
#       - w: a matrix where rows are gns and columns are average of all peaks stats that hit gene
create.gene.based.chip.scores.old <- function(pml.mat,gns){
  c.names.m <- c("gn_length", "kb_searched","num_peaks",
		"S_pval_mean","S_fc_mean","S_peak_mean","S_gene")
  # for simplicity
  w <- pml.mat
  # a place holder matrix in which i keep results for one experiment until moved to list of matrices
  m <- matrix(0,nr=length(gns),nc=length(c.names.m))
  colnames(m) <- c.names.m
  rownames(m) <- gns

  # get gene names for experiment i
  gn.nms.i <- unique(rownames(w))

  for (j in 1:length(gn.nms.i)){
    if((j%%500) == 0){cat(sep=""," ",round(j/length(gn.nms.i)*100,1),"% ")}
    ix <- which(rownames(w) == gn.nms.i[j])
    m[gn.nms.i[j],"gn_length"] <- w[ix[1],"gn_length"]
    # here i add 20kb to represent the effective region of a gene where peaks where searched
    m[gn.nms.i[j],"kb_searched"] <- m[gn.nms.i[j],"gn_length"] + upstream.gene.region + downstream.gene.region
    m[gn.nms.i[j],"num_peaks"] <- length(ix)
    m[gn.nms.i[j],"S_pval_mean"] <- mean(w[ix,"S_pval"])
    m[gn.nms.i[j],"S_fc_mean"] <- mean(w[ix,"S_fc"])
    m[gn.nms.i[j],"S_peak_mean"] <- mean(w[ix,"S_peak"])		
    # final score for gene in this chip experiment is (#peaks)*(mean peak score)/(kb_searched)
    m[gn.nms.i[j],"S_gene"] <- (m[gn.nms.i[j],"num_peaks"] * m[gn.nms.i[j],"S_peak_mean"]) / m[gn.nms.i[j],"kb_searched"]
  }
  cat("\n")
#  gml[[ names(pml)[i] ]] <- m
  return(m)
}


# calculate a median corrected zscore for each value on v
calc.median.corrected.zscore <- function(v) {
	return( (v-median(v))/sd(v) )
}


# assign a level: 1,2,...,10 for each value in v based on ranked numbers in r and store in o
# i.e. r is 10,40,90,210,300,400,500,620,800,900 and two values in v are 55 and 420 then the o is 3,7

assign.level.rank.increasing <- function(v,r=seq(0,1,by=.1)) {
	# validate that all values in v are smaller then largest value in r (other wise some values will not be ranked!)
	if(max(v) > max(r)) { stop("#####\ncan't assign ranks to v as max(v)>max(r). Bailing out...\n#####\n")}
	if(min(v) < min(r)) { stop("#####\ncan't assign ranks to v as min(v)<min(r). Bailing out...\n#####\n")}
	o <- numeric(length=length(v))
	cur.rank <- 1
	for (i in 2:length(r)){
		# find indices to entries in v that belong to current rank
		if(i==2){
			### first bin is[ r[1],r[2] ]###			
			cur.ix <- which( (v <= r[i]) & (v>=r[i-1]))	
		} else {
			### other bins are ( r[i-1],r[i] ], note the paren on the left hand side ###			
			cur.ix <- which( (v <= r[i]) & (v>r[i-1]))	
		}
		o[cur.ix] <- cur.rank
		cur.rank <- cur.rank + 1
	}
	# sanity check that all values in v have been assigned a value between 1-10
	if( length(which( o < 1 | o > (length(r)-1) )) > 0 ) { stop("#####\nsome values in o have not been assigned a correct level. Bailing out...\n#####\n")}
	return(o)
}

# assign a level: 1,2,...,10 for each value in v based on ranked numbers in r and store in o
# i.e. r is 10,40,90,210,300,400,500,620,800,900 and two values in v are 55 and 420 then the o is 3,7

assign.level.rank.decreasing <- function(v,r=seq(1,0,by=-.1)) {
	# validate that all values in v are smaller then largest value in r (other wise some values will not be ranked!)
	if(max(v) > max(r)) { stop("#####\ncan't assign ranks to v as max(v)>max(r). Bailing out...\n#####\n")}
	if(min(v) < min(r)) { stop("#####\ncan't assign ranks to v as min(v)<min(r). Bailing out...\n#####\n")}
	o <- numeric(length=length(v))
	cur.rank <- 1
	for (i in 2:length(r)){
		# find indices to entries in v that belong to current rank
		if(i==2){
			### first bin is[ r[1],r[2] ]###			
			cur.ix <- which( (v >= r[i]) & (v<=r[i-1]))	
		} else {
			### other bins are ( r[i-1],r[i] ], note the paren on the left hand side ###			
			cur.ix <- which( (v >= r[i]) & (v<r[i-1]))	
		}
		o[cur.ix] <- cur.rank
		cur.rank <- cur.rank + 1
	}
	# sanity check that all values in v have been assigned a value between 1-10
	if( length(which( o < 1 | o > (length(r)-1) )) > 0 ) { stop("#####\nsome values in o have not been assigned a correct level. Bailing out...\n#####\n")}
	return(o)
}


get.sequences.parallel <- function(ix, myseqs) {
  seqs <- getSeq(Mmusculus, names=myseqs$chr[ix], start=myseqs$start[ix], end=myseqs$end[ix],strand=myseqs$strand[ix])
  return(seqs)
}

deconvolve.peaks <- function(peaks,gn.nms) {
  pksList <- strsplit(peaks,";")
  pksList <- lapply(pksList,function(i) gsub("(","",i,fixed=TRUE))
  pksList <- lapply(pksList,function(i) gsub(")","",i,fixed=TRUE))
  
  l <- get.2recursive.length.of.list(pksList) # get number of elems in list of lists
  
  
  peak.mat <- matrix(0,nr=l,nc=10)
  row.nms.peak.mat <- character(l)
  colnames(peak.mat) <- c("peak_id","gn_length","TSS","TES","Summit","d_TSS","d_TES","pval","FC","FDR") # gene names will be row names
  cnt <- 0
#w=0
  for(i in 1:length(pksList)) {
    for(j in 1:length(pksList[[i]])) {
      cnt <- cnt+1
      # x has elems (1)peak_id,(2)summit,(3)d_TSS,(4)d_TES,(5)class,(6)pval,(7)fold_enrich,(8)FDR)
      # i keep (1)summit,(2)d_TSS,(3)d_TES,(4)pval,(5)fold_enrich,(6)FDR)
      x <- as.numeric(strsplit(pksList[[i]][j],",")[[1]][c(1,2,3,4,6,7,8)])
      row.nms.peak.mat[cnt] <- gn.nms[i]
	  #AM need to fix tss and tes calc for - strand (not used later on)
      # need to calculate tss and tes from summit,d_tss, and d_tes (this is buggy; only calc correctly for + strand)
      tss <- x[2]-x[3] #AM
      tes <- x[2]-x[4] #AM
      lnt <- abs(x[3]-x[4])/1000 # in kb
      peak.mat[cnt,] <- c(x[1],lnt,tss,tes,x[2],x[3],x[4],x[5],x[6],x[7])
    }
  }
  rownames(peak.mat) <- row.nms.peak.mat
  return(peak.mat)
}

get.2recursive.length.of.list <- function(l) {
  cnt <- 0
  for(i in 1:length(l)) {
    for(j in 1:length(l[[i]])) {
      cnt <- cnt+1
    }
  }
  return(cnt)
}

# for each peak find where it landed:
# if upstream (d_tss < 0) then in - bps of TSS (i.e. return d_TSS)
# if downstream (d_tes > 0) then in + bps of TES (i.e. return d_TES)
# if in gene (d_tss>0&d_tes<0) then return (summit-TSS)/abs(TSS-TES)
# return a vector of hits with the above values for each peak
get.norm.hit.locations <- function(peaks.mat) {
  hits <- numeric(dim(peaks.mat)[1])
  for(i in 1:dim(peaks.mat)[1]) {
    summit <- peaks.mat[i,"Summit"]
    tss <- peaks.mat[i,"TSS"]
    tes <- peaks.mat[i,"TES"]
    d_tss <- peaks.mat[i,"d_TSS"]
    d_tes <- peaks.mat[i,"d_TES"]
    # if peak is upstream of gene
    if(d_tss <= 0) {
      hits[i] <- peaks.mat[i,"d_TSS"]
    } else if(d_tss > 0 & d_tes < 0) { # if peak is inside gene
      hits[i] <- d_tss/abs(tss-tes)
    } else if(d_tes >= 0) { # if peak is downstream of gene
      hits[i] <- peaks.mat[i,"d_TES"]
    } else {
      stop("my logic is flowed in get.norm.hit.locations. Bailing out")
    }
  }
  return(hits)
}
# plots the peaks hits on normalized gene space -10kb upstream inside and +10kb downstream
plot.peaks.histogram <- function(peaks.mat,median.gn.length){
  hits <- get.norm.hit.locations(peaks.mat)
  ix.upstream <- which(hits<=0)
  ix.downstream <- which(hits>=1)
  ix.inside <- which(hits>0 & hits<1)
  X <- c(hits[ix.upstream]/10000,hits[ix.inside],1+hits[ix.downstream]/10000)
#  h <- hist(X,nclass=100,col="darkgreen",main=paste(f.nm," binding histogram",sep=""), xaxt="n",xlab="peak position (normalized in genes)", ylab = "counts")
  h <- hist(X,nclass=100,plot=FALSE)
  ix <- which(h$mids<=1 & h$mids >=0)
  h$counts[ix] <- h$counts[ix]/median.gn.length
  plot(h,col="darkgreen",main=paste(f.nm," binding histogram",sep=""), xaxt="n",xlab="peak position (normalized in genes)", ylab = "counts")
  axis(1,at=c(-1,0,1,2),labels=c(-10000,0,1,10000))
}

# plots the peaks hits on normalized gene space -10kb upstream inside and +10kb downstream
plot.peaks.on.gene <- function(peaks.mat.per.gene,gn.length,gn.nm,expt.nm){
  if (dim(peaks.mat.per.gene)[1]>0){
    hits <- get.norm.hit.locations(peaks.mat.per.gene)
    ix.upstream <- which(hits<=0)
    ix.downstream <- which(hits>=1)
    ix.inside <- which(hits>0 & hits<1)
    X <- c(hits[ix.upstream]/10000,hits[ix.inside],1+hits[ix.downstream]/10000)
#  h <- hist(X,nclass=100,col="darkgreen",main=paste(f.nm," binding histogram",sep=""), xaxt="n",xlab="peak position (normalized in genes)", ylab = "counts")
    h <- hist(X,nclass=1000,plot=FALSE)
    wr.hits.mids <- h$mids[which(h$counts>0)]
    plot(wr.hits.mids,rep(0,length(wr.hits.mids)),pch="|",cex=2,col="red",xlim=c(-1,2),ylim=c(0,1),yaxt="n",
         main=paste("Binding sites for tf ",expt.nm," on ",gn.nm,sep=""), xaxt="n",xlab=paste("peak position ([0,1]->[0,",round(gn.length/10^3,1),"kb]",sep=""),, ylab = "counts")
    axis(1,at=c(-1,0,1,2),labels=c(-10000,0,1,10000))
    legend("topright",paste(dim(peaks.mat.per.gene)[1],"hits"),pch="|",col="red")
  } else {
    plot(1,0,pch="",cex=2,col="red",xlim=c(-1,2),ylim=c(0,1),yaxt="n",
         main=paste("Binding sites for tf ",expt.nm," on ",gn.nm,sep=""), xaxt="n",xlab="peak position (NA)", ylab = "counts")
    axis(1,at=c(-1,0,1,2),labels=c(-10000,0,1,10000))
    legend("topright",paste(dim(peaks.mat.per.gene)[1],"hits"),pch="|",col="red")
  }
}

# a hack: takes from list of peak mats find if gene from gn.list have a peak, if so calc gene length
get.genes.length <- function(gn.list,peaks.mat.List){
  gn.lengths <- rep(1,length(gn.list))
  i <- 0
  seen.gns <- numeric(length(gn.list))
  while(any(seen.gns==0) & (i < length(peaks.mat.List))) {
    i <- i+1
    for(j in 1:length(gn.list)) {
      gn <- gn.list[j]
      ix <- which(toupper(rownames(peaks.mat.List[[i]])) %in% gn)
      if(length(ix)>0) {
        seen.gns[j] <- 1
        gn.lengths[j] <- abs(peaks.mat.List[[i]][ix[1],"TSS"]-peaks.mat.List[[i]][ix[1],"TES"])
      }
    }
  }
  return(gn.lengths)
}

# plots the peaks hits on normalized gene space -10kb upstream inside and +10kb downstream
plot.peaks.on.gene.multiple.tfs <- function(peaks.mat.per.gene.list,gn.length,gn.nm,cls,legnd.cex=1){
  num.tfs <- length(peaks.mat.per.gene.list)
  all.hits <- numeric()
  all.hits.tf.map <- numeric() # each tf is encoded by a number 1, 2, 3,...
  tfs.y.axis.loc <- seq(0,1,by=0.02)
  tfs.pch <- 1:num.tfs
  expt.nms <- names(peaks.mat.per.gene.list)
  for(i in 1:num.tfs) {
    num.peaks <- dim(peaks.mat.per.gene.list[[i]])[1]
    if (num.peaks>0){
      hits <- get.norm.hit.locations(peaks.mat.per.gene.list[[i]])
      all.hits <- c(all.hits,hits)
      all.hits.tf.map <- c(all.hits.tf.map,rep(i,length(hits)))
    }
  }
  ix.downstream <- which(all.hits<0)
  ix.upstream <- which(all.hits>1)
  if(length(ix.downstream)>0) {  all.hits[ix.downstream] <- all.hits[ix.downstream]/10^4   }
  if(length(ix.upstream)>0) {  all.hits[ix.upstream] <- 1+all.hits[ix.upstream]/10^4   }

  
  tfs.ix <- unique(all.hits.tf.map)
  for(j in 1:length(tfs.ix)) {
    ix <- which(all.hits.tf.map==tfs.ix[j])
    if(j==1){
      plot(all.hits[ix],rep(tfs.y.axis.loc[tfs.ix[j]],length(all.hits[ix])),pch="|",cex=2,cex.lab=1.3,col=cls[tfs.ix[j]],xlim=c(-1,2),ylim=c(0,1),yaxt="n",ylab="",
         main=paste("Binding sites for ",gn.nm,sep=""), xaxt="n",xlab=paste("peak position ([0,1]->[0,",round(gn.length/10^3,1),"kb])",sep=""))
      axis(1,at=c(-1,0,1,2),labels=c(-10000,0,1,10000),cex.axis=1.3)
     } else {
      points(all.hits[ix],rep(tfs.y.axis.loc[tfs.ix[j]],length(all.hits[ix])),pch="|",cex=2,col=cls[tfs.ix[j]],xlim=c(-1,2),ylim=c(0,1),yaxt="n", xaxt="n")
    }
  }
  # get number of hits per tf (per experiment)
  hits.vec <- numeric(num.tfs)
  for(i in 1:num.tfs){
    hits.vec[i] <- length(which(all.hits.tf.map == i))
  }
  legend("topright",paste(expt.nms," hits = ",hits.vec,sep=""),pch="|",col=cls,title = "TF map",cex=legnd.cex)
}



plot.hits.in.hist <- function(bg.scores,conf.vals,color="red",linewidth=2,main="",xLable="confidence score",whr.lgnd="topright",cx=1) {
  hist(bg.scores,nclass=1000,main="clr scores",xlab = "clr pseudo z score")
  zero.cnt <- 0
  for(j in 1:length(conf.vals)){
    lines(c(conf.vals[j],conf.vals[j]),c(0,1*10^6),lwd=2,col="red")
    if(conf.vals[j] == 0) {
      zero.cnt <- zero.cnt + 1
    }
  }
  ix <- sort(conf.vals,decreasing=TRUE,index.return=TRUE)$ix
  legend(whr.lgnd,paste(names(conf.vals[ix]),"=",round(conf.vals[ix],2)),lty=1,lwd=2,col="red",cex=1)
}










