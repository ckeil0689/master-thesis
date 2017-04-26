##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Apr 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '



## to get pvalues for clusters run some simulations
run.perms <- function(n.b=n.b,len.th17.per.chr=len.th17.per.chr,tf.th17.per.chr=tf.th17.per.chr,pval.th17.per.chr=pval.th17.per.chr,
         tf.th0.per.chr=tf.th0.per.chr,len.th0.per.chr=len.th0.per.chr,pval.th0.per.chr=pval.th0.per.chr,
          min.chr.loc=min.chr.loc,max.chr.loc=max.chr.loc,nms.th17.table,n.p){
  o <- mclapply(1:n.b,run.one.perm,len.th17.per.chr=len.th17.per.chr,tf.th17.per.chr=tf.th17.per.chr,pval.th17.per.chr=pval.th17.per.chr,
         tf.th0.per.chr=tf.th0.per.chr,len.th0.per.chr=len.th0.per.chr,pval.th0.per.chr=pval.th0.per.chr,
          min.chr.loc=min.chr.loc,max.chr.loc=max.chr.loc,nms.th17.table,mc.cores=n.p)
  th17.perm.table <- t(sapply(o,function(i) i[["th17.counts"]]))
  th0.perm.table <- t(sapply(o,function(i) i[["th0.counts"]]))
  th17.pval.table <- t(sapply(o,function(i) i[["th17.pval"]]))
  th0.pval.table <- t(sapply(o,function(i) i[["th0.pval"]]))
  return(list(th17.perm.table=th17.perm.table,
              th0.perm.table=th0.perm.table,
              th17.pval.table=th17.pval.table,
              th0.pval.table=th0.pval.table))
}
run.one.perm <- function(b,len.th17.per.chr=len.th17.per.chr,tf.th17.per.chr=tf.th17.per.chr,pval.th17.per.chr=pval.th17.per.chr,
         tf.th0.per.chr=tf.th0.per.chr,len.th0.per.chr=len.th0.per.chr,pval.th0.per.chr=pval.th0.per.chr,
          min.chr.loc=min.chr.loc,max.chr.loc=max.chr.loc,nms.th17.table){
  cat("working on bootstrap", b, "\n")
  ## make permuted rulers
  r.p.th17 <- list()
  r.p.th0 <- list()
  for(i in 1:length(chr)){
    r.p.th17[[i]] <- list()
    r.p.th0[[i]] <- list()
    r.p.th17[[i]][["s"]] <- sort(runif(len.th17.per.chr[i],min=min.chr.loc[i],max=max.chr.loc[i]))
    r.p.th17[[i]][["tf.to.s"]] <- sample(tf.th17.per.chr[[i]])
    r.p.th17[[i]][["pval"]] <- sample(pval.th17.per.chr[[i]])
    r.p.th0[[i]][["s"]] <- sort(runif(len.th0.per.chr[i],min=min.chr.loc[i],max=max.chr.loc[i]))
    r.p.th0[[i]][["tf.to.s"]] <- sample(tf.th0.per.chr[[i]])
    r.p.th0[[i]][["pval"]] <- sample(pval.th0.per.chr[[i]])    
  }
  ## calc perm rulers
  r.p.th17 <- assign.clusters.ids.to.peaks(r.p.th17)
  ll.p.th17 <- create.cluster.list.perm(r.p.th17,n.p)
  r.p.th0 <- assign.clusters.ids.to.peaks(r.p.th0)
  ll.p.th0 <- create.cluster.list.perm(r.p.th0)

  ## get talbe of each cluster occurances
  c.table.perm.th17 <- numeric(length=length(nms.th17.table))
  c.table.perm.th0 <- numeric(length=length(nms.th17.table))
  names(c.table.perm.th17) <- names(c.table.perm.th0) <- nms.th17.table
  
  clust.types.th17.perm <- sapply(ll.p.th17,function(i) paste( sort(unique(i[["tf.to.s"]])),collapse="_") )
  x <- table(clust.types.th17.perm)
  c.table.perm.th17[names(x)] <- x
  ## get pvalue medians for each clusters
  pval.list.th17 <- lapply(ll.p.th17,function(i) i[["pval"]] )
  names(pval.list.th17) <- clust.types.th17.perm
  pval.median.th17 <- unlist(sapply(1:length(nms.th17.table),function(i) if(length(which(names(pval.list.th17)==nms.th17.table[i])>0)){
              median(unlist(pval.list.th17[which(names(pval.list.th17)==nms.th17.table[i])]))}else{0} ))
  names(pval.median.th17) <- nms.th17.table
  
  clust.types.th0.perm <- sapply(ll.p.th0,function(i) paste( sort(unique(i[["tf.to.s"]])),collapse="_") )
  x <- table(clust.types.th0.perm)
  c.table.perm.th0[names(x)] <- x
  pval.list.th0 <- lapply(ll.p.th0,function(i) i[["pval"]] )
  names(pval.list.th0) <- clust.types.th0.perm
  pval.median.th0 <- unlist(sapply(1:length(nms.th17.table),function(i) if(length(which(names(pval.list.th0)==nms.th17.table[i])>0)){
              median(unlist(pval.list.th0[which(names(pval.list.th0)==nms.th17.table[i])]))}else{0} ))
  names(pval.median.th0) <- nms.th17.table
  
  return(list(th17.counts=c.table.perm.th17,
              th17.pval=pval.median.th17,
              th0.counts=c.table.perm.th0,
              th0.pval=pval.median.th0))
}


## print sungear gene names 
## input:
##     - ll: a cluster list
##     - f.nm: a file name (include path) to where you want files to print
## output
##     - a fasta file with > gene_nm
print.ll <- function(ll,f.nm="ll.xls"){
  cat(file=f.nm,c("loci","chromosome","tf.to.s","start","end","summit","trgts","d.tss"),sep="\t")
  cat(file=f.nm,"\n",append=TRUE)
  for(i in 1:length(ll)){
    loci <- ll[[i]][["l"]]
    chromosome <- ll[[i]][["chr"]][1]
    tf.to.s <- paste(sep="",sort(unique(unlist(ll[[i]]["tf.to.s"]))),collapse="_")
    start <- ll[[i]][["start"]]
    end <- ll[[i]][["end"]]
    summit <-  paste(sep="",ll[[i]][["s"]],collapse="_")
    trgts <- paste(sep="",unlist(ll[[i]]["trgts"]),collapse="_")
    d.tss <- paste(sep="",unlist(ll[[i]]["d.tss"]),collapse="_")
    cat(file=f.nm,c(loci,chromosome,tf.to.s,start,end,summit,trgts,d.tss),sep="\t",append=TRUE)
    cat(file=f.nm,"\n",append=TRUE)
  }
}


## pretty print (tab deim) to file requested elements out of chip cluster list, ll.
## input:
##     - ll: a cluster list
##     - f.nm: a file name (include path) to where you want files to print
##     - tfs: a list of the tfs we want to print the file for (the same as the tfs used for the peak clustering)
## output
##     - a tab delim file with clusers as rows and elems tab delim for each cluster
print.ll.verbose <- function(ll,f.nm="ll.xls",tfs){

  ## cat(file=f.nm,c("loci","chromosome","tf.to.s","start","end","summit","summit_names",
  ##       "pval_max_irf4","pval_max_batf","pval_max_maf","pval_max_stat3","pval_max_rorc","pval_max_p300","pval_mean","trgts"),sep="\t")
  cat(file=f.nm,c("loci","chromosome","tf.to.s","start","end","summit","summit_names",
        paste(sep="","pval_max_",tfs),"pval_mean","d.tss","trgts"),sep="\t")
  cat(file=f.nm,"\n",append=TRUE)
  for(i in 1:length(ll)){
    loci <- ll[[i]][["l"]]
    chromosome <- ll[[i]][["chr"]][1]
    tf.to.s <- unlist(ll[[i]]["tf.to.s"])
    tf.to.s.unique <- paste(sep="",sort(unique(tf.to.s)),collapse="_")
    start <- ll[[i]][["start"]]
    end <- ll[[i]][["end"]]
    summit <-  paste(sep="",ll[[i]][["s"]],collapse="_")
    summit.names <- paste(sep="",tf.to.s,collapse="_")
     # d.tss.min <- rep(NA,length(tfs))
    pval.max <- numeric(length(tfs))
    names(pval.max) <- tfs
    trgts <- paste(sep="",unlist(ll[[i]]["trgts"]),collapse="_")
	d.tss <- paste(sep="",unlist(ll[[i]]["d.tss"]),collapse="_")
    for(j in 1:length(tfs)){
      tf <- tfs[j]
      ix <- grep(tf,tf.to.s)
      if(length(ix)>=1){
        pval.max[tf] <- ll[[i]]$pval[ix][which.max(ll[[i]]$pval[ix])]
        # if(trgts[1]!=""){
        #   tss <- as.numeric(strsplit(ll[[i]]$d.tss[ix],"_")[[1]])
        #   if(length(tss)>0){ d.tss.min[tf] <- min(abs(tss)) }
        # }
      }
    }
    pval.mean <- sum(pval.max)/length(which(pval.max!=0))
    # ix.not.na <- which(!is.na(d.tss.min))
    # if(length(ix.not.na)>=1){
    #   d.tss.mean <- sum(d.tss.min[ix.not.na])/length(ix.not.na)
    # } else {
    #   d.tss.mean <- NA
    # }
    ## cat(file=f.nm,c(loci,chromosome,tf.to.s.unique,start,end,summit,summit.names,
    ##       pval.max.irf4,pval.max.batf,pval.max.maf,pval.max.stat3,pval.max.rorc,pval.max.p300,pval.mean,trgts),sep="\t",append=TRUE)
    cat(file=f.nm,c(loci,chromosome,tf.to.s.unique,start,end,summit,summit.names,
          pval.max,pval.mean,d.tss,trgts),sep="\t",append=TRUE)
    cat(file=f.nm,"\n",append=TRUE)
  }
}
## pretty print (tab deim) to file requested elements out of chip cluster list, ll.
## input:
##     - ll: a cluster list
##     - f.nm: a file name (include path) to where you want files to print
##     - tfs: a list of the tfs we want to print the file for (the same as the tfs used for the peak clustering)
## output
##     - a tab delim file with clusers as rows and elems tab delim for each cluster
print.ll.verbose.pioneer <- function(ll,f.nm="ll.xls",tfs){

  ## cat(file=f.nm,c("loci","chromosome","tf.to.s","start","end","summit","summit_names",
  ##       "pval_max_irf4","pval_max_batf","pval_max_maf","pval_max_stat3","pval_max_rorc","pval_max_p300","pval_mean","trgts"),sep="\t")
  cat(file=f.nm,c("loci","chromosome","tf.to.s","start","end","summit","summit_names",
        paste(sep="","pval_max_",tfs),"pval_mean","d.tss","trgts"),sep="\t")
  cat(file=f.nm,"\n",append=TRUE)
  for(i in 1:length(ll)){
    loci <- ll[[i]][["l"]]
    chromosome <- ll[[i]][["chr"]][1]
    tf.to.s <- unlist(ll[[i]]["tf.to.s"])
	tf.to.s.unique <- paste(sep="",sort(unique(tf.to.s)),collapse="_")
	lin.to.s <- unlist(ll[[i]]["expt"])
	lin.to.s[grep("Th17",lin.to.s,ignore.case=T)] <- "IRF4_Th17_SL2872_SL2876"
	lin.to.s <- sapply(lin.to.s,function(i) paste(sep="",strsplit(i,"_")[[1]][1:2],collapse="_"))
	lin.to.s.2 <- sapply(lin.to.s,function(i) strsplit(i,"_")[[1]][2])
	lin.to.s.unique <- paste(sep="",sort(unique(lin.to.s.2)),collapse="_")
    start <- ll[[i]][["start"]]
    end <- ll[[i]][["end"]]
    summit <-  paste(sep="",ll[[i]][["s"]],collapse="_")
    summit.names <- paste(sep="",lin.to.s.unique,collapse="_")
     # d.tss.min <- rep(NA,length(tfs))
    pval.max <- numeric(length(tfs))
    names(pval.max) <- tfs
    trgts <- paste(sep="",unlist(ll[[i]]["trgts"]),collapse="_")
	d.tss <- paste(sep="",unlist(ll[[i]]["d.tss"]),collapse="_")
    for(j in 1:length(lin.to.s)){
      tf <- lin.to.s[j]
      ix <- grep(paste(sep="",tf,"$"),lin.to.s,perl=T)
      if(length(ix)>=1){
        pval.max[tf] <- ll[[i]]$pval[ix][which.max(ll[[i]]$pval[ix])]
      }
    }
    pval.mean <- sum(pval.max)/length(which(pval.max!=0))
    # ix.not.na <- which(!is.na(d.tss.min))
    # if(length(ix.not.na)>=1){
    #   d.tss.mean <- sum(d.tss.min[ix.not.na])/length(ix.not.na)
    # } else {
    #   d.tss.mean <- NA
    # }
    ## cat(file=f.nm,c(loci,chromosome,tf.to.s.unique,start,end,summit,summit.names,
    ##       pval.max.irf4,pval.max.batf,pval.max.maf,pval.max.stat3,pval.max.rorc,pval.max.p300,pval.mean,trgts),sep="\t",append=TRUE)
    cat(file=f.nm,c(loci,chromosome,tf.to.s.unique,start,end,summit,summit.names,
          pval.max,pval.mean,d.tss,trgts),sep="\t",append=TRUE)
    cat(file=f.nm,"\n",append=TRUE)
  }
}
## pretty print (tab deim) to file requested elements out of chip cluster list, ll.
## input:
##     - ll: a cluster list
##     - f.nm: a file name (include path) to where you want files to print
##     - tfs: a list of the tfs we want to print the file for (the same as the tfs used for the peak clustering)
## output
##     - a tab delim file with clusers as rows and elems tab delim for each cluster
print.ll.verbose.all <- function(ll,f.nm="ll.xls"){
  options(digits=5)
  cat(file=f.nm,names(ll[[1]]),sep="\t")
  cat(file=f.nm,"\n",append=TRUE)
  for(i in 1:length(ll)){
	line <- ll[[i]][[1]] # put MTL number
	line <- paste(line,ll[[i]][[2]][1],sep="\t")
	for(j in 3:length(ll[[i]])){
		val <- ll[[i]][[j]]
		if(length(val) == 1){
			line <- paste(line,val,sep="\t")
		} else {
			line <- paste(line,paste(sep="",unlist(ll[[i]][j]),collapse="_"),sep="\t")
		}
	}
  	cat(file=f.nm,line,"\n",append=TRUE)	
  }
}


## ploting rutines:
plot.qc.2 <- function(c.table.th17,c.table.th0,cut.th17=cut.th17, ylim.pval=ylim.pval){
  labels <- names(c.table.th17)
  ix.th17 <- which(c.table.th17>cut.th17)
  mp <- barplot(c.table.th17[ix.th17],col="grey",main="Th-17",ylab="number of clusters",axes = FALSE, axisnames = FALSE )
  text(mp, par("usr")[3], labels = labels[ix.th17], srt = 45, adj = 1, xpd = TRUE,cex=.5)
  ## segments(mp, hh, mp, hh + 2*sqrt(1000*hh/100), col = mybarcol, lwd = 1.5) (error bars are too large)
  axis(2)
  ## plot bar plot of clusters numbers
  ##labels <- names(c.table.th0)
  mp <- barplot(c.table.th0[ix.th17],col="grey",main="Th-0",ylab="number of clusters",axes = FALSE,axisnames = FALSE)
  text(mp, par("usr")[3], labels = labels[ix.th17], srt = 45, adj = 1, xpd = TRUE,cex=.5)
  ## segments(mp, hh, mp, hh + 2*sqrt(1000*hh/100), col = mybarcol, lwd = 1.5) (error bars are too large)
  axis(2)
}
plot.qc.1 <- function(c.table.th17,
                      ctable.th0,
                      cut.th17=500,
                      cut.th0=500,
                      ylim.pval=c(0,2000),
                      ylim.pval.histones=c(0,3000),
                      plot.histone.marks=TRUE,
					  outliers=FALSE){
  ## plot bar plot of clusters numbers
  labels <- names(c.table.th17)
  ix.th17 <- which(c.table.th17>cut.th17)
  ix.th0.plot.boxplot <- which(c.table.th0>cut.th0)
  mp <- barplot(c.table.th17[ix.th17],col="grey",main="Th-17",ylab="number of clusters [10^3]",axes = FALSE, axisnames = FALSE,plot=F)
  ix <- seq(from=0,to=ceiling(round(max(c.table.th17[ix.th17])/1000,digits=0)))
  barplot(c.table.th17[ix.th17],col="grey",main="Th-17",ylab="number of clusters [10^3]",axes = FALSE, axisnames = FALSE,ylim=c(0,max(ix*1000)))
  text(mp, par("usr")[3], labels = labels[ix.th17], srt = 45, adj = c(1,2), xpd = TRUE,cex=.5)
  # axis(2,las=1,at=ix*1000,labels=paste(sep="",ix,"k"),cex.axis=0.8)
  axis(2,las=1,at=ix*1000,labels=ix,cex.axis=0.7)
  
  ## plot bar plot of clusters numbers
  ##labels <- names(c.table.th0)
  mp <- barplot(c.table.th0[ix.th17],col="grey",main="Th-0",ylab="number of clusters [10^3]",axes = FALSE,axisnames = FALSE,plot=F)
  ix <- seq(from=0,to=ceiling(round(max(c.table.th0[ix.th17[names( ix.th0.plot.boxplot)]])/1000,digits=0)))
  mp <- barplot(c.table.th0[ix.th17[names( ix.th0.plot.boxplot)]],col="grey",main="Th-0",ylab="number of clusters [10^3]",
		axes = FALSE,axisnames = FALSE,ylim=c(0,max(ix*1000)),xlim=c(0,length(ix.th17)))
  text(mp, par("usr")[3], labels = labels[ix.th17[names( ix.th0.plot.boxplot)]], srt = 45, adj = c(1,2), xpd = TRUE,cex=.5)
  axis(2,las=1,at=ix*1000,labels=ix,cex.axis=0.7)

  ## plot boxplot of clusters median pval 
  mp <- boxplot(pvals.th17[ix.th17],col="grey",main="Th-17",ylab="binding strength [-log10(pval)]",outline=outliers,
                axes = FALSE, axisnames = FALSE,pch=20,cex=.3,ylim=ylim.pval,xlim=c(0,length(ix.th17)),notch=T)
  text(1:length(mp$n), par("usr")[3], labels = labels[ix.th17], srt = 45, adj = 1, xpd = TRUE,cex=.5)
  axis(2,las=1,cex.axis=0.7)
 
  ylim.2 <- c(300,300,250,120,100)
  for(i in 1:length(tfs)){
	tf <- tfs[i]
	mp <- boxplot(pvals.th17.per.tf[[tf]][ix.th17],col="grey",main=paste(sep="",tf," in Th-17"),ylab="binding strength [-log10(pval)]",outline=outliers,
                axes = FALSE, axisnames = FALSE,pch=20,cex=.3,xlim=c(0,length(ix.th17)),plot=F)
	ylim.i <- c(ylim.pval[1],min(ylim.pval[2],max(mp$stats[which(!is.na(mp$stats))])))
	ylim.i[2] <- ylim.2[i]
	mp <- boxplot(pvals.th17.per.tf[[tf]][ix.th17],col="grey",main=paste(sep="",tf," in Th-17"),ylab="binding strength [-log10(pval)]",outline=outliers,
                axes = FALSE, axisnames = FALSE,pch=20,cex=.3,xlim=c(0,length(ix.th17)),ylim=ylim.i,notch=T)
	text(1:length(mp$n), par("usr")[3], labels = labels[ix.th17], srt = 45, adj = 1, xpd = TRUE,cex=.5)
	axis(2,las=1,cex.axis=0.7)
  }


  ## plot boxplot of clusters median pval
  mp <-  boxplot(pvals.th0[names(ix.th17[names( ix.th0.plot.boxplot)])],col="grey",main="Th-0",ylab="binding strength [-log10(pval)]",axes = FALSE,outline=outliers,
    axisnames = FALSE,pch=20,cex=.3,ylim=ylim.pval+.1,xlim=c(0,length(ix.th17[names( ix.th0.plot.boxplot)])),plot=F)
  ylim <-  c(ylim.pval[1],min(ylim.pval[2],max(mp$stats[which(!is.na(mp$stats))])))
  mp  <-  boxplot(pvals.th0[names(ix.th17[names( ix.th0.plot.boxplot)])],col="grey",main="Th-0",ylab="binding strength [-log10(pval)]",axes = FALSE,outline=outliers,
    axisnames = FALSE,pch=20,cex=.3,ylim=ylim,notch=T,xlim=c(0,length(ix.th17)))
  # mp  <-  boxplot(pvals.th0[names(ix.th17[names( ix.th0.plot.boxplot)][1])],col="grey",main="Th-0",ylab="binding strength [-log10(pval)]",axes = FALSE,outline=outliers,
  #   axisnames = FALSE,pch=20,cex=.3,ylim=ylim,xlim=c(0.7,length(ix.th17[names( ix.th0.plot.boxplot)])+0.1),notch=T)
  # for(i in 2:length(ix.th0.plot.boxplot)){
  #   ix <- ix.th0.plot.boxplot[i]
  #   boxplot(pvals.th0[names(ix)],at=which(names(ix.th17[names( ix.th0.plot.boxplot)])==names(ix)),col="grey",axes = FALSE,outline=outliers,
  #           axisnames = FALSE,pch=20,cex=.3,ylim=ylim.pval,xlim=c(1,15),add=T,notch=T)
  # }
  axis(2,las=1,cex.axis=0.7)
  text(1:length(ix.th17[names( ix.th0.plot.boxplot)]), par("usr")[3], labels = labels[ix.th17[names( ix.th0.plot.boxplot)]], srt = 45, adj = 1, xpd = TRUE,cex=.5)


  for(i in 1:length(tfs)){
	tf <- tfs[i]
	l.vec <- sapply(pvals.th0.per.tf[[tf]][names( ix.th0.plot.boxplot)],length)
	if(any(l.vec>0)){
		mp <- boxplot(pvals.th0.per.tf[[tf]][names( ix.th0.plot.boxplot)],col="grey",main=paste(sep="",tf," in Th-0"),notch=T,
				ylab="binding strength [-log10(pval)]",outline=outliers,
                axes = FALSE, axisnames = FALSE,pch=20,cex=.3,xlim=c(0,length(ix.th17))  )
		text(1:length(mp$n), par("usr")[3], labels = names( ix.th0.plot.boxplot), srt = 45, adj = 1, xpd = TRUE,cex=.5)
		axis(2,las=1,cex.axis=0.7)	
	}
  }

  ## boxplots of clusters dist from tss within 10kb up/down of TSS
  ## th17
  ylim <- c(-10^4,10^4)
  mp <- boxplot(dtss.th17[ix.th17],col="grey",main="Th-17",ylab="distance of TSS [kb]",outline=outliers,ylim=ylim,
                axes = FALSE, axisnames = FALSE,pch=20,cex=.3,xlim=c(0,length(ix.th17))  )
  text(1:length(mp$n), par("usr")[3], labels = labels[ix.th17], srt = 45, adj = 1, xpd = TRUE,cex=.5)
  ix.ylab <- seq(from=ylim[1],to=ylim[2],by=1000)
  # barplot(c.table.th17[ix.th17],col="grey",main="Th-17",ylab="number of clusters",axes = FALSE, axisnames = FALSE,ylim=c(0,max(ix*1000)))
  # text(mp, par("usr")[3], labels = labels[ix.th17], srt = 45, adj = c(1,2), xpd = TRUE,cex=.5)
  # axis(2,las=1,at=ix*1000,labels=paste(sep="",ix,"k"))
  axis(2,las=1,at=ix.ylab,labels=ix.ylab/1000,cex.axis=0.7)
  ## th0
  ylim <- c(-10^4,10^4)
  mp = boxplot(dtss.th0[names(ix.th17[names( ix.th0.plot.boxplot)])],col="grey",main="Th-0",ylab="distance of TSS [bp]",axes = FALSE,outline=outliers,ylim=ylim,
    axisnames = FALSE,pch=20,cex=.3,xlim=c(0,length(ix.th17)))
  # mp = boxplot(dtss.th0[names(ix.th17[names( ix.th0.plot.boxplot)][1])],col="grey",main="Th-0",ylab="distance of TSS [bp]",axes = FALSE,outline=outliers,ylim=ylim,
  #    axisnames = FALSE,pch=20,cex=.3,xlim=c(0,length(ix.th17)))
  #  for(i in 2:length(ix.th0.plot.boxplot)){
  #    ix <- ix.th0.plot.boxplot[i]
  #    boxplot(dtss.th0[names(ix)],at=which(names(ix.th17[names( ix.th0.plot.boxplot)])==names(ix)),col="grey",axes = FALSE,outline=outliers,ylim=ylim,
  #            axisnames = FALSE,pch=20,cex=.3,xlim=c(1,15),add=T)
  #  }
  axis(2,las=1,at=ix.ylab,labels=ix.ylab/1000,cex.axis=0.7)
  text(1:length(ix.th17[names( ix.th0.plot.boxplot)]), par("usr")[3], labels = labels[ix.th17[names( ix.th0.plot.boxplot)]], srt = 45, adj = 1, xpd = TRUE,cex=.5)

  if(plot.histone.marks==TRUE){
    ## now plots for histone marks
	  marks <- names(h.pvals.th17[[1]])
	  for(m in 1:length(marks)){
	    ## ####plot pvals####
	    ## for th17
	    x <- sapply(h.pvals.th17[ix.th17],function(i) eval( parse(text=paste(sep="","i$",marks[m]))) )
	    mp <- boxplot(x,col="grey",main=paste(sep="","Th-17 (",marks[m],")"),ylab="binding strength [-log10(pval)]",outline=outliers,
	                axes = FALSE, axisnames = FALSE,pch=20,cex=.3,ylim=ylim.pval,xlim=c(0,length(ix.th17))  )
	    text(1:length(mp$n), par("usr")[3], labels = labels[ix.th17], srt = 45, adj = 1, xpd = TRUE,cex=.5)
	    axis(2,las=1,cex.axis=0.7)
	    ## for th0
	    x <- sapply(h.pvals.th0[names(ix.th0.plot.boxplot)],function(i) eval( parse(text=paste(sep="","i$",marks[m]))) )
	    mp = boxplot(x[[1]],col="grey",main=paste(sep="","Th-0 (",marks[m],")"),ylab="binding strength [-log10(pval)]",axes = FALSE,outline=outliers,
	      axisnames = FALSE,pch=20,cex=.3,ylim=ylim.pval,xlim=c(0,length(ix.th17)))
	    for(i in 2:length(ix.th0.plot.boxplot)){
	      ix <- ix.th0.plot.boxplot[i]
	      boxplot(x[[names(ix)]] ,at=which(names(ix.th17)==names(ix)),col="grey",axes = FALSE,outline=outliers,
	              axisnames = FALSE,pch=20,cex=.3,ylim=ylim.pval,xlim=c(1,15),add=T)
	    }
	    axis(2,las=1,cex.axis=0.7)
	    text(1:length(ix.th17), par("usr")[3], labels = labels[ix.th17], srt = 45, adj = 1, xpd = TRUE,cex=.5)
	    ## ####plot fc####
	    ## for th17
	    x <- sapply(h.fc.th17[ix.th17],function(i) eval( parse(text=paste(sep="","i$",marks[m]))) )
	    mp <- boxplot(x,col="grey",main=paste(sep="","Th-17 (",marks[m],")"),ylab="Fold Change [log2]",outline=outliers,
	                axes = FALSE, axisnames = FALSE,pch=20,cex=.3,xlim=c(0,length(ix.th17))  )
	    text(1:length(mp$n), par("usr")[3], labels = labels[ix.th17], srt = 45, adj = 1, xpd = TRUE,cex=.5)
	    axis(2,las=1,cex.axis=0.7)
	    ## for th0
	    x <- sapply(h.fc.th0[names(ix.th0.plot.boxplot)],function(i) eval( parse(text=paste(sep="","i$",marks[m]))) )
	    mp = boxplot(x[[1]],col="grey",main=paste(sep="","Th-0 (",marks[m],")"),ylab="Fold Change [log2]",axes = FALSE,outline=outliers,
	      axisnames = FALSE,pch=20,cex=.3,xlim=c(0,length(ix.th17)))
	    for(i in 2:length(ix.th0.plot.boxplot)){
	      ix <- ix.th0.plot.boxplot[i]
	      boxplot(x[[names(ix)]] ,at=which(names(ix.th17)==names(ix)),col="grey",axes = FALSE,outline=outliers,
	              axisnames = FALSE,pch=20,cex=.3,xlim=c(1,15),add=T)
	    }
	    axis(2,las=1,cex.axis=0.7)
	    text(1:length(ix.th17), par("usr")[3], labels = labels[ix.th17], srt = 45, adj = 1, xpd = TRUE,cex=.5)   
	  } 
	}
}

## here we create a list of TF enriched loci (clusters)
## input:
##     - ruler: a line per chromosome with the locations of all tf binding sites (sorted from start of chromosome to end)
##     - n.p: number of processors we will use (if multicore)
## output:
##     -ll: a list of clusters ll where each cluster (elem in ll holds):
##        - l: the current loci (elem in ll)
##        - s: summits vector as before
##        - tf.to.s: tfs vector matching the summits
##        - pval: pvals vector matching the summits
##        - tfs: a vector of tfs matching s, the summits in loci l
##        - spans: a vector of spans matching s, the summits in loci l, where spans is the dist between start and end of a peak
##        - trgts: genes that are targeted by the cluster (10kb upstream, in gene, 10kb downstream)
create.cluster.list <- function(ruler,n.p=1){
  tmp <- list()
  ll <- list()
  for(j in 1:length(ruler)){
    r <- names(ruler)[j]
    cat("working on ",r,"\n")
    x <- ruler[[j]] # for short typing let x stand for ruler[[chr[j]]]
    l.vec <- unique(x[["l"]]) # the clusters ids on j chr (only unique)
    n <- length(l.vec) # iterate over n clusters
    tmp[[j]] <- lapply(1:n,get.cluster.params,x=x,l.vec=l.vec,r=r)
  }
  ## concatenate tmp into one list
  command <- paste(sep="","c(",paste(sep="","tmp[[",1:length(tmp),"]]",collapse=","),")")
  ll=eval(parse(text=command))
  return(ll)
}
get.cluster.params <- function(i,x,l.vec,r){
  ix <- which(x[["l"]]==l.vec[i])
  l <- l.vec[i]
  start <- min(x[["start"]][ix])
  end <- max(x[["end"]][ix])
  s <- x[["s"]][ix]
  tf.to.s <- x[["tf.to.s"]][ix]
  pval <- x[["pval"]][ix]
  span.tfs <- x[["end"]][ix]-x[["start"]][ix]
  span.l <- end-start
  peak.ids <- x[["peak.ids"]][ix]
  expt <- x[["expt"]][ix]
  trgts <- unique(x[["trgts"]][ix])
  d.tss <- unique(x[["d.tss"]][ix])
  chr <- rep(r,length(ix))
  return(list(l=l,chr=chr,start=start,end=end,s=s,tf.to.s=tf.to.s,pval=pval,span.tfs=span.tfs,span.l=span.l,peak.ids=peak.ids,expt=expt,trgts=trgts,d.tss=d.tss))
}

## here we create a list of TF enriched loci (clusters)
## input:
##     - ruler: a line per chromosome with the locations of all tf binding sites (sorted from start of chromosome to end)
## output:
##     -ll: a list of clusters ll where each cluster (elem in ll holds):
##        - l: the current loci (elem in ll)
##        - s: summits vector as before
##        - pval: pvals vector matching the summits
##        - tfs: a vector of tfs matching s, the summits in loci l
##        - spans: a vector of spans matching s, the summits in loci l, where spans is the dist between start and end of a peak
##        - trgts: genes that are targeted by the cluster (10kb upstream, in gene, 10kb downstream)
create.cluster.list.no.pml <- function(ruler){
  tmp <- list()
  ll <- list()
  for(j in 1:length(ruler)){
    r <- names(ruler)[j]
    # cat("working on ",r,"\n")
    x <- ruler[[j]] # for short typing let x stand for ruler[[chr[j]]]
    l.vec <- unique(x[["l"]]) # the clusters ids on j chr (only unique)
    n <- length(l.vec) # iterate over n clusters
    tmp[[j]] <- lapply(1:n,get.cluster.params.no.pml,x=x,l.vec=l.vec,r=r)
  }
  ## concatenate tmp into one list
  command <- paste(sep="","c(",paste(sep="","tmp[[",1:length(tmp),"]]",collapse=","),")")
  ll=eval(parse(text=command))
  return(ll)
}
get.cluster.params.no.pml <- function(i,x,l.vec,r){
  ix <- which(x[["l"]]==l.vec[i])
  l <- l.vec[i]
  start <- min(x[["start"]][ix])
  end <- max(x[["end"]][ix])
  s <- x[["s"]][ix]
  # tf.to.s <- x[["tf.to.s"]][ix]
  pval <- x[["pval"]][ix]
  pval.mean <- mean(pval)
  span.tfs <- x[["end"]][ix]-x[["start"]][ix]
  span.l <- end-start
  peak.ids <- x[["peak.ids"]][ix]
  expt <- x[["expt"]][ix]
  expt.alphanum.sorted <- sort(x[["expt"]][ix])
  trgt.prox <- unique(x[["trgt.prox"]][ix])
  trgt.dist <- unique(x[["trgt.dist"]][ix])
  dtss.prox <- unique(x[["dtss.prox"]][ix])
  dtss.dist <- unique(x[["dtss.dist"]][ix])
  chr <- rep(r,length(ix))
  return(list(l=l,chr=chr,expt=expt,expt.alphanum.sorted=expt.alphanum.sorted,start=start,
			  end=end,s=s,pval=pval,pval.mean=pval.mean,span.tfs=span.tfs,
			  span.l=span.l,peak.ids=peak.ids,trgt.prox=trgt.prox,trgt.dist=trgt.dist,
			  dtss.prox=dtss.prox,dtss.dist=dtss.dist))
}
## same as above two functions but for permutatoin analysis so much more lean on data
create.cluster.list.perm <- function(r.p,n.p=1){
  tmp <- list()
  ll <- list()
  for(j in 1:length(r.p)){
    cat(".")
    l.vec <- unique(r.p[[j]][["l"]]) # the clusters ids on j chr (only unique)
    n <- length(l.vec) # iterate over n clusters
    tmp[[j]] <- lapply(1:n,get.cluster.params.perm,x=r.p[[j]],l.vec=l.vec)
  }
  cat("\n")
  ## concatenate tmp into one list
  command <- paste(sep="","c(",paste(sep="","tmp[[",1:length(tmp),"]]",collapse=","),")")
  ll=eval(parse(text=command))
  return(ll)
}
get.cluster.params.perm <- function(i,x,l.vec){
  ix <- which(x[["l"]]==l.vec[i])
  l <- l.vec[i]
  s <- x[["s"]][ix]
  tf.to.s <- x[["tf.to.s"]][ix]
  pval <- x[["pval"]][ix]
  return(list(l=l,s=s,tf.to.s=tf.to.s,pval=pval))
}
## add cluster memberships based on ruler
## require no more than d.cut distance between tf summits
## cur.l is the current loci number (or cluster number)
assign.clusters.ids.to.peaks <- function(ruler,d.cut){
 cur.l <- 0
 for(j in 1:length(ruler)){
   s <- ruler[[j]][["s"]]
   l <- numeric(length=length(s))
   if(length(l)>0){
	   cur.l <- cur.l+1
   }
   if(length(l)==1){
   	l[1] <- cur.l
   } else if(length(l)>1) {
	l[1] <- cur.l
	   for(i in 2:length(l)){
	     d <- s[i]-s[i-1] # assumes s is sorted increasingly
	     if(d>d.cut){
	       cur.l <- cur.l+1
	     }
	     l[i] <- cur.l    
	   }
	}
   ruler[[j]][["l"]] <- l
 }
 return(ruler)
}

# ## add cluster memberships based on ruler
# ## require no more than d.cut distance between tfs
# ## cur.l is the current loci number (or cluster number)
# assign.clusters.ids.to.peaks <- function(ruler,d.cut){
#  cur.l <- 0
#  for(j in 1:length(ruler)){
#    cur.l <- cur.l+1
#    s <- ruler[[j]][["s"]]
#    l <- numeric(length=length(s))
#    if(length(l)==1){
#    l[1] <- cur.l
#    } else if(length(l)>1) {
# 	   for(i in 2:length(l)){
# 	     d <- s[i]-s[i-1] # assume s is sorted increasingly
# 	     if(d>d.cut){
# 	       cur.l <- cur.l+1
# 	     }
# 	     l[i] <- cur.l    
# 	   }
# 	}
#    ruler[[j]][["l"]] <- l
#  }
#  return(ruler)
# }

## make.ruler makes a ruler: a line per chromosome with the locations of all tf binding sites (sorted from start of chromosome to end)
# also match these summit locations with corresponding:
# pvals, tfs, peak start and peak end trgts
# input:
#     - chr: a list of chromosome names
#     - macs.list.per.chrom: a list of macs peaks for each chromosome
# output:
#     - o: a list each chormosome ruler as an element
make.ruler.no.pml <- function(chr,macs.list.per.chrom){
  x <- macs.list.per.chrom
  o <- list()
  for(j in 1:length(chr)){
    r <- chr[j] # chrmosome we go over
    s <- numeric()
    pval <- numeric()
    tf.to.s <- character()
    trgt.prox <- character()
    trgt.dist <- character()
    dtss.prox <- character() 
    dtss.dist <- character()
    start <- numeric()
    end <- numeric()
    trtmnts <- names(x) # the treatments name (expt names)
    ## debug parameters ###
    ## which experiment peaks come from
    expt <- character()
    ## what was the peak id in that experiment
    peak.ids <- numeric()
    ## this will allow us to always back track from a cluster to the actual peaks in it
    ## debug params end ###
    for(i in 1:length(trtmnts)){
      o[[r]] <- list()
      e <- trtmnts[i] #experiment name
      tf <- strsplit(e,"_")[[1]][1]
      s <- c(s,x[[e]][[r]][,"start"]+x[[e]][[r]][,"summit"])
      pval <- c(pval,x[[e]][[r]][,7]) # the name of 7th column is X.10.log10.pvalue. this is pval
      # tf.to.s <- c(tf.to.s,rep(tf,length(x[[e]][[r]][,7]))) # all summits belong to tf
      start <- c(start,x[[e]][[r]][,"start"])
      end <- c(end,x[[e]][[r]][,"end"])
      expt <- c(expt,rep(e,length(x[[e]][[r]][,"end"])))
      peak.ids <- c(peak.ids,x[[e]][[r]][,"peak.id"])    
      trgt.prox <- c(trgt.prox,x[[e]][[r]][,"trgt.prox"])
	  trgt.dist <- c(trgt.dist,x[[e]][[r]][,"trgt.dist"])
      dtss.prox <- c(dtss.prox,x[[e]][[r]][,"dtss.prox"])
      dtss.dist <- c(dtss.dist,x[[e]][[r]][,"dtss.dist"])
    }
    ix <- sort(s,index.return=TRUE)$ix
    o[[r]] <- list(s=s[ix],pval=pval[ix],expt=expt[ix],start=start[ix],end=end[ix],peak.ids=peak.ids[ix],
				   trgt.prox=trgt.prox[ix],trgt.dist=trgt.dist[ix],dtss.prox=dtss.prox[ix],dtss.dist=dtss.dist[ix])
  }
  return(o)
}
## make.ruler makes a ruler: a line per chromosome with the locations of all tf binding sites (sorted from start of chromosome to end)
# also match these summit locations with corresponding:
# pvals, tfs, peak start and peak end trgts
# input:
#     - chr: a list of chromosome names
#     - macs.list.per.chrom: a list of macs peaks for each chromosome
# output:
#     - o: a list each chormosome ruler as an element
make.ruler <- function(chr,macs.list.per.chrom){
  x <- macs.list.per.chrom
  o <- list()
  for(j in 1:length(chr)){
    r <- chr[j] # chrmosome we go over
    s <- numeric()
    pval <- numeric()
    tf.to.s <- character()
    trgts <- character()
    d.tss <- character() 
    start <- numeric()
    end <- numeric()
    trtmnts <- names(x) # the treatments name (expt names)
    ## debug parameters ###
    ## which experiment peaks come from
    expt <- character()
    ## what was the peak id in that experiment
    peak.ids <- numeric()
    ## this will allow us to always back track from a cluster to the actual peaks in it
    ## debug params end ###
    for(i in 1:length(trtmnts)){
      o[[r]] <- list()
      e <- trtmnts[i] #experiment name
      tf <- strsplit(e,"_")[[1]][1]
      s <- c(s,x[[e]][[r]][,"start"]+x[[e]][[r]][,"summit"])
      pval <- c(pval,x[[e]][[r]][,7]) # the name of 7th column is X.10.log10.pvalue. this is pval
      tf.to.s <- c(tf.to.s,rep(tf,length(x[[e]][[r]][,7]))) # all summits belong to tf
      start <- c(start,x[[e]][[r]][,"start"])
      end <- c(end,x[[e]][[r]][,"end"])
      expt <- c(expt,rep(e,length(x[[e]][[r]][,"end"])))
      peak.ids <- c(peak.ids,x[[e]][[r]][,"peak.id"])    
      trgts <- c(trgts,x[[e]][[r]][,"trgts"])
      d.tss <- c(d.tss,x[[e]][[r]][,"d.tss"])
    }
    ix <- sort(s,index.return=TRUE)$ix
    o[[r]] <- list(s=s[ix],pval=pval[ix],tf.to.s=tf.to.s[ix],start=start[ix],end=end[ix],expt=expt[ix],peak.ids=peak.ids[ix],trgts=trgts[ix],d.tss=d.tss[ix])
  }
  return(o)
}
## split.macs.list.to.chrmosomes
# input:
#      - macs.list: a list of macs expts: here are a few lines of one expt
## chr	start	end	length	summit	tags	#NAME?	fold_enrichment	FDR(%)
## chr1	4322210	4323069	860	494	55	158.95	6.03	0.05
## chr1	4797749	4798368	620	211	29	119.82	3.47	0.09
## chr1	4848182	4849113	932	494	46	105.42	2.9	0.09
#      - expts: a list of the expts names from macs.list that you want to process
#      - chr: chrmosomes names
# output:
#      - x: a list with as many elements as chr specified by input.
#      -x[[i]]:macs list per chr, with peak.id column added (these are the row numbers of MACS)
split.macs.list.to.chrmosomes.no.pml <- function(macs.list,expts,chr="chr1"){
  x <- list()
  n <- length(expts)
  for(i in 1:n){
    e <- expts[i] #experiment name
    cat("wroking on expt", e,"\n")
    x[[e]] <- lapply(chr,split.one.macs.expt.by.chromosome.no.pml,m=macs.list[[e]])
    names(x[[e]]) <- chr
  }
  return(x)
}
# a helper function for spliat.macs.list.to.chrmosomes, gives for one chromosome the macs rows for expt MACS matrix m
# input:
#     - r is chromosome
#     - m is macs matrix for expt e from above function
split.one.macs.expt.by.chromosome.no.pml <- function(r,m){
  ix.chr.i <- which(m[,"chr"]==r)
  # cat("working on",r,"\n")
  o <- list()
  o[[r]] <- m[ix.chr.i,]
  o[[r]]$peak.id <- ix.chr.i
  return(o[[r]])
}

## split.macs.list.to.chrmosomes
# input:
#      - macs.list: a list of macs expts: here are a few lines of one expt
## chr	start	end	length	summit	tags	#NAME?	fold_enrichment	FDR(%)
## chr1	4322210	4323069	860	494	55	158.95	6.03	0.05
## chr1	4797749	4798368	620	211	29	119.82	3.47	0.09
## chr1	4848182	4849113	932	494	46	105.42	2.9	0.09
#      - expts: a list of the expts names from macs.list that you want to process
#      - chr: chrmosomes names
#      - pml: a list of peaks mapped to their target genes (a peak can map to more than one gene)
##         peak_id gn_length       TSS       TES    Summit d_TSS d_TES   pval   FC  FDR
## BCL2A1D     777     8.563  88626792  88635355  88626740   -52 -8615 123.67 9.83  2.31
## GPR171      535     5.374  58905851  58911225  58905797   -54 -5428 118.65 9.50  2.50
## RNU11       594     0.108 131825989 131826097 131826045    56   -52 132.63 7.52  2.41
#      - n.p: number of processors you want to use if more than 1 assumes multicore library is loaded
#	   - mtbs.tss.dist: a minimum distance from TSS for which we include genes as targets of MTBS
#      - get.targets: a logical indicating if to get targets of each peak into output (time consuming)
# output:
#      - x: a list with as many elements as chr specified by input.
#      -x[[i]]:macs list per chr, with peak.id column added (these are the row numbers of MACS)
split.macs.list.to.chrmosomes <- function(macs.list,expts,chr="chr1",pml,n.p,mtbs.tss.dist,match.targets=FALSE){
  x <- list()
  n <- length(expts)
  for(i in 1:n){
    e <- expts[i] #experiment name
    # cat("wroking on expt", e,"\n")
    x[[e]] <- lapply(chr,split.one.macs.expt.by.chromosome,m=macs.list[[e]],p=pml[[e]],n.p=n.p,
	dtss.cut=mtbs.tss.dist,match.targets=match.targets)
    names(x[[e]]) <- chr
  }
  return(x)
}
# a helper function for spliat.macs.list.to.chrmosomes, gives for one chromosome the macs rows for expt MACS matrix m
# input:
#     - r is chromosome
#     - m is macs matrix for expt e from above function
#     - p is pml matrix for expt e from above function
#     - n.p is the number of processors we want to use from above function
split.one.macs.expt.by.chromosome <- function(r,m,p,n.p,dtss.cut,match.targets=FALSE){
  ix.chr.i <- which(m[,"chr"]==r)
  cat("working on",r,"\n")
  o <- list()
  o[[r]] <- m[ix.chr.i,]
  o[[r]]$peak.id <- ix.chr.i
  if(match.targets==TRUE){
    trgts <- character(length=length(ix.chr.i))
    d.tss <- character(length=length(ix.chr.i))
    if(length(ix.chr.i)>0){
      trgts <- unlist(mclapply(ix.chr.i,get.targets,p=p,dtss.cut=dtss.cut,mc.cores=n.p ))
      d.tss <- unlist(mclapply(ix.chr.i,get.dtss,p=p,dtss.cut=dtss.cut,mc.cores=n.p ))
    }
    o[[r]]$trgts <- trgts
    o[[r]]$d.tss <- d.tss
  } else {
	o[[r]]$trgts <- character(length=length(ix.chr.i))
    o[[r]]$d.tss <- character(length=length(ix.chr.i))
  }
  return(o[[r]])
}
get.dtss <- function(ix,p,dtss.cut){
  ix.trgts <- which(p[,"peak_id"]==ix)
  ix.trgts.good <- which(abs(p[ix.trgts,"d_TSS"])<dtss.cut)
  tmp <- paste(sep="",p[ix.trgts[ix.trgts.good],"d_TSS"],collapse="_")
  return(tmp)
}
get.targets <- function(ix,p,dtss.cut){
  ix.trgts <- which(p[,"peak_id"]==ix)
  ix.trgts.good <- which(abs(p[ix.trgts,"d_TSS"])<dtss.cut)
  tmp <- paste(sep="",unique(rownames(p)[ix.trgts[ix.trgts.good]]),collapse="_")
  return(tmp)
}
get.peaks.id <- function(ix,p){
  ix.trgts <- which(p[,"peak_id"]==ix)
  tmp <- paste(sep="",unique(rownames(p)[ix.trgts]),collapse="_")
  return(tmp)
}
# get.dtss <- function(ix,p){
#   ix.trgts <- which(p[,"peak_id"]==ix)
#   if(length(ix.trgts)>0){
#     dtss <- paste(sep="",p[ix.trgts,"d_TSS"],collapse="_")
#   } else {
#     dtss <- ""
#   }
#   return(dtss)
# }


