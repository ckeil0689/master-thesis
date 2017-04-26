##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## January 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# set paths
path.input.refseq <- "input/th17/used_for_paper/rawData/"
path.input.pre.macs <- paste(sep="","input/th17/used_for_paper/rawData/MACS_",date.macs.data.compiled,"/") # the files here not adhere to tab delim format
path.input.post.macs <- paste(sep="","input/th17/used_for_paper/rawData/MACS_tab_delim_",date.macs.data.compiled,"/") # the files here adhere to tab delim format
path.output <- paste(sep="","input/th17/used_for_paper/rawData/MACS_to_genes_",date.macs.data.compiled,"/")

# read refseq gene annotation table
## refseq <- read.delim(paste(sep="",path.input.refseq,"20101005_UCSC_mm9_kgXref.csv"))
refseq <- read.delim(paste(sep="",path.input.refseq,"UCSC_mm9_refseq_genes_",date.macs.data.compiled,".txt"))
cat("reading refseq table\n")
# refseq can have many transcripts for each gene
# here i make it have only one transcript for each gene (the longest one)
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
cat("making a unique (per gene name) refseq table (by taking the longest transcript of each gene)\n")
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
ix <- which(r.u$strand=="-")
tmp.tes <- r.u$txStart[ix]
tmp.tss <- r.u$txEnd[ix]
r.u[ix,c("txStart","txEnd")] <- cbind(tmp.tss,tmp.tes)

# read macs_id_to_expt_map.txt (so we can give meaningful names to each chipseq expt)
macs.id.to.expt.name <- read.delim(paste(sep="",path.input.refseq,"macs_id_to_expt_map.txt"))

# set file names
MACS.files <- list.files(path.input.pre.macs)
MACS.files.pos.peaks.pre <- paste(sep="",path.input.pre.macs,MACS.files[-grep("negative",MACS.files)])
MACS.files.pos.peaks.post <- paste(sep="",path.input.post.macs,MACS.files[-grep("negative",MACS.files)])

# make files tab delim proper!
system(paste(sep="","rm -rf ",path.input.post.macs))
system(paste(sep="","mkdir ",path.input.post.macs))
for(i in 1:length(MACS.files.pos.peaks.pre)){
  system(paste(sep=" ","tail -n+24",MACS.files.pos.peaks.pre[i],">",MACS.files.pos.peaks.post[i] )) # line 24 is header/ before is MACS params
}

# read MACS output files
macs.list <- list()
n <- length(MACS.files.pos.peaks.post)
# AM some experiments are trail macs runs! remove them Apr 14th
mistake.id.full <- c("SL2874_SL2875","SL2874_SL3315","SL2875_SL3316","SL3315_SL3316")
# AM some experiments did not clear quality control! remove them Sep 15th
mistake.id.full <- c(mistake.id.full,"SL1172_SL1171")
# AM some experiments are not for th17 project! remove them Sep 15th
mistake.id.full <- c(mistake.id.full,"SL6140_SL6142","SL6141_SL6142")

id.full <- character(length=n)
for(i in 1:n){
  id.full[i] <- gsub(".*(SL\\d+_SL\\d+).*","\\1",MACS.files.pos.peaks.post[i],perl=T) # e.g. "SL1040_SL972" (expt.id_ctrl.id)
  if(!any(mistake.id.full == id.full[i])){
    expt.nm <- as.character(macs.id.to.expt.name$expt[which(macs.id.to.expt.name$id == id.full[i])]) # map id to expt name
    macs.list[[ expt.nm ]] <- read.delim(file=MACS.files.pos.peaks.post[i])
  } else {
    cat("found a wrong MACS experiment",id.full[i],"\n")
  }
}

# make folder for gene centered macs files
system(paste(sep="","rm -rf ",path.output))
system(paste(sep="","mkdir ",path.output))
n <- length(macs.list)
expt.nms <- names(macs.list)
## gene.cntrd.pks <- list()
# for each expt (i goes over chip expts)
cat("writing gene centered macs files.  Current implementation is naive -> This may take a while...")
for(i in 1:n){
  # get the expt name
  nm <- expt.nms[i]
  # file name
  f.nm <- paste(sep="",path.output,nm,"_anno_",date.is,".txt")
  M <- macs.list[[i]]
  # get the number of genes
#  m <- dim(macs.list[[i]])[1]
  m <- length(gn.nms)
  # for each genes in the expt (j goes over genes)
  cat(file=f.nm,sep="\t","# Gene","n_peak","len","strand","(peak_id,summit,d_TSS,d_TES,class,pval,fold_enrich,FDR),(..)\n")
  peaks.summit <- (M[,"start"]+M[,"summit"])
  cat(sep="","working on ",nm,":\n")
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
    ix <- sort(union( c(which( abs(d.tss.all) < peak.dist.downstream ),which( abs(d.tes.all) < peak.dist.upstream )),
                 which( sign(d.tss.all) != sign(d.tes.all) )),decreasing=TRUE)
    ix <- ix[which(M[ix,"chr"]==chr)]    
    n.peaks <- length(ix)
    if(n.peaks > 0){
      # pring peaks for gene j    
      core.line <- paste(sep="\t",gn.nms[ j ],n.peaks,gn.lngth,strand)    
      # for each peak (l goes over peaks)
      peaks.line <- paste(sep="","(")
      for (l in 1:n.peaks){
        if(l>1){
          peaks.line <- paste(sep="",peaks.line,";(")
        }
        d.tss <- d.tss.all[ ix[l] ]
        d.tes <- d.tes.all[ ix[l] ]
        if(sign(d.tss)!=sign(d.tes)){
          class="intra"
        } else if(d.tss<=0){
          class="upstream"
        } else if(d.tss>0){
          class="downstream"
        }
        # (peak_id,summit,d_TSS,d_TES,class,pval,fold_enrich,FDR)
        peaks.line <- paste(sep="",peaks.line,paste(sep=",", ix[l], peaks.summit[ ix[l] ], d.tss, d.tes, class, M[ ix[l] ,7], M[ ix[l] ,8], M[ ix[l] ,9] ),")")
      }
      cat(file=f.nm,append=TRUE,sep="",core.line,"\t",peaks.line,"\n")
    }
  }
}

















