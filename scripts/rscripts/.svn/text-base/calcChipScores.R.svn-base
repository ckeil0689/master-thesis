##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## January 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

source("r_scripts/th17/used_for_paper/utilChipScores.R")
path.input <- "input/th17/used_for_paper/"
if(use.multicore==TRUE){
  library(multicore)
}

# set paths
## path.output <- "input/th17/used_for_paper/"
path.results <- "results/"
path.macs <- paste(sep="",path.input,"rawData/MACS_tab_delim_",date.macs.data.compiled,"/")
path.gene.cntric <- paste(sep="",path.input,"rawData/MACS_to_genes_",date.macs.data.compiled,"/")
path.genes.to.macs <- paste(sep="",path.input,"rawData/genes_to_MACS_",date.macs.data.compiled,"/")
path.input.post.macs <- paste(sep="","input/th17/used_for_paper/rawData/MACS_tab_delim_",date.macs.data.compiled,"/") # the files here adhere to tab delim format


# params
debug <- FALSE
getSequence <- FALSE # for motif finding this will add sequence information for clustered binding sites. (takes times)
plotGenesWithTFHits <- TRUE
# this is the amount we need to add up/down stream of each genes length to get the right kb region where we searched for peaks
upstream.gene.region <- peak.dist.upstream/1000
downstream.gene.region <- peak.dist.downstream/1000
######### set up ####################
all.files.gene.centric <- list.files(path.gene.cntric)
gene.centric.files <- all.files.gene.centric[grep("anno",all.files.gene.centric)]

# rearrange anno files to : (1) p300_th0 (2) irf4_th0 (3) batf_th0
#                           (4) maf_th0 (5) stat3_th0 (6) rorc_th0 (7) fosl2_th0
#                           (8) p300_th17 (9) irf4_th17 (10) batf_th17
#                           (11) maf_th17 (12) stat3_th17 (13) rorc_th17
#                           (14) fosl2_th17 
ix.p300.th0 <- grep("p300_th0",gene.centric.files,ignore.case = TRUE)
ix.irf4.th0 <- grep("irf4_th0",gene.centric.files,ignore.case = TRUE)
ix.batf.th0 <- grep("batf_th0",gene.centric.files,ignore.case = TRUE)
ix.maf.th0 <- grep("maf_th0",gene.centric.files,ignore.case = TRUE)
ix.stat3.th0 <- grep("stat3_th0",gene.centric.files,ignore.case = TRUE)
ix.rorc.th0 <- grep("rorc_th0",gene.centric.files,ignore.case = TRUE)
ix.fosl2.th0 <- grep("fosl2_th0",gene.centric.files,ignore.case = TRUE)

ix.p300.th17 <- c(grep("p300_th17_",gene.centric.files,ignore.case = TRUE),grep("p300_th17rorcwt_",gene.centric.files,ignore.case = TRUE))
ix.irf4.th17 <- c(grep("irf4_th17_",gene.centric.files,ignore.case = TRUE),grep("irf4_th17rorcwt_",gene.centric.files,ignore.case = TRUE))
ix.batf.th17 <- grep("batf_th17_",gene.centric.files,ignore.case = TRUE)
ix.maf.th17 <- grep("maf_th17_",gene.centric.files,ignore.case = TRUE)
ix.stat3.th17 <- c(grep("stat3_th17_",gene.centric.files,ignore.case = TRUE),grep("stat3_Th17rorcwt_",gene.centric.files,ignore.case = TRUE))
ix.rorc.th17 <- grep("rorc_th17_",gene.centric.files,ignore.case = TRUE)
ix.fosl2.th17 <- grep("fosl2_th17",gene.centric.files,ignore.case = TRUE)

gene.centric.files <- gene.centric.files[c(ix.p300.th0, ix.irf4.th0, ix.batf.th0,
                                           ix.maf.th0, ix.stat3.th0, ix.rorc.th0,ix.fosl2.th0,
                                           ix.p300.th17, ix.irf4.th17, ix.batf.th17,
                                           ix.maf.th17, ix.stat3.th17, ix.rorc.th17,ix.fosl2.th17)]

# keep all genes with a hit or more in here
gns.with.peaks <- character()
peaks.mat.List <- list()
######### read peaks and plot peaks for each tf on normalized gene space ####################
pdf(file="./results/chipSeqBindingSitesHist.pdf")
for(i in 1:length(gene.centric.files)) {
  # get chipseq data
  x <- read.table(paste(path.gene.cntric,gene.centric.files[i],sep=""))
  median.gn.length <-   median(x$V3)/10000 # for normalizing (in 10kb)
  tf.nm <- strsplit(gene.centric.files[i],"_")[[1]][1]
  expt.nm <- strsplit(gene.centric.files[i],"_")[[1]][2]
  expt.id <- paste(sep="_",strsplit(gene.centric.files[i],"_")[[1]][3],strsplit(gene.centric.files[i],"_")[[1]][4])
  f.nm <- paste(tf.nm,expt.nm,expt.id,sep="_") # label for hist
  cat(sep="","median gene length (in 10kb units) for tf-", tf.nm, " experiment-", expt.nm, " id-",expt.id , " is: ", median.gn.length, "\n")
  gsub("\"","",tf.nm)
  gn.nms <- toupper(as.character(x$V1))
  gns.with.peaks <- unique(c(gns.with.peaks,gn.nms))
  peaks <- as.character(x$V5) # get all peaks
  peaks.mat <- deconvolve.peaks(peaks,gn.nms)
  plot.peaks.histogram(peaks.mat,median.gn.length)
  # sort peaks by absolute distance to transcription start site
  ix <- sort(abs(peaks.mat[,"d_TSS"]),index.return=TRUE)$ix
  peaks.mat <- peaks.mat[ix,]
  peaks.mat.List[[f.nm]] <- peaks.mat
}
save(peaks.mat.List,file=paste(sep="",path.input,"peaks_mat_List.RData"))
dev.off()
# keep for each experiment duplicate the one experiment that is higher quality based on histogram and number of peaks
keepers <- c("P300_Th0_SL1948_SL1947","P300_Th17rorcwt_SL3594_SL3592",
             "IRF4_Th0_SL1235_SL1234","IRF4_Th17rorcwt_SL2872_SL2876",
             "BATF_Th0_SL3192_SL3190","BATF_Th17_SL3037_SL3036",
             "MAF_Th0_SL4424_SL4425","MAF_Th17_SL3032_SL2871",
             "STAT3_Th0_SL3780_SL3778","STAT3_Th17rorcwt_SL3315_SL3319",
             "RORC_Th0_SL3779_SL3778","RORC_Th17_SL2870_SL2871",
             "FOSL2_Th0_SL6500_SL6499","FOSL2_Th17_SL6498_SL6497")
# shortening name/keeping only 'best' experiments
pml <- peaks.mat.List[keepers]
for (i in 1:length(pml)){
  pml[[i]][,c("d_TES","d_TSS")] <- abs(pml[[i]][,c("d_TES","d_TSS")]/1000)
}

######### get number of peaks for each experiment ####################
cat("calculating the number of peaks per experiment\n")
# get peak distribution per experiment
pd <- list()
for (i in 1:length(pml)){
  x <- strsplit(names(pml)[i],"_")[[1]]
  e <- paste(sep="",x[[3]],"_",x[[4]])
  f.nm <- paste(sep="",e,"_peaks.xls")
  x <- read.table(paste(path.macs,f.nm,sep=""),header=T)
  # get pval distribution for experiment e
  pd[[i]] <- as.numeric(x[,7])
}
names(pd) <- names(pml)
rm(x)
######### create score for each ''''GENE'''' for each experiment ####################
cat("creating pval.peaks.num (from poisson model)\n") 
c.names.m <- c(	"gn_length", "kb_searched","num_peaks",
				"pval_macs_mean","pval_poisson")
# store results in gene mat list (gml)
gml <- list()
stop("AM")
if(use.multicore==TRUE){
  ## gml <- mclapply(mc.cores = 7,pml,create.gene.based.chip.scores,gns=gns.with.peaks)
  gml <- mclapply(mc.cores = n.pr,1:length(pml),create.gene.based.chip.scores,
                  pml=pml,
                  gns=gns.with.peaks,
                  peaks.dist.per.expt=pd,
                  mm9.effective.genome.size=mm9.effective.genome.size,
                  upstream.gene.region=upstream.gene.region,
                  downstream.gene.region=downstream.gene.region,
                  c.names.m=c.names.m)
  names(gml) <- names(pml)
} else {
  for (i in 1:length(pml)) {
    # for simplicity
    w <- pml[[i]]
    # a place holder matrix in which i keep results for one experiment until moved to list of matrices
    m <- matrix(0,nr=length(gns.with.peaks),nc=length(c.names.m))
    colnames(m) <- c.names.m
    rownames(m) <- gns.with.peaks

    # get gene names for experiment i
    gn.nms.i <- unique(rownames(w))
    for (j in 1:length(gn.nms.i)){
      ix <- which(rownames(w) == gn.nms.i[j])
      m[gn.nms.i[j],"gn_length"] <- w[ix[1],"gn_length"]
      ## here i add 20kb to represent the effective region of a gene where peaks where searched
      m[gn.nms.i[j],"kb_searched"] <- m[gn.nms.i[j],"gn_length"] + upstream.gene.region + downstream.gene.region
      ## lambda.per.kbp.searched <- m[gn.nms.i[j],"kb_searched"]*lambda.per.bp*1000
      m[gn.nms.i[j],"num_peaks"] <- length(ix)
      m[gn.nms.i[j],"pval_macs_mean"] <- mean(w[ix,"pval"])
      # m[gn.nms.i[j],"pval_macs_mean"] <- m[gn.nms.i[j],"pval_macs_mean"]/10*log(10) #ie. -10log10(pval)/10*ln(10)
      ## get lambda for poisson distribution (lambda = #peaks with at least as high mean pval over genome)/effective genome size* searched gene region (in bp)
      num.peaks.genome.wide <- length(which(pd[[i]]>=m[gn.nms.i[j],"pval_macs_mean"]))
      lambda <- num.peaks.genome.wide/mm9.effective.genome.size*m[gn.nms.i[j],"kb_searched"]*1000

      ## AM not sure why I have -1 here, what was i thinking.. need to get back to this!
      # m[gn.nms.i[j],"pval_poisson"] <- -log10(ppois(m[gn.nms.i[j],"num_peaks"]-1,lambda,lower.tail=FALSE))
	  m[gn.nms.i[j],"pval_poisson"] <- -log10(ppois(m[gn.nms.i[j],"num_peaks"],lambda,lower.tail=FALSE))
    }
    gml[[ names(pml)[i] ]] <- m
  }
}

save(gml,file=paste(sep="",path.input,"chipScoresPerGenePerTfList.RData"))

system(paste(sep="","rm -rf ",path.genes.to.macs))
system(paste(sep="","mkdir ",path.genes.to.macs))
# write the gmls to files
for(i in 1:length(gml)){
  f.nm <- names(gml)[i]
  write.table(gml[[i]],file=paste(path.genes.to.macs,f.nm,".xls",sep=""),sep="\t")
}

##############################################################
# create S_gene_th17_minus_th0 for each ''''GENE'''' for each TF
cat("create S_gene_th17_minus_th0_plusp300 for each ''''GENE'''' for each TF\n")

core.tfs <- c("P300",core.tfs)
# this are the columns to be substracted between th17 and th10
#substract.columns=c("num_peaks","S_pval_mean","S_fc_mean","S_dist_mean","S_peak_mean","pval_poisson")
substract.columns=c("num_peaks","pval_poisson")
# first keep pval_poisson score for p300; we will use it with each tf
p300.score <- (gml[["P300_Th17rorcwt_SL3594_SL3592"]][,"pval_poisson"]-gml[["P300_Th0_SL1948_SL1947"]][,"pval_poisson"])
# for BATF MAF IRF4 STAT3 RORC calculate Th17 scores
th17.list <- list()
for(i in 1:length(core.tfs)){
  tf <- core.tfs[i]
  ix.th17 <- grep(paste(sep="",tf,"_","Th17"),names(gml),ignore.case=TRUE)
  th17.list[[tf]] <- (gml[[ix.th17]][,"pval_poisson"])
}
# now calculate th0 scores
th0.list <- list()
for(i in 1:length(core.tfs)){
  tf <- core.tfs[i]
  ix.th0 <- grep(paste(sep="",tf,"_","Th0"),names(gml),ignore.case=TRUE)
  th0.list[[tf]] <- (gml[[ix.th0]][,"pval_poisson"])
}
# now calculate Th17 minus th0
th17.minus.th0.list <- list()
for(i in 1:length(core.tfs)){
  tf <- core.tfs[i]
  th17.minus.th0.list[[tf]] <- (th17.list[[tf]]-th0.list[[tf]])
}
# now add the p300 score to each TF score (i.e. add the lineage specific active transcription mark)
th17.minus.th0.plus.p300.list <- list()
for(i in 1:length(core.tfs)){
  tf <- core.tfs[i]
  th17.minus.th0.plus.p300.list[[tf]] <- (th17.minus.th0.list[[tf]]+p300.score)
}
# write output chipseq table
chip.scores <- matrix(0,
                      nr=length(th17.list[[1]]),
                      nc=length(th17.list)*4)
expts <- c("th17","th0","th17_minus_th0","th17_minus_th0_plus_p300")
cnt <- 0
c.names <- character()
for(e in 1:length(expts)){
    if(expts[e]=="th17"){
      for (i in 1:length(core.tfs)){
        tf <- core.tfs[i]
        c.names <- c(c.names,paste(sep="",tf,"_",expts[e]))              
        cnt <- cnt+1
        chip.scores[,cnt] <- th17.list[[tf]]
      }
    } else if(expts[e]=="th0") {
      for (i in 1:length(core.tfs)){
        tf <- core.tfs[i]
        c.names <- c(c.names,paste(sep="",tf,"_",expts[e]))                     
        cnt <- cnt+1        
        chip.scores[,cnt] <- th0.list[[tf]]
      }
    }else if(expts[e]=="th17_minus_th0") {
      for (i in 1:length(core.tfs)){
        tf <- core.tfs[i]
        c.names <- c(c.names,paste(sep="",tf,"_",expts[e]))                     
        cnt <- cnt+1        
        chip.scores[,cnt] <- th17.minus.th0.list[[tf]]
      }
    } else if(expts[e]=="th17_minus_th0_plus_p300"){
      for (i in 1:length(core.tfs)){
        tf <- core.tfs[i]
        c.names <- c(c.names,paste(sep="",tf,"_",expts[e]))              
        cnt <- cnt+1            
        chip.scores[,cnt] <- th17.minus.th0.plus.p300.list[[tf]]
      }
  }
}
colnames(chip.scores) <- c.names
rownames(chip.scores) <- names(th17.list[[1]])
cat("writing chipScores matrix!\n")
write.table(chip.scores,file=paste(path.input,"chip_scores_",date.macs.data.compiled,".xls",sep=""),sep="\t")
########## QC section for Chipscores ############

keepers <- c("P300_Th17rorcwt_SL3594_SL3592",
             "IRF4_Th17rorcwt_SL2872_SL2876",
             "BATF_Th17_SL3037_SL3036",
             "MAF_Th17_SL3032_SL2871",
             "STAT3_Th17rorcwt_SL3315_SL3319",
             "RORC_Th17_SL2870_SL2871")
qc.genes <- unique(c(core.tfs,known.genes))
qc.genes[1]="EP300"
pdf(file=paste(sep="",path.results, "qc_figs_chip_scores.pdf"))
boxplot(th17.list,col="green",main="th17 chip score")
boxplot(th17.minus.th0.list,col="green",main="th17_minus_th0 chip score")
boxplot(th17.minus.th0.plus.p300.list,col="green",main="th17_minus_th0_plus_p300(th17-th0)\n chip score")
for(i in 1:length(core.tfs)){
  tf <- core.tfs[i]
  x <- sort(th17.list[[tf]],decreasing=TRUE)
  th17 <- which(names(x) %in% qc.genes)
  names(th17) <- names(x)[which(names(x) %in% qc.genes)]
  x <- sort(th17.minus.th0.list[[tf]],decreasing=TRUE)
  th17.minus.th0 <- which(names(x) %in% qc.genes)
  names(th17.minus.th0) <- names(x)[which(names(x) %in% qc.genes)]
  x <- sort(th17.minus.th0.plus.p300.list[[tf]],decreasing=TRUE)
  th17.minus.th0.plus.p300 <- which(names(x) %in% qc.genes)
  names(th17.minus.th0.plus.p300) <- names(x)[which(names(x) %in% qc.genes)]
  x <- rep(1,length(qc.genes))
  y <- th17
  txt.pos <- 0.2
  colors <- topo.colors(length(qc.genes))
  names(colors) <- qc.genes
  cex.pts <- 0.8
  cex.txt <- 0.3
  plot(x,y,ylim=c(0,2000),xlim=c(0.5,5.5),cex=cex.pts,col=colors[names(y)],pch=20,xaxt = "n",
       xlab="ChIP analysis",ylab="rank",main=paste(sep="","QC for ",tf,"-ChIP analysis"))
  axis(1, 1:5, c("th17","th17-th0","th17-th0-p300",
                 "th17.pval","th17.fc"),cex.axis=.6)
  text(x+txt.pos,y,labels=names(y),cex=cex.txt)
  x=rep(2,length(qc.genes))
  y= th17.minus.th0
  points(x,y,cex=cex.pts,col=colors[names(y)],pch=20)
  text(x+txt.pos,y,labels=names(y),cex=cex.txt)
  x=rep(3,length(qc.genes))
  y= th17.minus.th0.plus.p300
  points(x,y,cex=cex.pts,col=colors[names(y)],pch=20)
  text(x+txt.pos,y,labels=names(y),cex=cex.txt)

  expt <- keepers[i]
  x <- sort(pml[[expt]][,"pval"],decreasing=TRUE)
  tmp <- which(names(x) %in% qc.genes)
  names(tmp) <- names(x)[which(names(x) %in% qc.genes)]
  gns <- unique(names(tmp))
  pval <- numeric()
  for(j in 1:length(gns)){pval <- c(pval,tmp[gns[j]][1])}

  x <- sort(pml[[expt]][,"FC"],decreasing=TRUE)
  tmp <- which(names(x) %in% qc.genes)
  names(tmp) <- names(x)[which(names(x) %in% qc.genes)]
  gns <- unique(names(tmp))
  fc <- numeric()
  for(j in 1:length(gns)){fc <- c(fc,tmp[gns[j]][1])}

  num.peaks <- length(pval)
  x <- rep(4,num.peaks)
  y <- pval
  colors <- topo.colors(length(qc.genes))
  names(colors) <- qc.genes
  cex.pts <- 0.8
  points(x,y,cex=cex.pts,col=colors[names(y)],pch=20)
  text(x+txt.pos,y,labels=names(y),cex=cex.txt)

  x=rep(5,num.peaks)
  y= fc
  points(x,y,cex=cex.pts,col=colors[names(y)],pch=20)
  text(x+txt.pos,y,labels=names(y),cex=cex.txt)
 
}
dev.off()
