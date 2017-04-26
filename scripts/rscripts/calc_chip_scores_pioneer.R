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
path.genes.to.macs <- paste(sep="",path.input,"rawData/genes_to_MACS_pioneer_",date.macs.data.compiled,"/")
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
ix.irf4 <- grep("irf4",gene.centric.files,ignore.case = TRUE)

gene.centric.files <- gene.centric.files[c(ix.irf4)]

# keep all genes with a hit or more in here
gns.with.peaks <- character()
peaks.mat.List <- list()
######### read peaks and plot peaks for each tf on normalized gene space ####################
pdf(file="./results/chipSeqBindingSitesHistIrf4AcrossLin.pdf")
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
save(peaks.mat.List,file=paste(sep="",path.input,"peaks_mat_List_pioneer.RData"))
dev.off()
# keep for each experiment duplicate the one experiment that is higher quality based on histogram and number of peaks
keepers <- c("IRF4_Th0_SL1235_SL1234","IRF4_Th1_SL3501_SL3500","IRF4_Th2_SL3194_SL3193",
			"IRF4_Th17rorcwt_SL2872_SL2876","IRF4_iTreg_SL3196_SL3195")
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
      m[gn.nms.i[j],"pval_macs_mean"] <- m[gn.nms.i[j],"pval_macs_mean"]/10*log(10) #ie. -10log10(pval)/10*ln(10)
      ## get lambda for poisson distribution (lambda = #peaks with at least as high mean pval over genome)/effective genome size* searched gene region (in bp)
      num.peaks.genome.wide <- length(which(peaks.dist.per.expt[[i]]>=m[gn.nms.i[j],"pval_macs_mean"]))
      lambda <- num.peaks.genome.wide/mm9.effective.genome.size*m[gn.nms.i[j],"kb_searched"]*1000

      ## AM not sure why I have -1 here, what was i thinking.. need to get back to this!
      m[gn.nms.i[j],"pval_poisson"] <- -log10(ppois(m[gn.nms.i[j],"num_peaks"]-1,lambda,lower.tail=FALSE))
    }
    gml[[ names(pml)[i] ]] <- m
  }
}
names(gml)[grep("IRF4_Th17rorcwt_SL2872_SL2876",names(gml))] <- "IRF4_Th17_SL2872_SL2876"
save(gml,file=paste(sep="",path.input,"chipScoresPerGenePerTfListPioneer.RData"))

system(paste(sep="","rm -rf ",path.genes.to.macs))
system(paste(sep="","mkdir ",path.genes.to.macs))
# write the gmls to files
for(i in 1:length(gml)){
  f.nm <- names(gml)[i]
  write.table(gml[[i]],file=paste(path.genes.to.macs,f.nm,".xls",sep=""),sep="\t")
}

##############################################################
# create S_gene_th17_minus_th0 for each ''''GENE'''' for each TF
cat("create S_gene_th17_minus_th0 for each ''''GENE'''' for each TF\n")

# this are the columns to be substracted between th17 and th10
# first calculate IRF4 Th17,Th1,Th2,iTreg scores
irf4.th.list <- list()
lineages <- c("Th0","Th1","Th2","Th17","iTreg")
for(i in 1:length(lineages)){
	ix <- grep(paste(sep="",lineages[i],"_"),names(gml),ignore.case=T)
	irf4.th.list[[lineages[i]]] <- (gml[[ix]][,"pval_poisson"])
}

# now calculate IRF4 Th17,Th1,Th2,iTreg scores - Th0
irf4.th.minus.th0.list <- list()
for(i in 1:length(lineages)){
	irf4.th.minus.th0.list[[lineages[i]]] <- (irf4.th.list[[lineages[i]]]) - (gml[["IRF4_Th0_SL1235_SL1234"]][,"pval_poisson"])
}

# write output chipseq table
n <- length(irf4.th.list)
chip.scores <- matrix(0,
                      nr=length(irf4.th.list[[1]]),
                      nc=n*2)
chip.scores <- cbind(sapply(irf4.th.list,function(i) i),sapply(irf4.th.minus.th0.list,function(i) i))
ix <- (n+1):(2*n) # change col names for lineages - th0
colnames(chip.scores)[ix] <- paste(sep="",colnames(chip.scores)[ix],"_minus_Th0")

cat("writing chipScores matrix!\n")
write.table(chip.scores,file=paste(path.input,"chip_scores_pioneer_",date.macs.data.compiled,".xls",sep=""),sep="\t")

# ## plot heatmap of results
# stop("AM")
# ix <- numeric()
# for(i in 1:n){
# 	ix <- c(ix,sort(chip.scores[,i],decreasing=T,index.return=T)$ix[1:1000])
# }
# ix.u <- unique(ix)
# heatmap(chip.scores[ix.u,1:5])












