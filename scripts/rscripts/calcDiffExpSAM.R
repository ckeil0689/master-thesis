##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
path.input <- "input/th17/used_for_paper/"
library(samr)  


## load entire rnaseq dataset
load(paste(path.input,"ranseqDatasetNoQuartNorm.RData",sep=""))
## th17.wt.ix <- grep("th17.*48.*wt",colnames(rnaseq.complete.matrix),value=T)
th17.wt.ix <- grep("th17_ss_48hr_.{1,6}_wt",colnames(rnaseq.complete.matrix),value=T,perl=T)
## th0.wt.ix <- grep("th0.*48.*wt",colnames(rnaseq.complete.matrix),value=T)
th0.wt.ix <-  grep("th0_ss_48hr_.{1,6}_wt",colnames(rnaseq.complete.matrix),value=T,perl=T)
## th17.wt.ix <- colnames(rnaseq.complete.matrix)[grep("th17.*48.*wt",colnames(rnaseq.complete.matrix))]
## th0.wt.ix <- colnames(rnaseq.complete.matrix)[grep("th0.*48.*wt",colnames(rnaseq.complete.matrix))]
## add time series 48 hr experiment to th17 and th0 ( "SL1858_th17_ts_48hr#1" and "SL2673_th0_ts_48hr#1")
th17.wt.ix <- c(th17.wt.ix,"SL1858_th17_ts_48hr#1")
th0.wt.ix <- c(th0.wt.ix,"SL2673_th0_ts_48hr#1")

# data set of only th17 and th0 conditions
th17.and.th0 <- cbind(rnaseq.complete.matrix[,th0.wt.ix],rnaseq.complete.matrix[,th17.wt.ix])

ps <- 1 # a pseudo count for fixing low expression values
x <- apply(th17.and.th0,1,sum)
# leave only genes with some expression level in them
## x.ix <- which(x>median(x))
x.ix <- which(x>5)
th17.and.th0.fltrd <- th17.and.th0[x.ix,]
th17.and.th0.fltrd.psdcnt <- th17.and.th0.fltrd + ps
# log2 transform
th17.and.th0.fltrd.log = log2(th17.and.th0.fltrd.psdcnt)
rownames(th17.and.th0.fltrd.log) = toupper(rownames(th17.and.th0.fltrd.log))

# run SAM
data=list(x=th17.and.th0.fltrd.log,y=c(rep(1,length(th0.wt.ix)),rep(2,length(th17.wt.ix))), geneid=as.character(1:dim(th17.and.th0.fltrd.log)[1]),genenames=rownames(th17.and.th0.fltrd.log), logged2=TRUE)
samr.obj <- samr(data,resp.type="Two class unpaired")
delta=0
samr.plot(samr.obj,delta)
     
delta.table <- samr.compute.delta.table(samr.obj)
     
siggenes.table <-samr.compute.siggenes.table(samr.obj,delta, data, delta.table)
dval.ranks.mat <- rbind(as.matrix(as.numeric(siggenes.table$genes.up[,"Score(d)"])),as.matrix(as.numeric(siggenes.table$genes.lo[,"Score(d)"])))
qval.ranks.mat <- rbind(as.matrix(as.numeric(siggenes.table$genes.up[,"q-value(%)"])),as.matrix(as.numeric(siggenes.table$genes.lo[,"q-value(%)"])))
final.mat <- cbind(dval.ranks.mat,qval.ranks.mat)
colnames(final.mat) <- c("Score_d","q-value")
rownames(final.mat) <- c(siggenes.table$genes.up[,"Gene ID"],siggenes.table$genes.lo[,"Gene ID"])

## sam does not output some of the genes with very small zscores add those with zscore = 0
miss.gns <- rownames(data$x)[which(! rownames(data$x) %in% rownames(final.mat))]
## calc zscores for this genes (ttest two sided unpaired)
miss.gns.zscore <- apply(data$x[miss.gns,],1, function(i) t.test(i[th17.wt.ix],i[th0.wt.ix])$statistic)
x <- cbind(miss.gns.zscore,50) # put 50% fdr on each of these missing genes zscores
final.mat <- rbind(final.mat,x)
ix <- sort(abs(final.mat[,"Score_d"]),decreasing=TRUE,index.return=TRUE)$ix
final.mat <- final.mat[ix,]

# also output rpkm median for each gene and how many 0 rpkm experiments are in the dataset
gns.median <- apply(rnaseq.complete.matrix,1,median)
gns.zero <- apply(rnaseq.complete.matrix,1,function(i) length(which(i==0)))
gns.median.rpkm.and.zero <- cbind(gns.median,gns.zero)
colnames(gns.median.rpkm.and.zero) <- c("median_rpkm","number_zero_rpkm_experiments")
gns.median.rpkm.and.zero <- gns.median.rpkm.and.zero[rownames(final.mat),]
final.mat <- cbind(final.mat,gns.median.rpkm.and.zero)
write.table(final.mat,file=paste(path.input,"samTh17VsTh0Zscores.xls",sep=""),sep="\t")
