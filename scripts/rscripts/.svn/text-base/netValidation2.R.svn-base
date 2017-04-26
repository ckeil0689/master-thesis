##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Oct 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# R code to create quality control plots for entire network pipeline

path.input <- "input/th17/used_for_paper/"
## path.output <- paste(sep="","results/validation_",date.is,"/")
path.output <- paste(sep="","results/validation/",date.is,"/")
path.input.sam <- "input/th17/used_for_paper/"

# define file names
i.f.nm <- paste(sep="",path.input,"infResults/immgen/savedData_numBoots_200.RData")
r.f.nm <- paste(sep="",path.input,"infResults/rnaseq/savedData_numBoots_200.RData")
th17.reg.inter.f.nm <- paste(sep="",path.input,"gold_standard_th17_reg_inter.txt")
th17.gns.f.nm <- paste(sep="",path.input,"gold_standard_th17_genes_",gold.stdrd.date,".txt")

# load files
load(i.f.nm)
res.i <- betaList[[1]]
gn.nms.i <- rownames(INPUT$general$dataset)
load(r.f.nm)
res.r <- betaList[[1]]
gn.nms.r <- rownames(INPUT$general$dataset)

gs.reg.inter <- toupper(read.delim(sep=" ",header=F,th17.reg.inter.f.nm, as.is=T))
gs.gns <- read.delim(sep="\t",header=T,th17.gns.f.nm, as.is=T)
gs.gns$Gene <- toupper(gs.gns$Gene)
# gs.gns.pos <- gs.gns$Gene[which(gs.gns$Effect=="positive")]
# gs.gns.neg <- gs.gns$Gene[which(gs.gns$Effect=="negative")]

x <- res.i[[1]]
w.i=numeric(length=max(x[,"trgt"]))
names(w.i) <- gn.nms.i
for(i in 1:max(x[,"trgt"])){w.i[i]=sum(x[which(x[,"trgt"]==i),"prd_xpln_var"])}

x <- res.r[[1]]
w.r=numeric(length=max(x[,"trgt"]))
names(w.r) <- gn.nms.r
for(i in 1:max(x[,"trgt"])){w.r[i]=sum(x[which(x[,"trgt"]==i),"prd_xpln_var"])}

gn.nms.both <- intersect(gn.nms.i,gn.nms.r)
# gs.gns.both.pos <- intersect(gs.gns.pos,gn.nms.both)
# gs.gns.both.neg <- intersect(gs.gns.neg,gn.nms.both)

f.nm <- paste(sep="",path.output,"predict_data_",date.is,".pdf")
pdf(f.nm)
plot(x=w.r[gn.nms.both],y=w.i[gn.nms.both],main="portion of explaind variance",
     xlab="portion of explaind variance (Inf on RNAseq)",
     ylab="portion of explaind variance (Inf on Immgen)",pch=20,cex=.3,col="gray")
lines(x=c(-1,2),y=c(0.5,0.5),lty="dashed")
lines(x=c(0.5,0.5),y=c(-1,2),lty="dashed")
# text(x=w.r[gs.gns.both.pos],y=w.i[gs.gns.both.pos],labels=gs.gns.both.pos,cex=.6,col="blue")
# text(x=w.r[gs.gns.both.neg],y=w.i[gs.gns.both.neg],labels=gs.gns.both.neg,cex=.6,col="darkgreen")
text(x=w.r[gs.gns$Gene],y=w.i[gs.gns$Gene],labels=gs.gns$Gene,cex=.7,col="black")
dev.off()
# legend("topleft",c("All","Gold Standard"),fill=c("gray","black"))

# legend("topleft",c("All","GS pos effect", "GS neg effect"),fill=c("black","blue","darkgreen"))
# 
# v=density(w.r[which(names(w.r) %in% gs.gns.pos)])
# plot(v,col="blue",main="predicting data (Inf on RNAseq)",xlab="portion of explaind variance")
# v=density(w.r)
# lines(v,col="black")
# v=density(w.r[which(names(w.r) %in% gs.gns.neg)])
# lines(v,col="darkgreen")
# legend("topleft",c("All","GS pos effect", "GS neg effect"),col=c("black","blue","darkgreen"),lty=rep(1,3))
# 
# v=density(w.i[which(names(w.i) %in% gs.gns.pos)])
# plot(v,col="blue",main="portion of explaind variance",xlab="genes")
# v=density(w.i)
# lines(v,col="black")
# v=density(w.i[which(names(w.i) %in% gs.gns.neg)])
# lines(v,col="darkgreen")
# legend("topleft",c("All","GS pos effect", "GS neg effect"),col=c("black","blue","darkgreen"),lty=rep(1,3))
# dev.off()
