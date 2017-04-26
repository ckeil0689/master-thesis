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

path.input <- "input/th17/used_for_paper/"
## path.output <- paste(sep="","results/validation_",date.is,"/")
path.output <- paste(sep="","results/validation/")
path.input.sam <- "input/th17/used_for_paper/"

# define file names
## combined analysis files
s.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_S_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
k.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_K_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
c.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_C_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
r.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_R_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
i.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_I_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
kc.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_KC_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
ri.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_RI_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
kcr.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_KCR_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
kcri.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_KCRI_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
## skcri.f.nm <- paste(sep="",path.input,"combinedAnalysis/","score_combine_SKCRI_matrix_zcut_",z.abs.cut,"_",type,"_",date.combine.data.run,".xls")
## sam.f.nm <-paste(sep="",path.input.sam,"samTh17VsTh0Zscores.xls")


## gold standard files
th17.reg.inter.f.nm <- paste(sep="",path.input,"gold_standard_th17_reg_inter.txt")
th17.gns.f.nm <- paste(sep="",path.input,"gold_standard_th17_genes_Sep_15_2011.txt")
## th17.gns.f.nm <- paste(sep="",path.input,"gold_standard_th17_genes.txt")

# load files
## combined analysis files
s <- as.matrix(read.table(s.f.nm,header=T,sep="\t"))
k <- as.matrix(read.table(k.f.nm,header=T,sep="\t"))
c <- as.matrix(read.table(c.f.nm,header=T,sep="\t"))
r <- as.matrix(read.table(r.f.nm,header=T,sep="\t"))
i <- as.matrix(read.table(i.f.nm,header=T,sep="\t"))
kc <- as.matrix(read.table(kc.f.nm,header=T,sep="\t"))
ri <- as.matrix(read.table(ri.f.nm,header=T,sep="\t"))
kcr <- as.matrix(read.table(kcr.f.nm,header=T,sep="\t"))
kcri <- as.matrix(read.table(kcri.f.nm,header=T,sep="\t"))
## skcri <- as.matrix(read.table(skcri.f.nm,header=T,sep="\t"))


## gold standard files
gs.reg.inter <- toupper(read.delim(sep=" ",header=F,th17.reg.inter.f.nm, as.is=T))
gs.gns <- read.delim(sep="\t",header=T,th17.gns.f.nm, as.is=T)
gs.gns$Gene <- toupper(gs.gns$Gene)

## get single expt type sum for the five core tfs
k.sum <- apply(k[,core.tfs],1,sum)
c.sum <- apply(c[,core.tfs],1,sum)
r.sum <- apply(r[,core.tfs],1,sum)
i.sum <- apply(i[,core.tfs],1,sum)

## get kc sum for the five core tfs
kc.sum <- apply(kc[,core.tfs],1,sum)

## get ri sum for the five core tfs
ri.sum <- apply(ri[,core.tfs],1,sum)
## get kcri sum for the five core tfs
## kcri.sum <- apply(cbind(kc.sum,ri.sum),1,sum)
kcr.sum <- apply(kcr[,core.tfs],1,sum)
kcri.sum <- apply(kcri[,core.tfs],1,sum)
## skcri.sum <- skcri[,"sum.score"]


## get sam diff exression zscores
## sam.score <- as.matrix(read.table(sam.f.nm,sep="\t"))
## sam.score <- sam.score[names(kcri.sum),"Score_d"]
sam.score <- as.numeric(s)
#####################
## get gene relative ranks
x <- cbind(sam.score,k.sum,c.sum,r.sum,i.sum,kc.sum,ri.sum,kcr.sum,kcri.sum)
gns <- rownames(x)
sam.rel.ranks <- gns[order(x[,"sam.score"],decreasing=T)]
k.rel.ranks <- gns[order(x[,"k.sum"],decreasing=T)]
c.rel.ranks <- gns[order(x[,"c.sum"],decreasing=T)]
r.rel.ranks <- gns[order(x[,"r.sum"],decreasing=T)]
i.rel.ranks <- gns[order(x[,"i.sum"],decreasing=T)]
ri.rel.ranks <- gns[order(x[,"ri.sum"],decreasing=T)]
kc.rel.ranks <- gns[order(x[,"kc.sum"],decreasing=T)]
kcr.rel.ranks <- gns[order(x[,"kcr.sum"],decreasing=T)]
kcri.rel.ranks <- gns[order(x[,"kcri.sum"],decreasing=T)]
## skcri.rel.ranks <- gns[order(x[,"skcri.sum"],decreasing=T)]

x.rel <- matrix(seq(from=0,to=1,length.out=dim(x)[1]),nr=dim(x)[1],nc=dim(x)[2])
x.rel.nms <- cbind(sam.rel.ranks,k.rel.ranks,c.rel.ranks,r.rel.ranks,i.rel.ranks,ri.rel.ranks,kc.rel.ranks,kcr.rel.ranks,kcri.rel.ranks)
colnames(x.rel) <- colnames(x.rel.nms)

## find the relative rank of gold standard Th17 known genes
w <- list()
w[["SAM"]] <- x.rel[which(x.rel.nms[,"sam.rel.ranks"]%in%gs.gns$Gene),"sam.rel.ranks"]
w[["k"]] <- x.rel[which(x.rel.nms[,"k.rel.ranks"]%in%gs.gns$Gene),"k.rel.ranks"]
w[["c"]] <- x.rel[which(x.rel.nms[,"c.rel.ranks"]%in%gs.gns$Gene),"c.rel.ranks"]
w[["r"]] <- x.rel[which(x.rel.nms[,"r.rel.ranks"]%in%gs.gns$Gene),"r.rel.ranks"]
w[["i"]] <- x.rel[which(x.rel.nms[,"i.rel.ranks"]%in%gs.gns$Gene),"i.rel.ranks"]
stop("AM")
w[["ri"]] <- x.rel[which(x.rel.nms[,"ri.rel.ranks"]%in%gs.gns$Gene),"ri.rel.ranks"]
w[["kc"]] <- x.rel[which(x.rel.nms[,"kc.rel.ranks"]%in%gs.gns$Gene),"kc.rel.ranks"]
## w[["kcr"]] <- x.rel[which(x.rel.nms[,"kcr.rel.ranks"]%in%gs.gns$Gene),"kcr.rel.ranks"]
w[["kcri"]] <- x.rel[which(x.rel.nms[,"kcri.rel.ranks"]%in%gs.gns$Gene),"kcri.rel.ranks"]
## w[["skcri"]] <- x.rel[which(x.rel.nms[,"skcri.rel.ranks"]%in%gs.gns$Gene),"skcri.rel.ranks"]

rel.ranks.all <- unlist(w)
tmp <- numeric()
for(i in 1:length(w)){tmp <- c(tmp,rep(i,length(w[["k"]])))}
names(rel.ranks.all) <- tmp

f.nm <- paste(sep="",path.output,"compare_net_performance_identify_th17_core_genes_",type,"_",date.is,".pdf")
pdf(f.nm)
for(i in length(w):1){
  if(i==length(w)){
    ## cls <- c("darkblue","darkgreen","darkred","black")
    cls <- c(rainbow(length(w)),"black")
    ## l.types <- c("dotted","dotted","dotted","dotted","twodash","dashed","solid","solid","solid")
    l.types <- c("dashed","dashed","dashed","dashed","dashed","solid","solid","solid","solid")
    v <- density(rel.ranks.all[which(names(rel.ranks.all)==i)],from=0, to=1)
    plot(v,xlim=c(0,1),col=cls[i],
        lty=l.types[i],xlab="relative rank of previously characterized th17 genes",lwd=2,cex.lab=1.5,main="",xaxt="n",yaxt="n",ylim=c(0,9))
    axis(1,at=seq(0,1,by=0.1),labels=seq(1,0,by=-0.1))
    axis(2,at=seq(0,9,by=1))
    v.slim <- v
    v.slim$x <- v$x[c(20)]
    v.slim$y <- v$y[c(20)]
    points(v.slim,pch=as.character(i))
  } else {
    v <- density(rel.ranks.all[which(names(rel.ranks.all)==i)],from=0, to=1)
    lines(v,xlim=c(0,1),col=cls[i],lty=l.types[i],lwd=2)
        v.slim <- v
    v.slim$x <- v$x[c(20)]
    v.slim$y <- v$y[c(20)]
    points(v.slim,pch=as.character(i))
  }  
}
v.slim$y <- rep(1,length(v.slim$y))
lines(x=c(0,1),y=c(1,1),lty=l.types[length(w)+1],lwd=2,col=cls[length(w)+1])
points(v.slim,xlim=c(0,1),col=cls[length(w)+1],lty=l.types[i],lwd=2,pch=as.character(0))
legend(x=c(0.8),y=c(7), c(rev(names(w)),"random"), col=c(rev(cls[-length(cls)]),cls[length(cls)]),
       lty=c(rev(l.types[-length(l.types)]),l.types[length(l.types)]),lwd=2,pch=as.character(length(w):0)) 
dev.off()

#####################

f.nm <- paste(sep="",path.output,"kc.ri.sum.xls")
cat(file=f.nm,sep="\t",c("gene_name",colnames(x)),"\n")
write.table(x,file=f.nm,sep="\t",append=T,col.names=F)


####################
## # get a plot of cluster number and composition as a function of tf.dist
## rm(list=ls()[-which(ls()=="GLOBAL")])

## date.is <- "Apr_24_2011"
## z.abs.cut <- 2.5
## path.input <- "results/cluster_analysis/"
## path.output <- paste(sep="","results/validation_",date.is,"/")

## # get all th17 ChIP clustering results
## fl.nms <- list.files(path.input)
## fl.nms <- fl.nms[grep("th17.*Aug_23_2011",fl.nms)]
## l.th17 <- list()
## tf.dist <- sort(sapply(fl.nms,function(i) as.numeric(strsplit(i,"_")[[1]][4])))
## for(i in 1:length(tf.dist)){
##   x <- read.table(paste(sep="",path.input,names(tf.dist)[i]),sep="\t",header=T)
##   l.th17[[as.character(tf.dist)[i]]] <- x$tf.to.s
## }

## len.l.th17 <- sapply(l.th17,length)
## plot(x=tf.dist,y=len.l.th17,ylab="cluster [#]",xlab="allowed tf summit distance [bp]")

## x <- table(l.th17[[1]])
## w <- matrix(0,nr=length(l.th17),nc=length(x) )
## colnames(w) <- names(x)
## rownames(w) <- names(l.th17)
## for(i in 1:length(l.th17)){
##   nm <- names(l.th17)[i]
##   x <- table(l.th17[[i]])
##   w[nm,names(x)] <- x
## }
## w.norm <- w/apply(w,1,sum)*100
## ix <- which(apply(w,2,mean)>500)
## w.n.s <- w.norm[,ix]
## w.s <- w[,ix]

## barplot(t(w.s),legend = colnames(w.s),col=rainbow(dim(w.s)[2]))
