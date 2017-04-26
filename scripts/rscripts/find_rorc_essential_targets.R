##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jun 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

#####################################################################
# 0. set up variables
#####################################################################
cat("getting essential targets for rorc\n")

path.input.deseq <- paste(sep="","input/th17/used_for_paper/DEseq/",date.deseq.run,"/")
file.nm.ko.pval <- paste(sep="","DEseq_pval_signed_",date.deseq.run,".xls")
file.nm.ko.fc <- paste(sep="","DEseq_log2fc_",date.deseq.run,".xls")

# get ko pvals for rorc ko/wt
f.nm <- paste(sep="",path.input.deseq,file.nm.ko.pval)
cat("loading KO scores file: ",f.nm,"\n")
k.pval <- as.matrix(read.table(f.nm,sep="\t"))
k.pval <- k.pval[,"Th17.rorc.wt.vs.Th17.rorc.ko"]
# get ko log2 fc for rorc ko/wt
f.nm <- paste(sep="",path.input.deseq,file.nm.ko.fc)
cat("loading KO scores file: ",f.nm,"\n")
k.fc <- as.matrix(read.table(f.nm,sep="\t"))
k.fc <- k.fc[,"Th17.rorc.wt.vs.Th17.rorc.ko"]

# get histone marks in response to rorc ko
x <- read.table(expts.file.nm,sep="\t",header=T,colClasses = "character")
p300.fc.ix <- grep("FC_SL3594",colnames(x),value=T)
p300.pval.ix <- grep("pval_SL3594",colnames(x),value=T)
k4me2.fc.ix <- grep("FC_SL4032",colnames(x),value=T)
k4me2.pval.ix <- grep("pval_SL4032",colnames(x),value=T)
k4me3.fc.ix <- grep("FC_SL4034",colnames(x),value=T)
k4me3.pval.ix <- grep("pval_SL4034",colnames(x),value=T)
ix.conds <- c(p300.fc.ix,p300.pval.ix,
				k4me2.fc.ix,k4me2.pval.ix,
				k4me3.fc.ix,k4me3.pval.ix)
# find genes that were targeted by rorc and had a |log2(FC)|>2 and -log10(pval)>10 based on histone marks and p300
ix.mtls.rorc <- grep("rorc",x[,"tf.to.s"],ignore.case=T)
ix.mtls.rorc.with.trgts <- ix.mtls.rorc[which(x[ix.mtls.rorc,"trgts"]!="")]
# here we only focus on rorc
m <- x[ix.mtls.rorc.with.trgts,]
for(i in 1:length(ix.conds)){
	m[,ix.conds[i]] <- as.numeric(m[,ix.conds[i]])	
}


abs.fc.cut <- 1
pval.cut <- 0
cond.1 <- numeric(nrow(m))
cond.2 <- numeric(nrow(m))
cond.3 <- numeric(nrow(m))
cond.4 <- numeric(nrow(m))
for(i in 1:nrow(m)){
	cond.1[i] <- abs(m[i,p300.fc.ix])>abs.fc.cut & m[i,p300.pval.ix]>pval.cut
	cond.2[i] <- abs(m[i,k4me2.fc.ix])>abs.fc.cut & m[i,k4me2.pval.ix]>pval.cut
	cond.3[i] <- abs(m[i,k4me3.fc.ix])>abs.fc.cut & m[i,k4me3.pval.ix]>pval.cut	
}
cond.marks <- cond.1+cond.2+cond.3

# get all genes that my cond.3
trgts <- unique(m[which(cond.3==1),"trgts"])
trgts.split <- unlist(sapply(trgts,function(i) strsplit(i,"_")[[1]]))
trgts.split <- trgts.split[-which(trgts.split=="")]
trgts <- sort(unique(trgts.split))


s=cbind(trgts,k.fc[trgts],k.pval[trgts])
f.nm <- "results/diff_chip_volcanoplots/Jul_27_2012/rorc_essential_targets.xls"
cat("gene_id","log2FC_rorc_ko","pval_rorc_ko","\n",file=f.nm,sep="\t")
write.table(file=f.nm,s,row.names=F,col.names=F,sep="\t",append=T)

# w <- matrix(0,nr=length(trgts),ncol=8)
# rownames(w) <- trgts
# colnames(w) <- c("p300-fc","k4me2-fc","k4me3-fc","ko-fc","p300-pval","k4me2-pval","k4me3-pval","ko-pval")
# for(i in 1:length(trgts)){
# 	trgt <- trgts[i]
# 	ix <- grep(trgt,m[,"trgts"])
# 	p300.fc <- m[ix,p300.fc.ix][which.max(abs(m[ix,p300.fc.ix]))]
# 	p300.pval <- m[ix,p300.pval.ix][which.max(abs(m[ix,p300.pval.ix]))]
# 	k4me2.fc <- m[ix,k4me2.fc.ix][which.max(abs(m[ix,k4me2.fc.ix]))]
# 	k4me2.pval <- m[ix,k4me2.pval.ix][which.max(abs(m[ix,k4me2.pval.ix]))]
# 	k4me3.fc <- m[ix,k4me3.fc.ix][which.max(abs(m[ix,k4me3.fc.ix]))]
# 	k4me3.pval <- m[ix,k4me3.pval.ix][which.max(abs(m[ix,k4me3.pval.ix]))]
# 	w[i,c("p300-fc","k4me2-fc","k4me3-fc","p300-pval","k4me2-pval","k4me3-pval")] <- 
# 			c(p300.fc,k4me2.fc,k4me3.fc,p300.pval,k4me2.pval,k4me3.pval)
# }
# 
# trgt <- m[i,"trgts"]
# cond.4[i] <- abs(k.fc[trgt])>abs.fc.cut & k.pval[trgt]>pval.cut	
# 
# 
# cond.all <- cond.marks+cond.4
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
