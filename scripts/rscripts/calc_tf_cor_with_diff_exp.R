##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Apr 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

path.output <- paste(sep="","results/validation/",date.is,"/")
system(paste(sep="","mkdir ", path.output))			
path.input.sam <- "input/th17/used_for_paper/"
path.input.combined <- paste(sep="","results/combinedAnalysis/",date.combine.data.run,"/")

add.str <- paste(sep="","_cut_prcnt_",prcnt.chng.cut,"_num_tfs_",num.tfs,"_sam_",z.abs.cut,"_deseq_cut_",deseq.pval.cut)

infl.nm <- paste(path.input.combined,"results_combine_data_exprsn",add.str,".Rdata",sep="")
# load input files
load(infl.nm)

# get kc whole results for each TF
tp <- comb.case # whole, abs, activation, or repression
c.tp <- c.type # th17, th0, th17 minus th0
cb.tp <- comb.type # k,c,r,i,kc,...
m <- matrix(0,nr=nrow(res[[tp]][[1]][[c.tp]]),nc=length(tfs))
rownames(m) <- rownames(res[[tp]][[1]][[c.tp]])
colnames(m) <- tfs
for(j in 1:length(tfs)){ # go over tfs
	tf <- tfs[j]
	scores.tf <- res[[tp]][[tf]][[c.tp]][,cb.tp]
	m[,tf] <- scores.tf
}

# get sam values for genes
if(diff.exp.by.sam.or.fc=="SAM"){
	sam.scores <- res[[tp]][[1]][[c.tp]][,"SAM"]
} else {
	sam.scores <- log2((res[[tp]][[1]][[c.tp]][,"th17_rpkm"]+1)/(res[[tp]][[1]][[c.tp]][,"th0_rpkm"]+1))
}

# get corelation for each TF with its own strong target genes differential expression
cor.tf.score.with.diff.exp <- numeric(length(tfs))
names(cor.tf.score.with.diff.exp) <- tfs
l <- list()
for(j in 1:length(tfs)){
	tf <- tfs[j]
	gns.regulated.by.tf.up <- names(which(m[,tf]>comb.type.cut))
	gns.regulated.by.tf.down <- names(which(m[,tf] < (-1*comb.type.cut) ))
	gns.regulated.by.tf <- names(which(abs(m[,tf])>comb.type.cut))
	scores.tf <- m[gns.regulated.by.tf,tf]
	sam.scores.tf.trgts <- sam.scores[gns.regulated.by.tf] 
	l[[paste(sep="",tf,"-","up.reg")]] <- sam.scores.tf.trgts[gns.regulated.by.tf.up]
	l[[paste(sep="",tf,"-","down.reg")]] <- sam.scores.tf.trgts[gns.regulated.by.tf.down]
	cor.tf.score.with.diff.exp[tf] <- cor(scores.tf,sam.scores.tf.trgts)
}



fl.nm <- paste(sep="",path.output,"scores_tf_vs_score_sam_prcnt_",prcnt.chng.cut,"_pval_",deseq.pval.cut,"_numtfs_",num.tfs,".pdf")
pdf(fl.nm)
mp <- boxplot(l,outline=F,xaxt = "n",ylab=diff.exp.by.sam.or.fc,col=c("green","red"))
lines(x=c(-10,100),y=c(0,0))
text(0.5:length(l), par("usr")[3], labels = paste(sep="",names(l),"\n (n=",mp$n,")"), srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.7)
dev.off()

# get corelation for each TF with its own strong target genes differential expression (only over genes that are regulated by more than 4 tfs)
w <- apply(m[,1:5],1,function(i) length(which( abs(i)> comb.type.cut )))
m2 <- m[which(w>=4),]
m2[m2==0] <- NA

cor.tf.score.with.diff.exp <- numeric(length(tfs))
names(cor.tf.score.with.diff.exp) <- tfs
l <- list()
for(j in 1:length(tfs)){
	tf <- tfs[j]
	gns.regulated.by.tf.up <- names(which(m2[,tf]>comb.type.cut))
	gns.regulated.by.tf.down <- names(which(m2[,tf] < (-1*comb.type.cut) ))
	gns.regulated.by.tf <- names(which(abs(m2[,tf])>comb.type.cut))
	scores.tf <- m2[gns.regulated.by.tf,tf]
	sam.scores.tf.trgts <- sam.scores[gns.regulated.by.tf] 
	l[[paste(sep="",tf,"-","up.reg")]] <- sam.scores.tf.trgts[gns.regulated.by.tf.up]
	l[[paste(sep="",tf,"-","down.reg")]] <- sam.scores.tf.trgts[gns.regulated.by.tf.down]
	cor.tf.score.with.diff.exp[tf] <- cor(scores.tf,sam.scores.tf.trgts)
}

fl.nm <- paste(sep="",path.output,"4-5way_scores_tf_vs_score_sam_prcnt_",prcnt.chng.cut,"_pval_",deseq.pval.cut,"_numtfs_",num.tfs,".pdf")
pdf(fl.nm)
mp <- boxplot(l,outline=F,xaxt = "n",ylab=diff.exp.by.sam.or.fc,col=c("green","red"))
lines(x=c(-10,100),y=c(0,0))
text(0.5:length(l), par("usr")[3], labels = paste(sep="",names(l),"\n (n=",mp$n,")"), srt = 45, adj = c(1,1.5), xpd = TRUE,cex=.7)
dev.off()






























