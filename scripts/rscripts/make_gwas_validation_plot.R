ix <- c("alzheimer","schizo","asthma","celiac","type2diab","type1diab","SLE","psoriasis","CD","MS","RA","UC","th17")

x.pr <- m.final.2[grep("AUCPR",rownames(m.final.2),ignore.case=T),]
x.pr.sum <- m.final["KC",grep("AUCPR",colnames(m.final),ignore.case=T)]
names(x.pr.sum) <- ix
rownames(x.pr) <- tfs
colnames(x.pr) <- ix

# remove th17 from comparison (otherwise it will dominate plot)
ix <- ix[-which(ix=="th17")]
x.pr <- x.pr[,-which(colnames(x.pr)=="th17")]
x.pr.sum <- x.pr.sum[-which(names(x.pr.sum)=="th17")]

# reorganize gold standard list (ctrl first ...)
ix <- c("alzheimer","schizo",
		"type2diab","asthma","SLE",
		"MS","UC","psoriasis","RA","type1diab","celiac","CD")
# cls.bars <- c("steelblue1", "steelblue1", 
# 			  "yellow", "yellow", "yellow",
# 			  "firebrick3", "firebrick3", "firebrick3", "firebrick3", "firebrick3", "firebrick3", "firebrick3")
cls.bars <- rep("gray",length(ix))


f.nm.res.pdf <- paste(sep="",path.output, "vldtn_GWAS_",comb.case,"_",
						c.type.all.in.one.plot,add.str,"_",gold.stdrd.date,".pdf")
pdf(f.nm.res.pdf)
ylim <- c(0,0.03)
# plot for precision recall
for(i in 1:length(ix)){ # go over datatypes
	d.type <- ix[i]
	for(k in 1:length(tfs)){ # go over tfs
		tf <- tfs[k]
		x <- x.pr[tf,d.type]
		x.sum <- x.pr.sum[d.type]
		# ylim <- c(0,max(c(x.pr,x.pr.sum)))
		xlim <- c(1,length(ix))
		if(k==1 & i==1){
			plot(1,main=paste(sep="---",type,c.tp),xlim=xlim,ylim=ylim,type="n",axes = FALSE,xlab="",ylab="Area Under Curve PR")
			axis(2,cex.axis=2)
			text(1:length(ix), par("usr")[3], labels = ix, srt = 60, adj = 1, xpd = TRUE,cex=1.5)
		}
		if(k==1){
			barplot(x.sum,space=i-0.5,width=1,col=cls.bars[i],add=T,names.arg="",axes= FALSE)
		}
		if(tf=="RORC"){
			points(y=x,x=i,pch=rorc.pch,cex=pt.cexs[k],col=cls[k],bg = cls[k])			
		} else {
			points(y=x,x=i,pch=pchs[k],cex=pt.cexs[k],col=cls[k])			
		}
	}
}
y=rand.pr
lines(x=c(-1,100),y=c(y,y),lwd=3,lty=2)
y <- sam.pr
lines(x=c(-1,100),y=c(y,y),lwd=3,lty=3)
legend(x=0.7,y=0.03,capwords(rev(c(tolower(tfs),"combined"))),col=rev(c(cls,"black")),pch=rev(c(pchs,rorc.pch,0)),pt.bg = cls[k],pt.cex=c(1,rev(pt.cexs)),cex=1.3)
# legend(x=6,y=0.04,capwords(rev(c("random","differential expression"))),col=rev(c("black","black")),lty=rev(c(2,3)),lwd=rev(c(3,3)),cex=1.3)
# legend("topright",capwords(rev(c(tolower(tfs),"combined"))),col=rev(c(cls,"black")),pch=rev(c(pchs,3,0)),cex=1.5)


ix <- c("alzheimer","schizo","asthma","type2diab","celiac","type1diab","SLE","psoriasis","CD","MS","RA","UC","th17")

x.pr <- m.final.2[grep("AUCROC",rownames(m.final.2),ignore.case=T),]
x.pr.sum <- m.final["KC",grep("AUCROC",colnames(m.final),ignore.case=T)]
names(x.pr.sum) <- ix
rownames(x.pr) <- tfs
colnames(x.pr) <- ix

# remove th17 from comparison (otherwise it will dominate plot)
ix <- ix[-which(ix=="th17")]
x.pr <- x.pr[,-which(colnames(x.pr)=="th17")]
x.pr.sum <- x.pr.sum[-which(names(x.pr.sum)=="th17")]

# reorganize gold standard list (ctrl first ...)
ix <- c("alzheimer","schizo",
		"type2diab","asthma","SLE",
		"MS","UC","psoriasis","RA","celiac","type1diab","CD")
# cls.bars <- c("steelblue1", "steelblue1", 
# 			  "yellow", "yellow", "yellow",
# 			  "firebrick3", "firebrick3", "firebrick3", "firebrick3", "firebrick3", "firebrick3", "firebrick3", "firebrick3")
cls.bars <- rep("gray",length(ix))


# plot for precision recall
for(i in 1:length(ix)){ # go over datatypes
	d.type <- ix[i]
	for(k in 1:length(tfs)){ # go over tfs
		tf <- tfs[k]
		x <- x.pr[tf,d.type]
		x.sum <- x.pr.sum[d.type]
		ylim <- c(0.4,1)
		xlim <- c(1,length(ix))
		if(k==1 & i==1){
			plot(1,main=paste(sep="---",type,c.tp),xlim=xlim,ylim=ylim,type="n",axes = FALSE,xlab="",ylab="Area Under Curve PR")
			axis(2,cex.axis=2)
			text(1:length(ix), par("usr")[3], labels = ix, srt = 60, adj = 1, xpd = TRUE,cex=1.5)
		}
		if(k==1){
			barplot(x.sum,space=i-0.5,width=1,col=cls.bars[i],add=T,names.arg="",axes= FALSE)
		}
		if(tf=="RORC"){
			points(y=x,x=i,pch=rorc.pch,cex=1.5,col=cls[k])			
		} else {
			points(y=x,x=i,pch=pchs[k],cex=1.2,col=cls[k])			
		}
	}
}
text(1:length(ix), par("usr")[3], labels = ix, srt = 60, adj = 1, xpd = TRUE,cex=1.5)
y=rand.roc
lines(x=c(-1,100),y=c(y,y),lwd=3,lty=2)
y <- sam.roc
lines(x=c(-1,100),y=c(y,y),lwd=3,lty=3)
legend(x=1,y=1,capwords(rev(c(tolower(tfs),"combined"))),col=rev(c(cls,"black")),pch=rev(c(pchs,rorc.pch,0)),pt.bg = cls[k],pt.cex=c(1,rev(pt.cexs)),cex=1.3)
legend(x=6,y=1,capwords(rev(c("random","differential expression"))),col=rev(c("black","black")),lty=rev(c(2,3)),lwd=rev(c(3,3)),cex=1.3)

# legend(x=1,y=1,capwords(rev(c(tolower(tfs),"combined"))),col=rev(c(cls,"black")),pch=rev(c(pchs,3,0)),cex=1.5)
dev.off()

# 
# f.nm.res.xls <- paste(sep="",path.output, "vldtn_perTF_AllInOne2_",comb.case,"_",
# 						c.type.all.in.one.plot,add.str,"_gs_data_","COMBINED","_",gold.stdrd.date,".xls")
# cat(sep="\t",file=f.nm.res.xls,"comb.type",colnames(m.final),"\n")						
# write.table(m.final,file=f.nm.res.xls,sep="\t",append=T,col.names=F,row.names=T)
# f.nm.res.xls <- paste(sep="",path.output, "vldtn_perTF_AllInOne2_",comb.case,"_",
# 						c.type.all.in.one.plot,add.str,"_gs_data_","SINGLE_KC","_",gold.stdrd.date,".xls")
# cat(sep="\t",file=f.nm.res.xls,"comb.type",colnames(m.final.2),"\n")						
# write.table(m.final.2,file=f.nm.res.xls,sep="\t",append=T,col.names=F,row.names=T)