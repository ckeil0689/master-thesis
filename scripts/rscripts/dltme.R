	pdf(f.nm)
	cls <- c("firebrick3",
				rgb(255, 106, 106,maxColorValue=255, 150),
				rgb(154, 205, 50,maxColorValue=255, 150),
				"darkgray")
	nf <- layout(rbind(1,2,3),heights=c(3,1,2))
	op <- par(mar=c(0,7,0.5,0.5),mgp=c(4, 1, 0))
	n <- length(RES[["KCRI"]])
	y <- RES[["KCRI"]]
	x <- 1:n
	plot(x=x,y=y,type="l",lwd=3,col=cls[4],cex.lab=2,ylim=c(0,1),
	ylab="Enrichment Score (ES)",xaxt="n",xlab="",cex.axis=2,las=1)
	y <- RES.1[["KCRI"]]
	lines(x=x,y=y,type="l",lwd=3,col=cls[1])
	legend(y=1,x=(n-8000),legend=c("K+C+R+I","Random"),lty=1,lwd=3,col=cls[c(1,4)],cex=2)

	plot(y=c(-0.05,0.05),x=c(1,1),col="gray89",type="l",lty=1,xlim=c(0,n),lwd=0.1,
	xlab="",ylab="",axes=F,ylim=c(-1,1))
	s <- IHIT[["KCRI"]]
	for(i in 1:n){
		if(s[i]==0){
			lines(y=c(-1,-0.1),x=c(i,i),col="gray89",lwd=0.1)
		} else {
			lines(y=c(-1,-0.1),x=c(i,i),col="darkred",lwd=2)				
		}
	}
	s <- IHIT.1[["KCRI"]]
	for(i in 1:n){
		if(s[i]==0){
			lines(y=c(0.1,1),x=c(i,i),col="gray89",lwd=0.1)
		} else {
			lines(y=c(0.1,1),x=c(i,i),col="darkred",lwd=2)				
		}
	}
	# add blue line where leading edge (enriched set end)
	lines(y=c(-1,-0.1),x=c(arg.ES[["KCRI"]],arg.ES[["KCRI"]]),col="darkblue",lwd=5)
	lines(y=c(0.1,1),x=c(arg.ES.1[["KCRI"]],arg.ES.1[["KCRI"]]),col="darkblue",lwd=5)
	text(y=c(0.5),x=c(arg.ES.1[["KCRI"]])+50,labels="ES max",cex=2,pos=4,col="darkblue")
	text(y=c(-0.5),x=c(arg.ES[["KCRI"]])+50,labels="ES max",cex=2,pos=4,col="darkblue")
	mtext("K+C+R+I",side=2,cex=1.5,las=1,line = -1,at=0.5)
	mtext("Random",side=2,cex=1.5,las=1,line = -1,at=-0.5)
		
	par(mar=c(4.5,7,0.5,0.5),mgp=c(3, 1, 0))
	x <- sort(o.sum.list.1[,"KCRI"],decreasing=T)
	plot(y=c(x[1]),x=c(1),col=cls[1],pch=20,lty=1,xlim=c(0,n),cex.axis=2,cex.lab=2,las=1,lwd=0.5,
	ylab="Network score",xlab="Ranked gene list",ylim=c(min(x),max(x)))
	legend(y=max(x)-1,x=(n-8000),legend=c("K+C+R+I","random"),lty=1,lwd=3,col=cls[c(1,4)],cex=2)
	lines(x=1:n,y=x,pch=20,col=cls[1],lwd=4)
	x <- sort(o.sum.list[,"KCRI"],decreasing=T)
	lines(x=1:n,y=x,pch=20,col=cls[4],lwd=4)			
	dev.off()