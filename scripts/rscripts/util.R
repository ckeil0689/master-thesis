##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## January 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# adapted to compare two prediction lists
GSEA.EnrichmentScore <- function(run.perm=FALSE,gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
#
# Computes the weighted GSEA score of gene.set in gene.list. 
# The weighted score type is the exponent of the correlation 
# weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
# necessary to input the correlation vector with the values in the same order as in the gene list.
#
# Inputs:
#   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
#   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
#  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
#
# Outputs:
#   ES: Enrichment score (real number between -1 and +1) 
#   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
#   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
#   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

	# if( !is.null(run.perm) | run.perm == T ){
	# 	gene.set <- sample(gene.list,length(gene.set))
	# }
	if( run.perm ){
		gene.set <- sample(gene.list,length(gene.set))
	}

   tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
   no.tag.indicator <- 1 - tag.indicator 
   N <- length(gene.list) 
   Nh <- length(gene.set) 
   Nm <-  N - Nh 
   if (weighted.score.type == 0) {
      correl.vector <- rep(1, N)
   }
   alpha <- weighted.score.type
   correl.vector <- abs(correl.vector**alpha)
   sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
   norm.tag    <- 1.0/sum.correl.tag
   norm.no.tag <- 1.0/Nm
   RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
   max.ES <- max(RES)
   min.ES <- min(RES)
   if (any(is.na(max.ES))){
	  RES <- numeric(length(RES))
      ES <- signif(0, digits = 5)
      arg.ES <- 1
   } else if (max.ES > - min.ES) {
#      ES <- max.ES
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
   } else {
#      ES <- min.ES
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
   }
   return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}

# validation script taken from DREAM3 and adapted to R by Aviv Madar
# The function evaluate evaluates the accuracy of a prediction compared to
# a gold standard
# 
# Usage: DreamEvaluationScript(testMat, goldStdrdMat)
#
# evaluate calcuates a precision and ROC values for a given prediction 
# (testMat) compared to a gold standard (goldStdrdMat).
# The input (testMat) holds the information on which goldstandard 
# should be used for comparison. fname is a file name for logging 
# purposes it is used to collect results from different executions of the script 
# in one file. 
# 
# All further information about the actual calculations can be found elsewhere.
#
#
# The algorithm was developed by Gustavo Stolovitzky and implemented by
# Bernd Jagla, Columbia University (baj2107_A_T_columbia.edu). All questions/suggestions should be
# directed to both, Gustavo and Bernd.
#
# Aviv Madar, New York University (am2654@nyu.edu), modified it and converted it from MATLAB to R.

evaluate.two.pred.lists <- function(testMat, goldStdrdMat) {
	#must modify this script when using not DREAM data, for example: in e.coli i,i interactions are allowed
	#because the regulators are a different set of genes than the tfs!
	maxN = 10000000;
	tstRnkdLst = sort(testMat,decreasing=TRUE,index.return=TRUE)[[2]]
# convert vectorized index to double array index
	tstRnkdLst1 = tstRnkdLst %% nrow(goldStdrdMat)
	tstRnkdLst1[which(tstRnkdLst1 == 0)] = nrow(goldStdrdMat)
	tstRnkdLst2 = tstRnkdLst/nrow(goldStdrdMat)
	tstRnkdLst2[which(! tstRnkdLst2%%1==0)] = as.integer(tstRnkdLst2[which(! tstRnkdLst2%%1==0)]+1)
# remove diagonal enteries we are not predicting auto reg
	rm = which(tstRnkdLst1 == tstRnkdLst2)
	tstRnkdLst1 = tstRnkdLst1[-rm]
	tstRnkdLst2 = tstRnkdLst2[-rm]
	
	
	# Analysis
	# initialization
	
	k=0
	Ak=0
	TPk=0
	FPk=0
	P = sum(goldStdrdMat)
	N = length(goldStdrdMat)-P # Diagonal is not considered
	T = P+N
	L = length(tstRnkdLst1)                                               
	rec = numeric(length=L)
	prec = numeric(length=L)
	tpr = numeric(length=L)
	fpr = numeric(length=L)
	
	if(length(tstRnkdLst1) > 0) {
		while(k < L) {
			k = k + 1
			# if ((k %% 1000) == 0) {	cat(k, "out of", L,"\n") } ##huh? what does this line mean??
			# if k pred is correct
			if(goldStdrdMat[tstRnkdLst1[k],tstRnkdLst2[k]] == 1) {
				TPk = TPk + 1
				if(k==1)	
					{ delta=1/P }  
				else		
					{ delta=(1-FPk*log(k/(k-1)))/P }
				Ak = Ak + delta
			} else if (goldStdrdMat[tstRnkdLst1[k],tstRnkdLst2[k]] == 0) {
				FPk = FPk + 1
			} else {
				cat("could not find ", tstRnkdLst1[k],tstRnkdLst2[k]);
			}
			rec[k] = TPk/P
			prec[k] = TPk/k
			tpr[k] = rec[k]
			fpr[k] = FPk/N
		}
	}

	TPL=TPk
	if (L < T) 	{ 
		rh = (P-TPL)/(T-L) 
	} else { 
		rh = 0 
	}
	
	if (L>0) {
		recL = rec[L]
	} else {
		recL = 0
	}
	while (TPk < P) {
		k = k + 1
		TPk = TPk + 1
		rec[k] = TPk/P
		if ( ((rec[k]-recL)*P + L * rh) != 0 ) {
			prec[k] = rh * P * rec[k]/((rec[k]-recL)*P + L * rh)
		} else {
			prec[k] = 0
		}
		tpr[k] = rec[k]
		FPk = TPk * (1-prec[k])/prec[k]
		fpr[k] = FPk/N
	}
	AL = Ak;
	if (!is.nan(rh) & rh != 0  & L != 0) {
		AUC = AL + rh * (1-recL) + rh * (recL - L * rh / P) * log((L * rh + P * (1-recL) )/(L *rh))
	} else if (L==0) {
		AUC = P/T
	} else {
		AUC = Ak
	}
	
	# Integrate area under ROC
	lc = fpr[1] * tpr[1] /2
	for (n in 1:(L+P-TPL-1)) {
		lc = lc + (fpr[n+1]+fpr[n]) * (tpr[n+1]-tpr[n]) / 2
	}
	
	AUROC = 1 - lc

	# specific precision values
	TrueP = c(1, 2, 5, 20, 100, 500) 
	prec0 = numeric(length=length(TrueP))
	names(prec0) = as.character(TrueP)
	rec0 = numeric(length=length(TrueP))
	for (i in 1:length(TrueP)) {
		if(TrueP[i]<=P) {
			rec0[i]=TrueP[i]/P
			j=which(rec == rec0[i])
			j=min(j) #In case there is more than 1 precision values for rec(i)
			prec0[i]=prec[j]
		}
	}

	# handle output
	object = list()
	object[["prec"]] = prec
	object[["rec"]] = rec
	object[["tpr"]] = tpr
	object[["fpr"]] = fpr
	object[["AUROC"]] = AUROC
	object[["AUPR"]] = AUC
	object[["specificPrecVals"]] = prec0
	return(object)

}
# same as above but for running random permutations of goldStdrdMat permutations
evaluate.two.pred.lists.perm <- function(dummy,testMat, goldStdrdMat) {
	#must modify this script when using not DREAM data, for example: in e.coli i,i interactions are allowed
	#because the regulators are a different set of genes than the tfs!
	goldStdrdMat=as.matrix(sample(goldStdrdMat))
	maxN = 10000000;
	tstRnkdLst = sort(testMat,decreasing=TRUE,index.return=TRUE)[[2]]
# convert vectorized index two double array index
	tstRnkdLst1 = tstRnkdLst %% nrow(goldStdrdMat)
	tstRnkdLst1[which(tstRnkdLst1 == 0)] = nrow(goldStdrdMat)
	tstRnkdLst2 = tstRnkdLst/nrow(goldStdrdMat)
	tstRnkdLst2[which(! tstRnkdLst2%%1==0)] = as.integer(tstRnkdLst2[which(! tstRnkdLst2%%1==0)]+1)
# remove diagonal enteries we are not predicting auto reg
	rm = which(tstRnkdLst1 == tstRnkdLst2)
	tstRnkdLst1 = tstRnkdLst1[-rm]
	tstRnkdLst2 = tstRnkdLst2[-rm]
	
	
	# Analysis
	# initialization
	
	k=0
	Ak=0
	TPk=0
	FPk=0
	P = sum(goldStdrdMat)
	N = length(goldStdrdMat)-P # Diagonal is not considered
	T = P+N
	L = length(tstRnkdLst1)                                               
	rec = numeric(length=L)
	prec = numeric(length=L)
	tpr = numeric(length=L)
	fpr = numeric(length=L)
	
	if(length(tstRnkdLst1) > 0) {
		while(k < L) {
			k = k + 1
			# if ((k %% 1000) == 0) {	cat(k, "out of", L,"\n") } ##huh? what does this line mean??
			# if k pred is correct
			if(goldStdrdMat[tstRnkdLst1[k],tstRnkdLst2[k]] == 1) {
				TPk = TPk + 1
				if(k==1)	
					{ delta=1/P }  
				else		
					{ delta=(1-FPk*log(k/(k-1)))/P }
				Ak = Ak + delta
			} else if (goldStdrdMat[tstRnkdLst1[k],tstRnkdLst2[k]] == 0) {
				FPk = FPk + 1
			} else {
				cat("could not find ", tstRnkdLst1[k],tstRnkdLst2[k]);
			}
			rec[k] = TPk/P
			prec[k] = TPk/k
			tpr[k] = rec[k]
			fpr[k] = FPk/N
		}
	}

	TPL=TPk
	if (L < T) 	{ 
		rh = (P-TPL)/(T-L) 
	} else { 
		rh = 0 
	}
	
	if (L>0) {
		recL = rec[L]
	} else {
		recL = 0
	}
	while (TPk < P) {
		k = k + 1
		TPk = TPk + 1
		rec[k] = TPk/P
		if ( ((rec[k]-recL)*P + L * rh) != 0 ) {
			prec[k] = rh * P * rec[k]/((rec[k]-recL)*P + L * rh)
		} else {
			prec[k] = 0
		}
		tpr[k] = rec[k]
		FPk = TPk * (1-prec[k])/prec[k]
		fpr[k] = FPk/N
	}
	AL = Ak;
	if (!is.nan(rh) & rh != 0  & L != 0) {
		AUC = AL + rh * (1-recL) + rh * (recL - L * rh / P) * log((L * rh + P * (1-recL) )/(L *rh))
	} else if (L==0) {
		AUC = P/T
	} else {
		AUC = Ak
	}
	
	# Integrate area under ROC
	lc = fpr[1] * tpr[1] /2
	for (n in 1:(L+P-TPL-1)) {
		lc = lc + (fpr[n+1]+fpr[n]) * (tpr[n+1]-tpr[n]) / 2
	}
	
	AUROC = 1 - lc

	# specific precision values
	TrueP = c(1, 2, 5, 20, 100, 500) 
	prec0 = numeric(length=length(TrueP))
	names(prec0) = as.character(TrueP)
	rec0 = numeric(length=length(TrueP))
	for (i in 1:length(TrueP)) {
		if(TrueP[i]<=P) {
			rec0[i]=TrueP[i]/P
			j=which(rec == rec0[i])
			j=min(j) #In case there is more than 1 precision values for rec(i)
			prec0[i]=prec[j]
		}
	}

	# handle output
	object = list()
	object[["prec"]] = prec
	object[["rec"]] = rec
	object[["tpr"]] = tpr
	object[["fpr"]] = fpr
	object[["AUROC"]] = AUROC
	object[["AUPR"]] = AUC
	object[["specificPrecVals"]] = prec0
	return(object)

}



# desciption: takes a matrix of TF-> target scores, and output a .sif and .eda files. 
# input:
#       - m: a matrix of scores for each TF(columns)->target(rows) interactions. Negative values correspond to repression edges, positive to activation.
#	- m.cut: a value that specify the interactions significance we will consider
#	- file.eda: a file name for the .eda file output (include path)
#	- file.sif: a file name for the .sif file output (include path)
#	- pos.edge: the annotation for the edge type for positive scores
#	- neg.edge: the annotation for the edge type for negative scores
#       - append: if append equals TRUE than keep appending to existing eda and sif files otherwise start a new file
# output:
#	- writes to file.eda and file.sif
write.cyto.net <- function(m,m.cut,file.eda="net.eda",file.sif="net.sif",pos.edge="positive",neg.edge="negative",append=FALSE){
  if(!append){
    cat( "confidence_score\n", file = file.eda)
	cat( "", file = file.sif)
  }
  tfs <- colnames(m)
  for(i in 1:length(tfs)){
    tf <- tfs[i]
    tf.targets.ix <- which( abs(m[,tf]) > m.cut )
    if(length(tf.targets.ix)>0){
      for(j in 1:length(tf.targets.ix)){
        trgt.gn <- names(tf.targets.ix)[j]
        if(m[trgt.gn,tf] > m.cut) {
          edge.type <- pos.edge
        } else if (m[trgt.gn,tf] < m.cut) {
          edge.type <- neg.edge
        } else {
          edge.type <- "neutral"
        }
        eda.line <- paste(tf, " (",edge.type,") ", trgt.gn, " = ", 
                          m[tf.targets.ix[j],tf], sep = "")
        eda.line <- paste(eda.line,"\n",sep="")
        sif.line <- paste( tf, edge.type, trgt.gn, sep = " ")
        sif.line <- paste(sif.line,"\n",sep="")
        cat( eda.line, file = file.eda, append = TRUE)
        cat( sif.line, file = file.sif, append = TRUE)
      }
    }
  }
}

# desciption: takes a matrix of TF-> target scores, and output a .sif and .eda files. 
# input:
#       - m: a matrix of scores for each TF(columns)->target(rows) interactions. Negative values correspond to repression edges, positive to activation.
#	- m.cut: a value that specify the interactions significance we will consider
#	- file.eda: a file name for the .eda file output (include path)
#	- file.sif: a file name for the .sif file output (include path)
#	- pos.edge: the annotation for the edge type for positive scores
#	- neg.edge: the annotation for the edge type for negative scores
#	- comb.data.vec: the annotation of what data type conbributed to the scores of each tf
#       - append: if append equals TRUE than keep appending to existing eda and sif files otherwise start a new file
# output:
#	- writes to file.eda and file.sif
write.cyto.net.extended <- function(m,m.cut,file.eda="net.eda",file.sif="net.sif",pos.edge="positive",neg.edge="negative",comb.data.vec,append=FALSE){
  if(!append){
    cat( "confidence_score\n", file = file.eda)
	cat( "", file = file.sif)
  }
  tfs <- colnames(m)
  for(i in 1:length(tfs)){
    tf <- tfs[i]
    tf.targets.ix <- which( abs(m[,tf]) > m.cut )
    if(length(tf.targets.ix)>0){
      for(j in 1:length(tf.targets.ix)){
        trgt.gn <- names(tf.targets.ix)[j]
        if(m[trgt.gn,tf] > m.cut) {
          edge.type <- paste(sep="",pos.edge,"_",comb.data.avilable[tf])
        } else if (m[trgt.gn,tf] < m.cut) {
          edge.type <- paste(sep="",neg.edge,"_",comb.data.avilable[tf])
        } else {
          edge.type <- paste(sep="","neutral","_",comb.data.avilable[tf])
        }
        eda.line <- paste(tf, " (",edge.type,") ", trgt.gn, " = ", 
                          m[tf.targets.ix[j],tf], sep = "")
        eda.line <- paste(eda.line,"\n",sep="")
        sif.line <- paste( tf, edge.type, trgt.gn, sep = " ")
        sif.line <- paste(sif.line,"\n",sep="")
        cat( eda.line, file = file.eda, append = TRUE)
        cat( sif.line, file = file.sif, append = TRUE)
      }
    }
  }
}


# desciption: takes a matrix of TF-> target scores, and output a .sif and .eda files. 
# input:
#       - m: a matrix of scores for each TF(columns)->target(rows) interactions. Negative values correspond to repression edges, positive to activation.
#	- m.cut: a value that specify the interactions significance we will consider
#	- file.txt: a file name for the .txt file output (include path)
#	- pos.edge: the annotation for the edge type for positive scores
#	- neg.edge: the annotation for the edge type for negative scores
#       - append: if append equals TRUE than keep appending to existing eda and sif files otherwise start a new file
# output:
#	- writes to file.eda and file.sif
write.mavisto.net <- function(m,m.cut,file.txt="net.txt",pos.edge="positive",neg.edge="negative",append=FALSE){
  if(!append){
	cat( "", file = file.txt)
  }
  tfs <- colnames(m)
  for(i in 1:length(tfs)){
    tf <- tfs[i]
    tf.targets.ix <- which( abs(m[,tf]) > m.cut )
    if(length(tf.targets.ix)>0){
      for(j in 1:length(tf.targets.ix)){
        trgt.gn <- names(tf.targets.ix)[j]
        if(m[trgt.gn,tf] > m.cut) {
          edge.type <- pos.edge
        } else if (m[trgt.gn,tf] < m.cut) {
          edge.type <- neg.edge
        } else {
          edge.type <- "neutral"
        }
        line <- paste("<",tf, "> <", trgt.gn, "> <", edge.type ,">", sep = "")
        line <- paste(line,"\n",sep="")
        cat( line, file = file.txt, append = TRUE)
      }
    }
  }
}

# desciption: takes a matrix of TF-> target scores, and output a .sif and .eda files. 
# input:
#       - m: a matrix of scores for each TF(columns)->target(rows) interactions. Negative values correspond to repression edges, positive to activation.
#	- m.cut: a value that specify the interactions significance we will consider
#	- file.txt: a file name for the .txt file output (include path)
#	- pos.edge: the annotation for the edge type for positive scores
#	- neg.edge: the annotation for the edge type for negative scores
#       - append: if append equals TRUE than keep appending to existing eda and sif files otherwise start a new file
# output:
#	- writes to file.eda and file.sif
write.fanmod.net <- function(m,m.cut,file.txt="net.txt",pos.edge=1,neg.edge=2,append=FALSE){
  if(!append){
	cat( "", file = file.txt)
  }
  tfs <- colnames(m)
  for(i in 1:length(tfs)){
	tf <- i
    tf.targets.ix <- which( abs(m[,tf]) > m.cut )
    if(length(tf.targets.ix)>0){
      for(j in 1:length(tf.targets.ix)){
		trgt.gn <- tf.targets.ix[j]
        if(m[trgt.gn,tf] > m.cut) {
          edge.type <- 1
        } else if (m[trgt.gn,tf] < m.cut) {
          edge.type <- 2
        } else {
          edge.type <- 0
        }
        line <- paste(tf," ",trgt.gn, " ", edge.type , sep = "")
        line <- paste(line,"\n",sep="")
        cat( line, file = file.txt, append = TRUE)
      }
    }
  }
}

# desciption: given a named vector of trgt.scores and a named vector of scores for the core th17 genes (such that s.trgts \in s.core),
#             get a running sum statistic vector p that describe how well the trgt scores vector matches (by rank and by score) the core scores vector
# input:
#       - l1: a named vector of the core scores to compare with (the gold standard)
#	- l2: a named vector of the trgts scores (the prediction)
#	- N: the number of permutation runs to perform in order to determine significance
# output:
#	- o: a list:
#	   - o$lb: a binned list of the run.sum.score achivable by best possible prediction
#	   - o$lw: a binned list of the run.sum.score achivable by worst possible prediction
#	   - o$lp: a binned list of the run.sum.score achived by actual prediction possible prediction
#	   - o$s: a binned list of the run.sum.stat achived by actual prediction possible prediction
#	   - o$sb: the best run sum stat achivable by prediction line (either enrichment, positive, or depletion, negative)
#	   - o$trgts: a list of targets that made it into the enrichment (a subset of total targets in l2)
#	   - o$z: a list of z scores for each bin
get.running.sum.statistic <- function(l1,l2,N,plot.it=TRUE,plot.name="TF"){
  # get best possible line for l2 predicting l1
  # 1) convert l1 distribution to be bin ranked the same as l2
  l1.mod <- l2
  names(l1.mod) <- names(l1)
  lb <- get.run.sum.score(l1,l1.mod)
  # get worst possible line
  lrev <- rev(l1.mod)
  names(lrev) <- names(l1)
  lw <- get.run.sum.score(l1,lrev )
  # get line for prediction as captured in l2
  lp <- get.run.sum.score(l1,l2 )
  # calculate distances for each bin
  d1 <- lb - lp
  d2 <- lp - lw
  s <- d2-d1
  # find best ix where tf.line is close to best line and far from worst line
  ix.min.max <- which.max(abs(s[-length(s)]))
  # keep best score
  sb <- s[ix.min.max]
  # score distribution in l2
  bins <- sort(unique(l2),decreasing=TRUE)
  # get a list of targets, up to ix.min.max (the targets that in l2 that are enriched in targets in l1)
  trgts <- names(l2)[which(l2>bins[ix.min.max])]
  # create N permutations of l2
  lprm.mat <- replicate(N,sample(l2))
  rownames(lprm.mat) <- names(l2)
  # get line for each permutation
  x <- apply(lprm.mat,2,get.run.sum.score, l1=l1)
  l.prm.mean <- apply(x,1,mean)
  l.prm.sd <- apply(x,1,sd)
  d1.mat <- lb-x
  d2.mat <- x-lw
  s.mat <- d2.mat-d1.mat
  # get permutation mean score for each bin
  s.mat.mean <- apply(s.mat,1,mean)
  # get permutation std for each bin
  s.mat.sd <- apply(s.mat,1,sd)
  # get zscore for each bin
  z <- (s-s.mat.mean)/s.mat.sd
  s.auc <- calc.auc(s)
  z.auc <- calc.auc(z)
  s.auc.bin3 <- calc.auc(s[1:3])
  z.auc.bin3 <- calc.auc(z[1:3])
  ix <- length(z)
  s[ix] <- z[ix] <- s.mat.mean[ix] <- s.mat.sd[ix] <- 0
  ## names(lb) <- names(lw) <- names(lp) <- names(s) <- names(trgts) <- names(z) <- names(s.mat.mean) <-
  ##                           names(s.mat.sd) <- names(d1) <- names(d2) <- names(l.prm.mean) <- names(l.prm.sd) <- bins
  names(lb) <- names(lw) <- names(lp) <- names(s) <- names(z) <- names(s.mat.mean) <-
                            names(s.mat.sd) <- names(d1) <- names(d2) <- names(l.prm.mean) <- names(l.prm.sd) <- bins  
  o <- list(lb=lb,lw=lw,lp=lp,s=s,sb=sb,trgts=trgts,z=z,s.perm.mean=s.mat.mean,s.perm.sd=s.mat.sd,d1=d1,d2=d2,
            l.prm.mean=l.prm.mean,l.prm.sd=l.prm.sd,bins=bins,s.auc=s.auc,z.auc=z.auc,s.auc.bin3=s.auc.bin3,z.auc.bin3=z.auc.bin3)
  if(plot.it==TRUE){
    make.qc.plot.enrichment(o,plot.name=plot.name)
  }
  return(o)
}

# desciption: a helper function for analysing (quality control) the results from get.running.sum.statistic
# input:
#       - lb: a named vector of the core scores (gold standard list, representing the best we can do)
#	- lw: a named vector of the worst scores (predicting the reverese of lb)
#	- lp: the actual predicted scored list
#	- lp.mat: a matrix of permutated predictions
# output:
#	- plot1: the best line, worst line, and predicted line in between
#	- plot2: difference plot showing how close/far the predicted line is from best and worst lines
make.qc.plot.enrichment <- function(o,plot.name="TF"){
  ylim <- c(0,1)
  xlim <- c(min(o$bins),max(o$bins))
  ix <- which.max(abs(o$s[-length(o$s)]))
  cex <- 1.5
  pch <- 20
  type <- "l"
  # plot best line
  plot(o$lb,xaxt="n",ylim=ylim,xlab="bins",ylab="portion of explianed scores",main=plot.name,cex=cex,type=type)
  axis(1,at=min(xlim):max(xlim)+1,labels=max(xlim):min(xlim))
  # plot worst line
  lines(o$lw,cex=cex,type=type)
  # plot predicted line
  lines(o$lp,cex=cex,type=type)
  # plot difference lines between best worst and predicted lines
  segments(x0=ix, y0=o$lp[ix], x1=ix, y1=o$lb[ix],col="red")
  segments(x0=ix, y0=o$lw[ix], x1=ix, y1=o$lp[ix],col="green")
  # put permutation line mean and +/-sd in plut
  m <- o$l.prm.mean
  std <- o$l.prm.sd
  segments(x0=(1:length(m)), y0=(m-std), x1=(1:length(m)), y1=(m+std),col="blue",cex=1.5)
  text(x=ix, y=mean(c(o$lb[ix],o$lp[ix])),round(o$d1[ix],3))
  text(x=ix, y=mean(c(o$lw[ix],o$lp[ix])),round(o$d2[ix],3))
  # plot run sum stat over bins
  plot(o$s,xaxt="n",xlab="bins",ylab="running sum statistic",main=plot.name,cex=cex,type=type)
  axis(1,at=min(xlim):max(xlim)+1,labels=max(xlim):min(xlim))
  segments(x0=ix,y0=0,x1=ix,y1=o$s[ix],col="red",lty=3,cex=1.5)
  text(x=ix, y=mean(c(0,o$s[ix])),paste(sep="","AUC = ",round(o$s.auc,3)))
  # plot zscores
  plot(o$z,xaxt="n",xlab="bins",ylab="zscore",main=plot.name,cex=cex,type=type)
  axis(1,at=min(xlim):max(xlim)+1,labels=max(xlim):min(xlim))
  segments(x0=ix,y0=0,x1=ix,y1=o$z[ix],col="red",lty=3,cex=1.5)
  text(x=ix, y=mean(c(0,o$z[ix])),paste(sep="","AUC = ",round(o$z.auc,3)))
}

# desciption: given a named vector of trgt.scores and a named vector of scores for the core th17 genes (such that s.trgts \in s.core),
#             get a running sum statistic vector p that describe how well the trgt scores vector matches (by rank and by score) the core scores vector
# input:
#       - l1: a named vector of the core scores to compare with
#	- l2: a named vector of the trgts scores
# output:
#	- p: a vector with increasing values p1<p2<...<pN
get.run.sum.score <- function(l1,l2){
  # get normalization factor
  norm.const <- sum(l1)
  # get sum statistic vector for tf
  levels <- sort(unique(l2),decreasing=T)
  p <- numeric(length(levels))
  for(j in 1:length(levels)){
    trgts.lvl.j <- names(l2)[which(l2==levels[j])]
    if(j==1) {
      p[j] <- (sum(l1[trgts.lvl.j]))/norm.const
    } else {
      p[j] <- p[j-1] + (sum(l1[trgts.lvl.j]))/norm.const
    }   
  }
  return(p)
}

# desciption: given a vector with increasing values starting with 0 and ending in 1, calculate the area under curve.
# input:
#	- p: a vector with increasing values p1<p2<...<pN
# output:
#	- auc: area under curve given by p
calc.auc <- function(p) {
  # define d.binned.rank.scores
  bin.length <- 1/length(p)
  bin.area <- p*bin.length
  auc <- sum(bin.area)
  return(auc)
}

# desciption: Uses score ranks (the higher the score the more confidence we have in an interaction) in order to transform confidence scores into quantile scores that range from 0 to 1.
#             1 is used for the highest confidence interactions, 0 for lowest cofidence, 0 for scores that had no confidence to begin with (i.e. 0 in beta matrix).
# input:
#	- d: a matrix with scores for each tf->target
# output:
#	- d.r.r : a matrix holding relative rank scores for each interaction.
convert.scores.to.relative.ranks.pos <- function(d) {
  if(is.matrix(d)==TRUE){
    ## define d.q: a matrix for the quantiles
    d.r.r <- matrix(0,nr=dim(d)[1],nc=dim(d)[2])
    dimnames(d.r.r) <- dimnames(d)
  } else {
    d.r.r <- numeric(length=length(d))
    names(d.r.r) <- names(d)
  }
  # find non-zero scores
  ix.non.zero <- which(d>0)
  # get scores as relative ranks
  ## l.total <- length(d)
  l.non.zero <- length(ix.non.zero)
  ## l.remainder <- l.total-l.non.zero
  score.remainder <- 0
  ix.non.zero.ranked.by.score <- sort(d[ix.non.zero],decreasing=T,index.return=T)$ix
  d.r.r[ix.non.zero[ ix.non.zero.ranked.by.score]] <- seq(1,0,length.out=length(ix.non.zero))
  return(d.r.r)
}

# desciption: Uses score ranks (the higher the score the more confidence we have in an interaction) in order to transform confidence scores into quantile scores that range from 0 to 1.
#             1 is used for the highest confidence interactions, 0 for lowest cofidence, 0 for scores that had no confidence to begin with (i.e. 0 in beta matrix).
# input:
#	- d: a matrix with scores for each tf->target
# output:
#	- d.binned.rank.scores : a matrix holding binned rank scores for each interaction.
convert.scores.to.relative.ranks <- function(d) {
  if(is.matrix(d)==TRUE){
    ## define d.q: a matrix for the quantiles
    d.r.r <- matrix(0,nr=dim(d)[1],nc=dim(d)[2])
    dimnames(d.r.r) <- dimnames(d)
  } else {
    d.r.r <- numeric(length=length(d))
    names(d.r.r) <- names(d)
  }
  # find non-zero scores
  ix.non.zero <- which(d!=0)
  # get scores as relative ranks
  ## l.total <- length(d)
  l.non.zero <- length(ix.non.zero)
  ## l.remainder <- l.total-l.non.zero
  score.remainder <- 0
  ix.non.zero.ranked.by.score <- sort(d[ix.non.zero],decreasing=T,index.return=T)$ix
  d.r.r[ix.non.zero[ ix.non.zero.ranked.by.score]] <- seq(1,0,length.out=length(ix.non.zero))
  return(d.r.r)
}

# desciption: Uses score ranks in order to transform confidence scores into integer scores that range from 0 to 10.
#             10 is used for the highest confidence interactions, 1 for lowest cofidence, 0 for scores that had no confidence to begin with (i.e. 0 in beta matrix).
# input:
#	- d: a matrix with scores for each tf->target
# output:
#	- d.binned.rank.scores : a matrix holding binned rank scores for each interaction.
convert.scores.to.rank.binned.scores <- function(d) {
  # define d.binned.rank.scores
  d.binned.rank.scores <- matrix(0,nr=dim(d)[1],nc=dim(d)[2])
  dimnames(d.binned.rank.scores) <- dimnames(d)
  # find non-zero scores
  ix.non.zero <- which(d!=0)
  # get scores as quantiles
  qunts <- seq(0,1,by=.1)
  scores.quantiles <- quantile(d[ix.non.zero],qunts)
  d.binned.rank.scores[ix.non.zero] <- assign.level.rank.increasing(d[ix.non.zero],scores.quantiles)
  return(d.binned.rank.scores)
}

# desciption: transforms confidence scores stored in 'v' into integer scores that range from 1 to length('r'-1),
#             where r is a vector with several quantile values for the vector 'v'.
#             length(r) is used for the highest confidence interactions, 1 for lowest cofidence.
# input:
#	- v: a vector with scores
# output:
#	- o : a vector of bin ranked scores
# Note: do not give this function zero confidence scores, only non zero scores.
assign.level.rank.increasing <- function(v,r) {
	# validate that all values in v are smaller then largest value in r (other wise some values will not be ranked!)
	if(max(v) > max(r)) { stop("#####\ncan't assign ranks to v as max(v)>max(r). Bailing out...\n#####\n")}
	if(min(v) < min(r)) { stop("#####\ncan't assign ranks to v as min(v)<min(r). Bailing out...\n#####\n")}
	o <- numeric(length=length(v))
	cur.rank <- 1
	for (i in 2:length(r)){
		# find indices to entries in v that belong to current rank
		if(i==2){
			### first bin is[ r[1],r[2] ]###			
			cur.ix <- which( (v <= r[i]) & (v>=r[i-1]))	
		} else {
			### other bins are ( r[i-1],r[i] ], note the paren on the left hand side ###			
			cur.ix <- which( (v <= r[i]) & (v>r[i-1]))	
		}
		o[cur.ix] <- cur.rank
		cur.rank <- cur.rank + 1
	}
	# sanity check that all values in v have been assigned a value between 1-10
	if( length(which( o < 1 | o > (length(r)-1) )) > 0 ) { stop("#####\nsome values in o have not been assigned a correct level. Bailing out...\n#####\n")}
	return(o)
}

# desciption: some genes may have more than one observation (row) take the row that has highest median value across the data (proxy for highest signal).
# input:
#	- d: a dataset where rownames may not be unique
# output:
#	- obj: a list with two objects:
#          -  obj$d: a smaller dataset where each row is unique
#          -  obj$numRepeatGenes: total number of repeats (e.g if gn x apeared 2 times and gene y apeared 3 times obj$numRepeatGenes=5)
take.median.gene.of.non.unique.genes <- function(d) {
	# average non-unique genes
	x <- sort(table(rownames(d)),decreasing=T)
	rpts.gns <- names(x[which(x>1)])
	if(length(rpts.gns) > 0){
		rpts.gns.avrg <- matrix(0,nrow=length(rpts.gns),ncol=ncol(d))
		rownames(rpts.gns.avrg) <- rpts.gns
		colnames(rpts.gns.avrg) <- colnames(d)
		# calc avrg of repeated genes from d
		for(i in 1:length(rpts.gns)) {
			rpts.ix <- which(rownames(d)==rpts.gns[i])
                        median.rpts <- apply(as.matrix(d[rpts.ix,]),1,median)
                        ix.max.median <- which.max(median.rpts)
                        rpts.gns.avrg[i,] <- d[rpts.ix[ ix.max.median ],]
		}
		# remove repeated genes from d
		ix.rm.gns <- numeric()
		for(i in 1:length(rpts.gns)) {
			ix.rm.gns <- c(ix.rm.gns,which(rownames(d)==rpts.gns[i]))
		}	
		d <- as.matrix(d[-ix.rm.gns,])
		# add repeated genes to the begining of d
		d <- rbind(rpts.gns.avrg,d)
		# return a list with uniqe data set d and number of repeated genes
		obj <- list()
		obj[["d"]] = d
		obj[["numRepeatGenes"]] = length(ix.rm.gns)
		return(obj)
	} else {
          obj <- list()
          obj[["d"]] = d
          obj[["numRepeatGenes"]] = 0       
          return(obj)
        }
}


# desciption: some genes may have more than one observation (row) average over those rows.
# input:
#	- d: a dataset where rownames may not be unique
# output:
#	- obj: a list with two objects:
#          -  obj$d: a smaller dataset where each row is unique
#          -  obj$numRepeatGenes: total number of repeats (e.g if gn x apeared 2 times and gene y apeared 3 times obj$numRepeatGenes=5)
average.non.unique.genes <- function(d) {
	# average non-unique genes
	x <- sort(table(rownames(d)),decreasing=T)
	rpts.gns <- names(x[which(x>1)])
	if(length(rpts.gns) > 0){
		rpts.gns.avrg <- matrix(0,nrow=length(rpts.gns),ncol=ncol(d))
		rownames(rpts.gns.avrg) <- rpts.gns
		colnames(rpts.gns.avrg) <- colnames(d)
		# calc avrg of repeated genes from d
		for(i in 1:length(rpts.gns)) {
			rpts.ix <- which(rownames(d)==rpts.gns[i])
			rpts.gns.avrg[i,] <- apply(as.matrix(d[rpts.ix,]),2,mean)
		}
		# remove repeated genes from d
		ix.rm.gns <- numeric()
		for(i in 1:length(rpts.gns)) {
			ix.rm.gns <- c(ix.rm.gns,which(rownames(d)==rpts.gns[i]))
		}	
		d <- as.matrix(d[-ix.rm.gns,])
		# add repeated genes to the begining of d
		d <- rbind(rpts.gns.avrg,d)
		# return a list with uniqe data set d and number of repeated genes
		obj <- list()
		obj[["d"]] = d
		obj[["numRepeatGenes"]] = length(ix.rm.gns)
		return(obj)
	} else {
          obj <- list()
          obj[["d"]] = d
          obj[["numRepeatGenes"]] = 0       
          return(obj)
        }
}


# desciption: Average dataset over repeats.  Requires naming of columns to be "expt.name" followed by "#" followed by repeat number.
#             e.g.  [1] "GSM399350_EA07068_52789_MoGene_preT.DN1.Th_.1.cel" [2] "GSM399351_EA07068_52790_MoGene_preT.DN1.Th_.2.cel"
#             [3] "GSM399352_EA07068_52791_MoGene_preT.DN1.Th_.3.cel" will create one experiment 
# input:
#	- d: a dataset with proper column names
# output:
#	- d.median: a smaller dataset where each experiment has now only one median observation
avg.over.reps.immgen <- function(d) {
  expt.nms <- colnames(d)
  # split expt names by _
  x=strsplit(expt.nms,"_")
  # rm GSM EA and other unique identifier numbers
  expt.matrix <- t(sapply(x,function(i) i[ (length(i)-1):length(i) ]))
  # replace . by _ in names
  expt.matrix <- gsub("\\.","_",expt.matrix)
  # fix names erroniously spliced names with MoGene where expt name should be (inconsistant naming scheme by immgen...)
  for(i in 1:dim(expt.matrix)[1]){
    if(expt.matrix[i,1]=="MoGene"){
      tmp = strsplit(expt.matrix[i,2],"_")[[1]]
      expt.matrix[i,1] <- paste(tmp[1:(length(tmp)-2)],collapse="_")
      expt.matrix[i,2] <- paste(tmp[ (length(tmp)-1): length(tmp)],collapse="_")
    }
  }
  # almost there ahhh...
  for(i in 1:dim(expt.matrix)[1]){
    tmp <- strsplit(expt.matrix[i,2],"_")[[1]]
    if(length(strsplit(tmp[1],"")[[1]])>1){
      expt.matrix[i,1] <- paste(expt.matrix[i,1],tmp[1],sep="_")
      expt.matrix[i,2] <- tmp[2]
    }
  }
  unique.expt <- unique(expt.matrix[,1])
  d.median <- matrix(0,nrow=dim(d)[1],ncol=length(unique.expt))
  colnames(d.median) <- unique.expt
  rownames(d.median) <- rownames(d)
  for(i in 1:length(unique.expt)){
    ix <- which(expt.matrix[,1]==unique.expt[i])
    if(length(ix)>1){
      d.median[,i] <- apply(d[,ix],1,median)
    } else {
      d.median[,i] <- d[,ix]
    }
  }
  return(d.median)
}


# desciption: calculate differential expression of any genes between two RNAseq conditions that have only one replicate each
#             calculation is done on two levels: 1) fold change and 2) absolute difference in rpkms.
#             each level is converted into zscores: 1) zscore fc and 2) zscore abs diff (based on fc and abs diff distributions repectively within each experiment)
#             zscores from fc and abs diff are combined using the stouffer method.
# input:
#	- wt: ranseq results (e.g. as rpkm) for wt (as a named vector with rpkm values and names of corresponding genes)
#	- ko: ranseq results (e.g. as rpkm) for ko (as a named vector with rpkm values and names of corresponding genes)
#	- c1: filter out genes with low rpkms in !!both!! experiments (defined by sum over the two experiments)
#	- ps.cnt: how many pseudo counts to add to wt and ko rpkms (this does not change rpkm.diff but does correct for small rpkm values with big fold change)
# output:
#	- stoufer corrected zscore of diff.rpkm.zscore and diff.fc.zscore for each gene
diff.exprs.analysis <- function(wt, ko, c1 = 5, ps.cnt = 1) {
  ix.good <- which((wt+ko) > c1)
  wt <- wt[ix.good]
  ko <- ko[ix.good]
  # calc diff in rpkm
  diff.rpkm <- wt - ko
  diff.fc <- log2((wt+ps.cnt)/(ko+ps.cnt))
  diff.zscore <- (diff.rpkm-median(diff.rpkm))/sd(diff.rpkm)
  fc.zscore <- (diff.fc-median(diff.fc))/sd(diff.fc)
  stoufer.zscore <- (diff.zscore+fc.zscore)/sqrt(2)
  return(stoufer.zscore)
}

# desciption: takes an n*m*b array of zscores and return an n*m matrix of median zscores (where median is taken over the b dimension of the array)
# input:
#	- a vector of zscores containing multiple entries for the same gene (by gene name)
# output:
#	- a vector of median zscores with unique entries for each gene

get.median.zscore <- function(zscores) {
  gn.nms <- unique(names(zscores))
  median.zscores <- numeric(length=length(gn.nms))
  names(median.zscores) <- gn.nms
  for(g in 1:length(gn.nms)) {
    median.zscores[g] <- median(zscores[ which(names(zscores)==gn.nms[g])])
  }
  return(median.zscores)
}

# desciption: create the inferelator data structure colMap (gives info about the different experiments in the dataset, e.g. timeseries vs. steadystate)
# input:
#	- ratios: the data matrix (rnaseq or microarray) m*n (m genes n conditions)
#	- is.ts: a vector with dimension n that for each experiments describe if it is part of time series or equalibrium
#	- del.t: a vector with dimension n that for each experiments describe what is the time difference between it and previous time point (if equalibrium time difference is set to zero)
# output:
#	- colMap (a list with each element being an experiment, e.g. one microarray experiment)

createColMap <- function(ratios,is.ts,del.t) {
	object <- list()
	# default is equilibrium
	object$isTs <- FALSE 
	object$is1stLast <- factor("e",levels=c("e","f","m","l"))
	object$prevCol <- as.logical("NA")
	object$del.t <- as.logical("NA")
	object$condName <- character()
	colMap <- list()
	expName <- colnames(ratios)

	for(i in 1:dim(ratios)[2]){
		colMap[[i]] <- object
                colMap[[i]]$isTs <- is.ts[i]
                colMap[[i]]$del.t <- del.t[i]
                if(is.ts[i] == TRUE) {
                  if(del.t[i] == 0) {
                    colMap[[i]]$is1stLast <- factor("f",levels=c("e","f","m","l"))
                  } else if(i == length(del.t) | (del.t[i+1]-del.t[i])<0) {
                    colMap[[i]]$is1stLast <- factor("l",levels=c("e","f","m","l"))
                  } else {
                    colMap[[i]]$is1stLast <- factor("m",levels=c("e","f","m","l"))
                  }
                }
                if(is.ts[i] == TRUE & colMap[[i]]$is1stLast != factor("f",levels=c("e","f","m","l"))) {
                  colMap[[i]]$prevCol <- colMap[[i-1]]$condName
                }
                colMap[[i]]$condName <- expName[i]         
	}
	return(colMap)
}


# create clusterStack (each cluster is a gene)
# load("input/nanoMed/ratios.RData")
createClusterStack <- function(ratios) {
	object = list()
	object$cols = colnames(ratios)
	object$ncols = length(colnames(ratios))
	object$rows = character()
	object$nrows = 1
	object$resid = as.logical("NA")
	object$k = integer()
	object$redExp = numeric()
	
	clusterStack = list()

	for (i in 1:dim(ratios)[1]){
		clusterStack[[i]] = object
		clusterStack[[i]]$rows = rownames(ratios)[i]
		clusterStack[[i]]$k = i
		clusterStack[[i]]$redExp = ratios[i,,drop=F] 
	}
	clusterStack[["K"]] = dim(ratios)[1]
        return(clusterStack)
#	save(clusterStack, file=paste(path,file,sep=""))
}

# Combine any number of four scores matrices (all have to have row.keepers in row names and col.keepers in colnames)
# note: if you want to combine less than four matrices, just give 0 as an input value for the matrix you don't want to combine.
# i.e. to get results only for k and c, knockout=rel.rank.ko,chip=rel.rank.chip,rnaseq=0,immgen=0
# input:(k-knockout scores, c-chip scores,r-rnaseq scores, and i-immgen scores, keepers is the row names that are to be kept)
# output: score combined matrix (n*m): s.c.m

combine_kcri <- function(knockout,chip,rnaseq,immgen,row.keepers,col.keepers){
                                        # define score.combined.matrix
  s.c.m <- matrix(0,nr=length(row.keepers),nc=length(col.keepers))
  rownames(s.c.m) <- row.keepers
  colnames(s.c.m) <- col.keepers

# fill score.combined.matrix by summing over the four int.score matrices
  for(i in 1:dim(s.c.m)[1]){
    # get scores for this gene
    gn <- rownames(s.c.m)[i]
    # initialize gns socres with all possible tfs to zero
    s.line <- numeric(length=dim(s.c.m)[2])
    names(s.line) <- colnames(s.c.m)
    i.l <- length(immgen);r.l <- length(rnaseq);k.l <- length(knockout);c.l <- length(chip)
    # for efficiency determine if each score matrix has non-zero values
    # add immgen scores
    if( (i.l>1) & (gn %in% rownames(immgen)) ){
      tfs <- colnames(immgen)
      s.line[tfs] <- s.line[tfs] + immgen[gn,tfs]
    }
    # add rnaseq scores
    if( (r.l>1) & (gn %in% rownames(rnaseq)) ){
      tfs <- colnames(rnaseq)
      s.line[tfs] <- s.line[tfs] + rnaseq[gn,tfs]
    }
    # add ko scores
    if( (k.l>1) & (gn %in% rownames(knockout)) ){
      tfs <- colnames(knockout)
      s.line[tfs] <- s.line[tfs] + knockout[gn,tfs]
    }
    # add chip scores
    if( (c.l>1) & (gn %in% rownames(chip)) ){
      tfs <- colnames(chip)
      s.line[tfs] <- s.line[tfs] + chip[gn,tfs]
    }
    s.c.m[i,] <- s.line
  }
  return(s.c.m)
}


# Combine two matrices with entries either >0 or 0
# map higest values in tmplt.mat (M1) with sml.mat (M2) that you want to combine to tmplt_mat

combine_mtrcs <- function(M1,M2){
	ix_m1 = sort(M1,decreasing = TRUE,index.return = TRUE)$ix
	ix_m2 = sort(M2,decreasing = TRUE,index.return = TRUE)$ix
	M = M1
	i=1
	while(M2[ix_m2[i]] > 0){
		M[ix_m2[i]] = sqrt(M[ix_m2[i]]^2+M[ix_m1[i]]^2)
		i = i+1
	}
	return(M)
}

# FisherEnrichedTFsToTargetGenes finds which tf's are enriched in reg inters with a target gene list compared to all other genes
# based on the fisher exact test.
# input:
# -m.tf is a matrix N*M matrix where N is the number of putative target genes and M is the num of tfs
#   each entry m_{i,j} holds the significance value of the interaction (confidence scores)
# - trgt.gn.lst is the target genes you want to find (enriched) tfs for
# output:
# - object containing 4 elements:
# 1. contngncy.matrices for each tf:
#                                       classifier.genes | non.classifier.genes
#        tf.regulated.genes                 a_{1,1}       |     a_{1,2}
#        non.tf.regulated.genes          a_{2,1}       |     a_{2,2}
# 2. fisher.summry the fisher test summary for each tf
# 3. p.vals.sorted is a sorted list of pvals for each tf
# 4. p.vals.ix.map is an index mapping btwn sorted p.vals list and the order of tfs in the original fisher.summry and cont.matrices lists

FisherEnrichedTFsToTargetGenes <- function(m.tf, trgt.gn.lst) {
  # get list of target genes regulated by each tf
  tfs.trgts.list <- apply(abs(m.tf),2,function(i) which(i>0))
  cg <- trgt.gn.lst
  non.cg <- rownames(m.tf)[which(!rownames(m.tf) %in% cg)]
  m.l <- list()
  m.l.f <- list()
  # here i keep correlations between tf and trgts over mRNA data (m.tf)
  cor.vec <- numeric(length=length(tfs.trgts.list))
  names(cor.vec) <- names(tfs.trgts.list)
  for (j in 1:length(tfs.trgts.list)) {
    tf <- names(tfs.trgts.list)[j]
    m <- matrix(0,nr=2,nc=2)
    rownames(m) <- c("tf.regulated.genes","non.tf.regulated.genes")
    colnames(m) <- c("classifier.genes","non.classifier.genes")
    m[1,1] <- length(which(cg %in% names(tfs.trgts.list[[j]])))
    m[1,2] <- length(which(non.cg %in% names(tfs.trgts.list[[j]])))
    m[2,1] <- length(cg) - m[1,1]
    m[2,2] <- length(non.cg) - m[1,2]
    m.l[[tf]] <- m
    m.l.f[[tf]] <- fisher.test(m,alternative = "greater")
  }
  p.vals <- sapply(m.l.f,function(i) i$p.value)
  p.vals.ix.map <- sort(p.vals,decreasing=F,index.return=T)$ix
  names(p.vals) <- names(p.vals.ix.map) <- names(tfs.trgts.list)
  
  obj <- list()
  obj[["contngncy.matrices"]] <- m.l
  obj[["fisher.summry"]] <- m.l.f
  obj[["p.vals"]] <- p.vals
  obj[["p.vals.ix.map"]] <- p.vals.ix.map
  return(obj)
}

# getCorTfWithTrgts finds the correlation vector between the tf and it's targets across d
# input:
# - d: a data matrix (typically mRNA data)
# - tf: the name of the tf of interest
# - m: a matrix of scores for all tfs with all possible trgts (only trgts with score > 0 are considered for each tf)
# - trgt.core: the names of the core trgt genes
# output:
# - cor.vec: a vector of correlation values between the tf and each of its trgts over d
getCorTfWithTrgts <- function(tf,d,m,trgt.core){
  trgts <- names(m[,tf][which(m[,tf]!=0)])
  trgts <- trgts[which(trgts %in% trgt.core)]
  if(length(trgts)>0){
    cor.vec <- cor(t(d[c(tf,trgts),]))[1,-1]
    cor.vec[which(cor.vec==1)] <- 0
  } else{
    cor.vec <- NA
  }
  return(cor.vec)
}
