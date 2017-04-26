##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

path.output.in.input <- paste(sep="","input/th17/used_for_paper/infResults/",GLOBAL[["type"]],"/")
# path.output.in.results <- paste(sep="","results/inf/",GLOBAL[["type"]],"/")
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 1- reads params, design and response matrices, found in PARAMS and INPUT list respectively
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

    source("./r_scripts/th17/used_for_paper/initRunNetInference.R")
	#making directory to save everything into
	system(paste("mkdir",PARAMS$general$saveToDir,sep=" ")) 
	#params for main (object PARAMS defined in init.R)
	b = 1 # this follow current iteration/bootstrap number
	N_b = PARAMS$"general"$"numBoots" # number of bootstraps
	btch_size = 10 # calculate this mumber of genes to all predictors MI scores in batches (to avoid running out of memory when calculating MI)
	percentCoverage <- PARAMS[["general"]][["percentCoverage"]] # (usually 100) percent of matrix that we want to resample
	lambda = PARAMS[["lars"]][["lambda"]] # set of l2 norm regularization weights to try in elastic net
	cleanUp <- FALSE # clear functions and other intermediate variables at end of run (leaves important variables more visible for end users)
	# response and design matrices for clr
	Y_clr = INPUT[["clr"]][["response_matrix"]]
        ###########
        # bug fix with naming for th17
        ## colnames(Y_clr) <- c(colnames(INPUT$general$dataset)[c(1,dim(INPUT$general$dataset)[2])],colnames(Y_clr)[-c(1,2)])
        ###########
	X_clr = INPUT[["clr"]][["design_matrix"]] # single predictors
	# response and design matrices for lars
	Y_lars = INPUT[["lars"]][["response_matrix"]]
        ###########
        # bug fix with naming for th17
        ## colnames(Y_lars) <- c(colnames(INPUT$general$dataset)[c(1,dim(INPUT$general$dataset)[2])],colnames(Y_lars)[-c(1,2)])
        ###########
	X_lars = INPUT[["lars"]][["design_matrix"]] # single predictors
	# store results (ODEs,Z scores, and error for each model for each bootstrap run, respectively)
	betaList = vector("list", N_b)
	modelErrorList = vector("list", N_b)
	#startTime <- date() #times how long a run takes
	allResults <- list() #list for storing all models

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 2- calculate Median corrected Zscores based on KO data
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 3- setup for bootstrap: create Pi-perm_vector/matrix, Y^pi,X^pi
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

while (b <= N_b) {
	#create permutation matrix
	cat("bootstrap #: ",b,"\n")
	if(b == 1){
		#here we want the original permutation, ie. getOrigPerm = TRUE (i.e. first bootstrap is exact dataset, no resampling)
		Pi_s_clr=createPermMatrix(cS=INPUT[["general"]][["clusterStack"]], allConds = colnames(Y_clr), getOrigPerm = TRUE, percentCoverage = percentCoverage)
		Pi_s_lars=Pi_s_clr
	} else {
		Pi_s_clr=createPermMatrix(cS=INPUT[["general"]][["clusterStack"]], allConds = colnames(Y_clr), getOrigPerm = FALSE, percentCoverage = percentCoverage)
		Pi_s_lars=Pi_s_clr
	}
	#create bicluster specific permutation matrix (ie. read from Pi_g, algorithm described in method comments)
	#this should be changed to be general for both cases where we have only single genes and cases where we havee biclusters
	Y_clr_p = permuteCols(Y_clr,Pi_s_clr)
	X_clr_p = permuteCols(X_clr,Pi_s_clr)

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 4- pass one: fill M - mutual information matrix or correlation matrix
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

	# dynamic MI scores stored here
	cat("calculating dynamic MI ")
	if(PARAMS[["general"]][["processorsNumber"]] > 1){
		Ms = calc_MI_one_by_one_parallel( Y_clr_p, X_clr_p, processorsNumber = PARAMS[["general"]][["processorsNumber"]], n.bins=PARAMS[["clr"]][["n.bins"]])
		# Ms = calc_MI_one_by_one_parallel( Y_clr, X_clr,Pi_s_clr, processorsNumber = PARAMS[["general"]][["processorsNumber"]], n.bins=PARAMS[["clr"]][["n.bins"]])
	} else {
		Ms = calc_MI_inBatces(Y_clr_p,X_clr_p,btch_size,n.bins=10)
	}
	diag(Ms) = 0
	cat("\n")
	# static MI scores stored here
	cat("calculating background MI ")
	if(PARAMS[["general"]][["processorsNumber"]] > 1){
		# Ms_bg = calc_MI_one_by_one_parallel( X_clr, X_clr,Pi_s_clr, processorsNumber = PARAMS[["general"]][["processorsNumber"]], n.bins=PARAMS[["clr"]][["n.bins"]])
		Ms_bg = calc_MI_one_by_one_parallel( X_clr_p, X_clr_p, processorsNumber = PARAMS[["general"]][["processorsNumber"]], n.bins=PARAMS[["clr"]][["n.bins"]])
	} else {	
		Ms_bg = calc_MI_inBatces(X_clr_p,X_clr_p,btch_size)
	}
	diag(Ms_bg) = 0
	cat("\n")

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 5- calculate mixed-CLR (or clr) matrix
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

	if(PARAMS[["general"]][["use_mixCLR"]]){
		cat("running mix-CLR ")
        Z_nt_fltrd = mixed_clr_parallel(Ms_bg,Ms,processorsNumber=PARAMS[["general"]][["processorsNumber"]])
	} else {
		cat("running CLR ")
		Z_nt_fltrd = clr(Ms)
	}
	cat("\n")
	colnames(Z_nt_fltrd) <- rownames(Z_nt_fltrd) <- rownames(X_clr)
	Z <- Z_nt_fltrd[,INPUT[["general"]][["tf_names"]]]

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 6- apply MCZ filter -- i.e. remove unlikely reg inters from further consideration by mixedCLR (and thus from Inf)
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
        if(PARAMS[["general"]][["use.priors"]]==TRUE){
	  # filter cutoff
          z.cut = PARAMS[["general"]][["rnaseq_ko_zscore_abs_filter"]]
          x <- INPUT$general$tf_knockOutConfList
          tfs <- names(x)
          for(i in 1:length(x)) {
            gns.rm <- toupper(names(which(abs(x[[ tfs[i] ]]) < z.cut)))
            gns.rm <- gns.rm[which(gns.rm %in% rownames(Z))]
            Z[gns.rm,tfs[i]] <- 0
          }
        }
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 7- run Inferelator (elastic net with ODE based modifications to response and design matrices) 
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
        cat("running elastic net ")
     	if(PARAMS[["general"]][["processorsNumber"]] > 1){
		x = calc_ode_model_weights_parallel(Xs = X_lars,Y = Y_lars, Pi = Pi_s_lars, M1 = Z, nS = PARAMS[["lars"]][["max_single_preds"]], 
								 lambda=lambda, processorsNumber = PARAMS[["general"]][["processorsNumber"]], plot.it = FALSE,
								 plot.file.name = paste(PARAMS$general$saveToDir,"/boot_",b,"_models.pdf",sep=""),verbose = FALSE)
	} else {
		x = calc_ode_model_weights(Xs = X_lars,Y = Y_lars, Pi = Pi_s_lars, M1 = Z, nS = PARAMS[["lars"]][["max_single_preds"]], lambda=lambda, 
								 plot.it = FALSE,plot.file.name = paste(PARAMS$general$saveToDir,"/boot_",b,"_models.pdf",sep=""),verbose = FALSE)
	}
	cat("\n")
	betaList[[b]] = x[[1]]
        ##AM fix model error
	## modelErrorList[[b]] = t(x[[5]])
        modelErrorList[[b]] = t(x[[5]])
        names(modelErrorList[[b]]) = rownames(X_clr)
	betaList[[b]]=add_weight_beta(bL=betaList[[b]],model_errors=modelErrorList[[1]],n=nrow(Y_lars),pS=nrow(X_lars),pD=0,col=4,col_name = "prd_xpln_var" )
	betaList[[b]]=add_zscore(bL=betaList[[b]],M1=Z,M2=NULL,col=5,col_name = "clr_zs")
	betaList[[b]]=add_bias_term(bL=betaList[[b]],bT=t(x[[3]]),col=6,col_name = "bias")
	# betaList holds (for each bootstrap, b) all reg inters founds together with a bunch of confidence values: specifically
	# betaList[[b]] is a matrix with rows representing each reg inter. col_1=trgt,col_2=regulator,col_3=beta weigth,col_4=portion of xplnd var,
	#	col_5=mixCLR score, col_7=bias term associated with each reg inter (the intercept for each trgt regression model)
	rm(x)
        
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 8- run heuristic to combine results from different methods (MCZ, mixCLR, and Inf)--- i.e. results from different pipelines
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

	# elastic net produces a model for each l2 norm weight (choose a model for each target from the l2 norm weight with minimum CV error)
	beta.mat = combine_l2_net_res(betaList[[b]],modelErrorList[[b]],col="beta")
	# beta list is a sparse matrix representation.  Turn it into a matrix
	beta.mat = unsparse(beta.mat ,matrix(0,dim(Z)[1],dim(Z)[2]) )
	# same as beta.mat only instead of having beta weight as values it has predictive value for each reg inter
	pred.mat.lnet = combine_l2_net_res(betaList[[b]],modelErrorList[[b]],col="prd_xpln_var")
	pred.mat.lnet = unsparse(pred.mat.lnet,matrix(0,dim(Z)[1],dim(Z)[2]) )
	# for each trgt get the bias term (needed to predict system's response to new perturbations)
	pred.mat.bias = combine_l2_net_res(betaList[[b]],modelErrorList[[b]],col="bias")
	# this is the heuristic described in DREAM3 and DREAM4 papers z = sqrt(z1^2+z2^2)^2
	#  first for DREAM3 pipeline (not additive with MCZ)
	pred.mat.lnet.mixCLR = combine_mtrcs(Z,pred.mat.lnet)
	#  second for DREAM4 pipeline (i.e. DREAM3 + MCZ)
	#pred.mat.lnet.mixCLR.zKo = combine_mtrcs(Z_ko,pred.mat.lnet.mixCLR)

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 9- for DREAM4 predict systems response to double ko's
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

	#if( chalName == "DREAM4" ){
	#	S = pred.mat.lnet.mixCLR.zKo
	#	bestVal <- testPredFilters( S, pred.mat.bias, beta.mat, pred.mat.lnet, X_ko,X_wt,plotIt = FALSE)
	#	pred <- makePredictions( S,pred.mat.bias,beta.mat,X_ko,X_wt, pred.mat.lnet, inCut <- 75,INPUT$general$dbl_ko,"combine")
	#	predRes <- calcDblKoRmsd(pred,INPUT$general$dblKo_goldStandard,INPUT$general$dbl_ko,X_ko,X_wt,"meanKo")
	#}

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 10- store current re-sampling results
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

	allResults[[b]] <- list()
	#allResults[[b]][["MCZ"]] <- X_ko_z
        allResults[[b]][["Z"]] <- Z
	allResults[[b]][["betaMat"]] <- beta.mat
	allResults[[b]][["bias"]] <- pred.mat.bias
       	allResults[[b]][["pred.mat.lnet"]] <- pred.mat.lnet
	allResults[[b]][["MixCLR.Inf"]] <- pred.mat.lnet.mixCLR
        allResults[[b]][["modelErrorList"]] <- modelErrorList
	# every saveInt runs save data as .RData file
	saveInt <- 10
	if( b %% saveInt == 0){
		  cat("saving! \n")
			#saveAndPlot(makePlots = FALSE, bootNum <- b, PARAMS$general$saveToDir )
			save(betaList,INPUT,PARAMS,N_b, b, allResults, file=paste(PARAMS$general$saveToDir,"/savedData_numBoots_",b,".RData",sep=""))
			prevBootNum <- b - saveInt
                  if(prevBootNum>0){
			#the line below is specific for linux based systems (ie. NOT WINDOWS)
			system(paste("rm ", paste(PARAMS$general$saveToDir,"/savedData_numBoots_",prevBootNum,".RData",sep=""), sep=""))
                      }
	}

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 11- while b<N_b increament b by 1 adn repeat steps 2-6
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

	b = b + 1
}
cat("writing mixCLR->Inferelator output files.\n")
save(modelErrorList,betaList,allResults,file=paste(PARAMS$general$saveToDir,"/all_results_",b-1,".RData",sep=""))

# get median conf scores
median.conf.scores <- getMedianNetworkFromBootstraps(allResults, "MixCLR.Inf")
dimnames(median.conf.scores) <- dimnames(allResults[[1]][["MixCLR.Inf"]])
# write.table(median.conf.scores,file=paste(path.output.in.results,"clrinf_cnfdnc","_zcut_",GLOBAL[["z.abs.cut"]],"_Nb_",b-1,"_",GLOBAL[["date.is"]],".xls",sep=""),sep="\t")
write.table(median.conf.scores,file=paste(path.output.in.input,"clrinf_cnfdnc","_zcut_",GLOBAL[["z.abs.cut"]],"_Nb_",b-1,"_",GLOBAL[["date.is"]],".xls",sep=""),sep="\t")
# get inferelator beta weights
median.beta.scores <- getMedianNetworkFromBootstraps(allResults, "betaMat")
dimnames(median.beta.scores) <- dimnames(allResults[[1]][["MixCLR.Inf"]])
# write.table(median.beta.scores,file=paste(path.output.in.results, "inf_betas","_zcut_",GLOBAL[["z.abs.cut"]],"_Nb_",b-1,"_",GLOBAL[["date.is"]],".xls",sep=""),sep="\t")
write.table(median.beta.scores,file=paste(path.output.in.input,"inf_betas","_zcut_",GLOBAL[["z.abs.cut"]],"_Nb_",b-1,"_",GLOBAL[["date.is"]],".xls",sep=""),sep="\t")
# get inferelator confidence scores
median.pred.mat.lnet.scores <- getMedianNetworkFromBootstraps(allResults, "pred.mat.lnet")
dimnames(median.pred.mat.lnet.scores) <- dimnames(allResults[[1]][["MixCLR.Inf"]])
# write.table(median.pred.mat.lnet.scores,file=paste(path.output.in.results,"inf_cnfdnc","_zcut_",GLOBAL[["z.abs.cut"]],"_Nb_",b-1,"_",GLOBAL[["date.is"]],".xls",sep=""),sep="\t")
write.table(median.pred.mat.lnet.scores,file=paste(path.output.in.input,"inf_cnfdnc","_zcut_",GLOBAL[["z.abs.cut"]],"_Nb_",b-1,"_",GLOBAL[["date.is"]],".xls",sep=""),sep="\t")

##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# 13- cleanup tmp variables functions
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

	if (cleanUp) {
		rm(numGenesInNet,make_final_design_and_response_matrix,add_bias_term,add_weight_beta,add_zscore,calc_MI_inBatces,calc_ode_model_weights,
			calcDblKoRmsd,calcFoldChange,calcZscores,create_Pi_g,create_Pi_s,create_Xpi,createPermMatrix,fName,get_all_perms,get_best_preds_idx,
			get_usr_chosen_dataset,get_usr_chosen_design_matrix,get_usr_chosen_response,let_usr_choose_dataset,let_usr_choose_design_matrix,
			let_usr_choose_response,load_gold_standard,load_predictions,make_sparse2,makePredictions,modelErrorList,percentCoverage,permuteCols,
			Pi_s_clr,Pi_s_lars,saveInt,splitDreamDataByType,btch_size)
	}
