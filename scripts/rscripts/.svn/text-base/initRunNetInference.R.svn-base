##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Jan 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# get required packages
# note: the following sign: #~ means that these lines are run from main
library(elasticnet)
library(mi)
# library(inline)

# source scripts
source("r_scripts/th17/used_for_paper/init_utilRunNetInference.R")
source("r_scripts/common/utils.R")
source("r_scripts/common/larsUtil.R")
source("r_scripts/common/clr.R")
source("r_scripts/common/bootstrapUtil.R")





# init PARAMS and INPUT
#~ PARAMS = list()
INPUT = list()

#~ PARAMS[["general"]] = list()
#~ PARAMS[["clr"]] = list()
#~ PARAMS[["lars"]] = list()
#~ PARAMS[["output"]] = list()

#~ ######################### most useful params ##################################
#~ # how many predictors (tf's) for elastic net to choose from?
#~ PARAMS[["lars"]][["max_single_preds"]] = 30
#~ # what l2 norm weights to use in elastic net? (lambda =0 same as LARS)
#~ PARAMS[["lars"]][["lambda"]] = c(0) # (0,1,100)
#~ # how many bootstrap runs? (one for each bootstraped dataset)
#~ PARAMS[["general"]][["numBoots"]] = 200
#~ # what is the maximum delta T to consider (in minutes)? (if delta T is bigger than it will be treated as steady state)
#~ PARAMS[["general"]][["delT_max"]] = 7*60
#~ # what is the time order (tau) of the reactinos (i.e. in what time do you expect reactions to happen in)?
#~ PARAMS[["general"]][["tau"]] = 3*60
#~ # for tfs for which we have ko data.  What zscore to consider for regulation (bigger than this number in abs terms)?
#~ PARAMS[["general"]][["rnaseq_ko_zscore_abs_filter"]] = 1
#~ # how many processors to use? if use more than one need to install the multicore package
#~ PARAMS[["general"]][["processorsNumber"]] = 7
#~ # What n fold cross validation for inferelator should we use?
#~ PARAMS[["lars"]][["nCv"]] = 10
#~ # How many bins to use to calculate mutual information?
#~ PARAMS[["clr"]][["n.bins"]] = 10
#~ ################################################################################
if(PARAMS[["general"]][["processorsNumber"]]>1){
	library(multicore)
}
PARAMS[["general"]][["inf.version"]] = "nwInf.1.2"
x = unlist(strsplit(date()," "))
PARAMS[["general"]][["date"]] = paste(x[2],x[3],x[5],x[4],sep="_")

#PARAMS[["general"]][["data"]] = let_usr_choose_dataset()

PARAMS[["general"]][["use_t0_as_steady_state"]] = FALSE
PARAMS[["general"]][["use_mixCLR"]] = TRUE #for DREAM4
PARAMS[["general"]][["use_delt_bigger_than_delT_max_as_steady_state"]] = TRUE
# making the directory name to save files into
PARAMS[["general"]][["saveToDir"]] = paste("results/distributions/z_abs_cut_",GLOBAL[["z.abs.cut"]],"_",gsub(":","-",PARAMS[["general"]]$date),"_nboots_",PARAMS$"general"$"numBoots",sep="")
PARAMS[["general"]][["percentCoverage"]] = 100

# are we running for 'time-sereies' only, 'steady state' only, or 'all' for lars and clr respectively
PARAMS[["clr"]][["what_final_design_response_matrix"]] = 'all' # choose here between ts, ss, or all
PARAMS[["lars"]][["what_final_design_response_matrix"]] = 'all' # choose here between ts, ss, or all

# x = let_usr_choose_response()
x = list("inf_1_all_intervals","inf_1_all_intervals")
PARAMS[["clr"]][["response_matrix"]] = x[[1]]
PARAMS[["lars"]][["response_matrix"]] = x[[2]]

# x = let_usr_choose_design_matrix()
x = list("time_delayed","time_delayed")
PARAMS[["clr"]][["design_matrix"]] = x[[1]]
PARAMS[["lars"]][["design_matrix"]] = x[[2]]

x = get_usr_chosen_dataset(PARAMS[["general"]][["data"]])
INPUT[["general"]][["dataset"]] = x[[1]]
INPUT[["general"]][["clusterStack"]] = x[[2]]
INPUT[["general"]][["colMap"]] = x[[3]]
INPUT[["general"]][["tf_names"]] = x[[4]]
INPUT[["general"]][["tf_knockOutConfList"]] = x[[5]]


#get clr design matrix
if (PARAMS[["clr"]][["response_matrix"]] == 'inf_1_all_intervals') { 
	x = 'all_intervals' 
} else { 
	x = 'consecutive' 
}

params = c(PARAMS[["general"]][["delT_max"]],PARAMS[["clr"]][["design_matrix"]], x,
			  PARAMS[["general"]][["use_t0_as_steady_state"]],PARAMS[["general"]][["use_delt_bigger_than_delT_max_as_steady_state"]])
# get clr design matrix: 1- steady_state, 2- time_series
x = get_usr_chosen_design_matrix(INPUT[["general"]][["colMap"]], INPUT[["general"]][["dataset"]],params)

INPUT[["clr"]][["design_matrix_steady_state"]] = x[[1]]
INPUT[["clr"]][["design_matrix_time_series"]] = x[[2]]

#get lars design matrix
if (PARAMS[["lars"]][["response_matrix"]] == 'inf_1_all_intervals') { 
	x = 'all_intervals' 
} else { 
	x = 'consecutive' 
}

params = c(PARAMS[["general"]][["delT_max"]],PARAMS[["lars"]][["design_matrix"]], x,
			  PARAMS[["general"]][["use_t0_as_steady_state"]],PARAMS[["general"]][["use_delt_bigger_than_delT_max_as_steady_state"]])
# get lars design matrix: 1- steady_state, 2- time_series
x = get_usr_chosen_design_matrix(INPUT[["general"]][["colMap"]], as.matrix(as.data.frame(INPUT[["general"]][["dataset"]])[INPUT[["general"]][["tf_names"]],]),params)
											
INPUT[["lars"]][["design_matrix_steady_state"]] = x[[1]]
INPUT[["lars"]][["design_matrix_time_series"]] = x[[2]]

#get clr response matrix
params =  c(PARAMS[["general"]][["delT_max"]],PARAMS[["clr"]][["response_matrix"]], PARAMS[["general"]][["tau"]],
				PARAMS[["general"]][["use_t0_as_steady_state"]],PARAMS[["general"]][["use_delt_bigger_than_delT_max_as_steady_state"]])
if( is.null(INPUT[["general"]][["redExp"]])){
	x = get_usr_chosen_response(INPUT[["general"]][["colMap"]], INPUT[["general"]][["dataset"]], params)
}else{
	cat("getting clr response for biclusters \n")
	x = get_usr_chosen_response(INPUT[["general"]][["colMap"]],INPUT[["general"]][["redExp"]], params)
}
						 
INPUT[["clr"]][["response_matrix_steady_state"]] = x[[1]]
INPUT[["clr"]][["response_matrix_time_series"]] = x[[2]]

#get lars response matrix
params =  c(PARAMS[["general"]][["delT_max"]],PARAMS[["lars"]][["response_matrix"]], PARAMS[["general"]][["tau"]],
				PARAMS[["general"]][["use_t0_as_steady_state"]],PARAMS[["general"]][["use_delt_bigger_than_delT_max_as_steady_state"]])
if( is.null(INPUT[["general"]][["redExp"]])){
	x = get_usr_chosen_response(INPUT[["general"]][["colMap"]], INPUT[["general"]][["dataset"]], params)
}else{
	cat("getting lars response for biclusters \n")
	x = get_usr_chosen_response(INPUT[["general"]][["colMap"]], INPUT[["general"]][["redExp"]], params)
}
						 
INPUT[["lars"]][["response_matrix_steady_state"]] = x[[1]]
INPUT[["lars"]][["response_matrix_time_series"]] = x[[2]]

# make final design/response matrices for clr
x = make_final_design_and_response_matrix(INPUT[["clr"]][["design_matrix_steady_state"]] , 
												  INPUT[["clr"]][["design_matrix_time_series"]] , 
												  INPUT[["clr"]][["response_matrix_steady_state"]], 
												  INPUT[["clr"]][["response_matrix_time_series"]], 
												  PARAMS[["clr"]][["what_final_design_response_matrix"]]) 
												  
INPUT[["clr"]][["response_matrix"]] = x[[1]]
INPUT[["clr"]][["design_matrix"]] = x[[2]]

# make final design/response matrices for lars
x = make_final_design_and_response_matrix(INPUT[["lars"]][["design_matrix_steady_state"]] , 
												  INPUT[["lars"]][["design_matrix_time_series"]] , 
												  INPUT[["lars"]][["response_matrix_steady_state"]], 
												  INPUT[["lars"]][["response_matrix_time_series"]], 
												  PARAMS[["lars"]][["what_final_design_response_matrix"]]) 

INPUT[["lars"]][["response_matrix"]] = x[[1]]
INPUT[["lars"]][["design_matrix"]] = x[[2]]

# remove helper variable
rm(x,params)








