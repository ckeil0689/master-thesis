##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Nov 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# create a summary table for results

# set paths
path.input <- "input/th17/used_for_paper/"
path.input.cytoscape <- paste(sep="","results/usedForPaper/cytoscape_",date.cytoscape.run,"/")
path.output <- "results/"

# set file names
f.nm.specificity <- paste(sep="",path.input.cytoscape,"specificity_mat_zcut_",z.abs.cut,"_",date.cytoscape.run,".txt") 
f.nm.diff.rorc <- paste(sep="",path.output,"diff_expression/Th17.rorc.wt.vs.Th17.rorc.ko_",date.diff.expsn.run,".xls")
f.nm.diff.rorc.rora <- paste(sep="",path.output,"diff_expression/Th17.rorawt.rorcwt.vs.Th17.rorako.rorcko_",date.diff.expsn.run,".xls")
f.nm.chip <- paste(sep="",path.input,"chip_scores_",date.chip.run,".xls")
f.nm.sam <- paste(sep="",path.input,"samTh17VsTh0Zscores.xls")
f.nm.kc <- paste(sep="",path.input,"combinedAnalysis/score_combine_KC_abs_matrix_Nov_17_2011.xls")


# read files
# get specificity data
specificity.mat <- read.delim(file=f.nm.specificity,sep="\t")
rownames(specificity.mat) <- specificity.mat[,"gene_id"]
specificity.mat <- specificity.mat[,-1]
min_specificity <- specificity.mat[,"min.specificity"]
names(min_specificity) <- rownames(specificity.mat)

# get differential expression data for rorc
diff.exprsn.rorc.mat <- read.delim(file=f.nm.diff.rorc,sep="\t")
rownames(diff.exprsn.rorc.mat) <- as.character(diff.exprsn.rorc.mat[,"id"])
diff.exprsn.rorc.mat <- diff.exprsn.rorc.mat[,-1]
diff.exprsn.rorc.mat <- diff.exprsn.rorc.mat[,c(-1,-7,-8,-9,-10,-11)]

# get differential expression data for rorc/rora double ko
diff.exprsn.rorc.rora.mat <- read.delim(file=f.nm.diff.rorc.rora,sep="\t")
rownames(diff.exprsn.rorc.rora.mat) <- as.character(diff.exprsn.rorc.rora.mat[,"id"])
diff.exprsn.rorc.rora.mat <- diff.exprsn.rorc.rora.mat[,-1]
diff.exprsn.rorc.rora.mat <- diff.exprsn.rorc.rora.mat[,c(-1,-7,-8,-9)]

# read chip scores
to.load <- c("BATF_th17_minus_th0","IRF4_th17_minus_th0", # have th0 ctrl
             "MAF_th17_minus_th0","STAT3_th17_minus_th0","RORC_th17_minus_th0") # don't have th0 ctrl
scores.chip <- as.matrix(read.table(f.nm.chip,sep="\t"))[,to.load]

# read sam scores
scores.sam <- as.matrix(read.table(f.nm.sam ,sep="\t"))
scores.sam <- scores.sam[,"Score_d"]

# get kc data
kc.mat <- read.delim(file=f.nm.kc,sep="\t")[,-7]
tfs <- colnames(kc.mat)[1:5]
colnames(kc.mat) <- c(paste(sep="","KC_",tfs),"KC_sum")

gn.nms.keepers <- intersect(unique(rownames(kc.mat),rownames(specificity.mat)),unique(rownames(diff.exprsn.rorc.mat),rownames(diff.exprsn.rorc.rora.mat)))
gn.nms.keepers <- intersect(gn.nms.keepers,rownames(scores.chip))
gn.nms.keepers <- intersect(gn.nms.keepers,names(scores.sam))
gn_id <- gn.nms.keepers

min_specificity <- min_specificity[gn.nms.keepers]
scores.sam <- scores.sam[gn.nms.keepers]

x <- cbind(gn_id,diff.exprsn.rorc.mat[gn.nms.keepers,],diff.exprsn.rorc.rora.mat[gn.nms.keepers,],scores.chip[gn.nms.keepers,],
	min_specificity,kc.mat[gn.nms.keepers,],scores.sam)

	
# write summary table focused on rora and rorc
fl.nm <- paste(path.output,"summary_table_rora_rorc_",date.is,".xls",sep="")
cat(colnames(x),sep="\t",append=FALSE,file=fl.nm)
cat("\n",append=TRUE,file=fl.nm)
write.table(sep="\t",x,quote=FALSE,append=TRUE,col.names=FALSE,row.names=FALSE,file=fl.nm)


