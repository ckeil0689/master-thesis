##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Apr 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# R code to create quality control plots
# output:
#   - chipseq steps improvement in ranks of known th17 specific genes

# paths
path.input <- "input/th17/used_for_paper/"
path.output <- paste(sep="","results/usedForPaper/qc_",date.is,"/")
system(paste(sep="","mkdir ",path.output))

# define file names
kc.f.nm <- paste(sep="",)

## score_combine_KC_matrix_zcut_2.5_Apr_24_2011.xls
# get chipScores after MACS for each experiment
C <- as.matrix(read.table(paste(path.input,"chipScores.xls",sep=""),header=T,sep="\t"))
gns <- rownames(C)

# calc quantiles of known genes to see how each step in chipseq compares to previous step
# p300
x <- sort(C[,"P300_Th17_1041_minus_p300_Th0_1948_S_gene"],decreasing=T)
ix.p300.th17.min.th0 <- which(names(x) %in% known.gns)
# batf
x <- sort(C[,"BATF_Th17_115_minus_BATF_Th0_130_S_gene"],decreasing=T)
ix.batf.th17.min.th0 <- which(names(x) %in% known.gns)
# batf with p300
x <- sort(C[,"BATF_Th17_115_minus_BATF_Th0_130_S_gene_plus_p300"],decreasing=T)
ix.batf.th17.min.th0.plus.p300 <- which(names(x) %in% known.gns)
# irf4
x <- sort(C[,"IRF4_Th17_971_minus_IRF4_Th0_1235_S_gene"],decreasing=T)
ix.irf4.th17.min.th0 <- which(names(x) %in% known.gns)
# irf4 with p300
x <- sort(C[,"IRF4_Th17_971_minus_IRF4_Th0_1235_S_gene_plus_p300"],decreasing=T)
ix.irf4.th17.min.th0.plus.p300 <- which(names(x) %in% known.gns)
# maf
x <- sort(C[,"Maf_Th17_102_S_gene"],decreasing=T)
ix.maf.th17.min.th0 <- which(names(x) %in% known.gns)
# maf with p300
x <- sort(C[,"Maf_Th17_102_S_gene_plus_p300"],decreasing=T)
ix.maf.th17.min.th0.plus.p300 <- which(names(x) %in% known.gns)
# stat3
x <- sort(C[,"STAT3_Th17_1040_S_gene"],decreasing=T)
ix.stat3.th17.min.th0 <- which(names(x) %in% known.gns)
# stat3 with p300
x <- sort(C[,"STAT3_Th17_1040_S_gene_plus_p300"],decreasing=T)
ix.stat3.th17.min.th0.plus.p300 <- which(names(x) %in% known.gns)
# rorc
x <- sort(C[,"RORc_Th17_93_S_gene"],decreasing=T)
ix.rorc.th17.min.th0 <- which(names(x) %in% known.gns)
# rorc with p300
x <- sort(C[,"RORc_Th17_93_S_gene_plus_p300"],decreasing=T)
ix.rorc.th17.min.th0.plus.p300 <- which(names(x) %in% known.gns)

































