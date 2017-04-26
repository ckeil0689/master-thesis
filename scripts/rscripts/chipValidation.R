##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Sep 2011 th17
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '


path.input <- "input/th17/used_for_paper/"
## path.output <- paste(sep="","results/validation_",date.is,"/")
path.output <- paste(sep="","results/validation/")
path.input.sam <- "input/th17/used_for_paper/"

# define file names
chip.f.nm <- paste(sep="",path.input,"chip_scores_",date.chip.scores.run,".xls")
sam.f.nm <- paste(sep="",path.input,"samTh17VsTh0Zscores.xls")

# load sam scores and chip scores
S <- as.matrix(read.table(sam.f.nm,header=T,sep="\t"))
C <- as.matrix(read.table(chip.f.nm,header=T,sep="\t"))
sam <- S[,"Score_d"]

# combine S and C over genes that are present in both and are differentially expressed
## gns <- names(sam)[which(abs(sam)>z.abs.cut)]
gns <- names(sam)
gns <- which(gns %in% rownames(C))
sam <- sam[gns]
x <- cbind(sam,C[gns,])
w <- cor(x)
heatmap(w, Rowv=NA, Colv=NA, scale="column", margins=c(5,10))
