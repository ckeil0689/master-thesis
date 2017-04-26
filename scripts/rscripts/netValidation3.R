##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Oct 2011 Th17 project (MCZ,tlCLR,Inferelator)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

# R code to create quality control plots for entire network pipeline

# here we validate chipseq processing

# set paths
path.input <- "input/th17/used_for_paper/"
path.output <- "results/"
# set file names
chip.f.nm <- paste(sep="","chip_scores_",date.chip.scores.run,".xls")
sam.f.nm <- "samTh17VsTh0Zscores.xls"

# load files
chip <- as.matrix(read.table(file=paste(sep="",path.input,chip.f.nm)))
sam <- as.matrix(read.table(file=paste(sep="",path.input,sam.f.nm)))

# limit to genes that are have a peak and SAM scores
gn.sam <- rownames(sam)
gn.chip <- rownames(chip)
gns.keep <- intersect(gn.sam,gn.chip)

# cbind chip and sam in one matrix
m <- cbind(chip[gns.keep,],sam[gns.keep,])

# set filter for genes with median rpkm (over th17/th0 data) > 5
ix.rpkm <- which(m[,"median_rpkm"] > 5)

# set filter for genes that have a diff expression score of 2/-2 or more
ix.diff <- which(abs(m[,"Score_d"]) > 2)

# get a filter for both diff and rpkm
ix.diff.and.rpkm <- intersect(ix.diff,ix.rpkm)

# set matrices for the filters
m.rpkm <- m[ix.rpkm,]
m.diff <- m[ix.diff,]
m.diff.and.rpkm <- m[ix.diff.and.rpkm,]

# calc correlation matrix for matrices
cor.m <- cor(m)
cor.m.rpkm <- cor(m.rpkm)
cor.m.diff <- cor(m.diff)
cor.m.diff.and.rpkm <- cor(m.diff.and.rpkm)

order.tfs <- c("P300","BATF","IRF4","MAF","STAT3","RORC")
# rm.expts <- c("Score_d","q.value","median_rpkm","number_zero_rpkm_experiments")
# row.rm <- which(colnames(cor.m) %in% rm.expts)
# row.keep <- (1:ncol(cor.m))[-row.rm]

pdf(paste(sep="",path.output,"validation3.pdf"))
# plot barplots for each dataset
for(i in 1:4){
	if(i==1){
		type="all genes"
		x=cor.m
	} else if (i==2){
		type="genes with median rpkm > 5"
		x=cor.m.rpkm
	} else if (i==3){
		type="genes with abs diff expression > 2"
		x=cor.m.diff
	} else if (i==4){
		type="genes with abs diff expression > 2 and rpkm > 5"
		x=cor.m.diff.and.rpkm
	}
	v <- x[,"Score_d"]
	v.mat <- cbind( v[grep(order.tfs[1],names(v))[c(2,1,3,4)]],
					v[grep(order.tfs[2],names(v))[c(2,1,3,4)]],
					v[grep(order.tfs[3],names(v))[c(2,1,3,4)]],
					v[grep(order.tfs[4],names(v))[c(2,1,3,4)]],
					v[grep(order.tfs[5],names(v))[c(2,1,3,4)]],
					v[grep(order.tfs[6],names(v))[c(2,1,3,4)]])
	rownames(v.mat) <- c("Th0","Th17","Th17-Th0","Th17-Th0+P300")
	colnames(v.mat) <- order.tfs
	v.mat[4,"P300"] <- 0
	barplot(v.mat, main=paste(sep="","cor with differential expression \n(", type,")"),ylab="correlation", col=rainbow(4),ylim=c(min(v.mat),0.4),
	  legend = rownames(v.mat), beside=TRUE)
}
dev.off()










