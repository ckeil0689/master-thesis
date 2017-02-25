# #!/usr/bin/env Rscript
# Master script for processing TF-target gene interaction data into network files for Cytoscape 
# (see Computational Methods paper, http://www.cell.com/action/showImagesData?pii=S0092-8674%2812%2901123-3 ) 
# 
# Input: 
#       C = customized MACS output (ChIP-seq) 
#       K = dESeq files for core TFs (RNA-seq KO) 
#       R = Inferelator matrix (RNA-seq compendium)
#       I = Inferelator matrix (Immgen microarray data)

# Input file directories
chipdir <- paste0(getwd(), "/../suppl/data/chipseq/")
deseqdir <- paste0(getwd(), "/../suppl/data/deseq/")

# Put all output into analysis folder
outpath <- paste0(getwd(), "/../suppl/data/analysis/")
if(!dir.exists(outpath)) {
  dir.create(outpath)
}

ALLOWED_COMBOS <- c("ri", "kc", "kcr", "kcri")
ALLOWED_CELLS <- c("Th17", "Th0")

# Get user input 
args <- commandArgs(trailingOnly=TRUE)
thx <- toupper(as.character(args[1]))
combo <- tolower(as.character(args[2]))
#thx <-'Th17'

# Check user input for valiitiy
print_usage <- function() {
  print(paste("Usage: ./script.r <combo> <cell-type>", "(k = rnaseq-ko, c = chipseq, r = rna-compendium, i = immgen) (Th17 | Th0)"))
}

if(length(args) != 2){
  print_usage
  stop("Incorrect number of arguments.")
}
   
if(combo == "" || !combo %in% ALLOWED_COMBOS) {
  err <- c("Problem with argument. Enter a valid combination: ", ALLOWED_COMBOS)
  print(paste(err, collapse = " "))
  print_usage
  stop("invalid argument")
}

if(thx == "" || !thx %in% ALLOWED_CELLS) {
  err <- c("Problem with argument. Enter a valid cell type: ", ALLOWED_CELLS)
  print(paste(err, collapse = " "))
  print_usage
  stop("invalid argument")
}

# 


# STEP 1: Read a ChIP-seq files from dir
chip_inpath <- paste0(getwd(), "/../suppl/data/chipseq/")
setwd(chip_inpath)

