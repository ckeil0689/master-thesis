# Read a single ChIP-seq file
setwd("/home/chrisk/Documents/uni/thesis/suppl/data/chipseq")

chipfiles <- list.files(getwd())

file=chipfiles[2]
cst=read.table(file, sep="\t", header=TRUE)

target=1

#gene size + 10kb on each side
Lg_len=cst$gn_length[target] + 20000
Lg_len

#returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\(|\\)$", "", x)

peaks=as.character(cst[target, 12])
peaks

peaks_sep <- unlist(strsplit(peaks, ";"))
peaks_sep

npeak=cst$genewide_n_peak[target]

# get pvals for gene
d_TSS=3
d_TES=4
pval_idx=6
nprox_peaks=0
gn_logpvals <- c()
gn_pvals <- c()
for(i in peaks_sep) {
  peak=unlist(strsplit(i, ","))
  
  #count prox_peaks (for validation)
  tss=abs(as.numeric(peak[d_TSS]))
  if(tss<=5000) {
    nprox_peaks = nprox_peaks+1
  }
  
  #get pvalues
  peak_pval=as.numeric(peak[pval_idx])
  peak_logpval=-log10(peak_pval)
  
  #append them to lists
  gn_logpvals=c(gn_logpvals, peak_logpval)
  gn_pvals=c(gn_pvals, peak_pval)
}

npeak
nprox_peaks

#genewide mean pval - maybe MACS -log10?
gn_pvals
gn_mean_pval=mean(gn_pvals)
gn_mean_pval

#genewide mean of -log10(pvals)
gn_logpvals
gn_mean_logpvals=mean(gn_logpvals)
gn_mean_logpvals

#length of genome (mm10)
genome_len=1.87e9

#count total number of peaks (genome-wide) with significance equal or greater to the average significance of peaks found in the gene
tf_strongpeaks=0
total_proxpeaks=0
macsout=cst$X.peak_id.summit.d_TSS.d_TES.class.pval.fold_enrich.FDR......
for(i in macsout) {
  peaklist <- unlist(strsplit(i, ";"))
  for(j in peaklist) {
    peak=unlist(strsplit(trim(j), ","))
    peak_pval=as.numeric(peak[pval_idx])
    #logpval=-log10(peak_pval)
    
    #count strong peaks
    #if(logpval >= gn_mean_logpvals) {
     # tf_strongpeaks=tf_strongpeaks + 1
    #}
    
    if(peak_pval >= gn_mean_pval) {
      tf_strongpeaks=tf_strongpeaks + 1
    }
    
    #count proximal peaks (only for validation)
    tss=abs(as.numeric(peak[d_TSS]))
    if(tss<=5000) {
      total_proxpeaks = total_proxpeaks+1
    }
  }
}

#ensure proximal peak count is valid
if(total_proxpeaks == sum(cst$prox_n_peak)) {
  print("Proximal peak count VALID.")
} else {
  print("Proximal peak count INVALID.")
}

total_npeak=sum(cst$genewide_n_peak)

#ratio of strong peaks 
tf_strongpeaks/total_npeak

#expected number of strong peaks per bp for TF 
px=tf_strongpeaks/genome_len
px

#expected number of peaks per bp for TF 
px2=total_npeak/genome_len
px2

#lambda for poisson dist
Lg_len
prox_lambda=px*Lg_len
prox_lambda

lambda=px2*Lg_len
lambda

#poisson dist
pois=ppois(npeak, lambda)
pois

log_pois=-log10(pois)
log_pois

#prox poisson dist
proxpois=ppois(nprox_peaks, prox_lambda)
proxpois

log_proxpois=-log10(proxpois)
log_proxpois

#rank the pois_model pvals (needs to be done for full matrix, not just this file)
scores=cst$genewide_pois_model_pval
r1=rank(abs(scores), na.last = "keep")
#r1

count=sum(scores!=0)
count

#get q-vals from ranks
q=(1-(r1/count))
sum(q>0.80)