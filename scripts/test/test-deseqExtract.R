#!/usr/bin/env Rscript

# Testing DESeq score matrix generation
source(paste0(getwd(), "/../" , "deseqExtract-fun.R"), chdir = TRUE)

# DEseq file directory relative to /scripts/
deseqdir.test <- paste0(getwd(), "/../../suppl/data/deseq/")

# Ensure we are in correct directory
if(!dir.exists(deseqdir.test)) stop("Cannot load DESeq-files (test) because the directory does not exist.")

# We consider all DESeq result files (as opposed to selective ChIP-seq loading)
deseqfiles.test <- list.files(deseqdir.test)
if(length(deseqfiles.test) == 0) stop(paste("No DESeq files found in: ", deseqdir, "Stopping (test)."))

context("TF extraction from DESeq files")
test_that("TF name is correctly extracted from DESeq file", {
  
  # Column name 2 and 3 as expected
  batf <- c("id",	"baseMean",	"Th17.batf.ko",	"Th17.batf.wt",	"foldChange",	"log2FoldChange",	
            "pval",	"padj",	"resVarA",	"resVarB",	"mean.rpkm.ko",	"mean.rpkm.wt",	
            "mean.rpkm.th0",	"mean.rpkm.th17",	"prcnt.chng")
  maf <- c("id",	"baseMean",	"Th17.maf.ko",	"Th17.maf.wt",	"foldChange",	"log2FoldChange",	
           "pval",	"padj",	"resVarA",	"resVarB",	"mean.rpkm.ko",	"mean.rpkm.wt",	
           "mean.rpkm.th0",	"mean.rpkm.th17",	"prcnt.chng")
  
  # Custom files
  maf.custom <- c("id",	"Th17.maf.ko", "log2FoldChange", "pval")
  
  # Incomplete
  stat <- c("id",	"baseMean",	"Th17.stat3.",	"Th17.stat3.",	"foldChange",	"log2FoldChange",	
            "pval",	"padj",	"resVarA",	"resVarB",	"mean.rpkm.ko",	"mean.rpkm.wt",	
            "mean.rpkm.th0",	"mean.rpkm.th17",	"prcnt.chng")
  rorc <- c("id",	"baseMean",	".rorc.ko",	".rorc.wt",	"foldChange",	"log2FoldChange",	
            "pval",	"padj",	"resVarA",	"resVarB",	"mean.rpkm.ko",	"mean.rpkm.wt",	
            "mean.rpkm.th0",	"mean.rpkm.th17",	"prcnt.chng")
  hifla <- c("id",	"baseMean",	".hif1a.",	".hif1a.",	"foldChange",	"log2FoldChange",	
             "pval",	"padj",	"resVarA",	"resVarB",	"mean.rpkm.ko",	"mean.rpkm.wt",	
             "mean.rpkm.th0",	"mean.rpkm.th17",	"prcnt.chng")
  
  # Faulty 
  fosl2 <- c("id",	"baseMean",	"Th17-fosl2-ko",	"Th17-fosl2-wt",	"foldChange",	"log2FoldChange",	
             "pval",	"padj",	"resVarA",	"resVarB",	"mean.rpkm.ko",	"mean.rpkm.wt",	
             "mean.rpkm.th0",	"mean.rpkm.th17",	"prcnt.chng")
  none1 <- c("id",	"baseMean",	"Th17..ko",	"Th17..wt",	"foldChange",	"log2FoldChange",	
             "pval",	"padj",	"resVarA",	"resVarB",	"mean.rpkm.ko",	"mean.rpkm.wt",	
             "mean.rpkm.th0",	"mean.rpkm.th17",	"prcnt.chng")
  none2 <- c("id",	"baseMean",	"Th17.ko",	"Th17.wt",	"foldChange",	"log2FoldChange",	
             "pval",	"padj",	"resVarA",	"resVarB",	"mean.rpkm.ko",	"mean.rpkm.wt",	
             "mean.rpkm.th0",	"mean.rpkm.th17",	"prcnt.chng")
  
  expect_that(extract.tf(batf, "batf-test"), is_identical_to("BATF"))
  expect_that(extract.tf(maf, "maf-test"), is_identical_to("MAF"))
  expect_that(extract.tf(maf.custom, "maf.custom-test"), is_identical_to("MAF"))
  expect_that(extract.tf(stat, "stat-test"), is_identical_to("STAT3"))
  expect_that(extract.tf(rorc, "rorc-test"), is_identical_to("RORC"))
  expect_that(extract.tf(hifla, "hif1a-test"), is_identical_to("HIF1A"))
  expect_that(tf <- extract.tf(fosl2, "fosl2-test"), gives_warning())
  expect_that(is.na(tf), is_true())
  expect_that(tf <- extract.tf(none1, "none1-test"), gives_warning())
  expect_that(is.na(tf), is_true())
  expect_that(tf <- extract.tf(none2, "none2-test"), gives_warning())
  expect_that(is.na(tf), is_true())
})

context("Setup of skeleton matrix for DESeq scores")
test_that("DESeq files are correctly parsed and skeleton matrix is setup as expected", {
  
  skel.mat <- get.skel.mat()
  expect_that(skel.mat, is_a("matrix"))
  expect_that(length(rownames(skel.mat)) > 0, is_true())
  expect_that(length(colnames(skel.mat)) > 0, is_true())
  
  # All values should be zero
  expect_that(length(skel.mat[skel.mat != 0]) == 0, is_true())
  
  # Column names should all be part of CORE_TFS
  expect_that(all(tolower(colnames(skel.mat)) %in% GLOBAL[["CORE_TFS"]]), is_true())
})

context("Data population of empty matrix from DESeq files")
test_that("DESeq-skeleton matrix is populated with p-val*log2 data for every TF-gene pair as expected.", {
  
  score.skel.mat <- get.skel.mat()
  deseq.scores <- populate.deseq.scores(score.skel.mat)
  
  # Check basic attributes
  expect_that(deseq.scores, is_a("matrix"))
  expect_that(length(deseq.scores[is.infinite(deseq.scores)]) == 0, is_true())
  expect_that(length(deseq.scores[is.na(deseq.scores)]) == 0, is_true())
  expect_that(length(rownames(deseq.scores)) > 0, is_true())
  expect_that(length(colnames(deseq.scores)) > 0, is_true())
  # Column names should all be part of CORE_TFS
  expect_that(all(tolower(colnames(deseq.scores)) %in% GLOBAL[["CORE_TFS"]]), is_true())
  
  # Compare known GEO samples to some hand-calculated results --> -log10(pval) * sign(log2foldchange)
  if("GSE40918_Th17.batf.wt.vs.Th17.batf.ko_Aug_2_2012.txt" %in% deseqfiles.test) {
    expect_that(deseq.scores["IKZF3", "BATF"], equals(0.1430428, tolerance = 1e-7))
    expect_that(deseq.scores["CROCC", "BATF"], equals(0.1019569, tolerance = 1e-7))
    expect_that(deseq.scores["ARL11", "BATF"], equals(0, tolerance = 1e-7))
  }
  
  if("GSE40918_Th17.maf.wt.vs.Th17.maf.ko_Aug_2_2012.txt" %in% deseqfiles.test) {
    expect_that(deseq.scores["ELMO1", "MAF"], equals(0.3227558, tolerance = 1e-7))
    expect_that(deseq.scores["INA", "MAF"], equals(0, tolerance = 1e-7))
    expect_that(deseq.scores["INPP4B", "MAF"], equals(-4.989882, tolerance = 1e-7))
  }
  
  if("GSE40918_Th17.rorc.wt.vs.Th17.rorc.ko_Aug_2_2012.txt" %in% deseqfiles.test) {
    expect_that(deseq.scores["DUSP6", "RORC"], equals(-3.361983, tolerance = 1e-7))
    expect_that(deseq.scores["GM5124", "RORC"], equals(0.02793168, tolerance = 1e-7))
    expect_that(deseq.scores["NOL6", "RORC"], equals(0.4853783, tolerance = 1e-7))
  }
})