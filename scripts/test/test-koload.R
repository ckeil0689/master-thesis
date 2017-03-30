#!/usr/bin/env Rscript

# Testing knockout score matrix generation
source(paste0(getwd(), "/../" , "deseqExtract-fun.R"), chdir = TRUE)

context("Testing TF extraction")

test_that("TF name is correctly extracted from DESeq file", {
  
  # Column name 2 and 3 as expected
  batf <- c("id",	"baseMean",	"Th17.batf.ko",	"Th17.batf.wt",	"foldChange",	"log2FoldChange",	
            "pval",	"padj",	"resVarA",	"resVarB",	"mean.rpkm.ko",	"mean.rpkm.wt",	
            "mean.rpkm.th0",	"mean.rpkm.th17",	"prcnt.chng")
  maf <- c("id",	"baseMean",	"Th17.maf.ko",	"Th17.maf.wt",	"foldChange",	"log2FoldChange",	
           "pval",	"padj",	"resVarA",	"resVarB",	"mean.rpkm.ko",	"mean.rpkm.wt",	
           "mean.rpkm.th0",	"mean.rpkm.th17",	"prcnt.chng")
  
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
  
  expect_that(extract.tf(batf), is_identical_to("BATF"))
  expect_that(extract.tf(maf), is_identical_to("MAF"))
  expect_that(extract.tf(stat), is_identical_to("STAT3"))
  expect_that(extract.tf(rorc), is_identical_to("RORC"))
  expect_that(extract.tf(hifla), is_identical_to("HIF1A"))
  expect_that(tf <- extract.tf(fosl2), gives_warning())
  expect_that(is.na(tf), is_true())
  expect_that(tf <- extract.tf(none1), gives_warning())
  expect_that(is.na(tf), is_true())
  expect_that(tf <- extract.tf(none2), gives_warning())
  expect_that(is.na(tf), is_true())
})

context("Testing setup of skeleton matrix")

test_that("DESeq files are correctly parsed and KO skeleton matrix is setup as expected", {
  
  expect_that(skel.mat <- get.skel.mat(), gives_warning()) # skips TF ikzf3 (not in CORE)
  expect_that(skel.mat, is_a("matrix"))
  expect_that(length(rownames(skel.mat)) > 0, is_true())
  expect_that(length(colnames(skel.mat)) > 0, is_true())
  
  # All values should be zero
  expect_that(length(skel.mat[skel.mat != 0]) == 0, is_true())
  
  # Column names should all be part of CORE_TFS
  expect_that(all(tolower(colnames(skel.mat)) %in% GLOBAL[["CORE_TFS"]]), is_true())
})

context("Testing data population of empty KO-matrix from DESeq files")

test_that("KO-skeleton matrix is populated with p-val*log2 data for every TF-gene pair as expected.", {
  
  score.skel.mat <- get.skel.mat()
  ko.scores <- populate.ko.scores(score.skel.mat)
  
  # Check basic attributes
  expect_that(ko.scores, is_a("matrix"))
  expect_that(length(ko.scores[is.infinite(ko.scores)]) == 0, is_true())
  expect_that(length(ko.scores[is.na(ko.scores)]) == 0, is_true())
  expect_that(length(rownames(ko.scores)) > 0, is_true())
  expect_that(length(colnames(ko.scores)) > 0, is_true())
  # Column names should all be part of CORE_TFS
  expect_that(all(tolower(colnames(ko.scores)) %in% GLOBAL[["CORE_TFS"]]), is_true())
  
  # Compare samples to some hand-calculated results.
  # expect_that(ko.scores[row, col], is_identical_to(num)) ...
})