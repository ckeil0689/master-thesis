#!/usr/bin/env Rscript

# Testing ChIP score matrix generation
source(paste0(getwd(), "/../" , "chipExtract-fun.R"), chdir = TRUE)

CORE_TFS <- c("batf", "irf4", "stat3", "maf", "rorc", "fosl2", "hif1a")

context("Fixing TF names")

test_that("TF names are fixed to match conventional output (matching Cell authors example matrices)", {
  expect_that(fix.tf.name("I_am_TF-GAMMA"), equals("iamtfgamma"))
  expect_that(fix.tf.name("RoRG"), equals("rorc"))
  expect_that(fix.tf.name("c-mAF"), equals("maf"))
  expect_that(fix.tf.name(""), equals(""))
  expect_that(fix.tf.name("----__-"), equals(""))
})

context("Extracting information from reference table")

test_that("TF names are extracted correctly from reference table", {
  # Files known as TF-Thx experiment
  rorc_th0 <- "GSM1004853_SL3779_SL3778_genes.txt"
  maf_th0 <- "GSM1004798_SL4424_SL4425_genes.txt"
  batf_th17 <- "GSM1004787_SL3037_SL3036_genes.txt"
  stat3_th17 <- "GSM1004865_SL3315_SL3319_genes.txt"
  p300_th0 <- "GSM1004842_SL1948_SL1947_genes.txt"
  p300_th17 <- "GSM1004851_SL3594_SL3592_genes.txt"
  
  non_existent <- "GSM0000000_SL00000_SL00000_genes.txt"
  gibberish <- "ceverwrgwegrfr.txt"
  empty <- ""
  
  expect_that(extract.tf.from.ref(rorc_th0, boost.p300 = TRUE, CORE_TFS), is_identical_to("rorc-th0"))
  expect_that(extract.tf.from.ref(maf_th0, boost.p300 = TRUE, CORE_TFS), is_identical_to("maf-th0"))
  expect_that(extract.tf.from.ref(batf_th17, boost.p300 = TRUE, CORE_TFS), is_identical_to("batf-th17"))
  expect_that(extract.tf.from.ref(stat3_th17, boost.p300 = TRUE, CORE_TFS), is_identical_to("stat3-th17"))
  expect_that(extract.tf.from.ref(p300_th0, boost.p300 = TRUE, CORE_TFS), is_identical_to("p300-th0"))
  expect_that(extract.tf.from.ref(p300_th17, boost.p300 = TRUE, CORE_TFS), is_identical_to("p300-th17"))
  expect_that(extract.tf.from.ref(non_existent, boost.p300 = TRUE, CORE_TFS), throws_error("No library match found."))
  expect_that(extract.tf.from.ref(gibberish, boost.p300 = TRUE, CORE_TFS), throws_error("No library match found."))
  expect_that(extract.tf.from.ref(empty, boost.p300 = TRUE, CORE_TFS), throws_error("No library match found."))
  expect_that(extract.tf.from.ref(NULL, boost.p300 = TRUE, CORE_TFS), throws_error("No experiment passed. Stopping."))
})

context("Generating the ChIP-seq confidence score matrix")

test_that("Skeleton matrix is created as expected", {
  # Files used as sample 
  all.chipfiles.boost <- c("GSM1004785_SL3192_SL3190_genes.txt", # Th0 BATF wt
                           "GSM1004824_SL1235_SL1234_genes.txt", # Th0 IRF4 wt
                           "GSM1004787_SL3037_SL3036_genes.txt", # Th17 BATF wt
                           "GSM1004833_SL2872_SL2876_genes.txt", # Th17 IRF4 rorc wt
                           "GSM1004842_SL1948_SL1947_genes.txt", # P300 Th0 wt
                           "GSM1004851_SL3594_SL3592_genes.txt") # P300 Th17 wt
  
  all.chipfiles.noboost <- c("GSM1004785_SL3192_SL3190_genes.txt", # Th0 BATF wt
                           "GSM1004824_SL1235_SL1234_genes.txt", # Th0 IRF4 wt
                           "GSM1004787_SL3037_SL3036_genes.txt", # Th17 BATF wt
                           "GSM1004833_SL2872_SL2876_genes.txt") # Th17 IRF4 rorc wt
  
  skel.matrix.boost <- get.skel.matrix(all.chipfiles.boost, TRUE, CORE_TFS)
  skel.matrix.noboost <- get.skel.matrix(all.chipfiles.noboost, FALSE, CORE_TFS)
  
  # Make sure the type is matrix
  expect_that(skel.matrix.boost, is_a("matrix"))
  expect_that(skel.matrix.noboost, is_a("matrix"))
  
  # Expected column names
  exp.colnames.boost <- c("batf-th0", "irf4-th0", "batf-th17", "irf4-th17", "p300-th0", "p300-th17")
  exp.colnames.noboost <- c("batf-th0", "irf4-th0", "batf-th17", "irf4-th17")
  
  # Make sure we have rownames defined and the column names exactly match the expected examples.
  expect_that(length(rownames(skel.matrix.boost)) > 0, is_true())
  expect_that(length(rownames(skel.matrix.noboost)) > 0, is_true())
  expect_that(length(colnames(skel.matrix.boost)) == 6, is_true())
  expect_that(length(colnames(skel.matrix.noboost)) == 4, is_true())
  expect_that(rownames(skel.matrix.boost), is_a("character"))
  expect_that(rownames(skel.matrix.noboost), is_a("character"))
  expect_that(colnames(skel.matrix.boost), is_a("character"))
  expect_that(colnames(skel.matrix.noboost), is_a("character"))
  expect_that(colnames(skel.matrix.boost), is_identical_to(exp.colnames.boost))
  expect_that(colnames(skel.matrix.noboost), is_identical_to(exp.colnames.noboost))
})

test_that("Skeleton matrix is filled with Poisson p-values as expected", {
  # Create a 3x3 skeleton dummy matrix
  skel.mat.boost <- matrix(0, nrow = 3, ncol = 3)
  colnames(skel.mat.boost) <- c("maf-th0", "batf-th17", "rorc-th0")
  rownames(skel.mat.boost) <- c("g1", "g2", "g3")
  
  # Create 3 tmp files from which we will read the Poisson values
  # File names, modified original file names by inserting '_tmp_' so library ID can still be extracted and TF name found in reference
  maf.tmpfile <- "GSM1004798_SL4424_SL4425_tmp_genes.txt"
  batf.tmpfile <- "GSM1004787_SL3037_SL3036_tmp_genes.txt"
  rorc.tmpfile <- "GSM1004853_SL3779_SL3778_tmp_genes.txt"
  tmp.chipfiles <- c(maf.tmpfile, # MAF Th0 wt
                     batf.tmpfile, #BATF Th17 wt
                     rorc.tmpfile) # RORC Th0 wt
                     
  # For MAF-th0
  maf.tmp <- matrix(0, nrow = 2, ncol = 2)
  colnames(maf.tmp) <- c("Gene_ID", "genewide_pois_model_pval")
  maf.tmp[,"Gene_ID"] <- c("g1", "g2")
  maf.tmp[,"genewide_pois_model_pval"] <- c(22.4198441743, 3.5005671754)
  write.table(maf.tmp, file = maf.tmpfile, sep = "\t", row.names = TRUE, col.names = NA)
  
  # For BATF-th17
  batf.tmp <- matrix(0, nrow = 2, ncol = 2)
  colnames(batf.tmp) <- c("Gene_ID", "genewide_pois_model_pval")
  batf.tmp[,"Gene_ID"] <- c("g1", "g3")
  batf.tmp[,"genewide_pois_model_pval"] <- c(1.9055811328, 5.3287439411)
  write.table(batf.tmp, file = batf.tmpfile, sep = "\t", row.names = TRUE, col.names = NA)
  
  # For RORC-th0
  rorc.tmp <- matrix(0, nrow = 3, ncol = 2)
  colnames(rorc.tmp) <- c("Gene_ID", "genewide_pois_model_pval")
  rorc.tmp[,"Gene_ID"] <- c("g1", "g2", "g3")
  rorc.tmp[,"genewide_pois_model_pval"] <- c(1.9957011261, 4.5485324703, 13.1669802036)
  write.table(rorc.tmp, file = rorc.tmpfile, sep = "\t", row.names = TRUE, col.names = NA)
  
  # 1) Normal test (skeleton matrix matches loaded ChIP-files)
  expected.result <- matrix(0, nrow = 3, ncol = 3)
  colnames(expected.result) <- c("maf-th0", "batf-th17", "rorc-th0")
  rownames(expected.result) <- c("g1", "g2", "g3")
  expected.result[,"maf-th0"] <- c(22.4198441743, 3.5005671754, 0)
  expected.result[,"batf-th17"] <- c(1.9055811328, 0, 5.3287439411)
  expected.result[,"rorc-th0"] <- c(1.9957011261, 4.5485324703, 13.1669802036)
  
  mat.pois.boost <- get.pois.vals(skel.mat.boost, tmp.chipfiles, TRUE, CORE_TFS)
  
  expect_that(mat.pois.boost, is_a("matrix"))
  expect_that(mat.pois.boost, is_identical_to(expected.result))
  
  # 2) Colnames in skeleton that do not occur in loaded files are skipped when filling pois.mat
  colnames(skel.mat.boost) <- c("stat3-th0", "batf-th17", "rorc-th0")
  colnames(expected.result)[[1]] <- "stat3-th0"
  expected.result[,"stat3-th0"] <- c(0, 0, 0) # no match -> values never filled
  
  mat.pois.boost <- get.pois.vals(skel.mat.boost, tmp.chipfiles, TRUE, CORE_TFS)
  
  expect_that(mat.pois.boost, is_a("matrix"))
  expect_that(mat.pois.boost, is_identical_to(expected.result))
  
  # 3) Broken/ bad (non-numerical) input is added as zero value
  rorc.tmp[,"genewide_pois_model_pval"] <- c("ghost", TRUE, 13.1669802036)
  write.table(rorc.tmp, file = rorc.tmpfile, sep = "\t", row.names = TRUE, col.names = NA)
  
  expected.result[,"rorc-th0"] <- c(0, 0, 13.1669802036)
  mat.pois.boost <- get.pois.vals(skel.mat.boost, tmp.chipfiles, TRUE, CORE_TFS)
  
  print(mat.pois.boost)
  expect_that(mat.pois.boost, is_a("matrix"))
  expect_that(mat.pois.boost, is_identical_to(expected.result))
  
  # Clean tmp files
  file.remove(maf.tmpfile)
  file.remove(batf.tmpfile)
  file.remove(rorc.tmpfile)
})

test_that("ChIP confidence scores are correctly calculated from Poisson p-value matrix", {
  
  genes.unique <- c("g1", "g2", "g3")
  tfs.list.unique <- c("batf", "stat3")
  
  # Sample Poisson value matrix as it would come from get.pois.vals()
  pois.mat <- matrix(0, nrow = 3, ncol= 6)
  rownames(pois.mat) <- genes.unique
  colnames(pois.mat) <- c("batf-th0", "batf-th17", "p300-th0", "p300-th17", "stat3-th0", "stat3-th17")
  
  # Fill with typicl data
  batf.th17 <- c(1.9957011261, 4.5485324703, 13.1669802036)
  batf.th0 <- c(0, 22.4198441743, 3.5005671754)
  p300.th17 <- c(1.9055811328, 0, 5.3287439411)
  p300.th0 <- c(0, 0, 13.1669802036)
  stat3.th17 <- c(3.1876487024, 2.6303732467, 4.5760732998)
  stat3.th0 <- c(7.5619903288, 2.4514891757, 1.1453294619)
  pois.mat[, "batf-th17"] <- batf.th17
  pois.mat[, "batf-th0"] <- batf.th0
  pois.mat[, "p300-th17"] <- p300.th17
  pois.mat[, "p300-th0"] <- p300.th0
  pois.mat[, "stat3-th17"] <- stat3.th17
  pois.mat[, "stat3-th0"] <- stat3.th0
  
  expected.scores.boost <- matrix(0, nrow = 3, ncol = 2)
  rownames(expected.scores.boost) <- toupper(genes.unique)
  colnames(expected.scores.boost) <- toupper(tfs.list.unique)
  
  boost.vals <- pois.mat[, "p300-th17"] - pois.mat[, "p300-th0"]
  expected.scores.boost[,"batf"] <- pois.mat[, "batf-th17"] - pois.mat[, "batf-th0"] + boost.vals
  expected.scores.boost[,"stat3"] <- pois.mat[, "stat3-th17"] - pois.mat[, "stat3-th0"] + boost.vals
  
  expected.scores.noboost <- matrix(0, nrow = 3, ncol = 2)
  rownames(expected.scores.noboost) <- toupper(genes.unique)
  colnames(expected.scores.noboost) <- toupper(tfs.list.unique)
  
  expected.scores.boost[,"batf"] <- pois.mat[, "batf-th17"] - pois.mat[, "batf-th0"]
  expected.scores.boost[,"stat3"] <- pois.mat[, "stat3-th17"] - pois.mat[, "stat3-th0"]
    
  scores.boost <- calc.chipscores(pois.mat, genes.unique, tfs.list.unique, TRUE)
  scores.noboost <- calc.chipscores(pois.mat, genes.unique, tfs.list.unique, FALSE)
  
  expect_that(scores.boost, is_identical_to(expected.scores.boost))
  expect_that(scores.noboost, is_identical_to(expected.scores.noboost))
})