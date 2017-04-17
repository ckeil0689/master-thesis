# !/usr/bin/env Rscript

# Testing the creation of an interaction list from the combined confidence score matrix.
source(paste0(getwd(), "/../" , "createInteractions-fun.R"), chdir = TRUE)

context("Generation of interaction list (edge list)")

test_that("Empty data table is created as expected", {
  total.edge.num <- 4
  empty.table <- create.empty.table(total.edge.num)
  
  expect_that(colnames(empty.table), is_identical_to(c("nodeA", "interaction", "nodeB", "confidence_score")))
  expect_that(dim(empty.table), equals(c(4, total.edge.num), tolerance=0.01))
  expect_that(all(is.na(empty.table[,"nodeA"])), is_true())
  expect_that(all(is.na(empty.table[,"nodeB"])), is_true())
  expect_that(unique(empty.table[,"interaction"]) == c("neutral"), is_true())
  cs.select <- empty.table[confidence_score == 0]
  expect_that(length(cs.select) == total.edge.num, is_true())
})

test_that("Edges are selected as expected", {
    tfs <- c("BATF", "MAF", "RORC")
    genes <- c("A", "B", "C", "D", "E")
    vals <- c(1.845, 0.502, -0.992, -1.762, -1.493, # BATF
             -1.998, -1.034, 0.934, 1.856, 0.995,  # MAF
              1.650, 1.651, 0.762, 0.000, 0.000)  # RORC
    cs.mat <- matrix(vals, nrow=length(genes), ncol=length(tfs), dimnames = list(genes, tfs))
    used.cut <- GLOBAL[["cs.abs.cut"]]
    total.edge.num <- length(cs.mat[abs(cs.mat) > used.cut])
    expected.filename <- paste0("kc_activator_", used.cut, "_cs-cut_", Sys.Date(), ".csv")
    empty.table <- create.empty.table(total.edge.num)
    
    cyt.table <- select.edges(cs.mat, empty.table, used.cut, "positive_KC", "negative_KC")
    
    expect_that(dim(cyt.table), equals(dim(empty.table)))
    expect_that(length(which(is.na(cyt.table[,"nodeA"]))) == 0, is_true())
    expect_that(length(which(is.na(cyt.table[,"nodeB"]))) == 0, is_true())
    
    # activator (pos.edge = positive_KC)
    expected.table <- create.empty.table(total.edge.num)
    expected.table[, "nodeA"] <- c("BATF", "BATF", "MAF", "MAF", "RORC")
    expected.table[, "interaction"] <- c("positive_KC", "negative_KC", "negative_KC", "positive_KC", "positive_KC")
    expected.table[, "nodeB"] <- c("A", "D", "A", "D", "B")
    expected.table[, "confidence_score"] <- c(1.845, -1.762, -1.998, 1.856, 1.651)
    
    print(expected.table)
    print(cyt.table)
    expect_that(cyt.table, is_identical_to(expected.table))
})

# test_that("Writing of the interaction list works as expected", {
#   outpath <- paste0(getwd(), "/")
#   genes <- c("A", "B", "C")
#   tfs <- c("BATF", "MAF", "RORC")
#   cs.vals <- c(1.845, 0.502, -0.992, # BATF
#                -1.998, -1.034, 0.934,# MAF
#                1.650, 1.500, 0.762)  # RORC
#   cs.mat <- matrix(cs.vals, nrow=length(genes), ncol=length(tfs), dimnames = list(genes, tfs))
#   used.cut <- 1.65
#   expected.filename <- paste0("kc_activator_", used.cut, "_cs-cut_", Sys.Date(), ".csv")
#   
#   # make sure we are not verifying the creation of a file that already exists
#   expect_that(file.exists(expected.filename), is_false())
#   # write
#   write.interactions(cs.mat, outpath, "kc", "activator", used.cut, FALSE)
#   # it should now exist
#   expect_that(file.exists(expected.filename), is_true())
#   
#   # load again (CSV file!)
#   reloaded.mat <- as.matrix(read.table(expected.filename, header = TRUE, sep = ",", row.names = 1))
#   print(reloaded.mat)
#   # should be identical with original
#   expect_that(reloaded.mat, is_identical_to(cs.mat))
#   
#   # now append the existing file
#   ext.genes <- c("D", "E")
#   ext.vals <- c(-0.762, -1.493,  # BATF
#                  0.856, 0.995,   # MAF
#                  0.000, 0.000)   # RORC
#   ext.mat <- matrix(ext.vals, nrow=length(ext.genes), ncol=length(tfs), dimnames = list(ext.genes, tfs))
#   
#   # same filename input --> append
#   write.interactions(ext.mat, outpath, "kc", "activator", used.cut, TRUE)
#   # it should exist
#   expect_that(file.exists(expected.filename), is_true())
#   
#   # load again (CSV file!)
#   reloaded.mat <- as.matrix(read.table(expected.filename, header = TRUE, sep = ",", row.names = 1))
#   
#   all.genes <- c("A", "B", "C", "D", "E")
#   all.vals <- c(1.845, 0.502, -0.992, -0.762, -1.493, # BATF
#                 -1.998, -1.034, 0.934, 0.856, 0.995,  # MAF
#                 1.650, 1.500, 0.762, 0.000, 0.000)  # RORC
#   
#   all.mat <-matrix(all.vals, nrow=length(all.genes), ncol=length(tfs), dimnames = list(all.genes, tfs))
#   
#   # appended file needs to be identical to all.mat
#   expect_that(all.mat, is_identical_to(reloaded.mat))
#   
#   # delete tmp files
#   file.remove(expected.filename)
#   
# })