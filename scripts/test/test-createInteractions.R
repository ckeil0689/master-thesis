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
  cs.select <- empty.table[confidence_score == 0] # selects subset of rows where this is true --> new data.table
  expect_that(dim(cs.select), equals(c(total.edge.num, 4)))
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
  empty.table <- create.empty.table(total.edge.num)
  
  # activator (pos.edge = positive_KC)
  cyt.table <- select.edges(cs.mat, empty.table, used.cut, "positive_KC", "negative_KC")
  
  expect_that(dim(cyt.table), equals(dim(empty.table)))
  expect_that(length(which(is.na(cyt.table[,"nodeA"]))) == 0, is_true())
  expect_that(length(which(is.na(cyt.table[,"nodeB"]))) == 0, is_true())
  
  expected.table <- create.empty.table(total.edge.num)
  expected.table[, "nodeA"] <- c("BATF", "BATF", "MAF", "MAF", "RORC")
  expected.table[, "interaction"] <- c("positive_KC", "negative_KC", "negative_KC", "positive_KC", "positive_KC")
  expected.table[, "nodeB"] <- c("A", "D", "A", "D", "B")
  expected.table[, "confidence_score"] <- c(1.845, -1.762, -1.998, 1.856, 1.651)
  
  # identical?
  expect_that(cyt.table, is_identical_to(expected.table))
  
  # repressor (pos.edge = negative_KC) --> interactions should be switched
  cyt.table <- select.edges(cs.mat, empty.table, used.cut, "negative_KC", "positive_KC")
  expected.table[, "interaction"] <- c("negative_KC", "positive_KC", "positive_KC", "negative_KC", "negative_KC")
  expect_that(cyt.table, is_identical_to(expected.table))
  
  # If all matrix vals are zero
  cs.mat.zero <- matrix(0, nrow=length(genes), ncol=length(tfs), dimnames = list(genes, tfs))
  total.edge.num <- length(cs.mat[abs(cs.mat.zero) > used.cut])
  empty.table <- create.empty.table(total.edge.num)
  expect_that(is.null(empty.table), is_false()) # dont want NULL! -- zero length table
  cyt.table <- select.edges(cs.mat.zero, empty.table, used.cut, "negative_KC", "positive_KC")
  expect_that(is.null(cyt.table), is_false())
  expect_that(dim(cyt.table), equals(c(0, 4)))
})

test_that("Writing of the interaction list works as expected", {
  outpath <- paste0(getwd(), "/")
  
  # test table for writing
  edges.cyt <- 5
  cyt.table <- create.empty.table(edges.cyt)
  cyt.table[, "nodeA"] <- c("BATF", "BATF", "MAF", "MAF", "RORC")
  cyt.table[, "interaction"] <- c("positive_KC", "negative_KC", "negative_KC", "positive_KC", "positive_KC")
  cyt.table[, "nodeB"] <- c("A", "D", "A", "D", "B")
  cyt.table[, "confidence_score"] <- c(1.845, -1.762, -1.998, 1.856, 1.651)
  
  used.cut <- GLOBAL[["cs.abs.cut"]]
  expected.filename <- paste0("kc_activator_", used.cut, "_cs-cut_", Sys.Date(), ".csv")
  
  # make sure we are not verifying the creation of a file that already exists
  expect_that(file.exists(expected.filename), is_false())
  # write
  write.interactions(cyt.table, outpath, "kc", "activator", used.cut, FALSE)
  # it should now exist
  expect_that(file.exists(expected.filename), is_true())
  
  # load again (CSV file!)
  reloaded.table <- as.data.table(read.table(expected.filename, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  # should be identical with original
  expect_that(reloaded.table, is_identical_to(cyt.table))
  
  # now append the existing file
  edges.ext <- 2
  ext.table <- create.empty.table(edges.ext)
  ext.table[, "nodeA"] <- c("MAF", "RORC")
  ext.table[, "interaction"] <- c("positive_KC", "positive_KC")
  ext.table[, "nodeB"] <- c("C", "F")
  ext.table[, "confidence_score"] <- c(1.945, 1.782)
  
  # same filename input --> append
  write.interactions(ext.table, outpath, "kc", "activator", used.cut, TRUE)
  # it should exist
  expect_that(file.exists(expected.filename), is_true())
  
  # load again (CSV file!)
  reloaded.table <- as.data.table(read.table(expected.filename, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  
  all.table <- create.empty.table((edges.cyt + edges.ext))
  all.table[, "nodeA"] <- c("BATF", "BATF", "MAF", "MAF", "RORC", "MAF", "RORC")
  all.table[, "interaction"] <- c("positive_KC", "negative_KC", "negative_KC", "positive_KC", 
                                  "positive_KC", "positive_KC", "positive_KC")
  all.table[, "nodeB"] <- c("A", "D", "A", "D", "B", "C", "F")
  all.table[, "confidence_score"] <- c(1.845, -1.762, -1.998, 1.856, 1.651, 1.945, 1.782)
  
  # appended file needs to be identical to all.mat
  expect_that(reloaded.table, is_identical_to(all.table))
  
  # delete tmp files
  file.remove(expected.filename)
})

test_that("Integrated generation of interactions list works as expected (main method)", {
  outpath <- paste0(getwd(), "/")
  tfs <- c("BATF", "MAF", "RORC")
  genes <- c("A", "B", "C", "D", "E")
  vals <- c(1.845, 0.502, -0.992, -1.762, -1.493, # BATF
            -1.998, -1.034, 0.934, 1.856, 0.995,  # MAF
            1.650, 1.651, 0.762, 0.000, 0.000)  # RORC
  cs.mat <- matrix(vals, nrow=length(genes), ncol=length(tfs), dimnames = list(genes, tfs))
  
  # activator run
  expected.filename <- paste0("kc_activator_", GLOBAL[["cs.abs.cut"]], "_cs-cut_", Sys.Date(), ".csv")
  expect_that(file.exists(expected.filename), is_false())
  create.interactions(cs.mat, outpath, "kc", "activator", pos.edge = "positive_KC", neg.edge = "negative_KC", append = FALSE)
  expect_that(file.exists(expected.filename), is_true())
  
  # repressor run (data remains the same in this example, but edge types are swapped)
  expected.filename.r <- paste0("kc_repressor_", GLOBAL[["cs.abs.cut"]], "_cs-cut_", Sys.Date(), ".csv")
  expect_that(file.exists(expected.filename.r), is_false())
  create.interactions(cs.mat, outpath, "kc", "repressor", pos.edge = "negative_KC", neg.edge = "positive_KC", append = FALSE)
  expect_that(file.exists(expected.filename.r), is_true())
  
  # append repressor run (full file)
  expected.filename.full <- paste0("kc_single_", GLOBAL[["cs.abs.cut"]], "_cs-cut_", Sys.Date(), ".csv")
  expect_that(file.exists(expected.filename.full), is_false())
  # create (for example activator -- data same in this case, but edge types swapped)
  create.interactions(cs.mat, outpath, "kc", "single", pos.edge = "positive_KC", neg.edge = "negative_KC", append = FALSE)
  # then append the other type (for example repressor)
  create.interactions(cs.mat, outpath, "kc", "single", pos.edge = "negative_KC", neg.edge = "positive_KC", append = TRUE)
  expect_that(file.exists(expected.filename.full), is_true())
  
  # Load written interaction table and check some attributes like dimension to confirm that the table is as expected (CSV file!)
  reloaded.table <- as.data.table(read.table(expected.filename.full, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  interactions.num.full <- length(which(abs(vals) > GLOBAL[["cs.abs.cut"]])) * 2 # doubled because cs.mat is used for activator and repressor
  expect_that(dim(reloaded.table), equals(c(interactions.num.full, 4))) # 4 columns are fixed by definition (nodeA, interaction, nodeB, cs)
  expect_that(all(is.na(reloaded.table[,"nodeA"])), is_false())
  expect_that(all(is.na(reloaded.table[,"nodeB"])), is_false())
  expect_that("neutral" %in% unique(reloaded.table[,"interaction"]), is_false())
  cs.select <- reloaded.table[confidence_score == 0] # selects subset of rows where this is true --> new data.table
  expect_that(dim(cs.select), equals(c(0, 4)))
  
  # delete tmp files
  file.remove(expected.filename)
  file.remove(expected.filename.r)
  file.remove(expected.filename.full)
})