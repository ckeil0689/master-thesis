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
  print(reloaded.table)
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
  
})