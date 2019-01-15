context("Test_Select.DNA")

data(Seq.DF4) # the clean table of sequences and metadata is called "CleanDataTable"
# and is the fifth object of the list "Seq.DF4".
# Run the function
Seq.DF5=Select.DNA(input = Seq.DF4$CleanDataTable, gene.list = c("co1", "16srrna"),
                   output = "Seq.DF4.dataTable")

test_that("Simple test that the output table in R envir has the good dimensions", {
  expect_equal(dim(Seq.DF5)[1], 10)
  expect_equal(dim(Seq.DF5)[2], 29)
  if(which(colnames(Seq.DF5) == "ProductGeneClean") == 29) {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(unique(Seq.DF5$ProductGeneClean)), 2)
})

Seq.DF5b = read.delim("Seq.DF4.dataTable.Select.DNA.txt", sep = "\t" , header = TRUE)

test_that("Simple test that the output table in WD has the good dimensions", {
  expect_equal(dim(Seq.DF5b)[1], 10)
  expect_equal(dim(Seq.DF5b)[2], 29)
  if(which(colnames(Seq.DF5b) == "ProductGeneClean") == 29) {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(unique(Seq.DF5b$ProductGeneClean)), 2)
})

### Create the artificial dataset
test = as.matrix(Seq.DF4$CleanDataTable)
TestTable = rbind(test, c("Titi_titi", rep("NA", 28), "co1"), c("Titi_titi", rep("NA", 28), "16srrna"),
                  c("Toto_toto", rep("NA", 28), "12srrna"),
                  c(test[11,c(1:3)], "Too_Long", test[11,c(5:30)]))
TestTable = as.data.frame(TestTable)

# run the function with the artificial dataset.
Seq.DF5c=Select.DNA(input = TestTable, gene.list = c("co1", "16srrna"),
                   output = "Seq.DF4c.dataTable")

test_that("Test the extract of a seqeunce using seqinr and the output in th R envir", {
  expect_equal(dim(Seq.DF5c)[1], 13)
  expect_equal(dim(Seq.DF5c)[2], 29)
  if(which(colnames(Seq.DF5c) == "ProductGeneClean") == 29) {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(unique(Seq.DF5c$ProductGeneClean)), 2)
  # Test the removal of the Too_long and replacement by the TRUE sequence
  if(length(grep("Too_Long", Seq.DF5c[,4])) == 0) {a = 1} else {a = 0}
  expect_equal(a, 1)
  # Test if the sequences with the accession number "JN640620" (must be preent twice in the output) is present twice.
  expect_equal(length(Seq.DF5c[grep(as.vector(test[11,3]), Seq.DF5c[,3]),4]), 2)
  # Test if the two seqeunces are identical
  expect_equal(length(unique(tolower(Seq.DF5c[grep(as.vector(test[11,3]), Seq.DF5c[,3]),4]))), 1)
})

Seq.DF5d = read.delim("Seq.DF4c.dataTable.Select.DNA.txt", sep = "\t" , header = TRUE)

test_that("Test the extract of a sequence using seqinr and the output in the WD", {
  expect_equal(dim(Seq.DF5d)[1], 13)
  expect_equal(dim(Seq.DF5d)[2], 29)
  if(which(colnames(Seq.DF5d) == "ProductGeneClean") == 29) {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(unique(Seq.DF5d$ProductGeneClean)), 2)
  # Test the removal of the Too_long and replacement by the TRUE sequence
  if(length(grep("Too_Long", Seq.DF5d[,4])) == 0) {a = 1} else {a = 0}
  expect_equal(a, 1)
  # Test if the sequences with the accession number "JN640620" (must be preent twice in the output) is present twice.
  expect_equal(length(Seq.DF5d[grep(as.vector(test[11,3]), Seq.DF5d[,3]),4]), 2)
  # Test if the two seqeunces are identical
  expect_equal(length(unique(tolower(Seq.DF5d[grep(as.vector(test[11,3]), Seq.DF5d[,3]),4]))), 1)
})

rm(Seq.DF5d, Seq.DF5b, Seq.DF5c, Seq.DF5, test, TestTable)

file.remove("Seq.DF4.dataTable.Select.DNA.txt", "Seq.DF4c.dataTable.Select.DNA.txt")
