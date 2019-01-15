context("Test_SpeciesGeneMat.Bl")

data(Seq.Diastocapen)
Seq.DF3 = Seq.Diastocapen$Seq.DF3

# Run the function
Seq.DF4=SpeciesGeneMat.Bl(input=Seq.DF3, output="Seq.DF4.")

test_that("Test the Species-by-gene matrix", {
  expect_equal(length(Seq.DF4), 4) # 4 objects exported in R envir
  expect_equal(dim(Seq.DF4[[1]])[1], 1)
  expect_equal(dim(Seq.DF4[[1]])[2], 4)
  if(Seq.DF4[[1]][1,1] == "Diastobranchus capensis") {a = 1} else {a = 0}
  expect_equal(a, 1)
})

test_that("Test the Species-by-gene matrix exported in the WD", {
  Seq.DF4WD = read.delim("Seq.DF4.SpDNA_Mat.txt", sep = "\t", header = TRUE)
  expect_equal(dim(Seq.DF4WD)[1], 1)
  expect_equal(dim(Seq.DF4WD)[2], 4)
  if(Seq.DF4WD[1,1] == "Diastobranchus capensis") {a = 1} else {a = 0}
  expect_equal(a, 1)
})

test_that("Test the DNA perspective table exported in the WD", {
  expect_equal(dim(Seq.DF4$Summary_DNA)[1], 3)
  expect_equal(dim(Seq.DF4$Summary_DNA)[2], 3)
  Seq.DF4DNA = read.delim("Seq.DF4.Summary_DNApers.txt", sep = "\t", header = TRUE)
  expect_equal(dim(Seq.DF4DNA)[1], 3)
  expect_equal(dim(Seq.DF4DNA)[2], 3)
  expect_equal(sum(Seq.DF4DNA[,2]), 3) ## test if the sum of the second column equal 3 (Species coverage of 3)
})

test_that("Test the Species perspective table exported in the WD", {
  expect_equal(dim(Seq.DF4$Summary_Species)[1], 1)
  expect_equal(dim(Seq.DF4$Summary_Species)[2], 3)
  Seq.DF4Sp = read.delim("Seq.DF4.Summary_SPECIESpers.txt", sep = "\t", header = TRUE)
  expect_equal(dim(Seq.DF4Sp)[1], 1)
  expect_equal(dim(Seq.DF4Sp)[2], 3)
  expect_equal(Seq.DF4Sp[,3], 11) ## test if 11 sequences are retrieved for Diasto. cap.
})


test_that("Test the Clean Data Table exported in the WD", {
  Seq.DF4cl = read.delim("Seq.DF4._CleanDataset.txt", sep = "\t", header = TRUE)
  expect_equal(dim(Seq.DF4cl)[1], 11) ### eleven sequences
  expect_equal(dim(Seq.DF4cl)[2], 30) ### 30 columns
  # test if the additional columns is called "ProductGeneClean"
  if(which(colnames(Seq.DF4cl) == "ProductGeneClean") == 30) {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(unique(Seq.DF4cl$ProductGeneClean)), 3)
})

file.remove("Seq.DF4._CleanDataset.txt")
file.remove("Seq.DF4.Summary_SPECIESpers.txt")
file.remove("Seq.DF4.Summary_DNApers.txt")
file.remove("Seq.DF4.SpDNA_Mat.txt")


Seq.DF4b=SpeciesGeneMat.Bl(input=Seq.DF3, output="Seq.DF4.Trash", NCBI.Trash = c("JX242943", "JN640620"),
                          BOLD.Trash = c("2754997"))

test_that("Test the NCBI and BOLD Blacklists", {
  expect_equal(length(Seq.DF4b), 4) # 4 objects exported in R envir
  expect_equal(dim(Seq.DF4b[[1]])[1], 1)
  expect_equal(dim(Seq.DF4b[[1]])[2], 3)
  if(Seq.DF4[[1]][1,1] == "Diastobranchus capensis") {a = 1} else {a = 0}
  expect_equal(a, 1)
})

test_that("Test the Clean data table exported in the WD using the NCBI and BOLD Blacklists", {
  Seq.DF4Cl2 = read.delim("Seq.DF4.Trash_CleanDataset.txt", sep = "\t", header = TRUE)
  dim(Seq.DF4Cl2)
  expect_equal(dim(Seq.DF4Cl2)[1], 8) #  only 8 sequences expected
  expect_equal(dim(Seq.DF4Cl2)[2], 30)
  expect_equal(length(unique(Seq.DF4Cl2$ProductGeneClean)), 2) # 2 type of genes regions are expected.
})


file.remove("Seq.DF4.Trash_CleanDataset.txt")
file.remove("Seq.DF4.TrashSummary_SPECIESpers.txt")
file.remove("Seq.DF4.TrashSummary_DNApers.txt")
file.remove("Seq.DF4.TrashSpDNA_Mat.txt")
