context("test-test_query-geome")

data(Seq.Diastocapen)
Seq.DF1 = Seq.Diastocapen$Seq.DF1

Seq.DF2 = Query.GeOMe.XY.R(input = Seq.DF1, Phylum = "Chordata", output = "Seq.DF2.txt")
dim(Seq.DF2)

test_that("Test the example, no sequence is present in GeOMe", {
  expect_equal(dim(Seq.DF2)[1], 11) # The original table is exported and includes 11 rows
  expect_equal(dim(Seq.DF2)[1], 26) # 26 columns
  expect_equal(length(grep("XY-GeOMe", Seq.DF2[,26])), 0)
})
file.remove("Seq.DF2.txt")
rm(Seq.DF2)

# Artificial dataset with 1 sequences present in GeoMe.
test.DF = matrix(c("Pomacentrus coelestis", NA, "KJ779898", rep(NA, 22), "NCBI"), ncol = 26)
colnames(test.DF) = colnames(Seq.DF1)
test.DF = as.data.frame(test.DF)

Seq.DF2 = Query.GeOMe.XY.R(input = test.DF, Phylum = "Chordata", output = "Seq.DF2.txt")

test_that("Test the example, 1 sequence is extracted from GeOMe", {
  expect_equal(dim(Seq.DF2)[1], 1) # The artificial table is exported and includes 1 row.
  expect_equal(dim(Seq.DF2)[1], 26) # 26 columns
  expect_equal(length(grep("XY-GeOMe", Seq.DF2[,26])), 1) # The info XY-GeOMe
  # is added in the column OriginDataBase
  expect_equal(Seq.DF2[,22], -26.64805)
  expect_equal(Seq.DF2[,23], 153.1827)
})

file.remove("Seq.DF2.txt")
rm(Seq.DF2)
