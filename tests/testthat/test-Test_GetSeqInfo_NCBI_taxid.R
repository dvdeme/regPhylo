context("Test_GetSeqInfo_NCBI_taxid")

test_that("Test if the function works with a positive control", {
  Splist=cbind(TaxID=c(443778), Species.Name=c("Diastobranchus capensis"))
  NCBI.output=GetSeqInfo_NCBI_taxid(splist=Splist, gene="ALL", filename="output.NCBI.txt")
  expect_equal(class(NCBI.output), "matrix")
  Output1=read.delim("output.NCBI.txt", sep="\t", h=TRUE)
  expect_equal(class(Output1), "data.frame")
  expect_equal(dim(Output1)[2], 23)
  expect_equal(dim(Output1)[1], as.vector(as.numeric(as.character(NCBI.output[1,2])))) ### test if the number of sequences retrieved
  ## is equal to the number of seqeunce printed in the summary table.
  file.remove(c("output.NCBI.txt", "data"))
})

test_that("Test if the function works with a negative control", {
  Splist=cbind(TaxID=c(0000), Species.Name=c("Titi toto"))
  NCBI.output2=GetSeqInfo_NCBI_taxid(splist=Splist, gene="ALL", filename="output2.NCBI.txt")
  expect_equal(class(NCBI.output2), "matrix")
  Output2=read.delim("output2.NCBI.txt", sep="\t", h=TRUE)
  expect_equal(class(Output2), "data.frame")
  expect_equal(dim(Output2)[2], 23)
  expect_equal(dim(Output2)[1], as.vector(as.numeric(as.character(NCBI.output2[1,2])))) ### test if the number of sequences retrieved
  ## is equal to the number of sequences printed in the summary table.
  expect_equal(dim(Output2)[1],0)
  file.remove(c("output2.NCBI.txt", "data"))
})

