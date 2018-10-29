context("Test_GetSeq_BOLD")

test_that("Test GetSeq_BOLD with positive control", {
  Splist=cbind(TaxID=c("Diastobranchus capensis"), Species.Name=c("Diastobranchus capensis"))
  BOLD.output=GetSeq_BOLD(splist=Splist, filename="output.BOLD.txt")
  expect_equal(class(BOLD.output), "list") ### test if the code of the function is not broken and export the proper object in the R environment
  expect_equal(length(BOLD.output), 2) ### similar than previous.
  Output1=read.delim("output.BOLD.txt", sep="\t", header=TRUE)
  expect_equal(dim(Output1)[2], 82) ### number of columns of the output in WD (test if the number of column remains
  # the same, or if bold changes the format of their output)
  if(dim(Output1)[1]<1){a=0} else {a=1} ### number of rows of the output in WD (test if the output in the WD contain the information)
  expect_equal(a, 1)
  file.remove("output.BOLD.txt")
})

test_that("Test GetSeq_BOLD with negative control", {
  Splist=cbind(TaxID=c("Titi toto"), Species.Name=c("Titi toto"))
  BOLD.output2=GetSeq_BOLD(splist=Splist, filename="output2.BOLD.txt")
  Output2=read.delim("output2.BOLD.txt", sep="\t", header=TRUE)
  if(dim(Output2)[1]<1){a=0} else {a=1}
  expect_equal(a, 0)
  expect_equal(dim(BOLD.output2[[2]])[1], 0)
  file.remove("output2.BOLD.txt")
})
