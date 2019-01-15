context("Test_AmMissData")

test_that("Simple example with 1 species and 100% complete edata matrix", {
  data(Seq.DF4) ## the first object of the list is the Species-by-gene matrix
  # Run the function
  res = AmMissData(input = Seq.DF4[[1]], gene.list = c("co1", "16srrna"))
  expect_equal(res, 0) # 0% of missing data expected
})

SpbyGeneMat = rbind(as.matrix(Seq.DF4[[1]]), c("Titi_titi",0, 2, 1), c("Toto_toto", 0, 0, 4))
row.names(SpbyGeneMat) = SpbyGeneMat[,1]
SpbyGeneMat = as.data.frame(SpbyGeneMat)

test_that("Simple example with 2 species and 2 genes", {
  # Run the function
  res = AmMissData(input = SpbyGeneMat[c(2:3),], gene.list = c("co1", "16srrna"))
  expect_equal(res, 25) # 25 percent of missing data expected
})

test_that("Simple example with 3 species and 2 genes", {
  # Run the function
  res = AmMissData(input = SpbyGeneMat, gene.list = c("12srrna", "16srrna"))
  expect_equal(res, 50) # 25 percent of missing data expected
})
