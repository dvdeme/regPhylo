context("Test_Matrix.Overlap")

test_that("General Test with 1 species and two genes", {
  data(Seq.DF4) ## the first object of the list is the Species-by-gene matrix
  # Run the function.
  res = Matrix.Overlap(input = Seq.DF4[[1]], gene.Sel = c("co1", "16srrna"))
  expect_equal(dim(res[[1]])[1], 2)
  expect_equal(dim(res[[2]])[1], 2)
  expect_equal(res[[2]][1,1], 100)
})


SpbyGeneMat=rbind(as.matrix(Seq.DF4[[1]]), c("Titi_titi",0, 2, 1), c("Toto_toto", 0, 0, 4))
row.names(SpbyGeneMat)=SpbyGeneMat[,1]
SpbyGeneMat=as.data.frame(SpbyGeneMat)

test_that("General Test with 3 species and 3 genes", {
  res = Matrix.Overlap(input = SpbyGeneMat, gene.Sel = c("co1", "16srrna", "12srrna"))
  expect_equal(dim(res[[1]])[1], 3)
  expect_equal(dim(res[[2]])[1], 3)
  expect_equal(sum(diag(res[[2]])), 300)
  expect_equal(sum(diag(res[[1]])), 6)
})
