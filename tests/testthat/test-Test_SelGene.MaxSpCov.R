context("Test_SelGene.MaxSpCov")

data(Seq.DF4) ## the first object of the list is the Species-by-gene matrix


test_that("Test simple function with original example", {
  res = SelGene.MaxSpCov(input = Seq.DF4[[1]])
  if(length(grep("There is only one species", res[[1]])) == 1){a= 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(dim(res[[2]])[2], 4)
})

SpbyGeneMat=rbind(as.matrix(Seq.DF4[[1]]), c("Titi_titi",0, 2, 1), c("Toto_toto", 0, 0, 4))
row.names(SpbyGeneMat)=SpbyGeneMat[,1]
SpbyGeneMat=as.data.frame(SpbyGeneMat)

test_that("Test simple function with 3 species example Gene name selection, 2genes", {
  res = SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = c("co1", "12srrna"))
  expect_equal(res[[1]], 100) # 100% species coverage
  expect_equal(dim(res$GeneContributionCoverage)[1], 2) #
  if(c(which(res[[4]] == "co1") + which(res[[4]] == "12srrna")) == 3) {a = 1} else {a = 0}
  expect_equal(a,1)
})

test_that("Test simple function with 3 species example, Gene name selection 1 gene", {
  res = SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = c("12srrna"))
  if(length(grep("12srrna covers", res[[1]])) == 1) {a = 1} else {a = 0}
  expect_equal(a,1)
  if(length(grep("33.33 %", res[[1]])) == 1) {b = 1} else {b = 0}
  expect_equal(b,1)
  expect_equal(dim(res[[2]])[1], 3)
})

test_that("Test simple function with 3 species example, Gene Number, 2", {
  res = SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = 2)
  expect_equal(res[[1]],2)
  expect_equal(res[[4]],3)
  if(length(setdiff(c("co1", "12srrna"), res[[2]])) == 0) {a = 1} else {a = 0}
})

test_that("Test simple function with 3 species example, Gene Number, 1", {
  res = SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = 1)
  if(length(grep("The co1 equally", res[[1]])) == 1) {a = 1} else {a = 0}
  expect_equal(a,1)
  if(length(grep("covers 100 % of species", res[[2]])) == 1) {b = 1} else {b = 0}
  expect_equal(b,1)
  expect_equal(sum(res[[3]][,1]), 3) ### one gene covers 100% of the 3 species.
})

file.remove("Rplots.pdf")
