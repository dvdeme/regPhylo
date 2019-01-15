context("Test_Detect.Outlier.Seq")

data(Example_16S_outlier)
Example_16S_outlier_align = Example_16S_outlier[[1]]

# Running the function with the 'Comb' option,
# a distance threshold of 0.6, and disabling the search for secondary outlier sequences.
S16_MisAlign0.6_1 = Detect.Outlier.Seq(inputal = Example_16S_outlier_align,
Strat.DistMat = "Comb", Dist.Th = 0.6, output = "Example_16S_outliers_1.txt",
Second.Outlier = "No")

test_that("Test simple example with primary outlier search only", {
  expect_equal(dim(S16_MisAlign0.6_1)[1], 36)
  expect_equal(dim(S16_MisAlign0.6_1)[2], 2)
})

test=read.delim("Example_16S_outliers_1.txt", sep = "\t", h =T)

test_that("Test simple example with primary oultier search, test output table in WD", {
  expect_equal(dim(test)[1], 36)
  expect_equal(dim(test)[2], 2)
})


# Running the function with the 'Comb' option, a distance threshold of 0.6,
# and allowing a search for secondary outlier sequences using local blast database.
S16_MisAlign0.6_2 = Detect.Outlier.Seq(inputal = Example_16S_outlier_align,
Strat.DistMat = "Comb", Dist.Th = 0.6, output = "Example_16S_outliers_2.txt",
Second.Outlier = "Yes", Bitsc.Th = 0.8)

test_that("Test simple example with secondary oultier search", {
  expect_equal(dim(S16_MisAlign0.6_2)[1], 85)
  expect_equal(dim(S16_MisAlign0.6_2)[2], 16)
  expect_equal(length(unique(S16_MisAlign0.6_2[,3])), 38) ##nb of unique potential outliers
})

test2=read.delim("Example_16S_outliers_2.txt", sep = "\t", h =T)

test_that("Test simple example with secondary oultier search, test output table in WD", {
  expect_equal(dim(test2)[1], 85)
  expect_equal(dim(test2)[2], 16)
  expect_equal(length(unique(test2[,3])), 38) ##nb of unique potential outliers
})


file.remove("Example_16S_outliers_1.txt", "Example_16S_outliers_2.txt")


