context("Test_Multi.Align")

# To run the example it might be better to copy the input alignment files
# provided by the package to a temporary directory created into the
# current working directory.
src.dir = system.file("extdata/multi.align", package = "regPhylo")
dir.create("TempDir.Multi.aligned")
# Set up the path to the TempDir folder.
dest.dir = paste(getwd(), "/TempDir.Multi.aligned", sep="")
file.names <- dir(src.dir)
# Copy all the files stored in regPhylo/extdata/multi.align"
# into a temporary folder.
sapply(file.names, function(x) {
file.copy(from = paste(src.dir, x, sep = "/"),
to = paste(dest.dir, x, sep = "/"),
overwrite = FALSE) })

output = "multi.alignedTest"
nthread = 3
methods = c("mafftfftnsi", "pasta")
input = "TempDir.Multi.aligned"
# run the function using two additional alignmenet program "mafftfftnsi" and "pasta"
Multi.Align(input = input, output = "multi.alignedTest", nthread = 3,
methods = c("mafftfftnsi", "pasta"))
test = list.files("multi.alignedTest")


test_that("multiplication works", {
  expect_equal(length(test), 12)
  expect_equal(length(unique(unlist(lapply(strsplit(test, "_"), function(x) x[1])))), 3)
})
