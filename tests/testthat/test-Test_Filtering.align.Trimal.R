context("Test_Filtering.align.Trimal")

src.dir = system.file("extdata/multi.align/multi.aligned", package = "regPhylo")
dir.create("TempDir.ToTrim")
# Set up the path to the TempDir folder.
dest.dir = paste(getwd(), "/TempDir.ToTrim", sep="")
file.names <- dir(src.dir)
# Copy all the files stored in regPhylo/extdata/multi.align"
# into a temporary folder.
sapply(file.names, function(x) {
file.copy(from = paste(src.dir, x, sep = "/"),
to = paste(dest.dir, x, sep = "/"),
overwrite = FALSE) })

input = "TempDir.ToTrim"
output = "Trimmed"

a = Filtering.align.Trimal(input = input, output = "Trimmed")

b = list.files("Trimmed")
length(b)

test_that("Test the function runs and exports the output", {
  # test the dimension of the output table
  expect_equal(dim(a)[1], 24)
  expect_equal(dim(a)[2], 5)
  # test the number of files generated
  expect_equal(length(b), 24)
  # test the length in pb of the Mafftfftns1_Alig_12srrna.fas after timming using trimAuto must be 122
  test1 = ape::read.dna(paste("Trimmed/", "trimAuto_Mafftfftns1_Alig_12srrna.fas", sep = ""), format = "fasta")
  expect_equal(dim(test1)[2], 122)
})

# removing temporary folder
unlink("TempDir.ToTrim", recursive = TRUE)
# removing temporary folder
unlink("Trimmed", recursive = TRUE)

