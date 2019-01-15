context("Test_First.Align.All")


# To run the example it might be better to copy the input alignment files
# provided by the package to a temporary directory created into the
# current working directory.
src.dir = system.file("extdata/FirstToAlign", package = "regPhylo")
dir.create("TempDir.FirstToAlign")
# Set up the path to the TempDir folder.
dest.dir = paste(getwd(), "/TempDir.FirstToAlign", sep="")
file.names <- dir(src.dir)
# Copy all the files stored in regPhylo/extdata/FirstToAlign"
# into a temporary folder.
sapply(file.names, function(x) {
  file.copy(from = paste(src.dir, x, sep = "/"),
            to = paste(dest.dir, x, sep = "/"),
            overwrite = FALSE) })

# Remove the empty file created for the folder name when
# importing the data file
file.remove("TempDir.FirstToAlign/FirstAligned")

input = "TempDir.FirstToAlign"
#dir.create("TempDir.FirstToAlign/FirstAligned")
output = "TempDir.FirstToAlign/FirstAligned"
nthread = 1
methods = c("mafftfftnsi", "pasta")
First.Align.All(input = input, output = output, nthread = 1, methods = c("mafftfftnsi"))


test_that("Test number output alignment", {
  expect_equal(length(list.files("TempDir.FirstToAlign/FirstAligned")), 6) # number of files
  a = list.files("TempDir.FirstToAlign/FirstAligned")
  i = 1
  for(i in 1:length(a)){
    a1 =ape::read.dna(paste("TempDir.FirstToAlign/FirstAligned/", a[i], sep=""), format = "fasta")
    if(classs(a1) == "DNAbin"){a = 1} else {a = 0}
    expect_equal(a,1)
  }

})


# To clean the file created while running the example do the following:
unlink("TempDir.FirstToAlign", recursive = TRUE)
