context("test-test_taxreport2sp")


src.dir = system.file("extdata/tax_export", package = "regPhylo")
dir.create("TempDir")
# Set up the path to the TempDir folder.
dest.dir = paste(getwd(), "/TempDir", sep="")
file.names <- dir(src.dir)
# Copy all the files stored in regPhylo/extdata/multi.align"
# into a temporary folder.
sapply(file.names, function(x) {
file.copy(from = paste(src.dir, x, sep = "/"),
to = paste(dest.dir, x, sep = "/"),
overwrite = FALSE) })

# Run the function using this example of tax_report.txt file
# for 30 species exported by NCBI taxonomic facility.
taxreport = Taxreport2Sp.List(input = "TempDir/tax_report.txt")
names(taxreport) # the first element of the list is the table
# for NCBI search, the second is for BOLD search.

## Table for NCBI search.
# head(taxreport$SpList.NCBI)
# dim(taxreport$SpList.NCBI) # one taxa did not have a taxid
# when the request was performed the 7/03/2019.

test_that("Test the output of Taxreport2Sp.List", {
  expect_equal(length(taxreport), 2)
  if(class(taxreport) == "list"){a=1} else {a=0}
  expect_equal(a, 1)
  expect_equal(dim(taxreport[[1]])[1], 29)
  expect_equal(dim(taxreport[[1]])[2], 2)
  if(names(taxreport)[1] == "SpList.NCBI"){b=1} else {b=0}
  expect_equal(b, 1)
  if(names(taxreport)[2] == "SpList.BOLD"){c=1} else {c=0}
  expect_equal(c, 1)
  expect_equal(dim(taxreport[[2]])[1], 32)
  expect_equal(dim(taxreport[[2]])[2], 2)
})

# Remove the Temporary folder
unlink("TempDir", recursive = TRUE)

