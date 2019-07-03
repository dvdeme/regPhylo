context("Test_Filtering.align.Gblocks")


# To run the example we have to copy the input alignment files
# provided by the package to a temporary directory created into the
# current working directory.
src.dir = system.file("extdata/multi.align/multi.aligned", package = "regPhylo")
dir.create("TempDir")
# Set up the path to the TempDir folder.
dest.dir = paste(getwd(), "/TempDir", sep="")
file.names <- dir(src.dir)
# Copy all the files stored in regPhylo/extdata/multi.align/multi.aligned"
# into a temporary folder.
sapply(file.names, function(x) {
file.copy(from = paste(src.dir, x, sep = "/"),
to = paste(dest.dir, x, sep = "/"),
overwrite = FALSE) })

# Run the function from the TempDir folder and store the outputs from
# Gblocks in the "Trimmed-Gblocks" folder.
a = Filtering.align.Gblocks(input = "TempDir", LessStringent = TRUE,
                            output = "TrimmedGblocks", Type = "d", remove.empty.align = TRUE)


b =  list.files("TrimmedGblocks")
filespace = vector()
ListDNAbin = list(1:12)
classAlign = vector()
i =1
for (i in 1:length(b)){
  filespace = c(filespace, file.info(paste("TrimmedGblocks/", b[i], sep = ""))$size)
  ListDNAbin[[i]] = ape::read.dna(paste("TrimmedGblocks/", b[i], sep = ""), format = "fasta")
  classAlign = c(classAlign, class(ListDNAbin[[i]]))
}



test_that("Test if the function run until the end, less stringent", {
  # test if all the intermediate files have been removed
  aa = list.files("TempDir")
  expect_equal(length(aa), 12)
  # Test if the output of the function in the R environment is a table with 12 rows
  expect_equal(dim(a)[1], 12)
  # Test if the output of the function in the R environment compute the sequence length.
  expect_equal(as.numeric(as.character(a[1, 5])), 463)
  expect_equal(as.numeric(as.character(a[12, 5])), 1086)
  # Test if the number of alignements exported in the output directory is correct must be 12
  expect_equal(length(b), 12)
  # Test if all the files have some content.
  expect_equal(length(which(filespace > 0)), 12)
  # Test that all the file can be read as DNAbin object
  expect_equal(length(which(classAlign == "DNAbin")), 12)
})

test_that("Test the length in pb of the first file, less stringent", {
  test = ape::read.dna(paste("TempDir/", "Mafftfftns1_Alig_12srrna.fas", sep = ""), format = "fasta")
  test1 = ape::read.dna(paste("TrimmedGblocks/", "Gblocksls_Mafftfftns1_Alig_12srrna.fas", sep = ""), format = "fasta")
  expect_equal(dim(test)[2], 483)
  expect_equal(dim(test1)[2], 463)
})

# Remove the folder with Gblocks outputs
unlink("TrimmedGblocks", recursive = TRUE)

rm(a, b, filespace, ListDNAbin, classAlign)

# Run the function from the TempDir folder and store the outputs from
# Gblocks in the "Trimmed-Gblocks" folder.
a = Filtering.align.Gblocks(input = "TempDir", LessStringent = FALSE,
                            output = "TrimmedGblocks", Type = "d", remove.empty.align = TRUE)


b =  list.files("TrimmedGblocks")
filespace = vector()
ListDNAbin = list(1:12)
classAlign = vector()
i =1
for (i in 1:length(b)){
  filespace = c(filespace, file.info(paste("TrimmedGblocks/", b[i], sep = ""))$size)
  ListDNAbin[[i]] = ape::read.dna(paste("TrimmedGblocks/", b[i], sep = ""), format = "fasta")
  classAlign = c(classAlign, class(ListDNAbin[[i]]))
}


test_that("Test if the function run until the end, more stringent", {
  # test if all the intermediate files have been removed
  aa = list.files("TempDir")
  expect_equal(length(aa), 12)
  # Test if the output of the function in the R environment is a table with 12 rows
  expect_equal(dim(a)[1], 12)
  # Test if the output of the function in the R environment compute the sequence length.
  expect_equal(as.numeric(as.character(a[1, 5])), 122)
  expect_equal(as.numeric(as.character(a[12, 5])), 615)
  # Test if the number of alignements exported in the output directory is correct must be 12
  expect_equal(length(b), 12)
  # Test if all the files have some content.
  expect_equal(length(which(filespace > 0)), 12)
  # Test that all the file can be read as DNAbin object
  expect_equal(length(which(classAlign == "DNAbin")), 12)
})

test_that("Test the length in pb of the first file, more stringent", {
  test = ape::read.dna(paste("TempDir/", "Mafftfftns1_Alig_12srrna.fas", sep = ""), format = "fasta")
  test1 = ape::read.dna(paste("TrimmedGblocks/", "Gblocksms_Mafftfftns1_Alig_12srrna.fas", sep = ""), format = "fasta")
  expect_equal(dim(test)[2], 483)
  expect_equal(dim(test1)[2], 122)
})

# Remove the folder with Gblocks outputs
unlink("TrimmedGblocks", recursive = TRUE)
rm(a, b, filespace, ListDNAbin, classAlign)

# Test the target.file option.
a = Filtering.align.Gblocks(input = "TempDir", target.file = "Mafftfftns1_Alig_co1.fas",
    LessStringent = TRUE, output = "TrimmedGblocks",
    Type = "d", remove.empty.align = TRUE)

b =  list.files("TrimmedGblocks")

test_that("Test if only one alignment is present in th output folder", {
  b =  list.files("TrimmedGblocks")
  expect_equal(length(b), 1)
})

test_that("Test the length in pb of the first file, more stringent", {
  test = ape::read.dna(paste("TempDir/", "Mafftfftns1_Alig_co1.fas", sep = ""), format = "fasta")
  test1 = ape::read.dna(paste("TrimmedGblocks/", "Gblocksls_Mafftfftns1_Alig_co1.fas", sep = ""), format = "fasta")
  expect_equal(dim(test)[2], 1551)
  expect_equal(dim(test1)[2], 650)
  expect_equal(as.numeric(as.character(a[5,1])), 650)
})


# Remove temporary folder
unlink("TempDir", recursive = TRUE)
# Remove the folder with Gblocks outputs
unlink("TrimmedGblocks", recursive = TRUE)

