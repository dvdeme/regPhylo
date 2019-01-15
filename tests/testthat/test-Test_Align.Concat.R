context("Test_Align.Concat")

# To run the example it might be better to copy the input alignment files
# provided by the package to a temporary directory created into the
# current working directory.
src.dir = system.file("extdata/multi.align/ForConcat", package = "regPhylo")
dir.create("TempDir.ForConcat")
# Set up the path of the TempDir folder.
dest.dir = paste(getwd(), "/TempDir.ForConcat", sep="")
file.names <- dir(src.dir)
# Copy all the files stored in regPhylo/extdata/multi.align/ForConcat"
# into a temporary folder.
sapply(file.names, function(x) {
file.copy(from = paste(src.dir, x, sep = "/"),
to = paste(dest.dir, x, sep = "/"),
overwrite = FALSE) })

# Run the function to build the supermatrix.
a = Align.Concat(input = "TempDir.ForConcat", Sp.List = NULL, outputConcat = NULL)

# Import the supermatrix in R
Supermatrix = ape::read.dna("TempDir.ForConcat/Concat.fas", format = "fasta")


test_that("Test the export of the output table and Supematrix aligmment in fasta and nexus format, No additional Species without DNA", {
  expect_equal(dim(a)[1], 4)
  expect_equal(dim(a)[2], 2)
  # test if the supermatrix is a DNAbin object
  if(class(Supermatrix) == "DNAbin"){a = 1} else {a = 0}
  expect_equal(a, 1)
  # test the size of the supermatrix in fasta
  expect_equal(dim(Supermatrix)[1], 16)
  expect_equal(dim(Supermatrix)[2], 3124)
  test1 = ape::read.nexus.data("TempDir.ForConcat/Concat.nex")
  test2 = ape::as.DNAbin(test1)
  # test if the nexus file could have been loaded and has been
  # properly been converted into a DNAbin object.
  if(class(test2) == "DNAbin"){b = 1} else {b = 0}
  expect_equal(b, 1)
  # test if the size of the supermatrix is identical between the fasta and nexus export
  expect_equal(length(labels(test2)), dim(Supermatrix)[1])
  expect_equal(length(test2[[1]]), dim(Supermatrix)[2])
  # test the presence partition file for RAxML
  part = readLines("TempDir.ForConcat/Partitions_Concat.txt")
  expect_equal(length(part), 4)

  # test the presence and the size of the convtab table.
  contab = read.delim("TempDir.ForConcat/convtab.txt", sep = "\t", header = TRUE)
  expect_equal(dim(contab)[1], 4)
  expect_equal(dim(contab)[2], 2)
  # test the number of files in the temporary directory
  a = list.files("TempDir.ForConcat")
  expect_equal(length(a), 12)

  # test if "Gblocksls_PASTA_Alig_12srrna.nex" alignement include 16 species
  test3 = readLines("TempDir.ForConcat/Gblocksls_PASTA_Alig_cytb.nex")
  # we use ReadeLines, because in R ape::read.nexus.data does not allow "|" in the seqeunce name.
  nbtaxa = as.numeric(strsplit(strsplit(test3[grep("DIMENSIONS NTAX=", test3)], " ")[[1]][2], "=")[[1]][2])
  expect_equal(nbtaxa, 16)
})

# Clean the environment before running thes econd example
rm(a, Supermatrix)


# Create another temporary file to build a supermatrix including species without DNA.
dir.create("TempDir.ForConcat2")
file.names <- dir("TempDir.ForConcat")
# select only the .fas alignment
file.names <- file.names[grep(".fas", file.names)]
 # remove the Concat.fas alignment just created above.
file.names <- file.names[-grep("Concat.fas", file.names)]
sapply(file.names, function(x) {
file.copy(from = paste("TempDir.ForConcat", x, sep = "/"),
to = paste("TempDir.ForConcat2", x, sep = "/"),
overwrite = FALSE) })


# Run the function to build a supermatrix including two species without DNA.
a = Align.Concat(input = "TempDir.ForConcat2",
Sp.List = c("Titi_titi", "Toto_toto"),
outputConcat = "TempDir.ForConcat2/Concat_2spNoDNA")


# Import the supermatrix in R
Supermatrix = ape::read.dna("TempDir.ForConcat2/Concat_2spNoDNA.fas", format = "fasta")


test_that("Test the export of the output table and Supematrix aligmment in fasta and nexus format, with 2 additional Species without DNA", {
  expect_equal(dim(a)[1], 4)
  expect_equal(dim(a)[2], 2)
  # test if the supermatrix is a DNAbin object
  if(class(Supermatrix) == "DNAbin"){a = 1} else {a = 0}
  expect_equal(a, 1)
  # test the size of the supermatrix in fasta
  expect_equal(dim(Supermatrix)[1], 18)
  expect_equal(dim(Supermatrix)[2], 3124)
  test1 = ape::read.nexus.data("TempDir.ForConcat2/Concat_2spNoDNA.nex")
  test2 = ape::as.DNAbin(test1)
  # test if the nexus file could have been loaded and has been
  # properly been converted into a DNAbin object.
  if(class(test2) == "DNAbin"){b = 1} else {b = 0}
  expect_equal(b, 1)
  # test if the size of the supermatrix is identical between the fasta and nexus export
  expect_equal(length(labels(test2)), dim(Supermatrix)[1])
  expect_equal(length(test2[[1]]), dim(Supermatrix)[2])
  # test the presence partition file for RAxML
  part = readLines("TempDir.ForConcat/Partitions_Concat.txt")
  expect_equal(length(part), 4)

  # test the presence and the size of the convtab table.
  contab = read.delim("TempDir.ForConcat2/convtab.txt", sep = "\t", header = TRUE)
  expect_equal(dim(contab)[1], 4)
  expect_equal(dim(contab)[2], 2)
  # test the number of files in the temporary directory
  a = list.files("TempDir.ForConcat2")
  expect_equal(length(a), 12)

  # test if "Gblocksls_PASTA_Alig_12srrna.nex" alignement include 16 species
  test3 = readLines("TempDir.ForConcat2/Gblocksls_PASTA_Alig_cytb.nex")
  # we use ReadeLines, because in R ape::read.nexus.data does not allow "|" in the seqeunce name.
  nbtaxa = as.numeric(strsplit(strsplit(test3[grep("DIMENSIONS NTAX=", test3)], " ")[[1]][2], "=")[[1]][2])
  # nbtaxa must be 18
  expect_equal(nbtaxa, 18)
})

# Clean the environment before running thes econd example
rm(a, Supermatrix)
unlink("TempDir.ForConcat", recursive = TRUE)
unlink("TempDir.ForConcat2", recursive = TRUE)


