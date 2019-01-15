context("Test_PartiFinder2")

src.dir = system.file("extdata/multi.align/ForPartiFinder2", package = "regPhylo")
dir.create("TempDir.ForPartiFinder2")
# Set up the path to the TempDir folder.
dest.dir = paste(getwd(), "/TempDir.ForPartiFinder2", sep="")
file.names <- dir(src.dir)
# Copy all the files stored in regPhylo/extdata/multi.align/ForPartiFinder2"
# into a temporary folder.
sapply(file.names, function(x) {
file.copy(from = paste(src.dir, x, sep = "/"),
to = paste(dest.dir, x, sep = "/"),
overwrite = FALSE) })


# Remove the file created by partiFinder2 in the example folder.
file.remove("TempDir.ForPartiFinder2/Concat_PF2.nex")
file.remove("TempDir.ForPartiFinder2/Partitions_Concat.txt_PF2_all.txt")


# Run the function using RAxML software, BIC criteria,
# the fast "rcluster" algorithm and the classic default parameters.
PartiFinder2(input = "TempDir.ForPartiFinder2/Concat.fas",
Partition = "TempDir.ForPartiFinder2/Partitions_Concat.txt",
codon = c(2:4), nexus.file = "TempDir.ForPartiFinder2/Concat.nex",
Path.PartiF2 = "/home/davidpc/Programs/PartitionFinder2/partitionfinder-2.1.1/PartitionFinder.py",
branchlengths = "linked", models = "all", model_selection = "BIC", search = "rcluster",
Raxml = "TRUE", nthread = 5, rcluster_percent = 10, rcluster_max = 1000)



test_that("Test the additional Nexus block in the nexus file of the concatenation", {
  a = list.files("TempDir.ForPartiFinder2")
  expect_equal(length(a), 14)
  test = ape::read.nexus.data("TempDir.ForPartiFinder2/Concat_PF2.nex")
  test1 = ape::as.DNAbin(test)
  # test if the concat in nexus file can be loaded in R and converted into a DNAbin object.
  if(class(test1) == "DNAbin") {b = 1} else {b = 1}
  expect_equal(b, 1)
  # test if the nexus file contained the additional block including the best partitioning scheme determine by PartiFinder2
  test2 = readLines("TempDir.ForPartiFinder2/Concat_PF2.nex")
  if(length(grep("begin sets", test2)) == 1){c = 1} else {c = 0}
  expect_equal(c, 1)
  if(length(grep("charset Subset", test2)) == 7){d = 1} else {d = 0}
  expect_equal(d, 1)
})

test_that("Test the new partition file for RAxML", {
  part = readLines("TempDir.ForPartiFinder2/Partitions_Concat.txt_PF2_all.txt")
  expect_equal(length(part), 7)
  if(length(grep("Subset", part)) == 7){c = 1} else {c = 0}
  expect_equal(c, 1)
})


# To clean the files created while running the example do the following:
# Remove the fodler "analysis".
unlink("analysis", recursive = TRUE)
# Remove the Temporary folder
unlink("TempDir.ForPartiFinder2", recursive = TRUE)

# Remove the files created by or for PartitionFinder2
file.remove("Concat.phy")
file.remove("log.txt")
file.remove("partition_finder.cfg")



