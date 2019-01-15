context("Test_Rm.OutSeq.Gap")

data(Example_16S_outlier)
Example_16S_outlier_align = Example_16S_outlier[[1]]
Table.Outlier.Seq.16S = Example_16S_outlier[[2]]

# Run the function. The name of all outliers sequences are stored
# in the third columns (i.e. "Hit_SeqName") of the table,
# becasue the same sequence might be present multiple time in the
# we use the function unique to remove the duplicated sequences (multiple hits by Blast).
# In the following example we also decided to remove another sequences
# ("_R_Polyprion_americanus|NA|AM158291|NA") which was not detected by
# the function Detect.Outlier.Seq.

New16SAlignment = Rm.OutSeq.Gap(input = c(as.character(unique(Table.Outlier.Seq.16S[,3])),
"_R_Polyprion_americanus|NA|AM158291|NA"), SeqInput = Example_16S_outlier_align,
AligoutputName = "16S_example_RMoutliers")

class(New16SAlignment) # it is a "DNAbin" object
New16SAlignment # 280 DNA sequences.

# To clean the file created while running the example do the following:
#' file.remove("16S_example_RMoutliers.fas")

test_that("Test the example, output in the r envir", {
  if(class(New16SAlignment) == "DNAbin") {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(labels(New16SAlignment)), 280)
  expect_equal(length(New16SAlignment), 1043280)
})

test = ape::read.dna("16S_example_RMoutliers.fas", format = "fasta")

test_that("Test the example, output in the WD", {
  if(class(test) == "DNAbin") {a = 1} else {a = 0}
  expect_equal(a, 1)
  expect_equal(length(labels(test)), 280)
  expect_equal(length(test), 1043280)
})

file.remove("16S_example_RMoutliers.fas")
