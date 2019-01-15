context("Test_ConstraintTaxo2newick")

data(TopoConstraints)
# The table storing the constraints include 22 topological constraints overall
# including constraints at the Family, Order, Series, Subdivision, Division,
# Subsection, Subcohort, Cohort, Supercohort, Infraclass, Subclass.
#
# The Classification table include 16 species from the New Zealand marine
# ray-finned fish species list.

# Create a Temporary folder to store the outputs of the function.
dir.create("TempDir.TopoConstraints")
# Run the function considering all the constraints
BackBoneTreeAll = ConstraintTaxo2newick(inputTaxo = TopoConstraints[[2]],
inputConst = TopoConstraints[[1]], outputNewick = "TempDir.TopoConstraints/BackboneTreeAll")


test_that("Test the newick tree", {
  if(class(BackBoneTreeAll[[2]]) == "phylo") {a = 1} else {a = 0}
  expect_equal(a, 1)
  # test the number of internal node must be 13
  expect_equal(BackBoneTreeAll[[2]]$Nnode, 13)
  # test the number of tips, must be 16
  expect_equal(length(BackBoneTreeAll[[2]]$tip.label), 16)
  # test if the tree is binary must be FALSE
  if(ape::is.binary(BackBoneTreeAll[[2]]) == FALSE) {b = 1} else {b = 0}
  expect_equal(b,1)
})

test_that("Test the newick tree exported in the WD", {
  tree = ape::read.tree("TempDir.TopoConstraints/BackboneTreeAll.txt")
  if(class(tree) == "phylo") {a = 1} else {a = 0}
  expect_equal(a, 1)
  # test the number of internal node must be 13
  expect_equal(tree$Nnode, 13)
  # test the number of tips, must be 16
  expect_equal(length(tree$tip.label), 16)
  # test if the tree is binary must be FALSE
  if(ape::is.binary(tree) == FALSE) {b = 1} else {b = 0}
  expect_equal(b,1)
  # test if the tree exported in the R envir and the WD are identical
  if(ape::all.equal.phylo(BackBoneTreeAll[[2]], tree) == TRUE) {c = 1} else {c = 0}
  expect_equal(c, 1)
})

test_that("Test classification table exported in the R envir", {
  expect_equal(dim(BackBoneTreeAll[[1]])[1], 16)
  expect_equal(dim(BackBoneTreeAll[[1]])[2], 11)
  # Test the proper completion of the table using few examples
  if(BackBoneTreeAll[[1]][5,1] == "Family_Bassanago_bulbiceps") {a = 1} else {a = 0}
  expect_equal(a, 1)
  if(BackBoneTreeAll[[1]][6,3] == "Series_Anguilliformes") {b = 1} else {b = 0}
  expect_equal(b, 1)
  if(BackBoneTreeAll[[1]][11,8] == "Cohort_Subcohort_Subsection_Division_Subdivision_Series_Stomiatiformes"){c = 1} else {c = 0}
  expect_equal(c, 1)
})



# Use only the constraint at the Family level
FamilyConst=TopoConstraints[[1]][TopoConstraints[[1]][,1]=="Family",]

# Run the function considering only the constraints at the family level.
BackBoneTreeFamily = ConstraintTaxo2newick(inputTaxo = TopoConstraints[[2]],
inputConst = FamilyConst, outputNewick = "TempDir.TopoConstraints/BackboneTreeFamily")


test_that("Test the newick tree, Family constraint only", {
  if(class(BackBoneTreeFamily[[2]]) == "phylo") {a = 1} else {a = 0}
  expect_equal(a, 1)
  # test the number of internal node must be 13
  expect_equal(BackBoneTreeFamily[[2]]$Nnode, 6)
  # test the number of tips, must be 10, because 4 taxa are not constrained at the family levels
  expect_equal(length(BackBoneTreeFamily[[2]]$tip.label), 10)
  # test if the tree is binary must be FALSE
  if(ape::is.binary(BackBoneTreeFamily[[2]]) == FALSE) {b = 1} else {b = 0}
  expect_equal(b,1)
})

test_that("Test the newick tree exported in the WD, Family constraint only", {
  tree = ape::read.tree("TempDir.TopoConstraints/BackboneTreeFamily.txt")
  if(class(tree) == "phylo") {a = 1} else {a = 0}
  expect_equal(a, 1)
  # test the number of internal node must be 13
  expect_equal(tree$Nnode, 6)
  # test the number of tips, must be 10, because 4 taxa are not constrained at the family levels
  expect_equal(length(tree$tip.label), 10)
  # test if the tree is binary must be FALSE
  if(ape::is.binary(tree) == FALSE) {b = 1} else {b = 0}
  expect_equal(b,1)
  # test if the tree exported in the R envir and the WD are identical
  if(ape::all.equal.phylo(BackBoneTreeFamily[[2]], tree) == TRUE) {c = 1} else {c = 0}
  expect_equal(c, 1)
})

test_that("Test classification table exported in the R envir, Family constraint only", {
  expect_equal(dim(BackBoneTreeFamily[[1]])[1], 10)
  expect_equal(dim(BackBoneTreeFamily[[1]])[2], 1)
  # Test the proper completion of the table using few examples
  if(length(unique(BackBoneTreeFamily[[1]][,1])) == 5) {a = 1} else {a = 0}
  expect_equal(a, 1)
  if(BackBoneTreeFamily[[1]][1,1] == "Aplodactylidae") {b = 1} else {b = 0}
  expect_equal(b, 1)
})

# To clean the files created.
unlink("TempDir.TopoConstraints", recursive = TRUE)
