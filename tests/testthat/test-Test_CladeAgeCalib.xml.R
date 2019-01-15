context("Test_CladeAgeCalib.xml")

src.dir = system.file("extdata/TopoConstraints", package = "regPhylo")
dir.create("TempDir.CladeAge")
# Set up the path of the TempDir folder.
dest.dir = paste(getwd(), "/TempDir.CladeAge", sep="")
file.names <- dir(src.dir)
# Copy all the files stored in regPhylo/extdata/TopoConstraints"
# into a temporary folder.
sapply(file.names, function(x) {
file.copy(from = paste(src.dir, x, sep = "/"),
to = paste(dest.dir, x, sep = "/"),
overwrite = FALSE) })

## Remove all the example files created by the function of interest from the example folder.
file.remove("TempDir.CladeAge/SimpleXml_ReadyForBEAST.xml", "TempDir.CladeAge/SimpleXml_2SpNoDNA_ReadyForBEAST.xml")

# We include 4 calibration constraints based on 4 clades (Elopomorpha,
# Anguilliformes, Stomiati, Perciformes)
# Import the table (i.e. "CalPointTable") listing the 4 clades constrained
# with the occurrence of the first fossil.
CalibrationTable4clades = read.delim("TempDir.CladeAge/Calib_CA_Fossil_4cl.csv",
sep="\t", header = TRUE)

#### Example restricted to taxa with at least a DNA sequence in the supermatrix.

# Load the classification table (the same that for the
# ConstraintTaxo2newick function), there are two way to do it:
# either through the .Rdata
data(TopoConstraints) # the second object of the list is the classification table
ClassifDF = TopoConstraints[[2]]

# Load the re-rooted tree (the same that for the
# ConstraintTaxo2newick function) in R (A rooted tree is available in the
# package and has been loaded in the temporary directory).

TreeRooted = ape::read.nexus("TempDir.CladeAge/RAxML_bipartitions.Concat_7GTR_Allconst_autoMRE_ReRooted")

# All diversification/turnover/sampling rates are from Matschiner et al. 2017.
CladeAgeCalib.xml(xml.input = "TempDir.CladeAge/SimpleXml_Wcont.xml", input.tree = TreeRooted,
output="TempDir.CladeAge/SimpleXml_ReadyForBEAST.xml", CalPointTable = CalibrationTable4clades,
MinDivRate = 0.041, MaxDivRate = 0.081, MinTurnoverRate = 0.0011,
MaxTurnoverRate = 0.37, MinSamplingRate = 0.0066, MaxSamplingRate = 0.01806,
xmltreename = "Subset1", inputTaxono = ClassifDF, Partitions = "TRUE")

test_that("multiplication works", {
  a = readLines("TempDir.CladeAge/SimpleXml_ReadyForBEAST.xml")
  # test the presence and the file of the xml file exported in the WD
  expect_equal(file.info("TempDir.CladeAge/SimpleXml_ReadyForBEAST.xml")$size, 88360)
  # test the presence of the CladeAge blocks, must be 4.
  expect_equal(length(grep("<fossilDistr id=\"FossilCalibration", a, fixed = TRUE)), 4)
  # test if the pop size of the random coalsecnt tree has been change to create the random tree.
  expect_equal(length(grep("<parameter id=\"randomPopSize.t:Subset7\" name=\"popSize\">100</parameter>", a, fixed = TRUE)), 1)
  # block of normal constraint without calibration point attached, must be 8
  expect_equal(length(grep("<distribution id=\"Node", a, fixed = TRUE)), 8)
})


#### Example including taxa without DNA in the supermatrix.

# Load the new classification table incluing the two additional taxa without DNA
# exported by the function MultiTopoConst.EditXML4BEAST2.
NewClassifDF = read.delim("TempDir.CladeAge/Classif18sp_2NoDNA.csv", sep = "\t", header = TRUE)

# Load the new rooted "RAxML" tree including the two additional taxa and also the
# bootstrap values for each node exported by the function MultiTopoConst.EditXML4BEAST2
NewTree = ape::read.tree("TempDir.CladeAge/BackboneTreeAll_2spNoDNA.txt")

# We load the calibartion table.
CalibrationTable4clades = read.delim("TempDir.CladeAge/Calib_CA_Fossil_4cl.csv",
sep="\t", header = TRUE)

# Run the function all the other settinga dn options remain unchanged
CladeAgeCalib.xml(xml.input = "TempDir.CladeAge/SimpleXml_2SpNoDNA_Wcont.xml",
input.tree = NewTree,
output = "TempDir.CladeAge/SimpleXml_2SpNoDNA_ReadyForBEAST.xml",
CalPointTable = CalibrationTable4clades,
MinDivRate = 0.041, MaxDivRate = 0.081, MinTurnoverRate = 0.0011,
MaxTurnoverRate = 0.37, MinSamplingRate = 0.0066, MaxSamplingRate = 0.01806,
xmltreename = "Subset1", inputTaxono = NewClassifDF, Partitions = "TRUE")



test_that("multiplication works", {
  a = readLines("TempDir.CladeAge/SimpleXml_2SpNoDNA_ReadyForBEAST.xml")
  # test the presence and the file of the xml file exported in the WD
  expect_equal(file.info("TempDir.CladeAge/SimpleXml_2SpNoDNA_ReadyForBEAST.xml")$size, 95712)
  # test the presence of the CladeAge blocks, must be 4.
  expect_equal(length(grep("<fossilDistr id=\"FossilCalibration", a, fixed = TRUE)), 4)
  # test if the pop size of the random coalsecnt tree has been change to create the random tree.
  expect_equal(length(grep("<parameter id=\"randomPopSize.t:Subset7\" name=\"popSize\">100</parameter>", a, fixed = TRUE)), 1)
  # block of normal constraint without calibration point attached, must be 9
  expect_equal(length(grep("<distribution id=\"Node", a, fixed = TRUE)), 9)
  # The name of th new species "Toto_toto" should appear 7 times
  expect_equal(length(grep("Toto_toto", a, fixed = TRUE)), 7)
  expect_equal(grep("<log idref=\"Perciformes_Node10.fossilprior\"/>", a, fixed = TRUE), 579)
  expect_equal(grep("<log idref=\"Anguilliformes_Node15.fossilprior\"/>", a, fixed = TRUE), 582)
  expect_equal(grep("<taxon idref=\"Titi_titi\" spec=\"Taxon\"/>", a, fixed = TRUE)[1], 435)
  expect_equal(grep("<taxon idref=\"Titi_titi\" spec=\"Taxon\"/>", a, fixed = TRUE)[2], 456)
})


rm(TreeRooted, CalibrationTable4clades, NewTree, NewClassifDF, ClassifDF)

# To clean the files created
unlink("TempDir.CladeAge", recursive = TRUE)

