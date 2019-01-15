
### Install all the packages used in regPhylo (all the one listed as "import")
?install.packages
install.packages("bold")
install.packages("seqinr")
install.packages("ape")
install.packages("geomedb")
install.packages("worrms")
install.packages("RJSONIO")
install.packages("stringr")
install.packages("fields")
install.packages("caper")
install.packages("phytools")

### Need to install separately all the regPhylo "import" R package 
install.packages("regPhylo_0.1.0.tar.gz" , repos=NULL, type="source", dependencies = TRUE)

#load the regPhylo R package
require(regPhylo)
# list all the object present in the regPhylo R package
ls("package:regPhylo")

# Test GetSeq_BOLD
Splist = cbind(TaxID = c("Diastobranchus capensis"), Species.Name = c("Diastobranchus capensis"))
BOLD.output = GetSeq_BOLD(splist = Splist, filename = "output.BOLD.txt")
dim(BOLD.output[[2]]) ### 33 by 81 normal include lot of record without a DNA sequence

data(Seq.Diastocapen)
BOLD.output = Seq.Diastocapen$Seq.BOLD
BOLD.output2=read.delim("output.BOLD.txt", sep="\t", h=T) 
dim(BOLD.output2) ### 8 by 82
dim(Seq.Diastocapen$Seq.BOLD) ### 8 by 82 Ok both toable output by the function in the working directory have the same number.


# Test the GetSeqInfo_NCBI_taxid
?GetSeqInfo_NCBI_taxid

Splist=cbind(TaxID=c(443778),
Species.Name=c("Diastobranchus capensis"))
# Run the function considering extracting all DNA sequences and associated metadata
## Not run: 
NCBI.output = GetSeqInfo_NCBI_taxid(splist = Splist, gene = "ALL",
filename = "output.NCBI.txt")
class(NCBI.output)
dim(NCBI.output)
NCBI.output2 = read.delim("output.NCBI.txt", sep="\t", h=T)
dim(NCBI.output2) ### 4 by 23 good its working


# Test the Congr.NCBI.BOLD.perReposit function
?Congr.NCBI.BOLD.perReposit

data(Seq.Diastocapen)
Seq.NCBI = Seq.Diastocapen$Seq.NCBI

# Load the table with BOLD data
Seq.BOLD = Seq.Diastocapen$Seq.BOLD

# Run the function with GenBank and BOLD data only.
AllSeqDF=Congr.NCBI.BOLD.perReposit(input.NCBI=Seq.NCBI,
input.BOLD=Seq.BOLD, output="AllSeq_NCBI_BOLD.txt")
dim(AllSeqDF)

# Load a personal repository
Seq.Perso = Seq.Diastocapen$Seq.Perso

# Run the function with GenBank, BOLD and the personal repository data.
AllSeqDF2=Congr.NCBI.BOLD.perReposit(input.NCBI=Seq.NCBI,
input.BOLD=Seq.BOLD, input.perReposit=Seq.Perso, perReposit="My.Rep",
output="AllSeq_NCBI_BOLD_perRep.txt")
dim(AllSeqDF2) # 11 by 25
head(AllSeqDF2)
tail(AllSeqDF2)


# Test the GeoCoord.WGS84
?GeoCoord.WGS84

data(Seq.Diastocapen)
Seq.DF = Seq.Diastocapen$Seq.DF

# Run the function
Seq.DF1=GeoCoord.WGS84(input=Seq.DF, output="Seq.DF.txt")
dim(Seq.DF1) # Ok good 11 by 26, add a new column
head(Seq.DF1) # looks all good

# Test the GeoCodeName
?GeoCodeName

data(Seq.Diastocapen)
Seq.DF1 = Seq.Diastocapen$Seq.DF1

# Run the function to retrieve geographic coordinates for sequences without coordinates.
Seq.DF3=GeoCodeName(input = Seq.DF1, output = "Seq.DF3.txt")
dim(Seq.DF3) ### ok good 11 by 28, 2 columns more "Location_used", and "Geo_accuracy"
colnames(Seq.DF3)
setdiff(colnames(Seq.DF3), colnames(Seq.DF1)) ### OK all good 
head(Seq.DF3)



# Test the SpeciesGeneMat.Bl
?SpeciesGeneMat.Bl

data(Seq.Diastocapen)
Seq.DF3 = Seq.Diastocapen$Seq.DF3

# Run the function
Seq.DF4=SpeciesGeneMat.Bl(input=Seq.DF3, output="Seq.DF4.")
class(Seq.DF4)
length(Seq.DF4)
dim(Seq.DF4[[1]]) # All good 1 by 4, 
Seq.DF4[[1]]
Seq.DF4[[2]]
Seq.DF4[[3]]
Seq.DF4[[4]]

# Test the SelGene.MaxSpCov
?SelGene.MaxSpCov

data(Seq.DF4) ## the first object of the list is the Species-by-gene matrix
# Run the function without a pre-selection of gene regions.
SelGene.MaxSpCov(input = Seq.DF4[[1]])
# Run the function with a pre-selection of gene regions.
SelGene.MaxSpCov(input = Seq.DF4[[1]], NBGene = c("co1", "12srrna"))

#build a fake species-by_gene  matrix
class(Seq.DF4[[1]])
SpbyGeneMat=rbind(as.matrix(Seq.DF4[[1]]), c("Titi_titi",0, 2, 1), c("Toto_toto", 0, 0, 4)) 

row.names(SpbyGeneMat)=SpbyGeneMat[,1]
SpbyGeneMat=as.data.frame(SpbyGeneMat)


SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = c("co1", "12srrna"))

SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = c("co1", "16srrna"))
SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = c("co1")) ## Problem see Word doc
SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = 1) ## Problem see Word doc

SelGene.MaxSpCov(input = SpbyGeneMat)


# Test the Matrix.Overlap
?Matrix.Overlap

Matrix.Overlap(input = Seq.DF4[[1]], gene.Sel = c("co1", "16srrna"))
Matrix.Overlap(input = SpbyGeneMat, gene.Sel = c("co1", "16srrna"))
Matrix.Overlap(input = SpbyGeneMat, gene.Sel = c("co1", "16srrna", "12srrna"))
# Ok all good


# Test AmMissData
?AmMissData

AmMissData(input = Seq.DF4[[1]], gene.list = c("co1", "16srrna"))
AmMissData(input = SpbyGeneMat, gene.list = c("co1", "16srrna"))
AmMissData(input = SpbyGeneMat, gene.list = c("co1", "16srrna", "12srrna"))
# Ok the function is working



# Test the Select.DNA
?Select.DNA

data(Seq.DF4) # the clean table of sequences and metadata is called "CleanDataTable"
# and is the fifth object of the list "Seq.DF4".
class(Seq.DF4)
length(Seq.DF4)

# Run the function
Seq.DF5=Select.DNA(input = Seq.DF4$CleanDataTable, gene.list = c("co1", "16srrna"),
output = "Seq.DF4.dataTable")
dim(Seq.DF5)
head(Seq.DF5)
Seq.DF5[,29]
# Ok the function ois working well.



# Test the SelBestSeq
?SelBestSeq

data(Seq.DF4) # the table is called "Seq.DF5" and constitute the
# sixth object of the list "Seq.DF4".
Seq.DF5 = Seq.DF4$Seq.DF5

# Run the function and export all the sequences in the alignment for
# each species and gene region.
## Not run: 
Seq.DF6=SelBestSeq(input = Seq.DF5, output = "Alig_Seq.DF5.All",
RefPoint = cbind(174.7976, -41.3355), perReposit = "My.Rep",
Alignment = T, MaxSeq = "ALL", gene.list = c("co1", "16srrna"),
SeqChoice = "Median")
dim(Seq.DF6)


# Run the function and export the best (the most proximal to the focal
area e.i. NZ) sequences in the alignment for each species and gene region,
using the sequence the with a median sequence length for each gene region and species.
Seq.DF7=SelBestSeq(input = Seq.DF5, output = "Alig_Seq.DF5.Best",
RefPoint = cbind(174.7976, -41.3355), perReposit = "My.Rep",
Alignment = T, MaxSeq =1, gene.list = c("co1", "16srrna"),
SeqChoice = "Median")
dim(Seq.DF7)


# Run the function and export the two most proximal sequences of the focal area (e.i. NZ)
in the alignment for each species and gene region, using the longuest sequence per
# gene region and species.
Seq.DF8=SelBestSeq(input = Seq.DF5, output = "Test/Alig_Seq.DF5.2Best",
RefPoint = cbind(174.7976, -41.3355), perReposit = "My.Rep",
Alignment = T, MaxSeq =2, gene.list = c("co1", "16srrna"),
SeqChoice = "Longest")

dim(Seq.DF8) # 3 by 33, the 16S got only one seqeunce anyway.

## OK the function is working well.


# Test Rm.OutSeq.Gap
?Rm.OutSeq.Gap

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
New16SAlignment # 280 DNA sequences. 3726 bp.

testalign=read.dna("16S_example_RMoutliers.fas", format="fasta")
testalign # 280 DNA sequences. 3726 bp.

### initial alignement size
as.DNAbin(Example_16S_outlier[[1]]) # 319 seqeunces, with a length of 4672 bp.

### Ok the function is working well


# Test the function Align.Concat
?Align.Concat

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
Align.Concat(input = "TempDir.ForConcat", Sp.List = NULL, outputConcat = NULL)

# Run the function to build a supermatrix including two species without DNA.
Align.Concat(input = "TempDir.ForConcat",
Sp.List = c("Titi_titi", "Toto_toto"),
outputConcat = "TempDir.ForConcat/Concat_2spNoDNA")

# To clean the files created while running the example do the following:
unlink("TempDir.ForConcat", recursive = TRUE)

### need to add one line of code to select only the .fas file in the 
# folder to run twice the function into the same folder.



# Test ConstraintTaxo2newick
?ConstraintTaxo2newick

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

# plot the constraining tree (the branch length does not matter, only the topology matters).
plot(BackBoneTreeAll[[2]], cex=0.8)

# Use only the constraint at the Family level
FamilyConst=TopoConstraints[[1]][TopoConstraints[[1]][,1]=="Family",]

# Run the function considering only the constraints at the family level.
BackBoneTreeFamily = ConstraintTaxo2newick(inputTaxo = TopoConstraints[[2]],
inputConst = FamilyConst, outputNewick = "TempDir.TopoConstraints/BackboneTreeFamily")

# plot the constraining tree (the branch length does not matter,
# only the topology matters), notice that only constrained taxa
# are present on the guiding tree, the unconstrained taxa will
# be positioned on the tree based on their molecular affinities.
plot(BackBoneTreeFamily[[2]], cex=0.8)

# Use only the constraint at the Family and Series levels.
FamilySeriesConst=TopoConstraints[[1]][c(which(TopoConstraints[[1]][,1] == "Family"),
which(TopoConstraints[[1]][,1] == "Series")),]

# Run the function considering only the constraints at the family and order levels.
BackBoneTreeFamilySeries = ConstraintTaxo2newick(inputTaxo = TopoConstraints[[2]],
inputConst = FamilySeriesConst, outputNewick = "TempDir.TopoConstraints/BackboneTreeFamilySeries")

# plot the constraining tree (the branch length does not matter,
# only the topology matters). notice that only constrained taxa
# are present on the guiding tree, the unconstrained taxa will
# be positioned on the tree based on their molecular affinities.
plot(BackBoneTreeFamilySeries[[2]], cex=0.8)

### OK the function is working well


# Test the function MultiTopoConst.EditXML4BEAST2

?MultiTopoConst.EditXML4BEAST2

src.dir = system.file("extdata/TopoConstraints", package = "regPhylo")
dir.create("TempDir.TopoConstraints2")
# Set up the path of the TempDir folder.
dest.dir = paste(getwd(), "/TempDir.TopoConstraints2", sep="")
file.names <- dir(src.dir)
# Copy all the files stored in regPhylo/extdata/TopoConstraints"
# into a temporary folder.
sapply(file.names, function(x) {
file.copy(from = paste(src.dir, x, sep = "/"),
to = paste(dest.dir, x, sep = "/"),
overwrite = FALSE) })

# Load the re-rooted tree in R (A rooted tree is available in the
# package and has been loaded in the temporary directory)
require(ape)
TreeRooted = read.nexus("TempDir.TopoConstraints2/RAxML_bipartitions.Concat_7GTR_Allconst_autoMRE_ReRooted")

#### Example restricted to taxa with at least a DNA sequence in the supermatrix.

# The input.xml file, called SimpleXml.xml, was prepare by Beauti, it is available
# in the folder "TempDir.TopoConstraints2/xmlfiles". This xml file includes
# the basic template of the xml file (including the supermatrix,
# the partions, the substitution and clock models, the priors,
# the MCMC parameters). However the hard constraints were not
# included in this input.xml file, but will be so by running the
# MultiTopoConst.EditXML4BEAST2 function below.

MultiTopoConst.EditXML4BEAST2(inputtree = TreeRooted,
output = "SimpleXml_Wcont.xml",
bootstrapTH = 100, xmltreename = "Subset1",
input.xml = "TempDir.TopoConstraints2/SimpleXml.xml",
Partitions = "TRUE")


# Load a table with the new taxa without DNA and the topological constraint
# This table will be used to fill the option "TaxaNoDNA".
TaxaNoDNA = cbind(c("Titi_titi", "Toto_toto"), c("Family", "Genus"), c("Congridae", "Scorpaena"))
colnames(TaxaNoDNA) = c("SpeciesName", "hier.level", "ConstraintName")
TaxaNoDNA = as.data.frame(TaxaNoDNA)

# Load the classification table (the same that for the
# ConstraintTaxo2newick function), there are two way to do it:
# either through the .Rdata
data(TopoConstraints) # the second object of the list is the classification table
dim(TopoConstraints[[2]]) # 16 by 23.
# or the classification table has been loaded into the temporary directory,
# and can be loaded into the R environment doing the following
ClassifDF = read.csv("TempDir.TopoConstraints2/Classif16sp.csv", header = TRUE)
dim(ClassifDF) # 16 by 23


# Run the function, but this time the input.xml file called "SimpleXml_2SpNoDNA.xml"
# generated by Beauti must be selected. This xml file template was built using
# an alignment (called "Concat_2spNoDNA_PF2.nex") including two empty sequences
# (i.e. "Titi_titi, "Toto_toto") built by the Align.Concat function.
# This alignment has been completed by the nexus block related to the best
# partition scheme defined by partitionFinder2 (here we used the same
# partition scheme that the one defined for the 16 taxa with at least a DNA sequence).

MultiTopoConst.EditXML4BEAST2(inputtree = TreeRooted,
output = "SimpleXml_2SpNoDNA_Wcont.xml",
bootstrapTH = 100, xmltreename = "Subset1",
input.xml = "TempDir.TopoConstraints2/SimpleXml_2SpNoDNA.xml",
Partitions = "TRUE", TaxaNoDNA = TaxaNoDNA,
TaxoTable = ClassifDF,
output.new.TaxoTable = "Classif18sp_2NoDNA.csv",
output.new.tree = "BackboneTreeAll_2spNoDNA.txt")

# Plot the newly constrained tree including the 2 taxa without DNA
NewTree = read.tree("TempDir.TopoConstraints2/BackboneTreeAll_2spNoDNA.txt")
plot(NewTree)

# To see the new classification table including the 2 taxa without DNA.
NewClassifDF = read.delim("TempDir.TopoConstraints2/Classif18sp_2NoDNA.csv",
sep = "\t", header = TRUE)
dim(NewClassifDF) # now there are 18 rows and still 23 columns.

### function working well once corrected the mistake in the example.



# Test the function CladeAgeCalib.xml
?CladeAgeCalib.xml

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
dim(TopoConstraints[[2]]) # 16 by 23.
# or the Classification table has been loaded into the temporary directory,
# and can be loaded into the r environment doing the following
ClassifDF = read.csv("TempDir.CladeAge/Classif16sp.csv", header = TRUE)
dim(ClassifDF) # 16 by 23

# Load the re-rooted tree (the same that for the
# ConstraintTaxo2newick function) in R (A rooted tree is available in the
# package and has been loaded in the temporary directory).
require(ape)
TreeRooted = read.nexus("TempDir.CladeAge/RAxML_bipartitions.Concat_7GTR_Allconst_autoMRE_ReRooted")

# All diversification/turnover/sampling rates are from Matschiner et al. 2017.
CladeAgeCalib.xml(xml.input = "TempDir.CladeAge/SimpleXml_Wcont.xml", input.tree = TreeRooted,
output="SimpleXml_ReadyForBEAST.xml", CalPointTable = CalibrationTable4clades,
MinDivRate = 0.041, MaxDivRate = 0.081, MinTurnoverRate = 0.0011,
MaxTurnoverRate = 0.37, MinSamplingRate = 0.0066, MaxSamplingRate = 0.01806,
xmltreename = "Subset1", inputTaxono = ClassifDF, Partitions = "TRUE")


#### Example including taxa without DNA in the supermatrix.

# Load the new classification table incluing the two additional taxa without DNA
# exported by the function MultiTopoConst.EditXML4BEAST2.
NewClassifDF = read.delim("TempDir.CladeAge/Classif18sp_2NoDNA.csv", sep = "\t", header = TRUE)

# Load the new rooted "RAxML" tree including the two additional taxa and also the
# bootstrap values for each node exported by the function MultiTopoConst.EditXML4BEAST2
require(ape)
NewTree = read.tree("TempDir.CladeAge/BackboneTreeAll_2spNoDNA.txt")

# We load the calibartion table.
CalibrationTable4clades = read.delim("TempDir.CladeAge/Calib_CA_Fossil_4cl.csv",
sep="\t", header = TRUE)

# Run the function all the other settinga dn options remain unchanged
CladeAgeCalib.xml(xml.input = "TempDir.CladeAge/SimpleXml_2SpNoDNA_Wcont.xml",
input.tree = NewTree,
output = "SimpleXml_2SpNoDNA_ReadyForBEAST.xml",
CalPointTable = CalibrationTable4clades,
MinDivRate = 0.041, MaxDivRate = 0.081, MinTurnoverRate = 0.0011,
MaxTurnoverRate = 0.37, MinSamplingRate = 0.0066, MaxSamplingRate = 0.01806,
xmltreename = "Subset1", inputTaxono = NewClassifDF, Partitions = "TRUE")

### The function is working well.









#################################################################################
#### Functuion depending on external software, see exemple using the phyloch R package.

require(phyloch)
require(ape)
data(woodmouse)
?mafft
mafft(woodmouse)
getwd()

testmafft=mafft(woodmouse, path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft")

testgblock=gblocks(testmafft, exec = "C:/Users/deme/Documents/Programs/Gblocks/Gblocks_Windows_0.91b/Gblocks_0.91b")
setwd(exec)

testprank=prank(testmafft, outfile = "TestPrank_woodmoose.fas", path="C:/Users/deme/Documents/Programs/Prank/prank.windows.140603/prank/bin/prank")

###





############################################################################################
##### Test the 29/10/2018 at 4:34pm for the the inclusion of Windows plateform for
### First.Align()
### 

### Need to install separately all the regPhylo "import" R package 
install.packages("regPhylo_0.1.0.tar.gz" , repos=NULL, type="source", dependencies = TRUE)


#load the regPhylo R package
require(regPhylo)
# list all the object present in the regPhylo R package
ls("package:regPhylo")



# Test First.Align
?First.Align.All


# To run the example it might be better to copy the input alignment files
# provided by the package to a temporary directory created into the
# current working directory.
src.dir = system.file("extdata/FirstToAlign", package = "regPhylo")
dir.create("TempDir.FirstToAlign")
# Set up the path of the TempDir folder.
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
output = "TempDir.FirstToAlign/FirstAligned"
nthread = 2
methods = "mafftfftnsi"
First.Align.All(input = input, output =
"TempDir.FirstToAlign/FirstAligned", nthread = 2,
methods = NULL, Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft")
Error in DFtempord[, 1] : incorrect number of dimensions

First.Align.All(input = input, output =
"TempDir.FirstToAlign/FirstAligned", nthread = 1,
methods = "mafftfftnsi", Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft")

Error in DFtempord[, 1] : incorrect number of dimensions




First.Align.All

function (input = NULL, output = NULL, nthread = NULL, methods = NULL, 
    Mafft.path = NULL) 
{
    AlignSelect = list.files(input)
    AlignSelect = AlignSelect[grep(".fas", AlignSelect)]
    AlignSelect2 = paste("revc_", AlignSelect, sep = "")
    dir.create(output)
    mafft = "mafft"
    os <- .Platform$OS
    if (os == "windows") {
        if (missing(Mafft.path)) {
            stop("The path to the mafft executable must be provided in Mafft.path")
        }
        mafft = Mafft.path
    }
    cl <- parallel::makeCluster(nthread)
    parallel::clusterExport(cl, varlist = c("output", "input", 
        "nthread")) #### HERE NEED TO ADD methods and 
    mafftfftns1.align = function(x) {
        a = paste(mafft, " --retree 1 --maxiterate 0 --adjustdirection ", 
            input, "/", x, " > ", input, "/", "revc_", x, sep = "")
        system(a)
    }
    parallel::parLapply(cl, AlignSelect, mafftfftns1.align)
    p = list.files(input)
    file.copy(from = paste(input, "/", p[grep("revc_", p, fixed = T)], 
        sep = ""), to = output)
    p1 = list.files(output)
    p2 = gsub("revc_", "Mafftfftns1_", p1, fixed = T)
    for (i in 1:length(p1)) {
        file.rename(from = paste(output, "/", p1[i], sep = ""), 
            to = paste(output, "/", p2[i], sep = ""))
    }
    parallel::stopCluster(cl)
    cl <- parallel::makeCluster(nthread)
    parallel::clusterExport(cl, varlist = c("output", "input", 
        "nthread"))
    if (length(which(methods == "mafftfftnsi")) == 1) {
        mafft = "fftnsi "
        os <- .Platform$OS
        if (os == "windows") {
            if (missing(Mafft.path)) {
                stop("The path to the mafft executable must be provided in Mafft.path")
            }
            mafft = paste(Mafft.path, " ", sep = "")
        }
        mafftfftnsi.align = function(x) {
            a = paste(mafft, input, "/", x, " > ", output, "/", 
                "Mafftfftnsi_", x, sep = "")
            system(a)
        }
        parallel::parLapply(cl, AlignSelect2, mafftfftnsi.align)
    }
    if (length(which(methods == "pasta")) == 1) {
        pasta.align = function(x) {
            a = paste("run_pasta.py -i ", input, "/", x, " -j PASTA_", 
                " --num-cpus 1", sep = "")
            system(a)
        }
        parallel::parLapply(cl, AlignSelect2, pasta.align)
        b = list.files(input)
        b1 = paste(input, "/", b[grep(".aln", b)], sep = "")
        bb = gsub(".aln", ".fas", gsub("(_[0-9]*.marker001.)", 
            "_", b1, perl = T), fixed = T)
        i = 1
        for (i in 1:length(bb)) {
            file.rename(from = b1[i], to = bb[i])
        }
        file.copy(from = bb, to = output)
        c = paste(input, "/", b[grep("PASTA_", b)], sep = "")
        file.remove(c)
    }
    parallel::stopCluster(cl)
    bn = list.files(input)
    file.remove(paste(input, "/", bn[grep("revc_", bn, fixed = T)], 
        sep = ""))
    b = list.files(output)
    bb = paste(output, "/", b[grep("_revc_", b, fixed = T)], 
        sep = "")
    bb2 = gsub("_revc", "", bb, fixed = T)
    for (i in 1:length(bb)) {
        file.rename(from = bb[i], to = bb2[i])
    }
    a = list.files(output)
    a1 = strsplit(a, "_")
    Nbprog = vector()
    Nbgene = vector()
    i = 1
    for (i in 1:length(a1)) {
        Nbprog = c(Nbprog, a1[[i]][1])
        Nbgene = c(Nbgene, a1[[i]][length(a1[[i]])])
    }
    Uniprog = unique(Nbprog)
    Unigene = unique(Nbgene)
    x = Unigene
    listgeneb = vector()
    listRevComp = matrix(NA, ncol = 2)[-1, ]
    i = 1
    for (i in 1:length(x)) {
        listgeneb = c(listgeneb, a[which(Nbgene == x[i])]) ## HERE A PROBLEM, don't need to 
# reinclude the previous gene so the line should be 
# listgeneb =  a[which(Nbgene == x[i])] 
# and the line above listgeneb = vector(), is not needed.
        j = 1
        for (j in 1:length(listgeneb)) {
            listAlig = tryCatch(seqinr::read.fasta(paste(output, 
                "/", listgeneb[j], sep = ""), as.string = T), 
                error = function(e) e)
            if (inherits(listAlig, "simpleError") == FALSE) {
                SeqName = labels(listAlig)
                SeqT = vector()
                k = 1
                for (k in 1:length(listAlig)) SeqT = c(SeqT, 
                  listAlig[[k]][1])
                DFtemp = cbind(SeqName, SeqT)
### Here need to add the following line of code to avoid 
# an error when the alignment get 1 sequence
# if(length(listAlig) == 1){
#DFtempord = t(as.matrix(DFtemp[order(DFtemp[, 1]), ]))
#} else {
#DFtempord = DFtemp[order(DFtemp[, 1]), ]
#}
                DFtempord = DFtemp[order(DFtemp[, 1]), ] # row need to be change by the above.
                Seq_Name = paste(">", paste(DFtempord[, 1], DFtempord[, 
                  2], sep = "|_|"), sep = "")
                AlignFasta = unlist(strsplit(Seq_Name, "|_|", 
                  Seq_Name, fixed = TRUE))
                write(AlignFasta, file = paste(output, "/", listgeneb[j], 
                  sep = ""))
            }
            else {
                warning(paste("The alignment ", listgeneb[j], 
                  " contains no sequence, either because the original alignment contains 1 sequence only, or because PASTA couldn't find a proper tree and crashed", 
                  sep = ""))
            }
        }
        listAli = tryCatch(seqinr::read.fasta(paste(output, "/", 
            a[grep(x[i], a)][1], sep = ""), as.string = T), error = function(e) e)
        if (inherits(listAlig, "simpleError") == FALSE) {
            SeqNa = labels(listAli)
            RevComp = SeqNa[grep("^_R_", SeqNa, perl = T)]
            if (length(RevComp) > 0) {
                listRevComp = rbind(listRevComp, cbind(rep(Unigene[i], 
                  length(RevComp)), RevComp))
            }
            else {
                listRevComp = rbind(listRevComp, cbind(rep(Unigene[i], 
                  1), NA))
            }
        }
    }
    colnames(listRevComp) = c("Gene", "SequenceName")
    utils::write.table(listRevComp, file = paste(input, "/ListSeq_RevCompl_FirstAlignAll.txt", 
        sep = ""), row.names = F, sep = "\t")
    return(listRevComp)
}
<environment: namespace:regPhylo>

###  CORRECTION DONE

# Test the Filtering.align.Gblocks
?Filtering.align.Gblocks


src.dir = system.file("extdata/multi.align/multi.aligned", package = "regPhylo")
dir.create("TempDir")
# Set up the path of the TempDir folder.
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
Filtering.align.Gblocks(input = "TempDir", LessStringent = TRUE,
output = "TrimmedGblocks", Type = "d", remove.empty.align = TRUE , 
Gblocks.path = "C:/Users/deme/Documents/Programs/Gblocks/Gblocks_Windows_0.91b/Gblocks_0.91b")

[1] TRUE
There were 13 warnings (use warnings() to see them)
> warnings()
Warning messages:
1: running command 'Gblocks TempDir/Mafftfftns1_Alig_12srrna.fas -t=d -b2=4 -b4=5 -b5=h -e=-gbls' had status 127
2: running command 'Gblocks TempDir/Mafftfftns1_Alig_co1.fas -t=d -b2=9 -b4=5 -b5=h -e=-gbls' had status 127
3: running command 'Gblocks TempDir/Mafftfftns1_Alig_cytb.fas -t=d -b2=3.5 -b4=5 -b5=h -e=-gbls' had status 127
4: running command 'Gblocks TempDir/Mafftfftns1_Alig_rag1.fas -t=d -b2=4.5 -b4=5 -b5=h -e=-gbls' had status 127
5: running command 'Gblocks TempDir/Mafftfftnsi_Alig_12srrna.fas -t=d -b2=4 -b4=5 -b5=h -e=-gbls' had status 127
6: running command 'Gblocks TempDir/Mafftfftnsi_Alig_co1.fas -t=d -b2=9 -b4=5 -b5=h -e=-gbls' had status 127
7: running command 'Gblocks TempDir/Mafftfftnsi_Alig_cytb.fas -t=d -b2=3.5 -b4=5 -b5=h -e=-gbls' had status 127
8: running command 'Gblocks TempDir/Mafftfftnsi_Alig_rag1.fas -t=d -b2=4.5 -b4=5 -b5=h -e=-gbls' had status 127
9: running command 'Gblocks TempDir/PASTA_Alig_12srrna.fas -t=d -b2=4 -b4=5 -b5=h -e=-gbls' had status 127
10: running command 'Gblocks TempDir/PASTA_Alig_co1.fas -t=d -b2=9 -b4=5 -b5=h -e=-gbls' had status 127
11: running command 'Gblocks TempDir/PASTA_Alig_cytb.fas -t=d -b2=3.5 -b4=5 -b5=h -e=-gbls' had status 127
12: running command 'Gblocks TempDir/PASTA_Alig_rag1.fas -t=d -b2=4.5 -b4=5 -b5=h -e=-gbls' had status 127
13: In file.remove(paste(input, "/", out[grep(".htm", out,  ... :
  cannot remove file 'TempDir/', reason 'Permission denied'




src.dir = paste(getwd(), "/", input, sep = "")
            file.names = list.files(src.dir)
            oldwd = getwd()
            i = 1
            for (i in 1:length(file.names)) {
### the line of code below require is the following
# file.copy(from = paste(src.dir, "/", file.names[i], 
# sep = ""), to = paste(Gblocks.path, "/", file.names[i], sep = ""))
                file.copy(from = paste(src.dir, "/", file.names[i], 
                  sep = ""), to = paste(Gblocks.path, "/", file.names[i], sep = "")) 
                setwd(Gblocks.path)
                Align = seqinr::read.fasta(file.names[i], as.string = TRUE)
                MinNbSeqFlankPos = (length(labels(Align))/2) + 
                  1
                if (LessStringent == TRUE) {
                  a = paste("./Gblocks ", file.names[i], " -t=", 
                    Type, " -b2=", MinNbSeqFlankPos, " -b4=5 -b5=h -e=-gbls", 
                    sep = "")
                  system(a)
                  out = list.files(Gblocks.path)
### need to remove the .htm files created 
file.remove(out[grep(".htm", out, fixed = TRUE)])
out2 = list.files(Gblocks.path)
### the following line need to be
                  outOri = out2[grep(".fas-gbls", out2, fixed = T)]
                  outRena = paste("Gblocksls_", gsub("-gbls", 
                    "", outOri, fixed = TRUE), sep = "")
                  file.rename(paste(Gblocks.path, "/", outOri, 
                    sep = ""), paste(oldwd, "/", output, "/", outRena, sep = ""))
### need to remove all the temporary files
out3 = lits.files(Gblocks.path)
file.remove(out3[grep(".fas", out3, fixed = TRUE)])
                }
                else {
                  a = paste("./Gblocks ", file.names[i], " -t=", 
                    Type, " -e=-gbms", sep = "")
                  system(a)
                  out = list.files(Gblocks.path)
                  ### need to remove the .htm files created 
file.remove(out[grep(".htm", out, fixed = TRUE)])
out2 = list.files(Gblocks.path)
### the following line need to be
                  outOri = out2[grep(".fas-gbls", out2, fixed = T)]
                  outRena = paste("Gblocksls_", gsub("-gbls", 
                    "", outOri, fixed = TRUE), sep = "")
                  file.rename(paste(Gblocks.path, "/", outOri, 
                    sep = ""), paste(oldwd, "/", output, "/", outRena, sep = ""))
### need to remove all the temporary files
out3 = lits.files(Gblocks.path)
file.remove(out3[grep(".fas", out3, fixed = TRUE)])

                }
                setwd(oldwd)











Filtering.align.Gblocks = function(input = NULL, LessStringent = NULL, Type = NULL, output = NULL, 
    remove.empty.align = NULL, Gblocks.path = NULL) 
{
    b = list.files(input)
    bb = b[grep(".fas", b, fixed = T)]
    Filespace = vector()
    i = 1
    for (i in 1:length(bb)) {
        Filespace = c(Filespace, file.info(paste(input, "/", 
            bb[i], sep = ""))[[1]])
    }
    pbfile = bb[which(Filespace == 0)]
    if (remove.empty.align) {
        if (length(which(Filespace == 0)) > 0) {
            bb = bb[-which(Filespace == 0)]
        }
    } else {
        warning(paste(pbfile, collapse = "\n"))
        stop(paste("Some alignments might be empty; empty alignments must be remove from the folder before running the function", 
            "\n", "Look at in priority the files targeted by the warning message", 
            "\n", sep = ""))
    }
    if (is.null(output)) {
        output = input
    } else {
        dir.create(output)
    }
    os <- .Platform$OS
    if (os == "windows") {
        if (missing(Gblocks.path)) {
            stop("The path to the Gblocks executable must be provided in Gblocks.path")
        } else {
            src.dir = paste(getwd(), "/", input, sep = "")
            file.names = list.files(src.dir)
            oldwd = getwd()
            i = 1
            for (i in 1:length(file.names)) {
                file.copy(from = paste(src.dir, "/", file.names[i], 
                  sep = ""), to = paste(Gblocks.path, "/", file.names[i], sep = ""))
                setwd(Gblocks.path)
                Align = seqinr::read.fasta(file.names[i], as.string = TRUE)
                MinNbSeqFlankPos = (length(labels(Align))/2) + 
                  1
                if (LessStringent == TRUE) {
                  a = paste("./Gblocks ", file.names[i], " -t=", 
                    Type, " -b2=", MinNbSeqFlankPos, " -b4=5 -b5=h -e=-gbls", 
                    sep = "")
                  system(a)
                  out = list.files(Gblocks.path)
                  ### need to remove the .htm files created 
file.remove(out[grep(".htm", out, fixed = TRUE)])
out2 = list.files(Gblocks.path)
### the following line need to be
                  outOri = out2[grep(".fas-gbls", out2, fixed = T)]
                  outRena = paste("Gblocksls_", gsub("-gbls", 
                    "", outOri, fixed = TRUE), sep = "")
                  file.rename(paste(Gblocks.path, "/", outOri, 
                    sep = ""), paste(oldwd, "/", output, "/", outRena, sep = ""))
### need to remove all the temporary files
out3 = list.files(Gblocks.path)
file.remove(out3[grep(".fas", out3, fixed = TRUE)])
                } else {
                  a = paste("./Gblocks ", file.names[i], " -t=", 
                    Type, " -e=-gbms", sep = "")
                  system(a)
                  out = list.files(Gblocks.path)
                  ### need to remove the .htm files created 
file.remove(out[grep(".htm", out, fixed = TRUE)])
out2 = list.files(Gblocks.path)
### the following line need to be
                  outOri = out2[grep(".fas-gbms", out2, fixed = T)]
                  outRena = paste("Gblocksms_", gsub("-gbms", 
                    "", outOri, fixed = TRUE), sep = "")
                  file.rename(paste(Gblocks.path, "/", outOri, 
                    sep = ""), paste(oldwd, "/", output, "/", outRena, sep = ""))
### need to remove all the temporary files
out3 = list.files(Gblocks.path)
file.remove(out3[grep(".fas", out3, fixed = TRUE)])
                }
                setwd(oldwd)
            }
        }
    } else {
        if (LessStringent == TRUE) {
            Gblocks.lessStr = function(x) {
                Align = seqinr::read.fasta(paste(input, "/", 
                  x, sep = ""), as.string = TRUE)
                MinNbSeqFlankPos = (length(labels(Align))/2) + 
                  1
                a = paste("Gblocks ", input, "/", x, " -t=", 
                  Type, " -b2=", MinNbSeqFlankPos, " -b4=5 -b5=h -e=-gbls", 
                  sep = "")
                system(a)
            }
            lapply(bb, Gblocks.lessStr)
            out = list.files(input)
            file.remove(paste(input, "/", out[grep(".htm", out, 
                fixed = TRUE)], sep = ""))
            out2 = list.files(input)
            outOri = out2[grep(".fas-gbls", out2, fixed = T)]
            outRena = paste("Gblocksls_", gsub("-gbls", "", outOri, 
                fixed = TRUE), sep = "")
            file.rename(paste(input, "/", outOri, sep = ""), 
                paste(output, "/", outRena, sep = ""))
        } else {
            Gblocks.def = function(x) {
                a = paste("Gblocks ", input, "/", x, " -t=", 
                  Type, " -e=-gbms", sep = "")
                system(a)
            }
            lapply(bb, Gblocks.def)
            out = list.files(input)
            file.remove(paste(input, "/", out[grep(".htm", out, 
                fixed = TRUE)], sep = ""))
            out2 = list.files(input)
            outOri = out2[grep(".fas-gbms", out2, fixed = TRUE)]
            outRena = paste("Gblocksms_", gsub("-gbms", "", outOri, 
                fixed = TRUE), sep = "")
            file.rename(paste(input, "/", outOri, sep = ""), 
                paste(output, "/", outRena, sep = ""))
        }
    }
}

### other test
Filtering.align.Gblocks(input = "TempDir", LessStringent = TRUE, output = "TrimmedGblocks", Type = "d", remove.empty.align = TRUE , Gblocks.path = "C:/Users/deme/Documents/Programs/Gblocks/Gblocks_Windows_0.91b/Gblocks_0.91b")

rm(input, output, out, out2, out3, outOri, outRena, MinNbSeqFlankPos, 
oldwd , Align, Type, Gblocks.path, remove.empty.align)


Need to remove the the empty folder "TrimmedGblocks" in the list.files(Gblocks.path)



### test how to start muscle 
a = paste("C:/Users/deme/Documents/Programs/Muscle/muscle3.8.31_i86win32.exe", 
" -in TempDir/Mafftfftns1_Alig_12srrna.fas", " -out Testmuscle", sep="")

system(a)

### This command line is working on windows from a remote WD.
system(paste("C:/Users/deme/Documents/Programs/Muscle/muscle3.8.31_i86win32.exe -in TempDir/Mafftfftns1_Alig_12srrna.fas -out Testmuscle", sep=""))


### Test with Prank

#it works
a = paste("C:/Users/deme/Documents/Programs/Prank/prank.windows.140603/prank/bin/prank.exe",
" -d=TempDir/Mafftfftns1_Alig_12srrna.fas -o=testPrnk.fas", sep="")
system(a)

# Form Multi.Align a = paste("prank -d=", input, "/", x, " -o=", output, 
                "/", "Prank_", x, sep = "")


# Test for TrimAl

a = paste("C:/Users/deme/Documents/Programs/TrimAl/trimal.v1.2rev59/trimAl/bin/trimal.exe -in TempDir/Mafftfftns1_Alig_12srrna.fas -out testrimal -automated1")
a
system(a)
# it works for TrimAl.

a = paste("trimal -in ", input, "/", x, " -out ", output, 
            "/", "trimAuto_", x, " -automated1", sep = "")
        system(a)


# Try with partitionFinder2
install python V2.7
set up Python27 in the PATH using the following
Exact steps for adding Python to the path on Windows 7+:
1. Computer -> System Properties (or Win+Break) -> Advanced System Settings
2. Click the Environment variables... button (in the Advanced tab)
3. Edit PATH and append ;C:\Python27 to the end (substitute your Python version)
4. Click OK. Note that changes to the PATH are only reflected in command prompts opened after the change took place.

system("python C:/Users/deme/Documents/Programs/PartitionFinder2/partitionfinder-2.1.1/partitionfinder-2.1.1/PartitionFinder.py --version")


# To run PartitionFidner2 we need to install the following package according PartitionFinder2 manual.
numpy
pandas
pytables
pyparsing
scipy
sklearn

To install those package do te following on windows

open the cmd 
cd C:\Python27\Scripts
pip.exe install numpy
pip.exe install pandas
pip.exe install pytables ### now this package is called "tables"
pip.exe install pyparsing
pip.exe install scipy
pip.exe install sklearn




### for NCBI BLAST+ tools
It works!
system("blastn -help")
system("blastn -help")








############################################################################################
##### Test the 01/11/2018 at 3:03pm for the the inclusion of Windows plateform for
### 

### Need to install separately all the regPhylo "import" R package 
install.packages("regPhylo_0.1.0.tar.gz" , repos=NULL, type="source", dependencies = TRUE)


#load the regPhylo R package
require(regPhylo)
# list all the object present in the regPhylo R package
ls("package:regPhylo")


# Test First.Align.All
?First.Align.All


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
output = "TempDir.FirstToAlign/FirstAligned"
nthread = 2
methods = c("mafftfftnsi", "pasta")
First.Align.All(input = input, output =
"TempDir.FirstToAlign/FirstAligned", nthread = 2,
methods = c("mafftfftnsi"), 
Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft")

Gene      SequenceName
[1,] "co1.fas" NA          
Warning message:
In First.Align.All(input = input, output = "TempDir.FirstToAlign/FirstAligned",  :
  The alignment Mafftfftnsi_Alig_Seq.DF5.All_16srrna.fas contains no sequence, either because the original alignment contains 1 sequence only, or because PASTA couldn't find a proper tree and crashed


First.Align.All(input = input, output =
"TempDir.FirstToAlign/FirstAligned", nthread = 2,
methods = NULL, 
Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft")

First.Align.All(input = input, output =
"TempDir.FirstToAlign/FirstAligned", nthread = 2, 
Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft")


First.Align.All(input = input, output =
"TempDir.FirstToAlign/FirstAligned", nthread = 2,
methods = c("mafftfftnsi", "pasta"), 
Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft")



if (length(which(methods == "pasta")) == 1) {
### Here need to to add an additionla ligne of code 
if(os == "windows") {
warning("Pasta cannot be used on a Windows plateform at the moment")
} else {
        pasta.align = function(x) {
            a = paste("run_pasta.py -i ", input, "/", x, " -j PASTA_", 
                " --num-cpus 1", sep = "")
            system(a)
        }
        parallel::parLapply(cl, AlignSelect2, pasta.align)
        b = list.files(input)
	  b1 = paste(input, "/", b[grep(".aln", b)], sep = "")
        bb = gsub(".aln", ".fas", gsub("(_[0-9]*.marker001.)", 
            "_", b1, perl = T), fixed = T)
        i = 1
        for (i in 1:length(bb)) {
            file.rename(from = b1[i], to = bb[i])
        }
        file.copy(from = bb, to = output)
        c = paste(input, "/", b[grep("PASTA_", b)], sep = "")
        file.remove(c)
    }
} # end if(os == "windows") { else


### overall the function works!!!




### Test the Multi.Align 
?Multi.Align


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
methods = c("mafftfftnsi")
input = "TempDir.Multi.aligned"
# run the function using two additional alignmenet program "mafftfftnsi" and "pasta"
Multi.Align(input = input, output = "multi.alignedTest", nthread = 3,
methods = c("mafftfftnsi", "muscle"), Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft",
Muscle.path = "C:/Users/deme/Documents/Programs/Muscle/muscle3.8.31_i86win32.exe")

# It works

# adding Prank 
Multi.Align(input = input, output = "multi.alignedTest", nthread = 3,
methods = c("mafftfftnsi", "muscle", "prank"), Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft",
Muscle.path = "C:/Users/deme/Documents/Programs/Muscle/muscle3.8.31_i86win32.exe",
Prank.path = "C:/Users/deme/Documents/Programs/Prank/prank.windows.140603/prank/bin/prank.exe")

## and mafftfftns2
Multi.Align(input = input, output = "multi.alignedTest", nthread = 3,
methods = c("mafftfftnsi", "mafftfftns2"), Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft",
Muscle.path = "C:/Users/deme/Documents/Programs/Muscle/muscle3.8.31_i86win32.exe",
Prank.path = "C:/Users/deme/Documents/Programs/Prank/prank.windows.140603/prank/bin/prank.exe")

# Everything work well.



### Test the "Filtering.align.Gblocks"
?Filtering.align.Gblocks

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
Filtering.align.Gblocks(input = "TempDir", LessStringent = TRUE,
output = "TrimmedGblocks", Type = "d", remove.empty.align = TRUE,
Gblocks.path = "C:/Users/deme/Documents/Programs/Gblocks/Gblocks_Windows_0.91b/Gblocks_0.91b")

### here need to precise that the name of the executable must not be included in the path of the software.
### need to add that in the man of the function!

Filtering.align.Gblocks(input = "TempDir", LessStringent = FALSE,
output = "TrimmedGblocks2", Type = "d", remove.empty.align = FALSE,
Gblocks.path = "C:/Users/deme/Documents/Programs/Gblocks/Gblocks_Windows_0.91b/Gblocks_0.91b")

### NEED To correct the following
if (remove.empty.align) {
        if (length(which(Filespace == 0)) > 0) {
            bb = bb[-which(Filespace == 0)]
        }
    } else {
	### Need to include the following
	if(length(which(Filespace == 0)) > 0){
        warning(paste(pbfile, collapse = "\n"))
        stop(paste("Some alignments might be empty; empty alignments must be remove from the folder before running the function", 
            "\n", "Look at in priority the files targeted by the warning message", 
            "\n", sep = ""))
	} ### need to include the closing brackets.
    }



Filtering.align.Gblocks(input = "TempDir", LessStringent = FALSE,
output = "TrimmedGblocks2", Type = "d", remove.empty.align = TRUE,
Gblocks.path = "C:/Users/deme/Documents/Programs/Gblocks/Gblocks_Windows_0.91b/Gblocks_0.91b")




### test the Filtering.align.trimal
?Filtering.align.Trimal

# To run the example it might be better to copy the input alignment files
# provided by the package to a temporary directory created into the
# current working directory.
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
Filtering.align.Trimal(input = input, output = "Trimmed", 
TrimAl.path = "C:/Users/deme/Documents/Programs/TrimAl/trimal.v1.2rev59/trimAl/bin/trimal.exe")

### it works well.



### Test PartiFinder2

?PartiFinder2

# To run the example it might be better to copy the input files
# provided by the package to a temporary directory created into the
# current working directory.
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

input = "TempDir.PartiFinder2/Concat.fas"
Partition = "TempDir.PartiFinder2/Partitions_Concat.txt"
# Open the convtab.txt document to know which genes need
# to be partitioned for the first , second, and third codon position.
read.delim("TempDir.ForPartiFinder2/convtab.txt", sep = "\t", header = TRUE)

# Run the function using RAxML software, BIC criteria,
# the fast "rcluster" algorithm and the classic default parameters.
PartiFinder2(input = "TempDir.ForPartiFinder2/Concat.fas",
Partition = "TempDir.ForPartiFinder2/Partitions_Concat.txt",
codon = c(2:4), nexus.file = "TempDir.ForPartiFinder2/Concat.nex",
Path.PartiF2 = "C:/Users/deme/Documents/Programs/PartitionFinder2/partitionfinder-2.1.1/partitionfinder-2.1.1/PartitionFinder.py",
branchlengths = "linked", models = "all", model_selection = "BIC", search = "rcluster",
Raxml = "TRUE", nthread = 5, rcluster_percent = 10, rcluster_max = 1000)


# To clean the files created while running the example do the following:
# Remove the fodler "analysis".
unlink("analysis", recursive = TRUE)
# Remove the Temporary folder
unlink("TempDir.ForPartiFinder2", recursive = TRUE)

# Remove the files created by or for PartitionFinder2
file.remove("Concat.phy")
file.remove("log.txt")
file.remove("partition_finder.cfg")

# Run the function using RAxML software, BIC criteria,
# the fast "rcluster" algorithm and the classic default parameters.
PartiFinder2(input = "TempDir.ForPartiFinder2/Concat.fas",
Partition = "TempDir.ForPartiFinder2/Partitions_Concat.txt",
codon = c(2:4), nexus.file = "TempDir.ForPartiFinder2/Concat.nex",
Path.PartiF2 = "C:/Users/deme/Documents/Programs/PartitionFinder2/partitionfinder-2.1.1/partitionfinder-2.1.1/PartitionFinder.py",
branchlengths = "linked", models = "HKY", model_selection = "BIC", search = "all",
Raxml = "FALSE", nthread = 5, rcluster_percent = 10, rcluster_max = 1000)


# To clean the files created while running the example do the following:
# Remove the fodler "analysis".
unlink("analysis", recursive = TRUE)
# Remove the Temporary folder
unlink("TempDir.ForPartiFinder2", recursive = TRUE)

# Remove the files created by or for PartitionFinder2
file.remove("Concat.phy")
file.remove("log.txt")
file.remove("partition_finder.cfg")




### Test the Detect.Outliers

?Detect.Outlier.Seq

### need to correct the man using the following sentence as description
The function helps to relatively quickly detect some outlier sequences (miss aligned sequences, which might be caused by different problems such as gene or species annotation problems, presence of paralogs sequences...) 


data(Example_16S_outlier)
Example_16S_outlier_align = Example_16S_outlier[[1]]
## Not run: 

# Running the function with the 'Comb' option,
# a distance threshold of 0.6, and disabling the search for secondary outlier sequences.
S16_MisAlign0.6_1 = Detect.Outlier.Seq(inputal = Example_16S_outlier_align,
Strat.DistMat = "Comb", Dist.Th = 0.6, output = "Example_16S_outliers_1.txt",
Second.Outlier = "No")
dim(S16_MisAlign0.6_1) ### 36 sequences detected as primary oultier sequences.
head(S16_MisAlign0.6_1)

# The output table is also present as an external data table provided in the regPhylo r package
and can be access by the following code:
# a = system.file("extdata/ExampleOutliers/Example_16S_outliers_1.txt", package = "regPhylo")
# Example_16S_outliers_1 = read.delim(a, sep="\t", header = TRUE)


# Running the function with the 'Comb' option, a distance threshold of 0.6,
# and allowing a search for secondary outlier sequences using local blast database.
S16_MisAlign0.6_2 = Detect.Outlier.Seq(inputal = Example_16S_outlier_align,
Strat.DistMat = "Comb", Dist.Th = 0.6, output = "Example_16S_outliers_2.txt",
Second.Outlier = "Yes", Bitsc.Th = 0.8)
length(unique(S16_MisAlign0.6_2[,3])) ### 38 primary and secondary outlier seqeunces detected.
head(S16_MisAlign0.6_2)
dim(S16_MisAlign0.6_2)


test=read.delim("Example_16S_outliers_2.txt", sep = "\t", h=T)
dim(test)

## the function work well on windows as well even including the BLAST option.




### Need to check for MUMSA











