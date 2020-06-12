
library(bold)
library(seqinr)
library(ape)
library(RJSONIO)
library(stringr)
library(fields)
library(parallel)
library(caper)
library(phytools)
library(geomedb)
library(regPhylo)
library(httr)


setwd("/home/iggy/BioinformatIg/Patiriella_Kermadecs/regPhylo_Tutorial/Tutorial/")
############ 1   ##########

# Load the species list and classification into R.
SpList.Classif = read.csv("SpeciesList_Classification.csv", h=TRUE)

# Extract the species list only.
Sp.List = SpList.Classif$SpeciesName

# Replace the "_" with a space between the genus and species name.
Sp.List = gsub("_", " ", Sp.List)
Sp.List[1:10] # Display the first 10 species.


# ############ 2 OLD WAY   ##########
# write.table(Sp.List, file="Sp.List_forNCBITaxo.txt", sep="\t", row.names = F,col.names = F, quote = F)
# #Submit the output to https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi and 
# #save output to tax_report.txt
# taxReport = read.delim("tax_report.txt",
#                        sep="\t",header = T,stringsAsFactors = FALSE)
# ############ 2 NEW WAY   ##########

taxReport <- GetTaxnum_NCBI(splist = Sp.List)
write.table(taxReport, file="tax_report.txt", sep= "\t", quote = FALSE, row.names = FALSE)

#Are any sp not in NCBI taxonomy database?
paste(taxReport$name[which(taxReport$code==3)])

# Run the function with the path to the file "tax_report.txt" exported by the 
# NCBI taxonomic web facility as input.
SpList.DF = Taxreport2Sp.List(input = "tax_report.txt")
SpList.DF$SpList.NCBI


# Extraction of the species list for NCBI search (first object of the list)
SpList.NCBI = SpList.DF$SpList.NCBI
head(SpList.NCBI, n = 3)
dim(SpList.NCBI)
SpList.BOLD = SpList.DF$SpList.BOLD
head(SpList.BOLD, n =3)
dim(SpList.BOLD) # 32 sp because 2 sp had other preferred name in ncbi

############ 3   ##########

#Fetch all seqs on Genbank
dir.create("Data_Extraction")
#This step takes a LONG TIME!
Seq.NCBI.info = regPhylo::GetSeqInfo_NCBI_taxid(splist = SpList.NCBI, gene = "ALL", 
                                    filename = "Data_Extraction/Seq.NCBI.txt")
  
Seq.NCBI.all = read.delim("Data_Extraction/Seq.NCBI.txt", sep = "\t", h = TRUE)
# #Fetch from Bold. Much faster compared to NCBI
Seq.BOLD.info = GetSeq_BOLD(splist = SpList.BOLD,
                            filename = "Data_Extraction/Seq.BOLD.txt")
Seq.BOLD=read.delim("Data_Extraction/Seq.BOLD.txt", sep="\t", h=T)
dim(Seq.BOLD) # 291 sequences are retrieved.

# Here we load a table with the same structure as the "Seq.NCBI.txt" 
# including 11 sequences with the associated metadata coming from a 
# personal repository.
# Seq.PerRep=read.delim("Data_Extraction/Seq.PerRep.txt", sep="\t", h=T)
# dim(Seq.PerRep) # 11 sequences.
#Concat all the results onto a single df
# AllSeqDF = Congr.NCBI.BOLD.perReposit(input.NCBI = Seq.NCBI.all, input.BOLD=Seq.BOLD,
#                                       output = "Data_Extraction/AllSeqDF.txt", 
#                                       input.perReposit = Seq.PerRep, perReposit = "PersonalRep")



AllSeqDF = Congr.NCBI.BOLD.perReposit(input.NCBI = Seq.NCBI.all, input.BOLD=Seq.BOLD,
                                      output = "Data_Extraction/AllSeqDF.txt")

dim(AllSeqDF) # 729 sequences in total after removing the duplicates
length(which(AllSeqDF$OriginDatabase == "NCBI")) # 427 sequences from NCBI
length(which(AllSeqDF$OriginDatabase == "BOLD")) # 176 sequences from BOLD
length(which(AllSeqDF$OriginDatabase == "NCBI-BOLD")) # 115 sequences duplicated 
# length(which(AllSeqDF$OriginDatabase == "PersonalRep")) # 115 sequences duplicated 

############ 4   ##########
#Improve geo metadata
dir.create("Geolocation")
AllSeqDF1 = GeoCoord.WGS84(input = AllSeqDF, output = "Geolocation/AllSeqDF_Geo1.txt")
dim(AllSeqDF1)
names(AllSeqDF1)
#get geo data from GeOMe. This example none found

#BROKEN?
# # ?Query.GeOMe.XY.R
# Phy="Echinodermata"
# # AllSeqDF2 = Query.GeOMe.XY.R(input = AllSeqDF1 , Phylum = "Chordata", 
#                              # output ="Geolocation/AllSeqDF_Geo2.txt")
# 
# 
# dftissue <- geomedb::queryMetadata('Tissue', query=paste('phylum=',Phy,sep=''), limit = 10000000)
# ID_info = unique(dftissue$Tissue[grep("KF012820", dftissue$Tissue$associatedSequences), c('materialSampleID', 'tissueID'),])
# dim(ID_info)[1] > 0
# #????
# # dfsamples <- geomedb::queryMetadata('Sample', query = 'phylum="Chordata"', limit = 10000000)
# # dfevent <- geomedb::queryMetadata('Event', query = 'phylum="Chordata"', limit = 10000000)
# # dftissue$Tissue$associatedSequences
# 
# # dftissue <- geomedb::queryMetadata('Tissue', query='phylum="Chordata"', limit = 10000000)
# AllSeqDF2 = read.delim("Geolocation/AllSeqDF_Geo2.txt", sep = "\t", h=T)
# length(grep("XY-GeOMe", AllSeqDF2$OriginDatabase))
#\BROKEN

#Refine vauge placename holders
# Run the function, without any correction for the place name.
AllSeqDF3a = GeoCodeName(input = AllSeqDF1, output = "Geolocation/AllSeqDF_Geo3a.txt")
# To detect the place name of the location that couldn't be found because 
# the location name requires some corrections (e.g "New Zealand; 100m 
# off-shore Muriwai beach" might be problematic and can easily be corrected
# to "New Zealand; Muriwai beach"), see the code below, and the first three examples.
as.character(unique(AllSeqDF3a$Location[AllSeqDF3a$Geo_accuracy == "NoLocationFound"]))[1:3]
# A two column table with the corrected place name can be loaded in R.
LocNameCorrected = read.delim("Geolocation/LocNameCorrected.csv", 
                              sep = "\t", h = T)
head(LocNameCorrected, n = 3)

# Run the function with the place correction.
AllSeqDF3 = GeoCodeName(input = AllSeqDF1, output = "Geolocation/AllSeqDF_Geo3.txt",
                        CorrTab = LocNameCorrected)
#Percent of DNA seqs with Geographic Coordinates
(table(AllSeqDF3$Geo_accuracy)/length(AllSeqDF3$Geo_accuracy))*100


############ 5   ##########
#Species by gene matrix and sequence cleaning
dir.create("CleanSeqPool")
# Prepare two lists with the accession number of the sequences that we would like 
# to remove from the pool of sequences.
# For BOLD 
BOLD.SeqTrash = c("3191061", "3214910")

# For GenBank
NCBI.SeqTrash = c("AB018233", "AF202547", "AF133061", 
                  "KP194660", "FJ896410", "EU366662")

# Run the SpeciesGeneMat.Bl function to get a clean pool of sequences and 
# an appropriate species-by-gene matrix.
Sp.DNAMat_cl=SpeciesGeneMat.Bl(input = AllSeqDF3, 
                               output = "CleanSeqPool/SpAll.DNA.Mat_cl_", 
                               NCBI.Trash = NCBI.SeqTrash, 
                               BOLD.Trash = BOLD.SeqTrash)
names(Sp.DNAMat_cl) # Names of the different elements of the list.
dim(Sp.DNAMat_cl$Species.Gene_matrix) # Dimension of the species-by-gene matrix.
# Display the first 5 columns (including the species name and the gene regions with the 
# best species coverage).
Sp.DNAMat_cl$Species.Gene_matrix[1:5, 1:5]
# Summary of the number of different gene regions ("NB_TypeDNA") and number of 
# sequences available for each species.
head(Sp.DNAMat_cl$Summary_Species, n = 3) 
# Summary of the number of different species for each gene region, and number of sequences.
head(Sp.DNAMat_cl$Summary_DNA, n = 5)

############ 6   ##########
#Select gene regions of interest
SelGene.Max = SelGene.MaxSpCov(input = Sp.DNAMat_cl$Species.Gene_matrix) 
#How many and which genes have full species coverage?
SelGene.Max$Minimum_Number_of_Gene_With_Full_Species_Coverage
SelGene.Max$List_Gene_Name_Full_Species_Coverage
# Based on the species-by-gene matrix, we select 2 mitochondrial 
# gene regions ("co1", "16srrna") and 3 nuclear gene regions 
# ("rag1", "myh6", "plagl2").
GeneSelection = c("co1", "16srrna", "rag1", "myh6", "plagl2")
Mat.overlap.GeneSelection = Matrix.Overlap(input = Sp.DNAMat_cl$Species.Gene_matrix, 
                                           gene.Sel = GeneSelection)
# The species overlap for each pairwise comparison of the selected gene regions is presented.
Mat.overlap.GeneSelection$NumberOfSpecies
# The average number of species that have a sequence and overlap among the selected gene regions.
diag(Mat.overlap.GeneSelection$NumberOfSpecies) = NA # Remove the diagonal.
mean(Mat.overlap.GeneSelection$NumberOfSpecies, na.rm = TRUE) 
#Missing data?
AmMissData(input = Sp.DNAMat_cl$Species.Gene_matrix, gene.list = GeneSelection)

# Load the pool of sequences exported by the SpeciesGeneMat.Bl function in r 
PoolSeq = read.delim("CleanSeqPool/SpAll.DNA.Mat_cl__CleanDataset.txt", sep="\t", h=T)
# Extract all the sequences (including from long annotated DNA sequences) and 
# associated information.
PoolAllSeq.5Genes = Select.DNA(input = PoolSeq, 
                               gene.list = GeneSelection, 
                               output = "CleanSeqPool/PoolAllSeq.5Genes", 
                               timeout = 15)
dim(PoolAllSeq.5Genes)

############ 7   ##########
#Export the best sequece per species and gene region based on seq len or geo criteria
dir.create("Alignments")
# Create the directory to store the alignments using geographic criteria.
dir.create("Alignments/OneBest.Geo")

# Run the function using sequence length and geographic proximity as the 
# criterion to select the sequence.
BestSeq.Geo.export = SelBestSeq(input = PoolAllSeq.5Genes, 
                                output = "Alignments/OneBest.Geo/Best.Geo", 
                                RefPoint = cbind(174.7976, -41.3355), perReposit = "PerRep", 
                                Alignment = T, MaxSeq = 1, gene.list = GeneSelection, 
                                SeqChoice = "Median")
dim(BestSeq.Geo.export) # 91 sequences have been exported.

# Create the directory to store the alignments, when the geographic criteria is disabled.
dir.create("Alignments/OneBest")
# Run the function disabling the geographic criteria to select the sequence.
BestSeq.export = SelBestSeq(input = PoolAllSeq.5Genes, 
                            output = "Alignments/OneBest/Best", perReposit = "PerRep", 
                            Alignment = T, MaxSeq = 1, gene.list = GeneSelection, 
                            SeqChoice = "Median")
dim(BestSeq.export) # 91 sequences have been exported.



############ 8   ##########
#Alignments
#Bug need to fix this. For now set parameters as objects in R environmnet.
output="Alignments/OneBest.Geo/MultiAlign"
input="Alignments/OneBest.Geo"
nthread=5
Multi.Align(input="Alignments/OneBest.Geo", output="Alignments/OneBest.Geo/MultiAlign", 
            nthread=5, methods = c("mafftfftns2", "mafftfftnsi", "muscle", 
                                   "prank"))

#Assess alignments with MUMSA
#Set Path in R. Rstudio does not read bash_profile $PATH
# Sys.setenv(PATH="/home/iggy/.local/bin:/home/iggy/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/home/iggy/.local/bin:/home/iggy/git_clones/mumsa-1.0/:/home/iggy/git_clones/Gblocks_0.91b")
MumsaRes = Mumsa.Comp(input = "Alignments/OneBest.Geo/MultiAlign", 
                      output = "Alignments/OneBest.Geo/Mumsa_Res_OneBest.Geo", 
                      remove.empty.align = TRUE)
MumsaRes # to display the table, with AOS and MOS score.
dir.create("Alignments/OneBest.Geo/ToTrim")
list.ali = list.files("Alignments/OneBest.Geo/MultiAlign")
list.mafft.ali = list.ali[grep("Mafftfftnsi", list.ali)]
file.copy(paste("Alignments/OneBest.Geo/MultiAlign/", list.mafft.ali, sep = ""), 
          "Alignments/OneBest.Geo/ToTrim")



############ 9   ##########
#Trim poorly Aligned regions
outtrimAl = Filtering.align.Trimal(input = "Alignments/OneBest.Geo/ToTrim", 
                                   output = "Alignments/OneBest.Geo/Trimmed_trimAl")
outtrimAl # to see the sequence length for the different alignments.



###########################WORK IN Progress after this###################################
# 
# 
# 
# # Four of the alignment are coding DNA and one is a non coding DNA (16srrna).
# list.ali = list.files("Alignments/OneBest.Geo/ToTrim")
# list.ali
# 
# # We prepare a vector with the Type of DNA for each alignment. 
# # The first alignment is the 16srrna, so it must 
# # be coded "d", all the others must be coded "c" (all coding DNA).
# 
# Type.ali = c("d", "c", "c", "c", "c")
# # We loop the function over the different alignments.
# ###########################WORK IN Progress after this
# # for(i in 1:length(list.ali)){
# #   outGblocks = Filtering.align.Gblocks(input="Alignments/OneBest.Geo/ToTrim", 
# #                                        target.file = list.ali[i],
# #                                        LessStringent="TRUE", Type=Type.ali[i], 
# #                                        output ="Alignments/OneBest.Geo/Trimmed_Gblocks", 
# #                                        remove.empty.align = TRUE)
# # }
# # outGblocks
# 
# 
# 
# ############ 10   ##########
# #Concat to single Supermatrix
# # Create the "ForConcat" folder.
# dir.create("Alignments/OneBest.Geo/ForConcat")
# list.ali = list.files("Alignments/OneBest.Geo/Trimmed_Gblocks")
# # Copy the files into the new folder. 
# file.copy(paste("Alignments/OneBest.Geo/Trimmed_Gblocks/", list.ali, sep = ""), 
#           "Alignments/OneBest.Geo/ForConcat")
# Align.Concat(input="Alignments/OneBest.Geo/ForConcat", 
#             Sp.List=NULL, 
#             outputConcat = "Alignments/OneBest.Geo/ForConcat/Concat")
# 
# 
# ############ 11   ##########
# #Run Partition Finder
# #problems running python in R. See the Tutorial for more info.
# 
# PartiFinder2(input = "Alignments/OneBest.Geo/ForConcat/Concat.fas", 
#              Partition = "Alignments/OneBest.Geo/ForConcat/Partitions_Concat.txt", 
#              codon = c(2:5), nexus.file = "Alignments/OneBest.Geo/ForConcat/Concat.nex", 
#              Path.PartiF2 = "/home/davidpc/Programs/PartitionFinder2/partitionfinder-2.1.1/PartitionFinder.py", 
#              branchlengths = "linked", models = "all", model_selection = "BIC", search = "greedy", 
#              Raxml = "TRUE", nthread = 5)
# bestpartition = readLines("Alignments/OneBest.Geo/ForConcat/Partitions_Concat.txt_PF2_all.txt")
# bestpartition # To display the 8 partitions that have been delineated.
# 
# ############ 12   ##########
# 
# ############ 13   ##########
# 
# ############ 14   ##########
# 
# ############ 15   ##########
# 
# ############ Appendix   ##########
