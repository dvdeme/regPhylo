library(rentrez)
library(seqinr)

library(regPhylo)


setwd("/home/iggy/BioinformatIg/Patiriella_Kermadecs/2019/regPhylo_Pat")

# Run the function with the path to the file "tax_report.txt" 
SpList.DF = Taxreport2Sp.List(input = "redo.txt")

# Extraction of the species list for NCBI search (first object of the list)
SpList.NCBI = SpList.DF$SpList.NCBI
SpList.NCBI.test = head(SpList.NCBI, n = 1)

#THIS is the old function and obviously works
#212 secs
Seq.NCBI.info.redo = GetSeqInfo_NCBI_taxid(splist = SpList.NCBI.test, gene = "ALL",
                                          filename = "~/Desktop/test1.txt", timeout = 20)

#3 secs
Seq.NCBI.info.redo2 = GetSeqInfo_NCBI_taxid_rentrez(splist = SpList.NCBI.test, gene = "ALL", filename = "~/Desktop/test2.txt")


#SANDBOX
#Set up seqinr
seqinr::choosebank("genbank", timeout = 25)
ee = seqinr::query("ee", paste("tid=", 60557, sep = ""))

##################TO REPLICATE######################
i=2
seqinr::autosocket()  # Automatically select the last ACNUC database.
ab = seqinr::getAnnot(ee$req[[i]], nbl = 5000)  # Extract all the annotations of the sequence.
## Edit the annotation text provided by the NCBI for each sequence.
abb = gsub("(^[ ]+)", "", ab, perl = TRUE)  # Remove the spaces at the start of each line.

ab
abb

# Extract the information from the LOCUS (row 1) to the FEATURES (but excluded)
# to split the information into two pieces of text (easier to deal with).  Define
# the boundaries where the information is stored for each section.
bound = grep("[A-Z]{3,}[ ]+", abb[c(1:grep("FEATURES", ab))], perl = TRUE)
boundinf = bound[c(1:length(bound) - 1)]
boundsup = bound[c(2:length(bound))]
m = vector()



####################################################




oo = rentrez::entrez_search(db="nucleotide",term = "txid60557 [Organism:exp]")

oo = tryCatch(rentrez::entrez_search(db="nucleotide",term = paste("txid", SpList.NCBI.test[1], " [Organism:exp]",sep = "")), error = function(e) e)
oo$count

# Fetch the sequences for all the records in a species as a tinyseq XML
oo.seqs = rentrez::entrez_fetch(db='nucleotide', id = oo$ids, rettype = 'fasta', retmode = 'xml', parsed = TRUE)
oo.seqs.xml = XML::xmlToList(oo.seqs)
imax = length(oo.seqs.xml)
# Fetch the annotation information for all the records in a species as genbank format text
oo.fetch = rentrez::entrez_fetch(db='nucleotide', id = oo$ids, rettype = 'gb', retmode = 'text')
oo.list = strsplit(oo.fetch, "//")[[1]]
oo.seqs.xml[i]$TSeq$TSeq_defline
as.numeric(oo.seqs.xml[i]$TSeq$TSeq_length)


obb
abb
identical('ab','ob')
identical(abb,obb)



oo.fetch2 = entrez_fetch(db='nucleotide', id=oo$ids,rettype = 'gb', retmode='xml',parsed=TRUE)
oo.xml = XML::xmlToList(oo.fetch2)
as.numeric(oo.xml[i]$GBSeq$GBSeq_length)
