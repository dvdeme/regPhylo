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
#Seq.NCBI.info.redo = GetSeqInfo_NCBI_taxid(splist = SpList.NCBI.test, gene = "ALL", 
 #                                          filename = "~/Desktop/test1.txt", timeout = 10)


Seq.NCBI.info.redo2 = GetSeqInfo_NCBI_taxid2(splist = SpList.NCBI.test, gene = "ALL", 
                                          filename = "~/Desktop/test2.txt", timeout = 10)


#SANDBOX
ee = seqinr::query("ee", paste("tid=", 60557, sep = ""))
oo = rentrez::entrez_search(db="nucleotide",term = "txid60557 [Organism:exp]")

oo = tryCatch(rentrez::entrez_search(db="nucleotide",term = paste("txid", SpList.NCBI.test[1], " [Organism:exp]",sep = "")), error = function(e) e)
oo$count
oo.fetch = entrez_fetch(db='nucleotide', id=oo$ids,rettype = 'gb', retmode='xml',parsed=TRUE)
oo.xml <- XML::xmlToList(oo.fetch)
i=2


