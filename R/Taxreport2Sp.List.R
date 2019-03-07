#' @title Convert a tax_report.txt file from the NCBI taxonomic facility to species lists used for GetSeqInfo_NCBI_taxid and GetSeq_BOLD functions.

#' @description This function converts the tax_report.txt file exported by the NCBI taxonomic facility
#' (https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi) into a list of two tables with
#' two columns with the species list appropriate to query DNA sequences in GenBank and associated databases
#' (through NCBI) using the function GetSeqInfo_NCBI_taxid, and in BOLD database using the function GetSeq_BOLD.

#' @param input path to the tax_report.txt file exported by the NCBI taxonomic facility (https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi)

#' @return The function reports a list of two tables: the first table $SpList.NCBI reports two columns the first one reports the NCBI taxid, and the second the binomial species names.
#' the second table $SpList.BOLD reports two columns the first one reports all binominal species names potentially including synonyms used in NCBI that will be used as query in BOLD, and the second table #' reports the accepted species name that is going to be reported in the output table of the GetSeq_BOLD.


Taxreport2Sp.List = function(input = NULL){
# open the tax_report.txt file
input = read.delim(input, sep = "\t", h=T)
input = input[,-c(2,4,6)]

### convert all columns as.character
input <- data.frame(lapply(input, as.character), stringsAsFactors=FALSE)

# Prepare the two column table for the NCBI.
SpList.NCBI = input[-which(is.na(input$taxid)==TRUE), c(4,2)]
colnames(SpList.NCBI) = c("taxid", "Sp.names")

# Prepare the two column table for BOLD.

## Detect when another preferred name is used.
P1 = seq(1, dim(input)[1])[-grep("^ $", input$preferred.name, perl = TRUE)]

if(length(P1) > 0){
SpList.BOLD = rbind(cbind(input$name, input$name), cbind(input$preferred.name[P1], input$name[P1]))

} else {
SpList.BOLD = cbind(input$name, input$name)
}

colnames(SpList.BOLD) = c("SpName.Bold.search", "Sp.names")
return(list(SpList.NCBI = SpList.NCBI, SpList.BOLD = SpList.BOLD))
}
