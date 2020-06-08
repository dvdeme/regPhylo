#' @title Convert a tax_report.txt file from the NCBI taxonomic facility to species lists used by GetSeqInfo_NCBI_taxid and GetSeq_BOLD functions.

#' @description This function converts the tax_report.txt file exported by the NCBI taxonomic facility
#' (\url{https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi}) into a list of two tables containing
#' the species lists to query DNA sequences in GenBank and associated databases
#' (through NCBI using the function \code{\link{GetSeqInfo_NCBI_taxid}}), and the BOLD database (using the
#' function \code{\link{GetSeq_BOLD}}).

#' @param input path to the "tax_report.txt" file exported by the NCBI taxonomic facility
#' (\url{https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi})

#' @return This function creates a list of two tables: the first table $SpList.NCBI has two columns,
#' the first one reports the NCBI taxid, and the second column contains the binomial species names.
#' The second table $SpList.BOLD also has two columns, the first column reports all binominal species names,
#' including synonyms used in NCBI taxon database that will be used to query BOLD,
#' and the second column reports the accepted species name that is going to be reported in the output table of \code{\link{GetSeq_BOLD}}.
#' Note: notice that in table $SpList.BOLD used to query BOLD, additional synonyms can be further added manually to those identified
#' in the NCBI taxon database.
#'
#' @export Taxreport2Sp.List
#'
#' @examples
#' \dontrun{
#' # To run the example, copy the input tax_report.txt file provided
#' # by the regPhylo package to a temporary directory created in the
#' # working directory.
#' src.dir = system.file("extdata/tax_export", package = "regPhylo")
#' dir.create("TempDir")
#' # Set up the path to the TempDir folder.
#' dest.dir = paste(getwd(), "/TempDir", sep="")
#' file.names <- dir(src.dir)
#' # Copy all the files stored in regPhylo/extdata/multi.align"
#' # into a temporary folder.
#' sapply(file.names, function(x) {
#' file.copy(from = paste(src.dir, x, sep = "/"),
#' to = paste(dest.dir, x, sep = "/"),
#' overwrite = FALSE) })
#'
#' # Run the function using this example of tax_report.txt file
#' # for 30 species exported by the NCBI taxonomic facility.
#' taxreport = Taxreport2Sp.List(input = "TempDir/tax_report.txt")
#' names(taxreport) # the first element of the list is the table
#' # for the NCBI search, the second is for the BOLD search.
#'
#' ## Table for NCBI search.
#' head(taxreport$SpList.NCBI)
#' dim(taxreport$SpList.NCBI) # one taxa did not have a taxid
#' # when the request was performed on the 7/03/2019.
#'
#' ## Table for BOLD search.
#' head(taxreport$SpList.BOLD)
#' dim(taxreport$SpList.BOLD) # two taxa have another preferred NCBI species name.
#' tail(taxreport$SpList.BOLD)
#' }

Taxreport2Sp.List = function(input = NULL){
# open the tax_report.txt file
input = read.delim(input, sep = "\t", header = TRUE)
input = input[,-c(2,4,6)]

### convert all columns as.character
input <- data.frame(lapply(input, as.character), stringsAsFactors=FALSE)

# Prepare the two column table for the NCBI.
nb.na = which(is.na(input$taxid)==TRUE)
if(length(nb.na) == 0){
  SpList.NCBI = input[, c(4, 2)]
} else {
SpList.NCBI = input[-nb.na, c(4,2)]
}
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
