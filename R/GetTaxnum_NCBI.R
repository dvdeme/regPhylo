  #' @title Retrieves NCBI taxa ID

#' @description This function retrieves the NCBI taxa ID from the NCBI taxonomy database
#' given a list of Gen sp values
#' @details This function requires httr
#' @return This function returns a table with the species names and the NCBI taxa ID
#' The output is retrieved from the following website: https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi



#' @param splist a character vector of species as "Genus species" to query the NCBI Taxonomy database.
#' 
#' 
#' @examples 
#' Create a species list and retrieve NCBI:txid codes
#' 
#' 
#' \dontrun{
#' Sp.List = c("Aplodactylus etheridgii", "Aplodactylus arctidens", "Asterorhombus filifer")
#' taxReport <- GetTaxnum_NCBI(splist = Sp.List)
#' write.table(taxReport, file="tax_report.txt", sep= "\t", quote = FALSE, row.names = FALSE)
#' 
#' }


#' @export GetTaxnum_NCBI

GetTaxnum_NCBI = function(splist = NULL) {
  #format list of names for POST using httr package
  string2post <- paste(splist, collapse="\n")
  r<-POST(url="https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi",
       encode="form",
       body=list(tax=string2post,
                 match=1,
                 button="Save in file"))
  #convert results of POST to a dataframe
  ncbiblock<- content(r, "text")[[1]]
  ncbitaxdf <- read.table(header=T, text=ncbiblock,sep='\t')
  return(ncbitaxdf)
} # End of GetTaxnum_NCBI

  


