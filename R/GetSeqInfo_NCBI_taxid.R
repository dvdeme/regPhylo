#' @title Extract DNA sequences and metadata from GenBank

#' @description This function extracts all the sequences, accession numbers, and related
#' information for all the gene regions or a given gene region, from a list of species using
#' their NCBI taxid.
#' @details This function extracts data from an external server and can be slow.
#' @return This function returns a table with the species names and the number of
#' occurrences retrieved from GenBank in the R environment, and exports a table
#' with the DNA sequences and the metadata into the working directory.
#'
#' @return The table exported in the working directory has the following
#'  headers:
#' \itemize{
#' \item 'TaxaName' (Binomial species name included in the second column of the NCBI
#' species list),
#' \item 'AccessNb' (NCBI accession number),
#' \item 'Sequence' (The sequence itself, or 'Too_Long' if more than 5000bp),
#' \item 'SeqLength' (Sequence length),
#' \item 'Definition' (NCBI definition field),
#' \item 'OrganismClassif' (Classification of the
#' organism),
#' \item 'Source' (NCBI source),
#' \item 'Title' (Title provided by NCBI),
#' \item 'Authors' (Authors provided by NCBI),
#' \item 'Journal' (Journal provided by NCBI),
#' \item 'Pubmed' (Pubmed references provided by NCBI),
#' \item 'Year' (Year provided by NCBI),
#' \item 'Organism' (Organism provided by NCBI),
#' \item 'Organelle' (Organelle provided by NCBI),
#' \item 'Mol_type' (Mol_type provided by NCBI),
#' \item 'Db_xref' (Db_xref for NCBI),
#' \item 'Product' (Product provided by NCBI),
#' \item 'Genes' (Gene provided by NCBI),
#' \item 'Location' (Field Country in NCBI),
#' \item 'isolation_source' (isolation_source provided by NCBI),
#' \item 'Lat_lon' (Lat_lon provided by NCBI),
#' \item 'Collection_date' (Collection_date provided by NCBI),
#' \item 'Date_Extract' (Date of the extraction for the data)
#' }


#' @param splist a two column table with the unique NCBI taxonomic ID (i.e. taxid) in
#' the first column and the accepted binomial species name in the second column
#' which will be retained as 'TaxaName' in the output file.
#' @param gene can be a
#' particular gene region or 'ALL' if all gene regions are
#' sought.
#' @param filename name of the output table (also needs to include the
#' extension, e.g. ".txt")
#' @param chunk_size number of records to be downloaded at a time from genbank. 
#' If connection drops try decreasing this value to 20. Default 50.

#' @examples # A table with two species and their unique NCBI taxa ID
#' Splist=cbind(TaxID=c(443778,189923),
#' Species.Name=c("Diastobranchus capensis", "Synaphobranchus affinis"))
#' # Run the function to extract all DNA sequences and associated metadata
#' \dontrun{
#' NCBI.output = GetSeqInfo_NCBI_taxid(splist = Splist, gene = "ALL",
#' filename = "output.NCBI.txt")
#'
#' # The output can be loaded by doing the following
#' data(Seq.Diastocapen)
#' NCBI.output = Seq.Diastocapen$Seq.NCBI
#'
#' # To remove the file created while running the example do the following:
#' file.remove("output.NCBI.txt")
#'
#' }


#' @export GetSeqInfo_NCBI_taxid

GetSeqInfo_NCBI_taxid = function(splist = NULL, gene = NULL, filename = NULL, chunk_size=50) {
    # Extract the character of string from the right.
    substrRight <- function(x, n) {
        substr(x, nchar(x) - n + 1, nchar(x))
    }
    
    print(paste("Starting processing:", date()))
    # Empty file to start sequence deposition.
    write(paste("Starting processing:", date(), filename))
    
    # Summary table providing the number of sequences for each taxa.
    TabSpGenSum = matrix("NA", ncol = 2)[-1, ]
    colnames(TabSpGenSum) = c("TaxaName", "NbSequences")
    
    # Create a table to store all the information for: 'TaxaName', 'AccessNb',
    # 'Sequence', 'SeqLength', 'Definition', 'OrganismClassif', 'Source', 'Title',
    # 'Authors', 'Journal', 'Pubmed', 'Year', 'Organism', 'Organelle', 'Mol_type',
    # 'Db_xref', 'Product', 'Genes', 'Location','Lon_lat', 'Collection_date',
    # 'Date_Extract_NCBI'.
    TabSpGenTot = matrix(NA, ncol = 23)[-1, ]
    colnames(TabSpGenTot) = c("TaxaName", "AccessNb", "Sequence", "SeqLength", "Definition",
                              "OrganismClassif", "Source", "Title", "Authors", "Journal", "Pubmed", "Year",
                              "Organism", "Organelle", "Mol_type", "Db_xref", "Product", "Genes", "Location",
                              "isolation_source", "Lat_lon", "Collection_date", "Date_Extract")
    utils::write.table(TabSpGenTot, file = filename, sep = "\t", quote = FALSE)
    
    
    # Loop over species.
    for (k in 1:dim(splist)[1]) {
        if (gene == "ALL") {
            # If gene='ALL' all the DNA sequences of that species are recorded.
            oo = rentrez::entrez_search(db="nucleotide",
                                        term = paste("txid", splist[k, 1]," [Organism:exp]",sep = ""),
                                        use_history=TRUE)
            
        } else {
            # # If the name of the gene is provided, the sequences of this gene of interest are recorded.
            oo = rentrez::entrez_search(db="nucleotide",
                                        term = paste("txid", splist[k, 1], " [Organism:exp]"," AND ",gene ,"[All Fields] ",sep = ""),
                                        use_history = TRUE)
        }
        # Filling up the summary table
        imax = as.numeric(oo$count)
        TabSpGenSum = rbind(TabSpGenSum, c(as.character(splist[k, 2]), imax))
        print(paste(as.character(splist[k, 2]), ": ", imax, "seq."))
        # Storing the information about the DNA sequences in a table.  Retrieve the
        # information only when at least 1 sequence is available.
        if (imax > 0){   
            for(seq_start in seq(0,imax,chunk_size)){ 
                #Loop to break up the genbank queries into chunks.
                # seq_start needs to start at 0 not 1
                # Fetch the sequences for records in a chunk as a Genbank XML
                oo.seqs = rentrez::entrez_fetch(db='nucleotide',
                                                web_history = oo$web_history, 
                                                rettype = 'gb', retmode = 'xml',
                                                retmax=chunk_size, retstart=seq_start,
                                                parsed = TRUE)
                oo.seqs.xml = XML::xmlToList(oo.seqs)
                # Fetch the annotation information for records in a chunk as genbank format text
                oo.fetch = rentrez::entrez_fetch(db='nucleotide', 
                                                 web_history = oo$web_history, 
                                                 rettype = 'gb', retmode = 'text',
                                                 retmax = chunk_size, retstart=seq_start)
                oo.list = strsplit(oo.fetch, "//\n\n")[[1]] # make a list of records
                for (i in 1:length(oo.list)) { #Loop through a chunk at a time
                    Infot = vector()
                    # Provide the accesion number provided by the list
                    Infot = c(Infot, as.character(splist[k, 2]), oo.seqs.xml[i]$GBSeq$`GBSeq_primary-accession`)
                    if (as.numeric(oo.seqs.xml[i]$GBSeq$GBSeq_length) > 5000) {
                        # If the sequence is longer than 5000bp, it is reported in the table as 'Too_Long'
                        Infot = c(Infot, "Too_Long")
                    } else {
                        # If the sequence length <5000bp, the sequence is directly reported in the table.
                        Infot = c(Infot, oo.seqs.xml[i]$GBSeq$GBSeq_sequence)
                    }
                    Infot = c(Infot, as.numeric(oo.seqs.xml[i]$GBSeq$GBSeq_length))
                    ab = strsplit(oo.list[i], "\n")[[1]] # Extract all the annotations of the sequence.
                    ## Edit the annotation text provided by the NCBI for each sequence.
                    abb = gsub("(^[ ]+)", "", ab, perl = TRUE)  # Remove the spaces at the start of each line.
                    # Extract the information from the LOCUS (row 1) to the FEATURES (but excluded)
                    # to split the information into two pieces of text (easier to deal with).  Define
                    # the boundaries where the information is stored for each section.
                    bound = grep("[A-Z]{3,}[ ]+", abb[c(1:grep("FEATURES", ab))], perl = TRUE)
                    boundinf = bound[c(1:length(bound) - 1)]
                    boundsup = bound[c(2:length(bound))]
                    m = vector()
                    for (z in 1:c(length(bound) - 1)) {
                        if ((boundsup[z] - boundinf[z]) > 1) {
                            m = c(m, paste(abb[boundinf[z]:c(boundsup[z] - 1)], sep = " ",
                                           collapse = " "))
                        } else {
                            m = c(m, abb[bound[z]])
                        }
                    }
                    
                    # Extract the: Definition OrganismClassif Source Authors Title Journal Pubmed
                    # year
                    
                    Infog = c("DEFINITION", "ORGANISM", "SOURCE", "TITLE", "AUTHORS",
                              "JOURNAL", "PUBMED")
                    for (j in 1:length(Infog)) {
                        if (length(grep(Infog[j], m)) > 0) {
                            Infot = c(Infot, gsub("(^[ ]+)", "", gsub(Infog[j], "", m[grep(Infog[j],
                                                                                           m, fixed = TRUE)][1], fixed = TRUE), perl = TRUE))
                        } else {
                            Infot = c(Infot, NA)
                        }
                    }
                    Infot = c(Infot, substrRight(m[grep("LOCUS", m)], 4))
                    
                    # Include the information of FEATURES.  Need to extract the information about:
                    # organism organelle mol_type db_xref product genes location isolation_source
                    # lon_lat collection_date
                    ftr = abb[c(grep("FEATURES", ab):length(abb))]
                    Infof = c("organism", "organelle", "mol_type")
                    for (j in 1:length(Infof)) {
                        if (length(ftr[grep(paste(Infof[j], "=", sep = ""), ftr, fixed = TRUE)]) >
                            0) {
                            Infot = c(Infot, gsub("\"", "", gsub(paste("/", Infof[j], "=\"",
                                                                       sep = ""), "", ftr[grep(paste(Infof[j], "=", sep = ""), ftr,
                                                                                               fixed = TRUE)], fixed = TRUE)[1], fixed = TRUE))
                        } else {
                            Infot = c(Infot, NA)
                        }
                    }
                    Infoh = c("db_xref", "product", "gene")
                    for (j in 1:length(Infoh)) {
                        if (length(ftr[grep(paste(Infoh[j], "=", sep = ""), ftr, fixed = TRUE)]) >
                            0) {
                            Infot = c(Infot, paste(unique(gsub("\"", "", gsub(paste("/",
                                                                                    Infoh[j], "=\"", sep = ""), "", ftr[grep(paste(Infoh[j],
                                                                                                                                   "=", sep = ""), ftr, fixed = TRUE)], fixed = TRUE), fixed = TRUE)),
                                                   sep = "", collapse = "; "))
                        } else {
                            Infot = c(Infot, NA)
                        }
                    }
                    # Test if the country is present in FEATURES, if yes the object Location will
                    # record it.
                    if (length(ftr[grep("country=", ftr, fixed = TRUE)]) > 0) {
                        bound2 = grep("^[/]", ftr, perl = TRUE)  # Check if the information about the country is present on more than one row.
                        # Need to test if the country is the last item provided in FEATURES.
                        last = bound2[length(bound2)] - grep("/country=\"", ftr, perl = TRUE)
                        if (last == 0) {
                            Infot = c(Infot, gsub("\"", "", gsub("/country=\"", "", paste(ftr[grep("/country=\"",
                                                                                                   ftr, perl = TRUE)]), fixed = TRUE), fixed = TRUE))
                        } else {
                            lim = (bound2[which(match(bound2, grep("/country=\"", ftr,
                                                                   perl = TRUE)) == 1) + 1] - grep("/country=\"", ftr, perl = TRUE))  # Check if the information about the country is present on more than one row.
                            if (lim > 1) {
                                countr = paste(ftr[grep("/country=\"", ftr, perl = TRUE):c(grep("/country=\"",
                                                                                                ftr, perl = TRUE) + lim - 1)], sep = " ", collapse = " ")
                                Infot = c(Infot, gsub("\"", "", gsub("/country=\"", "", countr,
                                                                     fixed = TRUE), fixed = TRUE))
                            } else {
                                Infot = c(Infot, gsub("\"", "", gsub("/country=\"", "", paste(ftr[grep("/country=\"",
                                                                                                       ftr, perl = TRUE)]), fixed = TRUE), fixed = TRUE))
                            }
                        }  # End if(last==0)
                    } else {
                        Infot = c(Infot, NA)
                    }  # End if(length(ftr[grep('country=', ftr, fixed=TRUE)])>0)
                    
                    # Test if isolation_source= is present in FEATURES
                    if (length(ftr[grep("isolation_source=", ftr, fixed = TRUE)]) >
                        0) {
                        bound2 = grep("^[/]", ftr, perl = TRUE)  # Check if the information about the isolation_source is present on more than one row.
                        # Need to test if the isolation_source is the last item provided in FEATURES.
                        last = bound2[length(bound2)] - grep("/isolation_source=\"",
                                                             ftr, perl = TRUE)
                        if (last == 0) {
                            Infot = c(Infot, gsub("\"", "", gsub("/isolation_source=\"",
                                                                 "", paste(ftr[grep("/isolation_source=\"", ftr, perl = TRUE)]),
                                                                 fixed = TRUE), fixed = TRUE))
                        } else {
                            lim = (bound2[which(match(bound2, grep("/isolation_source=\"",
                                                                   ftr, perl = TRUE)) == 1) + 1] - grep("/isolation_source=\"",
                                                                                                        ftr, perl = TRUE))  # Check if the information about the country is present on more than one row.
                            if (lim > 1) {
                                Isosour = paste(ftr[grep("/isolation_source=\"", ftr, perl = TRUE):c(grep("/isolation_source=\"",
                                                                                                          ftr, perl = TRUE) + lim - 1)], sep = " ", collapse = " ")
                                Infot = c(Infot, gsub("\"", "", gsub("/isolation_source=\"",
                                                                     "", Isosour, fixed = TRUE), fixed = TRUE))
                            } else {
                                Infot = c(Infot, gsub("\"", "", gsub("/isolation_source=\"",
                                                                     "", paste(ftr[grep("/isolation_source=\"", ftr, perl = TRUE)]),
                                                                     fixed = TRUE), fixed = TRUE))
                            }
                        }  # End if(last==0)
                    } else {
                        Infot = c(Infot, NA)
                    }  # End if(length(ftr[grep('isolation_source=', ftr, fixed=TRUE)])>0)
                    
                    # Test if the Geographic coordinate are present in FEATURES.  Record the latitude
                    # and longitude together (transformation into the appropriate coordinates later).
                    if (length(ftr[grep("lat_lon=", ftr, fixed = TRUE)]) > 0)
                        Infot = c(Infot, gsub("\"", "", gsub("/lat_lon=\"", "", ftr[grep("lat_lon=",
                                                                                         ftr, fixed = TRUE)], fixed = TRUE), fixed = TRUE)) else Infot = c(Infot, NA)
                    if (length(ftr[grep("collection_date=", ftr, fixed = TRUE)]) >
                        0)
                        Infot = c(Infot, gsub("\"", "", gsub("/collection_date=\"", "",
                                                             ftr[grep("collection_date=", ftr, fixed = TRUE)], fixed = TRUE),
                                              fixed = TRUE)) else Infot = c(Infot, NA)  # Date of collection.
                    Infot = c(Infot, date())
                    # Save the sequence and the information in the output file (if the connection
                    # breaks then we can re-start the function from where it stopped, and the
                    # information already extracted is saved in the file).
                    utils::write.table(t(as.matrix(Infot)), file = filename, sep = "\t", row.names = FALSE,
                                       col.names = FALSE, append = TRUE, quote = FALSE)
                }  # End for i (loop over the sequences).
            }  # End for if(is.null(ee$nelem)==FALSE).
        } # End for for( seq_start in seq(1,imax,50))
    }  # End for k (loop over the species ).
    
    print(paste("Finish processing:", date()))
    return(TabSpGenSum)
    
}  # End of the function.





