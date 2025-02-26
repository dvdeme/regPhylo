#' @title Merge outputs from GenBank, BOLD and a personal repository
#' (if available) into a single table and remove potential duplicate
#' sequences.

#' @description This function assembles a large common dataframe for all the data retrieved
#' from GenBank (through the NCBI) and BOLD, and also from a personal repository, if available.
#' The function checks for the presence of duplicated sequences and
#' selects the most relevant information (i.e. selects the longest sequence, and
#' maximise the metadata information across the different sources for the
#' location, geographic coordinates, and collection date).

#' @return The information contained in the output file includes:
#' \itemize{
#' \item 'TaxaName' (Binomial species name included in the second column of the NCBI
#' and BOLD species list),
#' \item 'AccessBold' (BOLD recordID which is also identical to
#' sequenceID in BOLD),
#' \item 'AccessNCBI' (NCBI sequence accession number),
#' \item 'Sequence' (The nucleotide sequence itself, or 'Too_Long' if more than 5000bp),
#' \item 'SeqLength' (Sequence length),
#' \item 'Definition' (NCBI definition field, or for BOLD a combination of
#' 'TaxaName', the string 'BOLD' and the BOLD 'processid'),
#' \item 'OrganismClassif' (Classification of the organism),
#' \item 'Source' (NCBI source or for BOLD combination of 'TaxaName' and 'sequenceID'),
#' \item 'Title' (Title provided by NCBI),
#' \item 'Authors' (Authors provided by NCBI, or for BOLD 'collectors'),
#' \item 'Journal' (Journal provided by NCBI),
#' \item 'Pubmed' (Pubmed references provided by NCBI),
#' \item 'Year' (Year provided by NCBI, or for BOLD year of the 'collectiondate_start'),
#' \item 'Organism' (Organism provided by NCBI or for BOLD 'TaxaName'),
#' \item 'Organelle' (Organelle provided by NCBI),
#' \item 'Mol_type' (Mol_type provided by NCBI or 'genomic DNA' for BOLD),
#' \item 'Db_xref' (Db_xref for NCBI, or for BOLD combination of 'processid',
#' 'sampleid', 'recordID', 'catalognum', 'fieldnum', 'institution_storing',
#' 'collection_code', separated by '; '),
#' \item 'Product' (Product provided by NCBI),
#' \item 'Genes' (Gene provided by NCBI or, for BOLD 'markercode'),
#' \item 'Location' (Field Country in NCBI, or for BOLD combination of 'country', 'province_state',
#' 'region', 'sector', 'exactsite'),
#' \item 'isolation_source' (isolation_source provided by NCBI),
#' \item 'Lat_lon' (Lat_lon provided by NCBI or for BOLD combination of 'lat' and 'lon'),
#' \item 'Collection_date' (Collection_date provided by NCBI, or for BOLD 'collectiondate_start'),
#' \item 'Date_Extract' (Date of the extraction of the data),
#' \item 'OriginDatabase' (Name of the Original database).
#' }

#' @param input.NCBI table coming from the function \code{\link{GetSeqInfo_NCBI_taxid}}.
#' @param input.BOLD table coming from the function \code{\link{GetSeq_BOLD}}.
#' @param output name of the output table in txt format.
#' @param input.perReposit a personal repository in a table format loaded
#' in the R environment with the following fields: 'TaxaName', 'AccessBold',
#' 'AccessNCBI', 'Sequence', 'SeqLength', 'Definition', 'OrganismClassif',
#' 'Source', 'Title', 'Authors', 'Journal', 'Pubmed', 'Year', 'Organism',
#' 'Organelle', 'Mol_type', 'Db_xref', 'Product', 'Genes', 'Location',
#' 'isolation_source', 'Lat_lon', 'Collection_date', 'Date_Extract'.  The fields
#' must follow the same definition as above, except that 'AccessBold',
#' 'AccessNCBI' should remain empty, and the internal unique id of the sequence
#' must be reported in the field 'Db_xref'. Any other information in this field
#' must be provided after the unique sequence id and separated by '; '.
#' @param perReposit if a name is provided then compare the NCBI and BOLD output to the
#' information contained in the personal repository, in case some sequences in the
#' personal repository have been already deposited in NCBI or BOLD.

#' @examples # Load the table with NCBI (GenBank) data
#' data(Seq.Diastocapen)
#' Seq.NCBI = Seq.Diastocapen$Seq.NCBI
#'
#' # Load the table with BOLD data
#' Seq.BOLD = Seq.Diastocapen$Seq.BOLD
#'
#' # Run the function with GenBank and BOLD data only.
#' AllSeqDF=Congr.NCBI.BOLD.perReposit(input.NCBI=Seq.NCBI,
#' input.BOLD=Seq.BOLD, output="AllSeq_NCBI_BOLD.txt")
#'
#' # Load a personal repository
#' Seq.Perso = Seq.Diastocapen$Seq.Perso
#'
#' # Run the function with GenBank, BOLD and the personal repository data.
#' AllSeqDF2=Congr.NCBI.BOLD.perReposit(input.NCBI=Seq.NCBI,
#' input.BOLD=Seq.BOLD, input.perReposit=Seq.Perso, perReposit="My.Rep",
#' output="AllSeq_NCBI_BOLD_perRep.txt")
#'
#' dim(AllSeqDF2)
#'


#' @export Congr.NCBI.BOLD.perReposit


Congr.NCBI.BOLD.perReposit = function(input.NCBI = NULL, input.BOLD = NULL, output = NULL,
    input.perReposit = NULL, perReposit = NULL) {

  if(is.null(input.BOLD)){
    input.BOLD2=matrix(NA, ncol=25)[-1,] # build an empty table.
  } else {
    # Remove the sequences without accession number in the BOLD ('processid)
    input.BOLD = input.BOLD[which(is.na(as.character(input.BOLD[, 2])) == "FALSE"),
                            ]
    # Change the BOLD table to match the NCBI table
    NoGapSeq = gsub("-","",as.character(input.BOLD[, 72]),fixed = TRUE) #Remove Gaps from BOLD seq
    SeqLength = nchar(NoGapSeq)
    Definition = paste(input.BOLD[, 1], "BOLD", input.BOLD[, 2], sep = " ")
    OrganismClassif = paste(input.BOLD[, 11], input.BOLD[, 13], input.BOLD[, 15],
                            input.BOLD[, 17], input.BOLD[, 19], input.BOLD[, 21],
                            input.BOLD[, 1], input.BOLD[,22], sep = "; ")
    Source = paste(input.BOLD[, 1], input.BOLD[, 69], sep = " ")
    Title = rep(NA, dim(input.BOLD)[1])
    Authors = as.character(input.BOLD[, 32])
    Journal = rep(NA, dim(input.BOLD)[1])
    Pubmed = rep(NA, dim(input.BOLD)[1])
    Year = substr(as.character(input.BOLD[, 33]), 1, 4)
    Organism = input.BOLD[, 1]
    Organelle = rep(NA, dim(input.BOLD)[1])
    Mol_type = rep("genomic DNA", dim(input.BOLD)[1])
    Db_xref = paste(input.BOLD[, 2], input.BOLD[, 3], input.BOLD[, 4],
                    input.BOLD[,5], input.BOLD[, 6], input.BOLD[, 7],
                    input.BOLD[, 8], sep = "; ")
    Product = rep(NA, dim(input.BOLD)[1])
    Genes = input.BOLD[, 70]
    Location = as.character(paste(input.BOLD[, 55], ": ",
                                  paste(input.BOLD[, 56], input.BOLD[, 57],
                                        input.BOLD[, 58], input.BOLD[, 59],
                                        sep = ", "), sep = ""))
    isolation_source = rep(NA, dim(input.BOLD)[1])
    Lat_lon = as.character(paste(input.BOLD[, 47], input.BOLD[, 48], sep = " "))
    Collection_date = input.BOLD[, 33]
    Date_Extract = input.BOLD[, 81]

    # Homogeneisation of the gene names in BOLD with NCBI.
    Genes = gsub("16S", "16srrna", Genes, fixed = TRUE)
    Genes = gsub("18S", "18srrna", Genes, fixed = TRUE)
    Genes = gsub("Rho", "rhod", Genes, fixed = TRUE)
    Genes = gsub("^ $", NA, Genes, perl = T)

    input.BOLD2 = cbind(input.BOLD[, c(1, 71)], NoGapSeq, SeqLength, Definition, OrganismClassif,
                        Source, Title, Authors, Journal, Pubmed, Year, Organism, Organelle, Mol_type,
                        Db_xref, Product, Genes, Location, isolation_source, Lat_lon, Collection_date,
                        Date_Extract)

    # Include the origin of the sequences, and add a column in each table specifying
    # the specific BOLD 'recordID' or 'sequenceID' which are the same.
    input.BOLD2 = cbind(input.BOLD2[, 1], AccessBold = as.character(input.BOLD[,
                                                                               69]),
                        input.BOLD2[, -1], OriginDatabase = rep("BOLD", dim(input.BOLD2)[1]))

  }

  if(is.null(input.NCBI)){
    input.NCBI=matrix(NA, ncol=25)[-1,] # build an empty table.
  } else {
    # Remove the sequences without accession number in the NCBI
    # databases, if any.
    input.NCBI = input.NCBI[which(is.na(as.character(input.NCBI[, 2])) == "FALSE"),
                            ]
    # Include the origin of the sequences, and add a column in each table specifying
    # the specific BOLD 'recordID' or 'sequenceID' which are the same.
    input.NCBI = cbind(input.NCBI[, 1], AccessBold = rep("NA", dim(input.NCBI)[1]),
                       input.NCBI[, -1], OriginDatabase = rep("NCBI", dim(input.NCBI)[1]))
  }

    NameCol = c("TaxaName", "AccessBold", "AccessNCBI", "Sequence", "SeqLength",
        "Definition", "OrganismClassif", "Source", "Title", "Authors", "Journal",
        "Pubmed", "Year", "Organism", "Organelle", "Mol_type", "Db_xref", "Product",
        "Genes", "Location", "isolation_source", "Lat_lon", "Collection_date", "Date_Extract",
        "OriginDatabase")

    colnames(input.NCBI) = NameCol
    colnames(input.BOLD2) = NameCol



    # If an input from a personal repository is provided, the function does the
    # following: removes duplicated sequences with BOLD and then merges into a single
    # table with the BOLD table, before further comparison with the NCBI table.
    if (is.null(perReposit) == "FALSE")
        {
            input.perReposit = cbind(input.perReposit, OriginDatabase = rep(perReposit,
                dim(input.perReposit)[1]))  ## add the origin of the database

            if(dim(input.BOLD2)[1]>0){ ### Check if the input.BOLD2 is an empty table

            # Test for the presence of duplicates sequences between the BOLD table and the
            # personnal repository using the BOLD 'processip' and the 'processip' of the
            # personal repository which has to be provided as the last element (using space
            # as separator) in the field 'Definition'.  Sub selection of the sequence coming
            # from BOLD.
            ee = strsplit(as.character(input.BOLD2[, 6]), " ", fixed = T)
            ef = strsplit(as.character(input.perReposit[, 6]), " ", fixed = T)

            processID_input.BOLD2 = vector()
            i = 1
            for (i in 1:length(ee)) {
                processID_input.BOLD2 = c(processID_input.BOLD2, ee[[i]][length(ee[[i]])])
            }
            processID_perReposit = vector()
            i = 1
            for (i in 1:length(ef)) {
                processID_perReposit = c(processID_perReposit, ef[[i]][length(ef[[i]])])
            }

            # Test the 'processID' in common in both tables.
            cc = intersect(processID_input.BOLD2, processID_perReposit)


            if (length(cc) > 0) {
                # if rows are duplicated, select the most informative one.
                input.BOLD2 = cbind(input.BOLD2, processID_input.BOLD2)
                input.perReposit = cbind(input.perReposit, processID_perReposit)
                # Select the duplicated rows in both tables.
                Dupli.input.BOLD2 = input.BOLD2[match(cc, input.BOLD2[, 26]), ]
                Dupli.input.perReposit = input.perReposit[match(cc, input.perReposit[,
                  26]), ]

                LocBOLD = gsub("NA: NA, NA, NA, NA", NA, Dupli.input.BOLD2[, 20],
                  fixed = TRUE)
                LatLong_perReposit = gsub("NA NA", NA, Dupli.input.perReposit[, 22],
                  fixed = TRUE)

                # In the case of duplicated rows we keep all the information from the public BOLD
                # database, except we check if additional information about the location can be
                # obtained and for any geographic coordinates.
                Location = vector()
                Lat_lon = vector()
                i = 1
                for (i in 1:length(cc)) {
                  # Test the best source of information for the location.
                  if (is.na(Dupli.input.perReposit[i, 20])) {
                    if (is.na(LocBOLD[i])) {
                      Location = c(Location, NA)
                    } else {
                      Location = c(Location, as.character(LocBOLD[i]))
                    }
                  } else {
                    Location = c(Location, as.character(Dupli.input.perReposit[i,
                      20]))
                  }
                  # Test the best source of information for the geographic coordinates.
                  if (is.na(LatLong_perReposit[i])) {
                    if (is.na(Dupli.input.BOLD2[i, 22])) {
                      Lat_lon = c(Lat_lon, NA)
                    } else {
                      Lat_lon = c(Lat_lon, as.character(Dupli.input.BOLD2[i, 22]))
                    }
                  } else {
                    Lat_lon = c(Lat_lon, LatLong_perReposit[i])
                  }
                }  # End for(i in 1:length(cc))

                # Assemble the new table.  Build the table for the duplicate sequences maximizing
                # the information.
                DupliDF = cbind(Dupli.input.BOLD2[, c(1:19)], Location, Dupli.input.BOLD2[,
                  21], Lat_lon, Dupli.input.BOLD2[, c(23:25)])
                colnames(DupliDF) = NameCol

                # Write the table after removing the duplicated sequences and extracting the
                # optimal information from the two sources of data.
                BOLDuniqDF = input.BOLD2[-match(cc, input.BOLD2[, 26]), ]  # Unique BOLD sequences
                perReposituniqDF = input.perReposit[-match(cc, input.perReposit[,
                  26]), ]  # Unique sequences from the personal repository.
                BOLDuniqDF = BOLDuniqDF[, -26]  # Remove the columns processID
                perReposituniqDF = perReposituniqDF[, -26]  # remove the columns processID

                DHTot = rbind(BOLDuniqDF, perReposituniqDF, DupliDF)  # Big table assembly.
                LocatClean = gsub("[ ]{3,}.+$", "", DHTot[, 20])  # Remove some code at the end of the line.
                input.BOLD2 = cbind(DHTot[, c(1:19)], Location = LocatClean, DHTot[,
                  c(21:25)])
            } else {
                input.BOLD2 = rbind(input.BOLD2, input.perReposit)  # Big table assembly.
            }  # end else length(cc)>0.
            } else {
              input.BOLD2 = rbind(input.BOLD2, input.perReposit)
            }
        }  # end if(perReposit=='TRUE').


    # Comparison between a single table merging BOLD and personal repository (if
    # provided) and NCBI data.

    # Select the duplicate sequences to be able to maximise the information for
    # Location, geographic coordinates and collection date.
    aa = intersect(as.character(input.NCBI[, 3]), as.character(input.BOLD2[, 3]))  # All the duplicated accession numbers among the two databases.
    if (length(aa) > 0) {
        DupliBOLD = input.BOLD2[match(aa, input.BOLD2[, 3]), ]
        DupliNCBI = input.NCBI[match(aa, input.NCBI[, 3]), ]
        LocBOLD = gsub("NA: NA, NA, NA, NA", NA, DupliBOLD[, 20], fixed = TRUE)
        LatLong_BOLD = gsub("NA NA", NA, DupliBOLD[, 22], fixed = TRUE)

        Sequence = vector()
        SeqLength = vector()
        Definition = vector()
        Product = vector()
        Genes = vector()
        Location = vector()
        Lat_lon = vector()
        Collection_date = vector()
        Date_Extract = vector()

         for (i in 1:length(aa)) {
           if (class(DupliNCBI[i,5])=="Factor") { #Need to change factor back to integer
             NCBI.len = as.integer(levels(DupliNCBI[i,5])[DupliNCBI[i,5]])
           } else{
             NCBI.len = DupliNCBI[i,5]
           }
           if (as.integer(NCBI.len) - DupliBOLD[i, 5] > 0) { 
                Sequence = c(Sequence, as.character(DupliNCBI[i, 4]))
                SeqLength = c(SeqLength, as.character(DupliNCBI[i, 5]))
                Definition = c(Definition, as.character(DupliNCBI[i, 6]))
                Product = c(Product, as.character(DupliNCBI[i, 18]))
                Genes = c(Genes, as.character(DupliNCBI[i, 19]))
                Date_Extract = c(Date_Extract, as.character(DupliNCBI[i, 24]))

            } else {
                Sequence = c(Sequence, as.character(DupliBOLD[i, 4]))
                SeqLength = c(SeqLength, as.character(DupliBOLD[i, 5]))
                Definition = c(Definition, as.character(DupliBOLD[i, 6]))
                Product = c(Product, as.character(DupliBOLD[i, 18]))
                Genes = c(Genes, as.character(DupliBOLD[i, 19]))
                Date_Extract = c(Date_Extract, as.character(DupliBOLD[i, 24]))

            }
            # Test the best source of information for the location.
            if (is.na(LocBOLD[i])) {
                if (is.na(DupliNCBI[i, 20])) {
                  Location = c(Location, NA)
                } else {
                  Location = c(Location, as.character(DupliNCBI[i, 20]))
                }
            } else {
                Location = c(Location, LocBOLD[i])
            }
            # Test the best source of information for the geographic coordinates.
            if (is.na(LatLong_BOLD[i])) {
                if (is.na(DupliNCBI[i, 22])) {
                  Lat_lon = c(Lat_lon, NA)
                } else {
                  Lat_lon = c(Lat_lon, as.character(DupliNCBI[i, 22]))
                }
            } else {
                Lat_lon = c(Lat_lon, LatLong_BOLD[i])
            }
            # Test the best source of information for the collection date.
            if (is.na(DupliBOLD[i, 23])) {
                if (is.na(DupliNCBI[i, 23])) {
                  Collection_date = c(Collection_date, NA)
                } else {
                  Collection_date = c(Collection_date, as.character(DupliNCBI[i,
                    23]))
                }
            } else {
                Collection_date = c(Collection_date, as.character(DupliBOLD[i, 23]))
            }

        }  # end for i


        OriginDatabase = rep("NCBI-BOLD", length(Sequence))

        # Build the table for the duplicate sequences maximizing the information.
        DupliDF = cbind(DupliBOLD[, c(1:3)], Sequence, SeqLength, Definition, DupliNCBI[,
            c(7:17)], Product, Genes, Location, DupliNCBI[, 21], Lat_lon, Collection_date,
            Date_Extract, OriginDatabase)
        colnames(DupliDF) = NameCol

        # Write the table after removing the duplicated sequences and extracting the
        # optimal information from the two sources of data.
        BOLDuniqDF = input.BOLD2[-match(aa, input.BOLD2[, 3]), ]  # unique BOLD/personal repository sequences.
        NCBIuniqDF = input.NCBI[-match(aa, input.NCBI[, 3]), ]  # Unique NCBI sequences.
        DHTot = rbind(NCBIuniqDF, BOLDuniqDF, DupliDF)  # Large table assembly.
        LocatClean = gsub("[ ]{3,}.+$", "", DHTot[, 20])  # Remove some code at the end of the line.
        DHTot = cbind(DHTot[, c(1:19)], Location = LocatClean, DHTot[, c(21:25)])
    } else {
        # No duplicated sequences.
        DHTot = rbind(input.NCBI, input.BOLD2)  # Large table assembly.
        LocatClean = gsub("[ ]{3,}.+$", "", DHTot[, 20])  # Remove some code at the end of the line.
        DHTot = cbind(DHTot[, c(1:19)], Location = LocatClean, DHTot[, c(21:25)])
    }  # End else if(length(aa)>0){
    DHTot[which(DHTot[, 3] == ""), 3] = NA

    utils::write.table(DHTot, file = output, sep = "\t", row.names = FALSE)
    return(DHTot)
}  # End of the function.



#'
#' @title Merge outputs from GenBank, BOLD (local queries) and a personal repository
#' (if available) into a single table and remove potential duplicate
#' sequences.

#' @description This function assembles a large common dataframe for all the data retrieved
#' from GenBank (through the NCBI) and BOLD (from a local database file 
#' that can be downloaded from https://boldsystems.org/data/data-packages/) ), and 
#' also from a personal repository, if available.
#' The function checks for the presence of duplicated sequences and
#' selects the most relevant information (i.e. selects the longest sequence, and
#' maximise the metadata information across the different sources for the
#' location, geographic coordinates, and collection date).

#' @return The information contained in the output file includes:
#' \itemize{
#' \item 'TaxaName' (Binomial species name included in the second column of the NCBI
#' and BOLD species list),
#' \item 'AccessBold' (concatenation of the BOLD processid and the marker_code = sequenceID in BOLDV5
#' web portal),
#' \item 'AccessNCBI' (NCBI sequence accession number),
#' \item 'Sequence' (The nucleotide sequence itself, or 'Too_Long' if more than 5000bp),
#' \item 'SeqLength' (Sequence length),
#' \item 'Definition' (NCBI definition field, or for BOLD a combination of
#' 'TaxaName', the string 'BOLD' and the BOLD 'processid'),
#' \item 'OrganismClassif' (Classification of the organism),
#' \item 'Source' (NCBI source or for BOLD combination of 'TaxaName', 'Voucher_type', 'processid', 
#' 'marker_code'),
#' \item 'Title' (Title provided by NCBI),
#' \item 'Authors' (Authors provided by NCBI, or for BOLD 'collectors'),
#' \item 'Journal' (Journal provided by NCBI),
#' \item 'Pubmed' (Pubmed references provided by NCBI),
#' \item 'Year' (Year provided by NCBI, or for BOLD year of the 'collection_date_start'),
#' \item 'Organism' (Organism provided by NCBI or for BOLD 'TaxaName'),
#' \item 'Organelle' (Organelle provided by NCBI),
#' \item 'Mol_type' (Mol_type provided by NCBI or 'genomic DNA' for BOLD),
#' \item 'Db_xref' (Db_xref for NCBI, or for BOLD combination of 'processid',
#'  'specimenid', sampleid', 'recordID', 'museumid', collection_code', 'inst', 'sovereign_inst' separated by '; '),
#' \item 'Product' (Product provided by NCBI),
#' \item 'Genes' (Gene provided by NCBI or, for BOLD 'marker_code'),
#' \item 'Location' (Field Country in NCBI, or for BOLD combination of 'country.ocean', 'province.state',
#' 'region', 'sector', 'site'),
#' \item 'isolation_source' (isolation_source provided by NCBI),
#' \item 'Lat_lon' (Lat_lon provided by NCBI or for 'coord' for BOLD),
#' \item 'Collection_date' (Collection_date provided by NCBI, or for BOLD 'collection_date_start'),
#' \item 'Date_Extract' (Date of the extraction of the data),
#' \item 'OriginDatabase' (Name of the Original database).
#' }

#' @param input.NCBI table coming from the function \code{\link{GetSeqInfo_NCBI_taxid}}.
#' @param input.BOLD table coming from the function \code{\link{GetSeq_BOLD}}.
#' @param output name of the output table in txt format.
#' @param input.perReposit a personal repository in a table format loaded
#' in the R environment with the following fields: 'TaxaName', 'AccessBold',
#' 'AccessNCBI', 'Sequence', 'SeqLength', 'Definition', 'OrganismClassif',
#' 'Source', 'Title', 'Authors', 'Journal', 'Pubmed', 'Year', 'Organism',
#' 'Organelle', 'Mol_type', 'Db_xref', 'Product', 'Genes', 'Location',
#' 'isolation_source', 'Lat_lon', 'Collection_date', 'Date_Extract'.  The fields
#' must follow the same definition as above, except that 'AccessBold',
#' 'AccessNCBI' should remain empty, and the internal unique id of the sequence
#' must be reported in the field 'Db_xref'. Any other information in this field
#' must be provided after the unique sequence id and separated by '; '.
#' @param perReposit if a name is provided then compare the NCBI and BOLD output to the
#' information contained in the personal repository, in case some sequences in the
#' personal repository have been already deposited in NCBI or BOLD.

#' @export Congr.NCBI.BOLD.perReposit.V2

Congr.NCBI.BOLD.perReposit.V2 = function(input.NCBI = NULL, input.BOLD = NULL, output = NULL,
                                         input.perReposit = NULL, perReposit = NULL) {
  
  if(is.null(input.BOLD)){
    input.BOLD2=matrix(NA, ncol=25)[-1,] # build an empty table.
  } else {
    # Remove the sequences without accession number in the BOLD ('processid)
    if(length(which(is.na(as.character(input.BOLD[, "processid"])))) > 0){
    input.BOLD = input.BOLD[which(is.na(as.character(input.BOLD[, "processid"])) == "FALSE"),]
    }
    
    # Change the BOLD table to match the NCBI table
    NoGapSeq = gsub("-","",as.character(input.BOLD[, "nuc"]), fixed = TRUE) #Remove Gaps from BOLD seq
    SeqLength = nchar(NoGapSeq)
    Definition = paste(input.BOLD[, "Species.names"], "BOLD", input.BOLD[, "processid"], sep = " ")
    OrganismClassif = paste(input.BOLD[, "phylum"], input.BOLD[, "class"], input.BOLD[, "order"],
                            input.BOLD[, "family"], input.BOLD[, "subfamily"], input.BOLD[, "genus"],
                            input.BOLD[, "species"], input.BOLD[,"taxid"], sep = "; ")
    Source = paste(input.BOLD[, "Species.names"], " ",  input.BOLD[, "voucher_type"], " ", input.BOLD[, "processid"], ".", input.BOLD[, "marker_code"], sep = "")
    Title = rep(NA, dim(input.BOLD)[1])
    Authors = as.character(input.BOLD[, "collectors"])
    Journal = rep(NA, dim(input.BOLD)[1])
    Pubmed = rep(NA, dim(input.BOLD)[1])
    Year = substr(as.character(input.BOLD[, "collection_date_start"]), 1, 4)
    Organism = input.BOLD[, "Species.names"]
    Organelle = rep(NA, dim(input.BOLD)[1])
    Mol_type = rep("genomic DNA", dim(input.BOLD)[1])
    Db_xref = paste(input.BOLD[, "processid"], input.BOLD[,"specimenid"], input.BOLD[, "sampleid"], 
                    input.BOLD[, "record_id"], input.BOLD[,"museumid"], input.BOLD[, "collection_code"], 
                    input.BOLD[, "inst"], input.BOLD[, "sovereign_inst"], sep = "; ")
    Product = rep(NA, dim(input.BOLD)[1])
    Genes = input.BOLD[, "marker_code"]
    Location = as.character(paste(input.BOLD[, "country.ocean"], ": ",
                                  paste(input.BOLD[, "province.state" ], input.BOLD[, "region"],
                                        input.BOLD[, "sector"], input.BOLD[, "site"],
                                        sep = ", "), sep = ""))
    isolation_source = rep(NA, dim(input.BOLD)[1])
    Lat_lon = gsub(",", "", gsub("]", "", gsub("[", "", input.BOLD[, "coord"], fixed = T), fixed = T), fixed = T)
    Collection_date = input.BOLD[, "collection_date_start"]
    Date_Extract = input.BOLD[, "Date_Extract"]
    
    # Homogeneisation of the gene names in BOLD with NCBI.
    Genes = gsub("16S", "16srrna", Genes, fixed = TRUE)
    Genes = gsub("18S", "18srrna", Genes, fixed = TRUE)
    Genes = gsub("Rho", "rhod", Genes, fixed = TRUE)
    Genes = gsub("^ $", NA, Genes, perl = T)
    
    input.BOLD2 = cbind(input.BOLD[, c("Species.names" , "insdc_acs")], NoGapSeq, SeqLength, Definition, OrganismClassif,
                        Source, Title, Authors, Journal, Pubmed, Year, Organism, Organelle, Mol_type,
                        Db_xref, Product, Genes, Location, isolation_source, Lat_lon, Collection_date,
                        Date_Extract)
    
    # Include the origin of the sequences, and add a column in each table specifying
    # the specific BOLD 'recordID' or 'sequenceID' which are the same.
    input.BOLD2 = cbind(input.BOLD2[, "Species.names"], AccessBold = paste(as.character(input.BOLD[,"processid"]), 
                        input.BOLD[,"marker_code"], sep = "."), input.BOLD2[, -1],
                        OriginDatabase = rep("BOLD", dim(input.BOLD2)[1]))
    
  }
  
  if(is.null(input.NCBI)){
    input.NCBI=matrix(NA, ncol=25)[-1,] # build an empty table.
  } else {
    # Remove the sequences without accession number in the NCBI
    # databases, if any.
    input.NCBI = input.NCBI[which(is.na(as.character(input.NCBI[, 2])) == "FALSE"),
    ]
    # Include the origin of the sequences, and add a column in each table specifying
    # the specific BOLD 'recordID' or 'sequenceID' which are the same.
    input.NCBI = cbind(input.NCBI[, 1], AccessBold = rep("NA", dim(input.NCBI)[1]),
                       input.NCBI[, -1], OriginDatabase = rep("NCBI", dim(input.NCBI)[1]))
  }
  
  NameCol = c("TaxaName", "AccessBold", "AccessNCBI", "Sequence", "SeqLength",
              "Definition", "OrganismClassif", "Source", "Title", "Authors", "Journal",
              "Pubmed", "Year", "Organism", "Organelle", "Mol_type", "Db_xref", "Product",
              "Genes", "Location", "isolation_source", "Lat_lon", "Collection_date", "Date_Extract",
              "OriginDatabase")
  
  colnames(input.NCBI) = NameCol
  colnames(input.BOLD2) = NameCol
  
  
  
  # If an input from a personal repository is provided, the function does the
  # following: removes duplicated sequences with BOLD and then merges into a single
  # table with the BOLD table, before further comparison with the NCBI table.
  if (is.null(perReposit) == "FALSE")
  {
    input.perReposit = cbind(input.perReposit, OriginDatabase = rep(perReposit,
                                                                    dim(input.perReposit)[1]))  ## add the origin of the database
    
    if(dim(input.BOLD2)[1]>0){ ### Check if the input.BOLD2 is an empty table
      
      # Test for the presence of duplicates sequences between the BOLD table and the
      # personnal repository using the BOLD 'processip' and the 'processip' of the
      # personal repository which has to be provided as the last element (using space
      # as separator) in the field 'Definition'.  Sub selection of the sequence coming
      # from BOLD.
      ee = strsplit(as.character(input.BOLD2[, 6]), " ", fixed = T)
      ef = strsplit(as.character(input.perReposit[, 6]), " ", fixed = T)
      
      processID_input.BOLD2 = vector()
      i = 1
      for (i in 1:length(ee)) {
        processID_input.BOLD2 = c(processID_input.BOLD2, ee[[i]][length(ee[[i]])])
      }
      processID_perReposit = vector()
      i = 1
      for (i in 1:length(ef)) {
        processID_perReposit = c(processID_perReposit, ef[[i]][length(ef[[i]])])
      }
      
      # Test the 'processID' in common in both tables.
      cc = intersect(processID_input.BOLD2, processID_perReposit)
      
      
      if (length(cc) > 0) {
        # if rows are duplicated, select the most informative one.
        input.BOLD2 = cbind(input.BOLD2, processID_input.BOLD2)
        input.perReposit = cbind(input.perReposit, processID_perReposit)
        # Select the duplicated rows in both tables.
        Dupli.input.BOLD2 = input.BOLD2[match(cc, input.BOLD2[, 26]), ]
        Dupli.input.perReposit = input.perReposit[match(cc, input.perReposit[,
                                                                             26]), ]
        
        LocBOLD = gsub("NA: NA, NA, NA, NA", NA, Dupli.input.BOLD2[, 20],
                       fixed = TRUE)
        LatLong_perReposit = gsub("NA NA", NA, Dupli.input.perReposit[, 22],
                                  fixed = TRUE)
        
        # In the case of duplicated rows we keep all the information from the public BOLD
        # database, except we check if additional information about the location can be
        # obtained and for any geographic coordinates.
        Location = vector()
        Lat_lon = vector()
        i = 1
        for (i in 1:length(cc)) {
          # Test the best source of information for the location.
          if (is.na(Dupli.input.perReposit[i, 20])) {
            if (is.na(LocBOLD[i])) {
              Location = c(Location, NA)
            } else {
              Location = c(Location, as.character(LocBOLD[i]))
            }
          } else {
            Location = c(Location, as.character(Dupli.input.perReposit[i,
                                                                       20]))
          }
          # Test the best source of information for the geographic coordinates.
          if (is.na(LatLong_perReposit[i])) {
            if (is.na(Dupli.input.BOLD2[i, 22])) {
              Lat_lon = c(Lat_lon, NA)
            } else {
              Lat_lon = c(Lat_lon, as.character(Dupli.input.BOLD2[i, 22]))
            }
          } else {
            Lat_lon = c(Lat_lon, LatLong_perReposit[i])
          }
        }  # End for(i in 1:length(cc))
        
        # Assemble the new table.  Build the table for the duplicate sequences maximizing
        # the information.
        DupliDF = cbind(Dupli.input.BOLD2[, c(1:19)], Location, Dupli.input.BOLD2[,
                                                                                  21], Lat_lon, Dupli.input.BOLD2[, c(23:25)])
        colnames(DupliDF) = NameCol
        
        # Write the table after removing the duplicated sequences and extracting the
        # optimal information from the two sources of data.
        BOLDuniqDF = input.BOLD2[-match(cc, input.BOLD2[, 26]), ]  # Unique BOLD sequences
        perReposituniqDF = input.perReposit[-match(cc, input.perReposit[,
                                                                        26]), ]  # Unique sequences from the personal repository.
        BOLDuniqDF = BOLDuniqDF[, -26]  # Remove the columns processID
        perReposituniqDF = perReposituniqDF[, -26]  # remove the columns processID
        
        DHTot = rbind(BOLDuniqDF, perReposituniqDF, DupliDF)  # Big table assembly.
        LocatClean = gsub("[ ]{3,}.+$", "", DHTot[, 20])  # Remove some code at the end of the line.
        input.BOLD2 = cbind(DHTot[, c(1:19)], Location = LocatClean, DHTot[,
                                                                           c(21:25)])
      } else {
        input.BOLD2 = rbind(input.BOLD2, input.perReposit)  # Big table assembly.
      }  # end else length(cc)>0.
    } else {
      input.BOLD2 = rbind(input.BOLD2, input.perReposit)
    }
  }  # end if(perReposit=='TRUE').
  
  
  # Comparison between a single table merging BOLD and personal repository (if
  # provided) and NCBI data.
  
  # Select the duplicate sequences to be able to maximise the information for
  # Location, geographic coordinates and collection date.
  aa = intersect(as.character(input.NCBI[, 3]), as.character(input.BOLD2[, 3]))  # All the duplicated accession numbers among the two databases.
  if (length(aa) > 0) {
    DupliBOLD = input.BOLD2[match(aa, input.BOLD2[, 3]), ]
    DupliNCBI = input.NCBI[match(aa, input.NCBI[, 3]), ]
    LocBOLD = gsub("NA: NA, NA, NA, NA", NA, DupliBOLD[, 20], fixed = TRUE)
    LatLong_BOLD = gsub("NA NA", NA, DupliBOLD[, 22], fixed = TRUE)
    
    Sequence = vector()
    SeqLength = vector()
    Definition = vector()
    Product = vector()
    Genes = vector()
    Location = vector()
    Lat_lon = vector()
    Collection_date = vector()
    Date_Extract = vector()
    
    for (i in 1:length(aa)) {
      if (class(DupliNCBI[i,5])=="Factor") { #Need to change factor back to integer
        NCBI.len = as.integer(levels(DupliNCBI[i,5])[DupliNCBI[i,5]])
      } else{
        NCBI.len = DupliNCBI[i,5]
      }
      if (as.integer(NCBI.len) - DupliBOLD[i, 5] > 0) { 
        Sequence = c(Sequence, as.character(DupliNCBI[i, 4]))
        SeqLength = c(SeqLength, as.character(DupliNCBI[i, 5]))
        Definition = c(Definition, as.character(DupliNCBI[i, 6]))
        Product = c(Product, as.character(DupliNCBI[i, 18]))
        Genes = c(Genes, as.character(DupliNCBI[i, 19]))
        Date_Extract = c(Date_Extract, as.character(DupliNCBI[i, 24]))
        
      } else {
        Sequence = c(Sequence, as.character(DupliBOLD[i, 4]))
        SeqLength = c(SeqLength, as.character(DupliBOLD[i, 5]))
        Definition = c(Definition, as.character(DupliBOLD[i, 6]))
        Product = c(Product, as.character(DupliBOLD[i, 18]))
        Genes = c(Genes, as.character(DupliBOLD[i, 19]))
        Date_Extract = c(Date_Extract, as.character(DupliBOLD[i, 24]))
        
      }
      # Test the best source of information for the location.
      if (is.na(LocBOLD[i])) {
        if (is.na(DupliNCBI[i, 20])) {
          Location = c(Location, NA)
        } else {
          Location = c(Location, as.character(DupliNCBI[i, 20]))
        }
      } else {
        Location = c(Location, LocBOLD[i])
      }
      # Test the best source of information for the geographic coordinates.
      if (is.na(LatLong_BOLD[i])) {
        if (is.na(DupliNCBI[i, 22])) {
          Lat_lon = c(Lat_lon, NA)
        } else {
          Lat_lon = c(Lat_lon, as.character(DupliNCBI[i, 22]))
        }
      } else {
        Lat_lon = c(Lat_lon, LatLong_BOLD[i])
      }
      # Test the best source of information for the collection date.
      if (is.na(DupliBOLD[i, 23])) {
        if (is.na(DupliNCBI[i, 23])) {
          Collection_date = c(Collection_date, NA)
        } else {
          Collection_date = c(Collection_date, as.character(DupliNCBI[i,
                                                                      23]))
        }
      } else {
        Collection_date = c(Collection_date, as.character(DupliBOLD[i, 23]))
      }
      
    }  # end for i
    
    
    OriginDatabase = rep("NCBI-BOLD", length(Sequence))
    
    # Build the table for the duplicate sequences maximizing the information.
    DupliDF = cbind(DupliBOLD[, c(1:3)], Sequence, SeqLength, Definition, DupliNCBI[,
                                                                                    c(7:17)], Product, Genes, Location, DupliNCBI[, 21], Lat_lon, Collection_date,
                    Date_Extract, OriginDatabase)
    colnames(DupliDF) = NameCol
    
    # Write the table after removing the duplicated sequences and extracting the
    # optimal information from the two sources of data.
    BOLDuniqDF = input.BOLD2[-match(aa, input.BOLD2[, 3]), ]  # unique BOLD/personal repository sequences.
    NCBIuniqDF = input.NCBI[-match(aa, input.NCBI[, 3]), ]  # Unique NCBI sequences.
    DHTot = rbind(NCBIuniqDF, BOLDuniqDF, DupliDF)  # Large table assembly.
    LocatClean = gsub("[ ]{3,}.+$", "", DHTot[, 20])  # Remove some code at the end of the line.
    DHTot = cbind(DHTot[, c(1:19)], Location = LocatClean, DHTot[, c(21:25)])
  } else {
    # No duplicated sequences.
    DHTot = rbind(input.NCBI, input.BOLD2)  # Large table assembly.
    LocatClean = gsub("[ ]{3,}.+$", "", DHTot[, 20])  # Remove some code at the end of the line.
    DHTot = cbind(DHTot[, c(1:19)], Location = LocatClean, DHTot[, c(21:25)])
  }  # End else if(length(aa)>0){
  DHTot[which(DHTot[, 3] == ""), 3] = NA
  
  utils::write.table(DHTot, file = output, sep = "\t", row.names = FALSE)
  return(DHTot)
}  # End of the function.



