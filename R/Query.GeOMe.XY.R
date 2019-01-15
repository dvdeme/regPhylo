#' @title Retrieve accurrate geographic coordinates from GeOMe database.
#'
#' @description  This function queries the GeOMe database to retrieve location
#' metadata (i.e. geographic coordinates) for sequences lacking geographic coordinates
#' (\url{https://www.geome-db.org/} (Deck et al. 2017).

#' @param input is a table coming from the function \code{\link{GeoCoord.WGS84}}.
#' @param Phylum is the name of the phylum
# (e.g. 'Chordata') that should be queried in GeOMe.
#' @param output name of the output table exported in the working directory.
#'
#' @details The function checks for species name synonymies using Worms
#' (\url{http://www.marinespecies.org}) through the worrms R package
#' (Chamberlain 2018, \url{https://CRAN.R-project.org/package=worrms}).
#' @details The names of the accepted phyla names in GeOMe are available
#' at the following web page: (\url{https://www.geome-db.org/template}).
#' Select any project, and refer to 'DEF' next to 'Phylum'.

#' @return Two identical tables are exported, one in the R environment
#' and one in the working directory with the output name.
#' The output is similar to the input table except that the fields
#' Lat Long have been updated with proper geographic coordinates if they are
#'  found in GeOMe, finally the field 'OriginDatabase' is amended by
#' '_XY-GeOMe' for the relevant sequences.
#'
#' @references Deck et al. 2017, DOI: 10.1371/journal.pbio.2002925

#' @examples # Import Seq.DF1
#' #' data(Seq.Diastocapen)
#' Seq.DF1 = Seq.Diastocapen$Seq.DF1
#' \dontrun{
#' # Run the function
#' Seq.DF2 = Query.GeOMe.XY.R(input = Seq.DF1, Phylum = "Chordata", output = "Seq.DF2.txt")
#' }
#'
#' @export Query.GeOMe.XY.R

Query.GeOMe.XY.R = function(input = NULL, Phylum = NULL, output=NULL) {

    # Select all the species without lat lon coordinates.  Add a temporary column
    # recording the order of the sequences.
    inputa = cbind(input, TempSeqOrder = seq(1, dim(input)[1], 1))
    inputb = inputa[which(is.na(inputa[, 22]) == "TRUE" & is.na(inputa[, 3]) == "FALSE"),
        c(1, 2, 3, 22, 23, 27)]

    # Extract all the metadata associated with 'Chordata' in GeOMe.
    testOme = geomedb::queryMetadata(expeditions = list(), query = Phylum, names = NULL)

    # Create a field combining genus and species name, that will be queried.
    SpNameGeOME = paste(testOme[, 12], " ", testOme[, 13], sep = "")

    ## A new table with the full species name.
    testOme2 = cbind(SpecieName = SpNameGeOME, testOme)

    # Build a table with the species name used in GeOMe, then look for the Worms
    # unique ID, and all the accepted synonymised names.
    SpNameGeOMEuni = unique(SpNameGeOME)

    options(show.error.messages = FALSE)  # Remove the errors message if a species does not have any synonyms.
    ResTab = do.call("rbind", lapply(seq(1, length(SpNameGeOMEuni), 1), function(x) {
        a = try(worrms::wm_name2id(name = as.character(SpNameGeOMEuni[x])))  # Retrieve the unique worms id.
        if (inherits(a, "try-error")) {
            a2 = c(as.character(SpNameGeOMEuni[x]))
            cbind(SpNameGeoME = rep(as.character(SpNameGeOMEuni[x]), length(a2)),
                wormsid = NA, SynonymisedName = NA)
        } else {
            a1 = try(worrms::wm_synonyms_(id = a))  # Retrieve the list of synonymised names for the species.
            if (inherits(a1, "try-error")) {
                a2 = as.character(SpNameGeOMEuni[x])
            } else {
                a2 = c(as.character(SpNameGeOMEuni[x]), a1$scientificname)
            }  ## End try a1.
            cbind(SpNameGeoME = rep(as.character(SpNameGeOMEuni[x]), length(a2)),
                wormsid = rep(a, length(a2)), SynonimisedName = a2)
        }  # End try a.
    }))
    options(show.error.messages = TRUE)  # Disable the suppression of error messages.


    # For each species without XY coordinates.
    SptoQuery = as.character(unique(inputb[, 1]))
    # Extract the species matching between the query and the GeOMe species list.
    SpTotarget = ResTab[stats::na.omit(match(SptoQuery, ResTab[, 3])), 3]

    if (length(SpTotarget) > 0) {
        # Extract all the accession numbers for the target species.
        DFTemp = do.call("rbind", lapply(seq(1, length(SpTotarget), 1), function(x) {
            inputb[which(inputb[, 1] == SpTotarget[x]), ]
        }))

        # For the target species, looks for each NCBI accession number without geographic
        # coordinates if there is a match in the GeOMe db.
        DFresult = do.call("rbind", lapply(seq(1, length(DFTemp[, 3]), 1), function(x) {
            # a=grep(as.character(DFTemp[x,3]), testOme2[, 45]) # Old version, GeOMe
            # increases the number of columns, so replace the 23/01/2018, by name of the
            # column see next row.
            a = grep(as.character(DFTemp[x, 3]), testOme2$associatedSequences)
            if (length(a) > 0) {
                cbind(testOme2[a, c(1, 4, 7, 6)], AccNB = as.character(DFTemp[x,
                  3]), GeomeOrder = x, OriginalSeqOrder = DFTemp[x, 6])
            }
        }))

        # Remove the duplicated sequences retrieved from the GeOMe database.
        DFresult = unique(DFresult)

        # Attribute the X and Y to the proper sequence in the original file.
        Lat = inputa[, 22]
        Long = inputa[, 23]
        Lat[DFresult[, 7]] = DFresult[, 4]
        Long[DFresult[, 7]] = DFresult[, 3]

        # Indicate in the field 'OriginDatabase', if the XY coordinates are coming from
        # GeOMe.
        OriginDatabase = as.character(inputa[, 26])
        OriginDatabase[DFresult[, 7]] = paste(OriginDatabase[DFresult[, 7]], "_XY-GeOMe",
            sep = "")
        inputc = cbind(input[, c(1:21)], Lat, Long, input[, c(24:25)], OriginDatabase)
    } else {
        inputc = inputa[, -dim(inputa)[2]]
    }
    utils::write.table(inputc, file = output, sep = "\t", row.names = FALSE)
    return(inputc)
}
