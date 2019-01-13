#' @title Extract DNA sequences and metadata from BOLD (Barcode Of Life Database)
#'

#' @description This function extracts sequences and associated information from the Barcode Of
#' Life Database (BOLD) for all the species in the species list
#' @param splist a two column table, the first column contains all the binomial species name including
#' synonymes and the second column contains only the names that will be reported
#' in the output table as 'species_name'.
#' @param filename name of the output table. The format of the output follows the
#' format provided by the bold R package.
#'
#' @return The function returns 3 tables: a first table with the raw output from bold package plus
#' an additional column with the date of extraction is exported into the working directory
#' (this table reports only metadata when the DNA sequence is present), and a
#' list of two tables exported in the R environment, a summary table reporting the number of
#' occurrences retrieved for each species, and the second table reporting all
#' the data retrieved from bold (including metadata in absence of DNA sequence,
#' as such the data present in this table and the table exported in the working directory might differ)

#' @examples # a table with 1 species
#' Splist = cbind(TaxID = c("Diastobranchus capensis"), Species.Name = c("Diastobranchus capensis"))
#' \dontrun{
#' # Run the function considering extracting all DNA sequences and associated metadata
#' BOLD.output = GetSeq_BOLD(splist = Splist, filename = "output.BOLD.txt")
#'
#' # The output can be loaded doing the following
#' data(Seq.Diastocapen)
#' BOLD.output = Seq.Diastocapen$Seq.BOLD
#' }

#' @export GetSeq_BOLD

GetSeq_BOLD = function(splist = NULL, filename = NULL) {
    NameTab = c("species_name", "processid", "sampleid", "recordID", "catalognum",
        "fieldnum", "institution_storing", "collection_code", "bin_uri", "phylum_taxID",
        "phylum_name", "class_taxID", "class_name", "order_taxID", "order_name",
        "family_taxID", "family_name", "subfamily_taxID", "subfamily_name", "genus_taxID",
        "genus_name", "species_taxID", "subspecies_taxID", "subspecies_name", "identification_provided_by",
        "identification_method", "identification_reference", "tax_note", "voucher_status",
        "tissue_type", "collection_event_id", "collectors", "collectiondate_start",
        "collectiondate_end", "collectiontime", "collection_note", "site_code", "sampling_protocol",
        "lifestage", "sex", "reproduction", "habitat", "associated_specimens", "associated_taxa",
        "extrainfo", "notes", "lat", "lon", "coord_source", "coord_accuracy", "elev",
        "depth", "elev_accuracy", "depth_accuracy", "country", "province_state",
        "region", "sector", "exactsite", "image_ids", "image_urls", "media_descriptors",
        "captions", "copyright_holders", "copyright_years", "copyright_licenses",
        "copyright_institutions", "photographers", "sequenceID", "markercode", "genbank_accession",
        "nucleotides", "trace_ids", "trace_names", "trace_links", "run_dates", "sequencing_centers",
        "directions", "seq_primers", "marker_codes", "Date_Extract")

    TabTot = matrix(NA, ncol = 81)[-1, ]
    colnames(TabTot) = NameTab

    TabSum = matrix(NA, ncol = 2)[-1, ]

    cat(t(as.matrix(NameTab)), sep = "\t", "\n", file = filename)
    i = 1
    for (i in 1:dim(splist)[1]) {
        tab = bold::bold_seqspec(taxon = as.character(splist[i, 1]))
        if (is.null(dim(tab)[1])) {
            TabSum = rbind(TabSum, c(as.character(splist[i, 1]), 0))
        } else {
            if (dim(tab)[1] == 0) {
                TabSum = rbind(TabSum, c(as.character(splist[i, 1]), 0))
            } else {
                tab = cbind(rep(as.character(splist[i, 2]), dim(tab)[1]), tab[, -c(22)],
                  rep(date(), dim(tab)[1]))
                colnames(tab) = NameTab
                utils::write.table(tab, file = filename, sep = "\t", row.names = FALSE,
                  col.names = FALSE, append = T)
                TabTot = rbind(TabTot, tab)
                TabSum = rbind(TabSum, c(as.character(splist[i, 1]), dim(tab)[1]))
            }
        }
    }

    # Remove the rows without DNA sequences.
    TabAll = utils::read.delim(filename, sep = "\t", header = TRUE)
    NBChar = apply(TabAll, 1, function(x) {
        # Count the number of nucleotides for each sequence.
        nchar(as.character(x[72]))
    })
    if (length(which(NBChar == 0)) > 0) {
        ## When the sequence is empty, the row is removed from the table.
        TabAll2 = TabAll[-which(NBChar == 0), ]
    } else {
        TabAll2 = TabAll
    }
    if (length(which(is.na(TabAll2[, 72]) == TRUE)) > 0) {
        # remove the rows without sequences or an NA.
        TabAll2 = TabAll2[-which(is.na(TabAll2[, 72]) == TRUE), ]
    } else {
        TabAll2 = TabAll2
    }
    utils::write.table(TabAll2, file = filename, sep = "\t", row.names = FALSE)  # export the table and overwrite the previous one.

    colnames(TabSum) = c("Species_Name", "Nb_Occ")
    return(list(TabSum, TabTot))  # Report a table with all the species listed in the list and the number of occurrences found in BOLD (Note that some of occurrences are not associated with a sequence and so the total number of occurrences may be different from the number of occurrences retained in the output table (Only those with DNA sequences are retained)).
}  # End of the function.
