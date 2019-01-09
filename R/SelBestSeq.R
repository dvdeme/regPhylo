#' @title Select the sequences based on geographic proximity
#' and/or sequence length, and export sequences (fasta format) and metadata for the selected gene regions
#'
#' @description This function selects and exports the sequences for a list of gene
#' regions for all the species present in the input data table. The selection of
#' the sequence for a given species and gene region can be optimised considering the most geographically proximate
#' sequence to the RefPoint. It can additionally, or alternatively, select a
#' sequence according to the length of the seqeunce.
#'
#' @details The choice based on the sequence length can be made using the longest sequence (option
#' SeqChoice = 'Longest') or the sequence with a length closest to the median
#' of the sequence length distribution for a given species and a given gene region
#' (option SeqChoice = 'Median') (See Antonelli et al. 2017, Systematic Biology DOI:
#' 10.1093/sysbio/syw066 for a similar approach).  The sequence length is computed
#' after discarding the potential indels ('-') and the ambiguous nucleotides.
#'
#' @return The function returns a table with all the information about the selected sequences
#' including 3 additional columns: 'SeqLengthNoIndelsNoAmbig' = Sequence length
#' after removing the indels and ambiguous nucleotides; 'Dist2RefPoint' = Great circle
#' distance in kilometers from the RefPoint; 'dist2med' = distance to the median
#' of the sequence length distribution. The alignment of the gene region are exported in
#' fasta format.

#' @param input a data table coming from the function \code{\link{Select.DNA}} or
#' \code{\link{SpeciesGeneMat.Bl}} (the data table should have 29 columns).
#' @param output the name of the output table (including the path if necessary,
#' but the folder must be created first).
#' @param RefPoint a vector with geographic coordinates with
#' longitude and latitude in decimal degrees WGS84 of the geographic reference
#' point of the focal area. (e.g. the reference point for New Zealand, is 174.7976
#' degrees longitude, and -41.3355 degrees latitude). If RefPoint is NULL then the selection
#' based on geographic proximity is disabled, and only the selection based on the length
#' of the sequence is performed.
#' @param perReposit if the name of the repository is provided (i.e. the same name as used in the function
#' 'Congr.NCBI.BOLD.perReposit.R' option 'perReposit') it implies that several
#' sequences are coming from a personal repository and the function will use the
#' ID of the sequence provided by the field 'Db_xref' as a unique ID for the
#' sequence. The function uses the first element (using ';' as separator) present
#' in the column 'Db_xref'.
#' @param Alignment if Alignment = 'T' then fasta files of the sequences for the selected gene regions are
#' exported (requiring the seqinr R package).
#' @param MaxSeq Maximum number of sequences to be exported for a species
#' when multiple sequences are available for a gene region
#' (sequences maximising the spatial optimisation and/or length criteria are exported first).
#' If MaxSeq = 'ALL', then all the sequences available for that gene region are exported.
#' @param gene.list a vector of selected gene regions (gene names have to
#' be consistent with the header of the table with the suffix '_SpDNA_Mat.txt'
#' exported by the function \code{\link{SpeciesGeneMat_Bl}}).  For example Cytochrome c
#' oxydase subunit 1, should have to be written 'co1' and not 'COI' or 'COX1'.
#' @param SeqChoice select the longest sequence after removing ambiguous nucleotides
#' (option SeqChoice = 'Longest') or the sequence with the length closest to the
#' median sequence length (option SeqChoice = 'Median').

#' @examples # Load the data table exported by the function Select.DNA
#' data(Seq.DF4) # The table is called "Seq.DF5" and constitutes the
#' # sixth object of the list "Seq.DF4".
#' Seq.DF5 = Seq.DF4$Seq.DF5
#'
#' # Run the function and export all the sequences in the alignment for
#' # each species and gene region.
#' \dontrun{
#' Seq.DF6=SelBestSeq(input = Seq.DF5, output = "Alig_Seq.DF5.All",
#' RefPoint = cbind(174.7976, -41.3355), perReposit = "My.Rep",
#' Alignment = T, MaxSeq = "ALL", gene.list = c("co1", "16srrna"),
#' SeqChoice = "Median")
#'
#' dim(Seq.DF6) # 10 sequences are reported.
#'
#'
#' # Run the function and export the best (the most proximal to the focal
#' # area i.e. NZ) sequences in the alignment for each species and gene region,
#' # selecting the sequence with a median sequence length.
#' Seq.DF7=SelBestSeq(input = Seq.DF5, output = "Alig_Seq.DF5.Best",
#' RefPoint = cbind(174.7976, -41.3355), perReposit = "My.Rep",
#' Alignment = T, MaxSeq =1, gene.list = c("co1", "16srrna"),
#' SeqChoice = "Median")
#'
#' dim(Seq.DF7) # 2 sequences are reported.
#'
#'
#' # Run the function and export the two most proximal sequences of the focal area (i.e. NZ)
#' # in the alignment for each species and gene region, selecting the longest sequence.
#' Seq.DF8=SelBestSeq(input = Seq.DF5, output = "Alig_Seq.DF5.Best",
#' RefPoint = cbind(174.7976, -41.3355), perReposit = "My.Rep",
#' Alignment = T, MaxSeq =2, gene.list = c("co1", "16srrna"),
#' SeqChoice = "Longest")
#'
#' dim(Seq.DF8) # 3 sequences are reported: 2 CO1 sequences, and the only 16srrna
#' # sequence available.
#'
#'
#' # Run the function to select the two longest sequences per gene region and species,
#' # regardless of geographic proximity.
#' Seq.DF9=SelBestSeq(input = Seq.DF5, output = "Alig_Seq.DF5.Best",
#' perReposit = "My.Rep", Alignment = T, MaxSeq = 2,
#' gene.list = c("co1", "16srrna"), SeqChoice = "Longest")
#'
#' dim(Seq.DF9) # 3 sequences are reported: 2 CO1 sequences and the only 16srrna
#' # sequence available.
#' }
#'
#' @references Antonelli et al. 2017, DOI: 10.1093/sysbio/syw066
#'
#'
#' @export SelBestSeq



SelBestSeq = function(input = NULL, output = NULL, RefPoint = NULL, perReposit = NULL,
    Alignment = NULL, MaxSeq = NULL, gene.list = NULL, SeqChoice = NULL) {
    DF = input
    perRepositID = rep(NA, dim(DF)[1])

    if (is.null(perReposit) == FALSE)
        {
            # If a personal repository is provided the function uses the first element (using
            # ';' as separator) present in the column 'Db_xref'.
            pos.Tepapa = grep(perReposit, DF[, 28])
            a = strsplit(as.character(DF[pos.Tepapa, ][, 17]), ";", fixed = T)
            ab = vector()
            i = 1
            for (i in 1:length(a)) {
                ab = c(ab, a[[i]][1])
            }
            perRepositID[pos.Tepapa] <- ab
        }  # End if
    DF = cbind(DF, perRepositID)


    Table_InfoTot = matrix(NA, ncol = 32)[-1, ]
    colnames(Table_InfoTot) = c(names(DF), "SeqLengthNoIndelsNoAmbig", "GreatCircleDistanceKm")

    k = 1
    for (k in 1:length(gene.list)) {
        # Loop over multiple genes.
        DF1 = as.data.frame(DF[grep(gene.list[k], DF[, 29], fixed = TRUE), ])
        SpUniq = unique(DF1[, 1])  # List of species.

        # Table info.
        Table_Info = matrix(NA, ncol = 32)[-1, ]
        colnames(Table_Info) = c(names(DF1), "SeqLengthNoIndelsNoAmbig", "GreatCircleDistanceKm")

        j = 1
        for (j in 1:length(SpUniq)) {
            # Loop for each species.
            DFtemp = DF1[which(DF1[, 1] == SpUniq[j]), ]
            SeqLengthNoIndelsNoAmbig = vector()
            SeqLengthNoIndels = vector()
            i = 1
            for (i in 1:dim(DFtemp)[1]) {
                # For each sequence per species.
                nbam = nchar(gsub("-", "", as.character(DFtemp[i, 4]))) - sum(stringr::str_count(tolower(as.character(DFtemp[i,
                  4])), c("a", "c", "t", "g")))  # Count the number of ambiguous nucleotides.
                SeqLengthNoIndelsNoAmbig = c(SeqLengthNoIndelsNoAmbig, nchar(gsub("-",
                  "", as.character(DFtemp[i, 4]))) - nbam)
                SeqLengthNoIndels = c(SeqLengthNoIndels, nchar(gsub("-", "", as.character(DFtemp[i,
                  4]))))
            }  # End for i.

            # Allowing the possibility to disable the selection based on geographic criterion.
            if(is.null(RefPoint)){
              Dist2NZ = rep(NA, dim(DFtemp)[1])
            } else {
            Dist2NZ = t(fields::rdist.earth(RefPoint, cbind(as.numeric(as.character(DFtemp[,
                24])), as.numeric(as.character(DFtemp[, 23]))), miles = FALSE))  # Great circle distance between the RefPoint and the sampling site of the sequence.
            }

            DFtemp[, 5] = SeqLengthNoIndels  # Sequence length excluding the indels.
            dist2med = abs(SeqLengthNoIndelsNoAmbig - stats::median(SeqLengthNoIndelsNoAmbig))  # Distance expressed as an absolute number of nucleotide differences between each sequence and the median sequence length.
            DFtemp1 = cbind(DFtemp, SeqLengthNoIndelsNoAmbig, Dist2RefPoint, dist2med)


            # Define the order of priorities to select the best sequence.  If none of the
            # sequence has geographic coordinates, just base the selection on SeqChoice
            # criteria, which can be the longest sequence or the one with a length closest to
            # the median length after excluding indels and ambiguous positions.  If the
            # sequence also has the geographic coordinates, then the sequence sampled from a
            # location closest to the RefPoint will be prioritised, and then the final choice
            # will be made according to the SeqChoice.
            if (length(stats::na.omit(DFtemp1[, 32])) == 0) {
                if (SeqChoice == "Longest") {
                  DFtempord = DFtemp1[order(as.numeric(as.character(DFtemp1[, 31])),
                    decreasing = TRUE), ]
                }
                if (SeqChoice == "Median") {
                  DFtempord = DFtemp1[order(as.numeric(as.character(DFtemp1[, 33])),
                    decreasing = FALSE), ]
                }
            } else {
                if (SeqChoice == "Longest")
                  {
                    DFtempord = DFtemp1[order(as.numeric(as.character(DFtemp1[, 32])),
                      -as.numeric(as.character(DFtemp1[, 31])), na.last = TRUE),
                      ]
                  }  # Order by geographic proximity and then by sequence length (the longest first and then decreasing length), and any NAs (attributed to the geographic proximity) are put at the end.
                if (SeqChoice == "Median")
                  {
                    DFtempord = DFtemp1[order(as.numeric(as.character(DFtemp1[, 32])),
                      as.numeric(as.character(DFtemp1[, 33])), na.last = TRUE), ]
                  }  # Order by geographic proximity then by sequence length (sequence with the closest distance to the median sequence length and then increasing distance), and any NAs (attributed to the geographic proximity) are put at the end.
            }


            if (MaxSeq == "All") {
                # Keep all the sequences.
                Table_Info = rbind(Table_Info, DFtempord)
            } else {
                if (dim(DFtemp)[1] < MaxSeq)
                  imax = dim(DFtemp)[1] else imax = MaxSeq
                Table_Info = rbind(Table_Info, DFtempord[1:imax, ])
            }
        }  # End for j.

        # Build the alignments if needed.
        if (Alignment == T)
            {
                Seq = as.character(Table_Info[, 4])
                Seq[which(nchar(Seq) <= 2)] <- seqinr::c2s(rep("-", 100))  # All individuals without a sequence (here signalled by a sequence of length '2' because NA has two characters) are replaced by a blank sequence (repetition of -), this should not happened.
                if (is.null(perReposit) == FALSE) {
                  Seq_Name = paste(">", paste(gsub(" ", "_", Table_Info[, 1], fixed = TRUE),
                    "|", Table_Info[, 2], "|", Table_Info[, 3], "|", Table_Info[,
                      30], "_%_", Seq, sep = ""), sep = "")  # Name of the sequence and the sequence in one string.
                } else {
                  Seq_Name = paste(">", paste(gsub(" ", "_", Table_Info[, 1], fixed = TRUE),
                    "|", Table_Info[, 2], "|", Table_Info[, 3], "_%_", Seq, sep = ""),
                    sep = "")  # Name of the sequence and the sequence in one string
                }  # End else.
                AlignFasta = unlist(strsplit(Seq_Name, "_%_", Seq_Name, fixed = TRUE))
                write(AlignFasta, file = paste(output, "_", gene.list[k], ".fas",
                  sep = ""))  # Write the alignments.
            }  # End if Alignment==T.

        # Feed the large table with the selected sequences and their associated
        # information for all the species and each gene.
        Table_InfoTot = rbind(Table_InfoTot, Table_Info)
    }  # End for k (multiple genes).

    # Write the table of selected sequences with all the information.
    if (is.null(perReposit) == FALSE) {
      utils::write.table(Table_InfoTot, file = paste(output, "_InfoTabSelSeq.txt", sep = ""),
            quote = FALSE, sep = "\t", row.names = FALSE)
    } else {
        Table_InfoTot = Table_InfoTot[, -30]
        utils::write.table(Table_InfoTot, file = paste(output, "_InfoTabSelSeq.txt", sep = ""),
            quote = FALSE, sep = "\t", row.names = FALSE)
    }
    return(Table_InfoTot)
}  # End of the function
