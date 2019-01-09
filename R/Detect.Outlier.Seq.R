#' @title Help to detect outlier sequences in an alignment

#' @description The function helps to relatively quickly detect some outlier sequences (miss
#' aligned sequences, which might be caused by different problems such as gene or
#' species annotation problems, presence of paralogs sequences...) that should be
#' potentially removed from the pool of sequences per species and gene regions
#' before selecting the best sequence.  The function is far from being accurate
#' and is only provided as a tool to help detecting potential outlier sequences.
#' Great care must be taken to the results that may erroneously include outlier
#' sequences or on the contrary omit outlier sequences. Eyes inspection and
#' checking must always be carry on to ensure a satisfaying results.


#' @details The function follows the following general strategy to detect outlier
#' sequences.  First a distance matrix is computed among all sequences (See
#' different strategies available below), second a BIONJs tree (i.e. an improved
#' Neighbhor joining allowing missing values, Criscuolo & Gascuel (2008)
#' using the disance matrix is build, and the root to tip distance is computed.
#' All the sequence above a certain distance threshold (we
#' use 50-60%, see Chen et al. 2015 for similar thresholds) of the maximun root
#' to tip distances are extracted as potential 'primary' outlier sequences.  In
#' option, we can blast the primary outlier sequences on a local blast database to
#' look for secondary outlier sequences that might be similar to the primary
#' sequences but just below the distance threshold. All the sequences with a
#' bitscore above a bitscore threshold (We usually use 80%) obtained when blasting
#' the query sequence against itself are considered as potential 'secondary'
#' outliers. We used the bitscore because it is not influenced by the size of the
#' local database contrary to the E.value
#' (https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html).

#' @return The function return a table for all the retained sequences with the following headers:
#' 'Query_SeqName', 'SeqLen_Query', 'Hit_SeqName', 'Hit_SeqLen', 'evalue',
#' 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
#' 'bitscore', 'qcovs' Name of the query sequence, length of the query sequence,
#' name of the hit sequence, length of the hit sequence, E-value, means Percentage
# 'of identical matches, means Alignment length, means Number of mismatches, means
#' Number of gap openings, means start of alignment in query, means End of
#' alignment in query, means Start of alignment in subject, means End of alignment
#' in subject, means Bit score, means Query Coverage Per Subject (for all HSPs)
#' (See https://www.ncbi.nlm.nih.gov/books/NBK279675/ for more details and
#' options)

#' @details Three options are available to compute the distance matrix among sequences.
#' The first otpion 'MisAli' uses the number of indelblocks to compute the
#' distance matrix in order to preferentially detect misaligned sequences
#' generating gaps in the alignment.  The second option 'DivSeq' uses the TN93
#' substitution model to compute the distance matrix in order to detect very
#' divergent sequences. The third option 'Comb' combined the two previous
#' distance matrices to detect both misaligned and very divergent sequence. The
#' two distance matrix are first re-scaled as percentage of the their highest
#' distances respectively and then summed.
#' @details The function requires BLAST installed and in the PATH.

#' @param inputal it can be an object of class "alignment" (seqinr R pcakage),
#' or "DNAbin" (ape R package) or the name (including the path if necessary) of the input fasta file.
#' @param Strat.DistMat can be 'MisAli', 'DivSeq' or 'Comb', see details.
#' @param Dist.Th distance threshold to detect the primary set of
#' outlier sequences (between [0-1]).
#' @param output name of the output table exported in the working directory.
#' @param Second.Outlier if 'TRUE', blasts the primary outlier seqeunces
#' against a local database based on the sequencs alignment to detect secondary outliers.
#' @param Bitsc.Th bitscore similarity threshold to detect 'secondary'
#' outlier sequences, should be a value between [0-1].

#' @examples # Load the alignment file (class "alignment"), this is the first
#' # object of the list called Example_16S_outlier.
#' data(Example_16S_outlier)
#' Example_16S_outlier_align = Example_16S_outlier[[1]]
#' \dontrun{
#'
#' # Running the function with the 'Comb' option,
#' # a distance threshold of 0.6, and disabling the search for secondary outlier sequences.
#' S16_MisAlign0.6_1 = Detect.Outlier.Seq(inputal = Example_16S_outlier_align,
#' Strat.DistMat = "Comb", Dist.Th = 0.6, output = "Example_16S_outliers_1.txt",
#' Second.Outlier = "No")
#' dim(S16_MisAlign0.6_1) ### 36 sequences detected as primary oultier sequences.
#' head(S16_MisAlign0.6_1)
#'
#' # The output table is also present as an external data table provided in the regPhylo r package
#' and can be access by the following code:
#' # a = system.file("extdata/ExampleOutliers/Example_16S_outliers_1.txt", package = "regPhylo")
#' # Example_16S_outliers_1 = read.delim(a, sep="\t", header = TRUE)
#'
#'
#' # Running the function with the 'Comb' option, a distance threshold of 0.6,
#' # and allowing a search for secondary outlier sequences using local blast database.
#' S16_MisAlign0.6_2 = Detect.Outlier.Seq(inputal = Example_16S_outlier_align,
#' Strat.DistMat = "Comb", Dist.Th = 0.6, output = "Example_16S_outliers_2.txt",
#' Second.Outlier = "Yes", Bitsc.Th = 0.8)
#' length(unique(S16_MisAlign0.6_2[,3])) ### 38 primary and secondary outlier seqeunces detected.
#' head(S16_MisAlign0.6_2)
#'
#' # The output table is also present as an external data table provided in the regPhylo r package
#' and can be accessed by the following code:
#' # a = system.file("extdata/ExampleOutliers/Example_16S_outliers_2.txt", package = "regPhylo")
#' # Example_16S_outliers_2 = read.delim(a, sep="\t", header = TRUE)
#'
#' # To clean the file created while running the example do the following:
#' file.remove("16S_example_RMoutliers.fas")
#'
#' }
#'
#' @references Criscuolo & Gascuel 2008, DOI: 10.1186/1471-21.5-9-166
#' @references Chen et al. 2015, DOI: 10.1093/sysbio/syv059

#' @export Detect.Outlier.Seq


Detect.Outlier.Seq = function(inputal = NULL, Strat.DistMat = NULL, Dist.Th = NULL,
    output = NULL, Second.Outlier = NULL, Bitsc.Th = NULL) {
  # if the class of the inputal object is "alignment" object recognised by ape and seqinr r packages.
  if(class(inputal)=="alignment"){
    cytb2=inputal
    cytb2a = ape::as.DNAbin(inputal) # Convert in a DNAbin object.
  }
  # if the path to the alignment file in fasta format is provided as inputal
  if(class(inputal)=="character"){
    cytb2 = seqinr::read.alignment(inputal, format = "fasta")
    cytb2a = ape::as.DNAbin(cytb2) # Convert in a DNAbin object.
  }
  # if the class of the inputal object is a "DNAbin" object.
  if(class(inputal)=="DNAbin"){
    cytb2a=inputal
    cytb2=ape::as.alignment(inputal) # convert to an "alignment" object.
  }

    if (Strat.DistMat == "DivSeq") {
        # Compute the distance between all sequences using the TN93 substitution model,
        # enabling pairwise.deletion option.
        a = ape::dist.dna(cytb2a, model = "TN93", pairwise.deletion = TRUE)
        if (is.na(mean(a, na.rm = T))) {
            aprop = a
            aprop[is.na(aprop)] = 0
        } else {
            aprop = a/max(a, na.rm = T)  # Re-scale the matrix in 0-1 ratio (percentage) based on its highest value.
            aprop[is.na(aprop)] = 1  # if NA then the max distance value is return to 1.
        }
        matdis = aprop
    }

    if (Strat.DistMat == "MisAli") {
        # Compute the distance between all sequences counting the number of indel blocks.
        b = ape::dist.dna(cytb2a, model = "indelblock")
        if (is.na(mean(b, na.rm = T))) {
            bprop = b
            bprop[is.na(aprop)] = 0
        } else {
            bprop = b/max(b, na.rm = T)  # Re-scale the matrix in 0-1 ratio (percentage) based on its highest value.
            bprop[is.na(bprop)] = 1  # if NA then the max distance value is return to 1.
        }
        matdis = bprop
    }

    if (Strat.DistMat == "Comb") {
        # Compute the distance between all sequences using the TN93 substitution model
        a = ape::dist.dna(cytb2a, model = "TN93", pairwise.deletion = TRUE)
        if (is.na(mean(a, na.rm = T))) {
            aprop = a
            aprop[is.na(aprop)] = 0
        } else {
            aprop = a/max(a, na.rm = T)  # Re-scale the matrix in 0-1 ratio (percentage) based on its highest value.
            aprop[is.na(aprop)] = 1  # if NA then the max distance value is return to 1.
        }
        # Compute the distance between all sequences counting the number of indel blocks.
        b = ape::dist.dna(cytb2a, model = "indelblock")
        if (is.na(mean(b, na.rm = T))) {
            bprop = b
            bprop[is.na(aprop)] = 0
        } else {
            bprop = b/max(b, na.rm = T)  # Re-scale the matrix in 0-1 ratio (percentage) based on its highest value.
            bprop[is.na(bprop)] = 1  # if NA then the max distance value is return to 1.
        }

        # Sum the two re-scaled distance matrices.
        matdis = aprop + bprop

    }


    # Build a BIONJs tree.
    anj2 = ape::bionjs(matdis)
    anj2$edge.length = abs(anj2$edge.length)  # To avoid very short negative branch length.

    # Compute distance root to tips for all the tips.
    distroot = diag(ape::vcv.phylo(anj2))


    # Detect all the tips with branch length > Dist.Th.
    distrootProp = distroot/max(distroot)
    Outliertips = labels(which(distrootProp > Dist.Th))

    if (!Second.Outlier == "Yes") {

        ## Warning if the number of 'primary' outliers represents more than 30% of the
        ## overall sequences in the alignments.
        if ((length(Outliertips)/length(anj2$tip.label)) > 0.3) {
            warning("Very high proportion of primary outlier sequences, ", round((length(Outliertips)/length(anj2$tip.label)) *
                100, 2), "%.", "\n", "Be cautious, the outlier sequences might not be true outlier sequences.",
                sep = "")
        }

        # Export the primary outlier sequences
        Outliertips = gsub(" no comment", "", Outliertips, fixed = T)  ## remove potential ' no comment' on the sequence name.
        DF = cbind(Outliertips, as.vector(distrootProp[which(distrootProp > Dist.Th)]))
        colnames(DF) = c("PrimaryOulierSeq", paste("Dist.Prop.", Strat.DistMat, sep = ""))
        utils::write.table(DF, file = output, sep = "\t", row.names = F)
        return(DF)

    } else {
        # if(Second.Outlier=='Yes')


        ## Stop if the number of 'primary' outliers represents more than 30% of the
        ## overall sequences in the alignments.
        if ((length(Outliertips)/length(anj2$tip.label)) > 0.3) {
            # Export the primary outlier sequences
            Outliertips = gsub(" no comment", "", Outliertips, fixed = T)  ## remove potential ' no comment' on the sequence name.
            DF = cbind(Outliertips, as.vector(distrootProp[which(distrootProp > Dist.Th)]))
            colnames(DF) = c("PrimaryOulierSeq", paste("Dist.Prop.", Strat.DistMat,
                sep = ""))
            utils::write.table(DF, file = output, sep = "\t", row.names = F)
            warning("Very high proportion of primary outlier sequences, be cautious, the outlier sequences might not be true outlier sequences.",
                "\n", "Search for secondary outlier sequences aborted.", sep = "")
            return(DF)
        }

        # Create a fasta file with only the outlier sequences and performed a local all
        # vs all BLASTN approach to reveal other potential 'secondary' outlier sequences.
        # Desaligned the sequences (remove all the indels) for BLAST.
        cytb2_desaligned = lapply(cytb2$seq, function(x) {
            gsub("-", "", x)
        })
        cytb2_des = cytb2
        cytb2_des$seq = cytb2_desaligned


        # Export the desaligned fasta sequences to prepare a local BLAST database
        output1 = paste(output, "ForLocalBlast_allseq.fas", sep = "")
        i = 1
        for (i in 1:length(cytb2_des$seq)) {
            cat(file = output1, append = T, paste(">", cytb2_des$nam[i], sep = ""),
                "\n", sep = "")
            cat(file = output1, append = T, cytb2_des$seq[[i]], "\n", sep = "")
        }

        # Export the desaligned outlier sequences in fasta format (=query sequences).
        Seq.outlier = cytb2_des$seq[match(Outliertips, cytb2_des$nam)]
        output2 = paste(output, "Seq.outlier.fas", sep = "")
        i = 1
        for (i in 1:length(Seq.outlier)) {
            cat(file = output2, append = T, paste(">", Outliertips[i], sep = ""),
                "\n", sep = "")
            cat(file = output2, append = T, Seq.outlier[[i]], "\n", sep = "")
        }

        # Build a local BLAST database.
        aa = paste("makeblastdb -in ", output1, " -dbtype nucl -out ", paste(gsub(".fas",
            "", output1, fixed = T), "_BlastDB", sep = ""), sep = "")
        system(aa)

        # All vs all Blastn using the outlier sequences as query sequences on the local
        # BLAST database.
        aa = paste("blastn -db ", paste(gsub(".fas", "", output1, fixed = TRUE), "_BlastDB",
            sep = ""), " -query ", output2, " -out ", paste(gsub(".fas", "", output2),
            "_OutPutBLAST.txt", sep = ""), " -outfmt \"7 qacc qlen sacc slen evalue pident length mismatch gapopen qstart qend sstart send bitscore qcovs\"",
            sep = "")
        system(aa)

        # Read the local Blast output file.
        outputBlast = readLines(paste(gsub(".fas", "", output2), "_OutPutBLAST.txt",
            sep = ""))

        # Extract the info from the Blastn output.
        outputBlast2 = outputBlast[-grep("# ", outputBlast, fixed = TRUE)]
        outputBlast3 = strsplit(outputBlast2, "\t", fixed = TRUE)
        outputBlast4 = do.call("rbind", outputBlast3)
        colnames(outputBlast4) = c("Query_SeqName", "SeqLen_Query", "Hit_SeqName",
            "Hit_SeqLen", "evalue", "pident", "length", "mismatch", "gapopen", "qstart",
            "qend", "sstart", "send", "bitscore", "qcovs")

        utils::write.table(outputBlast4, file = paste(output, "temp4.txt", sep = ""), sep = "\t",
            row.names = FALSE)
        outputBlast4 = read.delim(paste(output, "temp4.txt", sep = ""), sep = "\t",
            header = TRUE)
        file.remove(paste(output, "temp4.txt", sep = ""))

        # Filtering other potential 'secondary' outlier sequences using the bit score.
        # if the bit score of the hit is above the bitscore threshold defined as a
        # proportion of the highest bitscore possible when blasting the query against
        # itself then the sequence is considered as a 'secondary' outlier sequence.
        Outliertips = gsub(" no comment", "", Outliertips, fixed = T)
        bitscoreProp = vector()
        i = 1
        for (i in 1:length(Outliertips)) {
            DFtemp = outputBlast4[which(outputBlast4[, 1] == Outliertips[i]), ]
            bitscoreProp = c(bitscoreProp, DFtemp[, 14]/DFtemp[1, 14])
        }
        BitSco_outlier = outputBlast4[which(bitscoreProp > Bitsc.Th), ]

        # Export a table of the results and then the user can decide to select what is
        # the best criteria to refine the search for 'secondary' outlier sequences.
        outputBlast5 = cbind(BitSco_outlier, BitScore.Prop = bitscoreProp[which(bitscoreProp >
            Bitsc.Th)])

        utils::write.table(outputBlast5, file = output, sep = "\t", row.names = T)


        file.remove(output1, output2, paste(gsub(".fas", "", output1, fixed = T),
            "_BlastDB.nhr", sep = ""), paste(gsub(".fas", "", output1, fixed = T),
            "_BlastDB.nin", sep = ""), paste(gsub(".fas", "", output1, fixed = T),
            "_BlastDB.nsq", sep = ""), paste(gsub(".fas", "", output2), "_OutPutBLAST.txt",
            sep = ""))  # Remove the temporary files.

        return(outputBlast5)
    }
}
