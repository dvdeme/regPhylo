#' @title Help to detect outlier sequences in an alignment

#' @description This function helps to detect some potential outlier sequences
#' (i.e. mis-aligned sequences, which might be caused by different problems such
#' as gene or species annotation problems, presence of paralogous sequences, etc)
#' that should be removed from the pool of sequences (per species and gene regions)
#' before selecting the best sequence.  This function is only provided as a tool
#' to HELP detect potential outlier sequences, great care must be taken to
#' ensure sequences are not erroneously omitted as potential outliers,
#' or remain included undetected outlier sequences. We recommend checking
#' alignments by eye to ensure your confidence in the retained sequences.


#' @details This function uses the following general strategy to detect outlier
#' sequences.  First, a distance matrix is computed among all sequences (see
#' different strategies available below), second a BIONJs tree (i.e. an improved
#' Neighbor joining method that allows missing values, Criscuolo & Gascuel 2008)
#' is built based on the distance matrix, and the tip to root distance is computed.
#' All the sequences above a certain distance threshold (we
#' use tip to root distances at least 50-60\%, see Chen et al. 2015 for similar thresholds)
#' are extracted as potential 'primary' outlier sequences.  In
#' addition, we can blast these primary outlier sequences back against the original alignment
#' converted into a local blast database (the alignement is automatically converted
#' by the function into a local database) to look for secondary outlier sequences that might
#' be similar to the primary sequences but just below the tip to root distance threshold.
#' All the sequences above a bitscore threshold (we use 80\%) are considered as potential 'secondary'
#' outliers. A bitscore threshold of 100\% corresponds to blasting the query sequence against itself.
#' We used the bitscore because it is not influenced by the size of the local database contrary to the E.value
#' (\url{https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html}).

#' @return The function return a table for all the retained sequences with the following headers:
#' \itemize{
#' \item 'Query_SeqName' Name of the query sequence,
#' \item 'SeqLen_Query' length of the query sequence,
#' \item 'Hit_SeqName' name of the hit sequence,
#' \item 'Hit_SeqLen' length of the hit sequence,
#' \item 'evalue' E-value,
#' \item 'pident' means percentage of identical matches,
#' \item 'length' means alignment length,
#' \item 'mismatch' means number of mismatches,
#' \item 'gapopen' means number of gap openings,
#' \item 'qstart' means start of alignment in query,
#' \item 'qend' means end of alignment in query,
#' \item 'sstart' means start of alignment in subject,
#' \item 'send' means end of alignment in subject,
#' \item 'bitscore' means bitscore,
#' \item 'qcovs' means query coverage per subject,
#' }
#' (See \url{https://www.ncbi.nlm.nih.gov/books/NBK279675/} for more details and
#' options)

#' @details Three options are available to compute the distance matrix among sequences.
#' The first option 'MisAli' uses the number of indelblocks to compute the
#' distance matrix in order to preferentially detect misaligned sequences
#' generating gaps in the alignment.  The second option 'DivSeq' uses the TN93
#' substitution model to compute the distance matrix in order to detect very
#' divergent sequences. The third option 'Comb' combines the two
#' distance matrices described above to detect both misaligned and very divergent sequence.
#' To do this, the two distance matrices are first re-scaled as percentages of the their highest
#' distances respectively and then summed. For a additional details see Eme et al. (2019) Appendix 4.
#' @details This function requires BLAST+ installed and in the PATH in order to detect 'secondary outliers'.
#' To download and install BLAST+ software locally go
#' to \url{https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download},
#' specific installation instructions for the different operating systems (OS) can be find
#' at \url{https://www.ncbi.nlm.nih.gov/books/NBK279671/}.

#' @param inputal an object of class "alignment" (seqinr R package),
#' or "DNAbin" (ape R package) or the name (including the path if necessary) of the input fasta file.
#' @param Strat.DistMat can be 'MisAli', 'DivSeq' or 'Comb', see details.
#' @param Dist.Th distance threshold to detect the primary set of
#' outlier sequences (between [0-1]).
#' @param output name of the output table exported in the working directory.
#' @param Second.Outlier if 'TRUE', blasts the primary outlier sequences
#' against the local database of retrieved and aligned sequences
#' (i.e. the alignment provided to 'inputal') in order to detect secondary outliers.
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
#' dim(S16_MisAlign0.6_1) ### 36 sequences detected as primary outlier sequences.
#' head(S16_MisAlign0.6_1)
#'
#' # The output table is also present as an external data table provided in the regPhylo r package
#' and can be accessed by the following code:
#' # a = system.file("extdata/ExampleOutliers/Example_16S_outliers_1.txt", package = "regPhylo")
#' # Example_16S_outliers_1 = read.delim(a, sep="\t", header = TRUE)
#'
#'
#' # Running the function with the 'Comb' option, a distance threshold of 0.6,
#' # and allowing a search for secondary outlier sequences using local blast database.
#' S16_MisAlign0.6_2 = Detect.Outlier.Seq(inputal = Example_16S_outlier_align,
#' Strat.DistMat = "Comb", Dist.Th = 0.6, output = "Example_16S_outliers_2.txt",
#' Second.Outlier = "Yes", Bitsc.Th = 0.8)
#' length(unique(S16_MisAlign0.6_2[,3])) ### 38 primary and secondary outlier sequences detected.
#' head(S16_MisAlign0.6_2)
#'
#' # The output table is also present as an external data table provided in the regPhylo r package
#' and can be accessed by the following code:
#' # a = system.file("extdata/ExampleOutliers/Example_16S_outliers_2.txt", package = "regPhylo")
#' # Example_16S_outliers_2 = read.delim(a, sep="\t", header = TRUE)
#'
#' # To remove the file created while running the example do the following:
#' file.remove(c( "Example_16S_outliers_1.txt", "Example_16S_outliers_2.txt"))
#'
#' }
#'
#' @references Criscuolo & Gascuel 2008, DOI: 10.1186/1471-21.5-9-166
#' @references Chen et al. 2015, DOI: 10.1093/sysbio/syv059
#' @references Eme et al. (2019). An integrated pathway for building regional
#' phylogenies for ecological studies. \emph{Global Ecology and Biogeography}, Accepted.

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

        # Test if BLAST+ tools have been installed and setup in the PATH.
        # Test with makeblastdb tools to build the local blast database.
        zz = try(system("makeblastdb -version", intern = TRUE))
        if(class(zz) == "try-error") {
          stop("BLAST+ software cannot be found in the PATH ('makeblastdb' cannot be found). For secondary outlier detection BLAST+ software have to be installed locally and in the PATH.", "\n",
                "To download and install BLAST+ software locally go to https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download" )
        }

        # Test with blastn tools to blast the nucleotide sequences.
        zz = try(system("blastn -version", intern = TRUE))
        if(class(zz) == "try-error") {
          stop("BLAST+ software cannot be found in the PATH ('blastn' cannot be found). For secondary outlier detection BLAST+ software have to be installed locally and in the PATH.", "\n",
               "To download and install BLAST+ software locally go to https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download" )
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
