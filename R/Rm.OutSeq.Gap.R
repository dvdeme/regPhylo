#' @title Remove outlier sequences from an alignment and delete the gap sites only

#' @description This function removes the outlier sequences from the alignment and deletes the
#' gap sites only.

#' @param input vector with the name of the sequences to remove.
#' @param SeqInput the original alignment as an "alignment" (seqinr) or a "DNAbin" (ape) object or
#' the name (can include the path) of the original sequence alignment in fasta format to import into R.
#' @param AligoutputName name of the output (output will be in fasta format with .fas extension).

#' @examples # Load the original alignment including the outlier sequences and the
#' # table of potential outlier sequences detected by Detect.Outlier.Seq function.
#' # Example_16S_outlier is an .Rdata object providing a list of: 1) the original alignment; and
#' # 2) the table of potential outlier sequences.
#' \dontrun{
#' data(Example_16S_outlier)
#' Example_16S_outlier_align = Example_16S_outlier[[1]]
#' Table.Outlier.Seq.16S = Example_16S_outlier[[2]]
#'
#' # Run the function. The name of all outlier sequences are stored
#' # in the third columns (i.e. "Hit_SeqName") of the table,
#' # because the same sequence might be present multiple times,
#' # we use the function unique to remove the duplicated sequences
#' # (i.e. multiple hits by Blast).
#' # In the following example we also decided to remove another sequences
#' # (e.g. "_R_Polyprion_americanus|NA|AM158291|NA") which was not detected by
#' # the function Detect.Outlier.Seq, but that we would like to remove.
#'
#' New16SAlignment = Rm.OutSeq.Gap(input = c(as.character(unique(Table.Outlier.Seq.16S[,3])),
#' "_R_Polyprion_americanus|NA|AM158291|NA"), SeqInput = Example_16S_outlier_align,
#' AligoutputName = "16S_example_RMoutliers")
#'
#' class(New16SAlignment) # it is a "DNAbin" object
#' New16SAlignment # 280 DNA sequences.
#'
#' # To remove the file created while running the example do the following:
#' file.remove("16S_example_RMoutliers.fas")
#'
#' # The output files (and another example of outlier table "Example_16S_outliers_1.txt")
#' # are also present as an external data in the regPhylo package
#' # and can be accessed by the following code:
#' # a = system.file("extdata/ExampleOutliers", package = "regPhylo")
#' # list.files(a)
#'
#'
#' }

#' @return Two alignments are exported: one "DNAbin" alignment in the R environment;
#'  and one alignment exported in fasta format with the .fas extension to
#'  the location defined in the AligoutputName option.

#' @export Rm.OutSeq.Gap

Rm.OutSeq.Gap = function(input = NULL, SeqInput = NULL, AligoutputName = NULL) {

  # if the class of the SeqInput object is "alignment" object recognised by ape and seqinr r packages.
  if(class(SeqInput)=="alignment"){
    cytb2S=SeqInput
  }
  # if the path to the alignment file in fasta format is provided as SeqInput
  if(class(SeqInput)=="character"){
    cytb2S = seqinr::read.alignment(SeqInput, format = "fasta")
  }
  # if the class of the inputal object is a "DNAbin" object.
  if(class(SeqInput)=="DNAbin"){
    cytb2S=ape::as.alignment(SeqInput) # convert to an "alignment" object.
  }

    cytb2S$nam = gsub(" no comment", "", cytb2S$nam, fixed = T)  # remove the ' no comment' string after the sequence name if necessary.
    # remove the outlier sequences
    cytb2S$seq = cytb2S$seq[-match(as.character(input), cytb2S$nam)]
    cytb2S$nam = cytb2S$nam[-match(as.character(input), cytb2S$nam)]
    cytb2S$nb = length(cytb2S$nam)
    cytb2Smat = seqinr::as.matrix.alignment(cytb2S)  # convert in a matrix

    # remove the deletion gap site only.
    nbindels = apply(as.matrix(seq(1, dim(cytb2Smat)[2], 1)), 1, function(x) {
        sum(stringr::str_count(cytb2Smat[, x], "-"))
    })
    NewAlign.Noindels = cytb2Smat[, -which(nbindels == dim(cytb2Smat)[1])]

    # convert the matrix in string for all the sequences
    NewAlign.Noindels2 = vector()
    i = 1
    for (i in 1:dim(cytb2Smat)[1]) {
        NewAlign.Noindels2 = c(NewAlign.Noindels2, seqinr::c2s(NewAlign.Noindels[i, ]))
    }
    cytb2S$seq = NewAlign.Noindels2
    cytb2S$nam = unlist(lapply(strsplit(cytb2S$nam, " "), function(x) {
        x[1]
    }))  # remove the comments in the sequence name that might be automatically added by Seaview (Gouy et al. 2010, DOI: 10.1093/molbev/msp259)
    cytb2Sa = ape::as.DNAbin(cytb2S)  # convert to DNAbin object and export as fasta format.
    ape::write.dna(cytb2Sa, file = paste(AligoutputName, ".fas", sep = ""), format = "fasta",
        nbcol = -1, colsep = "")
    return(cytb2Sa)
}
