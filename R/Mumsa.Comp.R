#' @title Compare congruence among multiple alignments: a wrapper for MUMSA software

#' @description This function compares multiple alignments using the MUMSA software (Lassmann &
#' Sonnhammer 2005) and reports the MOS (Multiple
#' Overlap Score) and the AOS (Average Overlap Score) statistics.

#' @details The MOS is provided for each alignment enabling us to check the support (more exactly the
#' consensuality) of those alignments. The alignment with the highest MOS is
#' supposed to be the most consensual alignment (Lassmann & Sonnhammer 2005).
#' Alignments with equal MOS are equally supported. The AOS gives an idea as to
#' whether the sequences are too divergent to be aligned. A AOS less than 0.5
#' indicates that the sequences are too divergent to be properly aligned
#' (potentially caused by saturation effects).
#'
#' @param input path to the folder storing the alignments. The names of the alignments
#' must follow the following template: name of the alignment program, the name of
# interest, and the name of the gene region just before the extension, and using
# '_' as separator.  The alignments have to be in fasta format (extension
# '.fas').
#'
#' @param output path and name of the file where you want to store the
#' results of MUMSA analysis. The output table reports the AOS for each gene
#' region and the MOS for each alignment and each gene region.
#' @param remove.empty.align If TRUE, the empty alignments are excluded from the computation.
#'
#' @details The function requires that the MUMSA software is installed and in the PATH.
#'
#' @return A tables with the MOS for each alignement and programs and the AOS for each
#' alignment exported as output txt document.
#'
#' @references Lassmann & Sonnhammer 2005, DOI: 10.1093/nar/gki1020
#'
#' @examples # Run the function without removing empty alignments
#' # Here the toy example consist of the alignment of four gene regions
#' # (co1, 12srrna, cytb, rag1) for a subset of 16 species.
#' # Each gene region was aligned using three programs mafftfftns1,
#' # mafftfftnsi, and pasta.
#' \dontrun{
#' input = system.file("extdata/multi.align/multi.aligned", package = "regPhylo")
#' Mumsa.Comp(input = input, output = "Mumsa_output",
#' remove.empty.align = FALSE)
#' # The function report an error message if some alignment are empty.
#' # Run the function without removing empty alignments
#' input = system.file("extdata/multi.align/multi.aligned", package = "regPhylo")
#' Mumsa.Comp(input = input, output = "Mumsa_output",
#' remove.empty.align = TRUE)
#'
#' # The output file is also present as an external data table in the regPhylo package
#' # and can be access by the following code:
#' # a = system.file("extdata", "Mumsa_output.txt", package = "regPhylo")
#' # a = read.delim(a, sep="\t", header = TRUE)
#'
#' # To clean the file created while running the example do the following:
#' file.remove("Mumsa_output.txt")
#'
#' }
#'
#' @export Mumsa.Comp
#'

Mumsa.Comp = function(input = NULL, output = NULL, remove.empty.align = NULL) {
    # Automatically detect the number of alignment programs used and the number of
    # gene regions introduced.
    a = list.files(input)

    # Check for empty alignment checking the disk space for each file.
    Filespace=vector()
    i=1
    for(i in 1:length(a)){
      Filespace=c(Filespace, file.info(paste(input, "/", a[i], sep=""))[[1]])
    }
    pbfile=a[which(Filespace==0)]
    # If remove.empty.align = TRUE, then the empty alignments are eliminated.
    if(remove.empty.align){
      if(length(which(Filespace==0))>0){
      a = a[-which(Filespace==0)]
      }
    } else {
      if(length(which(Filespace==0))>0){
    warning(paste(pbfile, collapse = "\n"))
    stop(paste("Some alignments might be empty; empty alignments must be remove from the folder before running the function", "\n",
                 "Look at in priority the files targeted by the warning message", "\n",  sep = ""))
      }
    }

    a1 = strsplit(a, "_")

    Nbprog = vector()
    Nbgene = vector()
    i = 1
    for (i in 1:length(a1)) {
        Nbprog = c(Nbprog, a1[[i]][1])
        Nbgene = c(Nbgene, a1[[i]][length(a1[[i]])])
    }
    Uniprog = unique(Nbprog)
    Unigene = unique(Nbgene)
    # For each gene prepare the MUMSA input file with the list of alignments for each
    # gene region.
    listgene = vector()
    i = 1
    for (i in 1:length(Unigene)) {
        listgene = c(listgene, paste(paste(input, "/", a[grep(Unigene[i], a)], sep = ""),
            collapse = " "))
    }

    # Run the MUMSA software
    i = 1
    for (i in 1:length(listgene)) {
        b = paste("mumsa -s ", listgene[i], " >> ", paste(output, ".txt", sep = ""),
            sep = "")
        system(b)
    }
    # Build a table with the AOS for each gene region and the MOS for each gene
    # region and each alignment program.
    bb = utils::read.table(paste(output, ".txt", sep = ""))
    colnames(bb) = c("AOS", paste(Uniprog, "_MOS", sep = ""))
    rownames(bb) = Unigene

    # Make a list of the most consensual alignments (higher MOS) for each gene
    # region, and also provide the alignments with the lowest score (sometimes the
    # most different alignments are the one making the most sens, for instance, it is
    # possible for one software to give an outstanding result compared to the
    # majority and consensus poor alignments!).
    tt = apply(bb[, -1], 1, max)
    ttmin = apply(bb[, -1], 1, min)

    Best.Align = vector()
    MinMos = vector()
    i = 1
    for (i in 1:length(tt)) {
        tt1 = match(tt[i], bb[i, -1])
        tt2 = match(ttmin[i], bb[i, -1])
        if (length(tt1) > 0)
            tt1 = tt1[1]
        Best.Align = c(Best.Align, gsub("_MOS", "", names(bb[, -1])[tt1], fixed = T))
        MinMos = c(MinMos, gsub("_MOS", "", names(bb[, -1])[tt2], fixed = T))
    }  # End for i


    bb = cbind(rownames(bb), bb, Best.Align, MinMos)
    colnames(bb) = c("Gene", "AOS", paste(Uniprog, "_MOS", sep = ""), "Most.Consensual.Align",
        "Most.Divergent.Align")
    utils::write.table(bb, file = paste(output, ".txt", sep = ""), row.names = F, sep = "\t")
    return(bb)
}  # End of the function.
