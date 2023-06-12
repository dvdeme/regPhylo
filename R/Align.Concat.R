#' @title Concatenate alignments from different gene regions into a supermatrix at the species level

#' @description This function concatenates the alignments from different gene regions into a
#' single supermatrix in nexus and fasta formats, at the species level.  The function also allows the
#' inclusion of species without DNA sequences, if necessary (for instance, to then use BEAST to resolve
#' polytomies).

#' @return This function returns: 1) the alignments from different gene regions
#' in nexus format '.nex', including taxa that do not have DNA sequence
#' information (nucleotides replace by '-') (all these files can be loaded separately
#' into BEAST); 2) a concatenation (a supermatrix) of all the sequences of
#' the different gene regions in nexus and fasta format; 3) a partition file in
#' txt format 'Partitions_Concat.txt' including the partitions of the different gene regions in the
#' concatenated file (this file is designed to be RAxML compatible, see RAxML manual
#' v8.2, https://sco.h-its.org/exelixis/resource/download/NewManual.pdf);
#' 4) a conversion table 'convtab.txt' betwen the gene region names used in the partition file and the
#' true name of the gene region.

#' @param input the path to the folder storing the alignments (alignments have to be in
#' fasta format with the '.fas' extension)
#' @param Sp.List.NoDNA an optional  vector of the species without DNA sequences that should be included
#' in the alignment, or the option can be NULL, in which case the function automatically creates a complete
#'  species list of all the species present in the different alignments.
#' @param outputConcat Name of the supermatrix (can include the path as well).
#' @param split the split between the information in the sequence name, 
#' by default the separator is "_".
#' 
#'@param chunk.names.to.Keep the number of chunk of information (splitted by split) within the 
#'sequence name to retain to build the super-matrix, this name must be common to the different 
#'alignment files. By default the first two chunks of information splitted by the split (e.g. "_").
#'
#'
#' @examples # Run the function to build a supermatrix
#' \dontrun{
#'
#' # To run the example, copy the input alignment files
#' # provided by the package to a temporary directory created in the
#' # current working directory.
#' src.dir = system.file("extdata/multi.align/ForConcat", package = "regPhylo")
#' dir.create("TempDir.ForConcat")
#' # Set up the path of the TempDir folder.
#' dest.dir = paste(getwd(), "/TempDir.ForConcat", sep="")
#' file.names <- dir(src.dir)
#' # Copy all the files stored in regPhylo/extdata/multi.align/ForConcat"
#' # into a temporary folder.
#' sapply(file.names, function(x) {
#' file.copy(from = paste(src.dir, x, sep = "/"),
#' to = paste(dest.dir, x, sep = "/"),
#' overwrite = FALSE) })
#'
#' # Run the function to build the supermatrix.
#' Align.Concat(input = "TempDir.ForConcat", Sp.List = NULL, outputConcat = NULL)
#'
#' #' # Import the supermatrix in R
#' require(ape)
#' Supermatrix = read.dna("TempDir.ForConcat/Concat.fas", format = "fasta")
#'
#' # Create another temporary file to build a supermatrix including species without DNA.
#' dir.create("TempDir.ForConcat2")
#' file.names <- dir("TempDir.ForConcat")
#' # select only the .fas alignment
#' file.names <- file.names[grep(".fas", file.names)]
#' # remove the Concat.fas alignment just created above.
#' file.names <- file.names[-grep("Concat.fas", file.names)]
#' sapply(file.names, function(x) {
#' file.copy(from = paste("TempDir.ForConcat", x, sep = "/"),
#' to = paste("TempDir.ForConcat2", x, sep = "/"),
#' overwrite = FALSE) })
#'
#'
#' # Run the function to build a supermatrix including two species without DNA.
#' Align.Concat(input = "TempDir.ForConcat2",
#' Sp.List = c("Titi_titi", "Toto_toto"),
#' outputConcat = "TempDir.ForConcat2/Concat_2spNoDNA")
#'
#' # Import the supermatrix into R.
#' Supermatrix2SpNoDNA = read.dna("TempDir.ForConcat2/Concat_2spNoDNA.fas",
#' format = "fasta")
#'
#'
#' # To remove the files created while running the example do the following:
#' unlink("TempDir.ForConcat", recursive = TRUE)
#' unlink("TempDir.ForConcat2", recursive = TRUE)
#'
#' }
#'
#' @export Align.Concat

Align.Concat = function(input = NULL, Sp.List.NoDNA = NULL, outputConcat = NULL, 
                        split = "_", chunk.names.to.Keep = 2) {
    listAli = paste(input, list.files(input), sep = "/")
    listAl = listAli[grep(".fas", listAli)]
    if(length(listAl) < length(listAli)) {
      stop(paste("Other files in non fasta format using '.fas' extension are present
                   in the input folder, only '.fas' alignment files must be present."))
    }


    AlignList = list(1:length(listAl))
    SeqName = vector()
    DFmat = matrix(NA, ncol = 3)[-1, ]
    i = 1
    for (i in 1:length(listAl)) {
        AlignList[[i]] = seqinr::read.fasta(listAl[i], as.string = TRUE)  # Store the alignments in a list.
        SeqT = vector()
        k = 1
        for (k in 1:length(AlignList[[i]])) SeqT = c(SeqT, gsub(" ", "", AlignList[[i]][[k]][1],
            fixed = TRUE))
        DFmat = rbind(DFmat, cbind(Alignment = rep(i, length(SeqT)), Seq.Name = labels(AlignList[[i]]),
            Sequences = SeqT))  # Extract the sequence names, the sequences and the numbers of the alignment.
    }  # End for i.
    Seq.Name.cor = gsub("_R_", "", DFmat[, 2], fixed = TRUE)  # Remove the '_R_' pattern when the sequences have been reversed complemented.
    Seq.Name.cor = gsub(".", "", Seq.Name.cor, fixed = TRUE)  # Remove the '.', in the sequence name.
    Seq.Name.cor = gsub("?", "", Seq.Name.cor, fixed = TRUE)  # Remove the '?', in the sequence name.
    Seq.Name.cor = gsub("-", "", Seq.Name.cor, fixed = TRUE)  # Remove the '-', in the sequence name.
    Seq.Name.cor = gsub("_sp_", "_sp", Seq.Name.cor, fixed = TRUE)  # Remove the '_' between sp and the letter or number defining a species not yet assigned a binomial species name.
    Seq.Name.cor = gsub("_nsp_", "_nsp", Seq.Name.cor, fixed = TRUE)  # Remove the '_' between nsp and the letter or number defining a new species not yet assigned a binomial species name.
    a = strsplit(Seq.Name.cor, split, fixed = T)  # Split the sequence name using '_' to extract the genus and species name.
    #a1 = lapply(a, function(x) x[1])
    #a2 = lapply(a, function(x) unlist(strsplit(x[2], "|", fixed = T))[1])
    #Sp.Name = unlist(lapply(seq(1, length(a1)), function(x) paste(a1[x], "_", a2[x],
    #   sep = "")))  # Extract the species name as the first two elements of each item in the list.
    
    Sp.Name = unlist(lapply(a, function(x) paste(x[1:chunk.names.to.Keep], collapse = "_")))  # Extract the species name as the first two elements of each item in the list.
    
    Sp.Name.list = unique(Sp.Name)  # The species list present in the different alignments

    # Include the option to also provide additional species without DNA.
    if (is.null(Sp.List.NoDNA)) {
        Sp.DF = as.data.frame(cbind(Sp.Name = Sp.Name.list, PresenceOverall = rep(1,
            length(Sp.Name.list))))
    } else {
        Sp.List.NoDNA = gsub(" ", "_", Sp.List.NoDNA, fixed = TRUE)  # Remove the spaces in the species names, for additional species without DNA.
        Sp.List.NoDNA = gsub(".", "", Sp.List.NoDNA, fixed = TRUE)  # Same syntax correction for the sequence name.
        Sp.List.NoDNA = gsub("?", "", Sp.List.NoDNA, fixed = TRUE)
        Sp.List.NoDNA = gsub("-", "", Sp.List.NoDNA, fixed = TRUE)
        Sp.List.NoDNA = gsub("_sp_", "_sp", Sp.List.NoDNA, fixed = TRUE)
        Sp.List.NoDNA = gsub("_nsp_", "_nsp", Sp.List.NoDNA, fixed = TRUE)
        Sp.DF = as.data.frame(cbind(Sp.Name = c(Sp.Name.list, Sp.List.NoDNA), PresenceOverall = rep(1,
            length(c(Sp.Name.list, Sp.List.NoDNA)))))
    }

    # Large table storing all the sequences for all the alignments of interest.
    DFmat2 = as.data.frame(cbind(Alignment = DFmat[, 1], Sp.Name = Sp.Name, Seq.Name.cor,
        Sequences = DFmat[, 3]))  # include the species name for each sequence.

    # Prepare the nexus extension of the file names.
    listAl.nexus = gsub(".fas", ".nex", listAl, fixed = TRUE)


    # Prepare the supermatrix.
    SuperMat = matrix(NA, ncol = 1)[-1, ]

    # Create an alignment with all the species when considering all the alignments
    # together.
    i = 1
    for (i in 1:length(unique(DFmat2[, 1]))) {
        DFtemp = DFmat2[which(DFmat2[, 1] == i), ]  # Select the alignment
        AlignTemp = merge(Sp.DF, DFtemp, by.x = 1, by.y = 2, all.x = TRUE)
        AlignTemp = as.matrix(AlignTemp)  # Convert into a matrix
        # Test is mutliple sequences per species are present in the alignment
        if(dim(AlignTemp)[1]>dim(Sp.DF)[1]) {
          stop(paste("The alignment ", listAl[i], " certainly contains multiple sequences for the same species:
                     to be concatenated at the species level, only one sequence per species per alignment must be provided!"))
        }
        AlignTemp[which(is.na(AlignTemp[, 4]) == "TRUE"), 5] <- paste(rep("-", nchar(as.character(DFtemp[1,
            4]))), collapse = "")  # Replace the empty sequence by a long string of '----'
        AlignTemp[which(is.na(AlignTemp[, 4]) == "TRUE"), 4] = as.character(AlignTemp[which(is.na(AlignTemp[,
            4]) == "TRUE"), 1])  # Replace the sequence name of the empty sequence by the species name.

        # Feed the supermatrix.
        SuperMat = cbind(SuperMat, AlignTemp[, 5])

        # Create a nexus file including those empty sequences.
        NBChar = nchar(as.character(AlignTemp[1, 5]))

        cat(file = listAl.nexus[i], "#NEXUS", "\n", "\n", "BEGIN DATA;", "\n", "\t",
            paste("DIMENSIONS NTAX=", dim(AlignTemp)[1], sep = ""), paste(" NCHAR=",
                NBChar, ";", sep = ""), sep = "", append = TRUE)
        cat(file = listAl.nexus[i], "\n", "\t", "FORMAT DATATYPE=DNA GAP=-;", "\n",
            "MATRIX", "\n", sep = "", append = T)
        utils::write.table(AlignTemp[, c(4, 5)], file = listAl.nexus[i], sep = "\t", append = TRUE,
            col.names = FALSE, row.names = FALSE, quote = FALSE)
        cat(file = listAl.nexus[i], "\t", ";", "\n", "END;", sep = "", append = TRUE)
    }  ## End for i

    # Create a large supermatrix.
    concat = vector()
    i = 1
    for (i in 1:dim(SuperMat)[1]) {
        concat = c(concat, paste(SuperMat[i, ], collapse = ""))  # Concatenate the different alignments in one long sequence.
    }

    SuperMatDF = cbind(sort(as.character(Sp.DF[, 1])), concat)  # Add the species name for all the sequences, and order the sequences alphabetically to match the species names in the supermatrix.

    # Create the nexus file for the supermatrix.
    NBChar = nchar(as.character(concat[1]))
    b = unlist(strsplit(listAl.nexus[1], "/", fixed = TRUE))
    # ConcatName = paste(paste(b[-length(b)], collapse = "/"), "Concat.nex", sep = "/")

    # If the option outputConcat is null, the name of the concat will be "Concat"
    if(is.null(outputConcat)){
      ConcatName = paste(paste(b[-length(b)], collapse = "/"), "Concat.nex", sep = "/")
    } else {
      ConcatName = paste(outputConcat, ".nex", sep="")
    }
    cat(file = ConcatName, "#NEXUS", "\n", "\n", "BEGIN DATA;", "\n", "\t", paste("DIMENSIONS NTAX=",
        dim(SuperMatDF)[1], sep = ""), paste(" NCHAR=", NBChar, ";", sep = ""), sep = "",
        append = TRUE)
    cat(file = ConcatName, "\n", "\t", "FORMAT DATATYPE=DNA GAP=-;", "\n", "MATRIX",
        "\n", sep = "", append = TRUE)
    utils::write.table(SuperMatDF, file = ConcatName, sep = "\t", append = TRUE, col.names = FALSE,
        row.names = FALSE, quote = FALSE)
    cat(file = ConcatName, "\t", ";", "\n", "END;", sep = "", append = TRUE)

    # Create a fasta file for the supermatrix.
    Seq.name.seq = paste(paste(">", SuperMatDF[, 1], sep = ""), SuperMatDF[, 2],
        sep = "+++")
    FastaAlign = unlist(strsplit(Seq.name.seq, "+++", fixed = TRUE))
    if(is.null(outputConcat)){
      ConcatName = paste(paste(b[-length(b)], collapse = "/"), "Concat.fas", sep = "/")
    } else {
      ConcatName = paste(outputConcat, ".fas", sep="")
    }
    write(FastaAlign, file = ConcatName)  # write the alignemnt


    GeneName=gsub(".fas", "", unlist(lapply(strsplit(listAl, "_"), function(x) x[length(x)])), fixed=TRUE)

    # Create a partition file for RAxML identifying the beginning and the end of each
    # gene region.
    genelength = vector()
    i = 1
    for (i in 1:dim(SuperMat)[2]) {
        genelength = c(genelength, nchar(as.character(SuperMat[1, i])))
    }
    LimSup = cumsum(genelength)  # Upper limit of the gene region.
    liminf = c(1, LimSup[-length(LimSup)] + 1)  # Lower limit of the gene region.
    # Print the partition file (compatible with RAxML), and a conversion table providing information between the code of the gene region used by
    # PartitionFinder2 and the true name of the gene region.
    convtab=matrix(NA, ncol=2)[-1,]
    i = 1
    for (i in 1:dim(SuperMat)[2]) {
        cat(file = paste(paste(b[-length(b)], collapse = "/"), "Partitions_Concat.txt",
            sep = "/"), "DNA, gene", i, " = ", liminf[i], "-", LimSup[i], "\n", sep = "",
            append = TRUE)
      convtab = rbind(convtab, c(paste("gene", i, sep = ""), GeneName[i]))
    }
    colnames(convtab) = c("Name.PartitionFinder2", "Common.Gene.Name")
    utils::write.table(convtab, file=paste(input, "/convtab.txt", sep=""), sep="\t", row.names=FALSE)
return(convtab)
}  # End of the function.
