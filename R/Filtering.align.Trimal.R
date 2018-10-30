#' @title Remove poorly aligned nucleotides position using trimAl software.

#' @description This function uses the software trimAl (Capella-Gutierrez et al. 2009, DOI:
#' 10.1093/bioinformatics/btp348) and in particular two approaches to filter out
#' the data using the option -gappyout or -automated1 see (Capella-Gutierrez et
#' al.  2009).

#' @param input path to the folder storing the alignments in fasta format (with the
# extension '.fas').

#' @param output path to the folder storing the trimmed alignments.
#' The output folder is created automatically. The suffix 'trimAuto_' or
#' 'trimGapy_' are added to the alignment file name according to the chosen
#' heuristic.

#' @return An ouptut folder is created with the trimmed alignments. In the R environment
#' the function returns a table with the length of the different alignments for
#' each gene region and alignment programs.

#' @details The option "-gappyout" is very conservative method, it keeps the maximum nucleotide information and remove
#' the most gappy positions). This option is one of the best according Tan et al.  2015, DOI: 10.1093/sysbio/syv033).
#' The option "-automated1" automatically adjusts between -strict -strictplus and -gappyout heuristics which all provide a more
#' stringent selection potentially than the -gappyout option (see Capella-Gutierrez et al. 2009, DOI:
#' 10.1093/bioinformatics/btp348, for more information)

#' @details The function requires trimAl to be installed and set up in the PATH.

#' @examples # Run the function
#' \dontrun{
#' # To run the example it might be better to copy the input alignment files
#' # provided by the package to a temporary directory created into the
#' # current working directory.
#' src.dir = system.file("extdata/multi.align/multi.aligned", package = "regPhylo")
#' dir.create("TempDir.ToTrim")
#' # Set up the path of the TempDir folder.
#' dest.dir = paste(getwd(), "/TempDir.ToTrim", sep="")
#' file.names <- dir(src.dir)
#' # Copy all the files stored in regPhylo/extdata/multi.align"
#' # into a temporary folder.
#' sapply(file.names, function(x) {
#' file.copy(from = paste(src.dir, x, sep = "/"),
#' to = paste(dest.dir, x, sep = "/"),
#' overwrite = FALSE) })
#'
#' input = "TempDir.ToTrim"
#' output = "Trimmed"
#' Filtering.align.Trimal(input = input, output = "Trimmed")
#'
#' # To clean the file created while running the example do the following:
#' # Remove the temporary folder
#' unlink("TempDir.ToTrim", recursive = TRUE)
#' # Remove the folder with Gblocks outputs
#' unlink("Trimmed", recursive = TRUE)
#' }
#'
#' @export Filtering.align.Trimal
#'
#'
#' @references Capella-Gutierrez et al. 2009, DOI:
#' 10.1093/bioinformatics/btp348

Filtering.align.Trimal = function(input = NULL, output = NULL) {
    b = list.files(input)
    bb = b[grep(".fas", b, fixed = T)]

    # Very conservative method (keep the maximum nucleotide information and remove
    # the most gappy positions) -gappyout option is one of the best according Tan et
    # al.  2015, DOI: 10.1093/sysbio/syv033).
    dir.create(output)
    trimal.gappy = function(x) {
        a = paste("trimal -in ", input, "/", x, " -out ", output, "/", "trimGapy_",
            x, " -gappyout", sep = "")
        system(a)
    }
    tryCatch(lapply(bb, trimal.gappy), error = function(e) e)

    # Automatically adjust between -strict -strictplus and -gappyout heuristics (more
    # stringent selection potentially than the -gappyout option).
    trimal.automated = function(x) {
        a = paste("trimal -in ", input, "/", x, " -out ", output, "/", "trimAuto_",
            x, " -automated1", sep = "")
        system(a)
    }
    tryCatch(lapply(bb, trimal.automated), error = function(e) e)

    a0 = list.files(output)
    a1 = strsplit(a0, "_")

    Nbprog = vector()
    Nbgene = vector()
    Nbtrim = vector()
    i = 1
    for (i in 1:length(a1)) {
        Nbtrim = c(Nbtrim, a1[[i]][1])
        Nbprog = c(Nbprog, a1[[i]][2])
        Nbgene = c(Nbgene, a1[[i]][length(a1[[i]])])
    }
    Uniprog = unique(Nbprog)
    Unigene = unique(Nbgene)
    Unitrim = unique(Nbtrim)

    # Estimate the length of each alignment.
    resDF = matrix(NA, ncol = 5)[-1, ]
    k = 1
    for (k in 1:length(Unitrim)) {
        Atemp = a0[grep(Unitrim[k], a0)]
        j = 1
        for (j in 1:length(Unigene)) {
            Atemp1 = Atemp[grep(Unigene[j], Atemp)]
            i = 1
            for (i in 1:length(Uniprog)) {
                Align = seqinr::read.fasta(paste(output, "/", Atemp1[grep(Uniprog[i],
                  Atemp1)], sep = ""), as.string = T)
                resDF = rbind(resDF, c(Atemp1[grep(Uniprog[i], Atemp1)], Unitrim[k],
                  Unigene[j], Uniprog[i], nchar(Align[1])))
            }  ## End for i
        }  # End for j
    }  # End for k
    colnames(resDF) = c("Align.Name", "Trim.Method", "Gene.Name", "Program", "SeqLength")
    # Order the table according to the gene name.
    resDFord = resDF[order(resDF[, 3]), ]
    resDFord = as.data.frame(resDFord)
    return(resDFord)
}  # End of the function.

