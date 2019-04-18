#' @title Remove poorly aligned nucleotide positions using trimAl software

#' @description This function uses the software trimAl (Capella-Gutierrez et al. 2009)
#' and in particular two approaches to filter out
#' poorly aligned nucleotide positions using the option -gappyout or -automated1 (see Capella-Gutierrez et
#' al. 2009).

#' @param input path to the folder storing the alignments in fasta format (with the
#' extension '.fas').

#' @param output path to the folder storing the trimmed alignments.
#' The output folder is created automatically. The suffix 'trimAuto_' or
#' 'trimGapy_' is added to the alignment file name according to the chosen
#' heuristic.
#'
#' @param TrimAl.path for the Windows platform, a character string which provides the path
#' to the trimAl executable
#' (e.g. "C:/Users/deme/Documents/Programs/TrimAl/trimal.v1.2rev59/trimAl/bin/trimal.exe").
#' For Linux the trimAl software must be in the $PATH.

#' @return An output folder is created with the trimmed alignments. In the R environment
#' the function returns a table with the length of the different alignments for
#' each gene region and alignment programs.

#' @details The option "-gappyout" is a very conservative method, it keeps the maximum nucleotide information and remove
#' the most gappy positions). This option is one of the best according Tan et al.  2015).
#' The option "-automated1" automatically adjusts between the '-strict', '-strictplus' and '-gappyout' heuristics which all provide a more
#' stringent selection potentially than the '-gappyout' option (for more information see Capella-Gutierrez et al. 2009).

#' @details The function requires trimAl to be installed and set up in the PATH for Linux platform.
#' Online documentation, including download page, are available at \url{http://trimal.cgenomics.org/downloads}.

#' @examples # Run the function
#' \dontrun{
#' # To run the example copy the input alignment files
#' # provided by the package to a temporary directory created into the
#' # current working directory.
#' src.dir = system.file("extdata/multi.align/multi.aligned", package = "regPhylo")
#' dir.create("TempDir.ToTrim")
#' # Set up the path to the TempDir folder.
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
#'
#' # Run the function. For windows users, remember to
#' # additionally specify the TrimAl.path.
#' Filtering.align.Trimal(input = input, output = "Trimmed")
#'
#' # To remove the file created while running the example do the following:
#' # Remove the temporary folder
#' unlink("TempDir.ToTrim", recursive = TRUE)
#' # Remove the folder with TrimAl outputs
#' unlink("Trimmed", recursive = TRUE)
#' }
#'
#' @export Filtering.align.Trimal
#'
#' @references Capella-Gutierrez et al. 2009, DOI:
#' 10.1093/bioinformatics/btp348
#' Tan et al.  2015, DOI: 10.1093/sysbio/syv033

Filtering.align.Trimal = function(input = NULL, output = NULL, TrimAl.path = NULL) {
    b = list.files(input)
    bb = b[grep(".fas", b, fixed = T)]

    # Very conservative method (keep the maximum nucleotide information and remove
    # the most gappy positions) -gappyout option is one of the best according Tan et
    # al.  2015, DOI: 10.1093/sysbio/syv033).
    dir.create(output)


    trimal = "trimal"
    os <- .Platform$OS
    if(os == "windows"){
      if(missing(TrimAl.path)){
        stop("The path to the TrimAl executable must be provided in TrimAl.path")
      }
      trimal = TrimAl.path
    }

    trimal.gappy = function(x) {
        a = paste(trimal, " -in ", input, "/", x, " -out ", output, "/", "trimGapy_",
            x, " -gappyout", sep = "")
        system(a)
    }
    tryCatch(lapply(bb, trimal.gappy), error = function(e) e)

    # Automatically adjust between -strict -strictplus and -gappyout heuristics (more
    # stringent selection potentially than the -gappyout option).
    trimal.automated = function(x) {
        a = paste(trimal, " -in ", input, "/", x, " -out ", output, "/", "trimAuto_",
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

    # Test if multiple alignment belonging to the same gene occur in the input file.
    #if(length(which(table(Nbgene)>2))>0){
    #  warning("There are multiple alignments of the same gene region in the input folder,
    # the function coudn't export the table with the sequences length of the alignment,
    #          per program, gene region, and trimming algorithm")
    #}

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
                  Atemp1)], sep = ""), as.string = TRUE)
                resDF = rbind(resDF, c(Atemp1[grep(Uniprog[i], Atemp1)], Unitrim[k],
                  Unigene[j], Uniprog[i], nchar(gsub("[ ]?", "", Align[1], perl = TRUE))))
            }  ## End for i
        }  # End for j
    }  # End for k
    colnames(resDF) = c("Align.Name", "Trim.Method", "Gene.Name", "Program", "SeqLength")
    # Order the table according to the gene name.
    resDFord = resDF[order(resDF[, 3]), ]
    resDFord = as.data.frame(resDFord)
    return(resDFord)
}  # End of the function.

