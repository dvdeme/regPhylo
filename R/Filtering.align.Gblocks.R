#' @title Remove poorly aligned nucleotide positions using Gblocks software

#' @description This function uses the software Gblocks V 0.91b (Castresana 2000)
#' to filter out poorly aligned positions and divergent regions in an alignment.
#' It allows the user to specify if the sequences are proteins, DNA-non coding, or DNA coding regions.
#' The function offers the possibility to use the default parameter (more stringent selection)
#' or the less stringent selection approach (parameters equivalent to that on the Gblocks server,
#' \url{http://molevol.cmima.csic.es/castresana/Gblocks_server.html}).

#' @param input name of the folder storing the alignments in fasta format (with the
#' extension '.fas').
#' @param target.file the name of the file(s) in the input folder that will be used.
#' By default the function uses all the files with a .fas extension in the input folder.
#' @param LessStringent if 'TRUE' opts for the less
#' stringent selection (equivalent to ticking the following three options in Gblocks server,
#' 'Allow smaller final blocks', 'Allow gap positions within the final blocks',
#' 'Allow less strict flanking positions', the output has the prefix Gblocksls),
#' otherwise the default parameters of Gblocks will be used providing a more stringent selection
#' (output with the prefix Gblocksms).
#' @param Type Type of sequences can be Protein (p), DNA (d), or Codons (c).
#' @param output path to the folder storing the trimmed alignments.
#' The output folder is created automatically is necessary.
#' @param remove.empty.align If TRUE, the empty alignments are excluded from the computation.
#' @param Gblocks.path for the Windows plateform, a character string which provides the path
#' to the Gblocks executable but without the name of the executable
#' (e.g. "C:/Users/deme/Documents/Programs/Gblocks/Gblocks_Windows_0.91b/Gblocks_0.91b").
#' For Linux the Gblocks software must be in the $PATH.
#'
#' @details This function requires Gblocks to be installed and set up in the PATH.
#' Online documentation, including instruction for installing the program is available at
#' \url{http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html},
#' and the download page is available at
#' \url{http://molevol.cmima.csic.es/castresana/Gblocks.html}.

#' @return An ouptut folder is created with the trimmed alignments. In the R environment
#' the function returns a table with the length of the different alignments for
#' each gene region and the alignment programs.
#'
#' @examples # Run the function for DNA, using the less stringent
#' # selection heuristic, and remove potentially empty alignments.
#' \dontrun{
#' # To run the example copy the input alignment files
#' # provided by the package to a temporary directory created in the
#' # current working directory.
#' src.dir = system.file("extdata/multi.align/multi.aligned", package = "regPhylo")
#' dir.create("TempDir")
#' # Set up the path to the TempDir folder.
#' dest.dir = paste(getwd(), "/TempDir", sep="")
#' file.names <- dir(src.dir)
#' # Copy all the files stored in regPhylo/extdata/multi.align/multi.aligned"
#' # into a temporary folder.
#' sapply(file.names, function(x) {
#' file.copy(from = paste(src.dir, x, sep = "/"),
#' to = paste(dest.dir, x, sep = "/"),
#' overwrite = FALSE) })
#'
#' # Run the function from the TempDir folder and store the outputs from
#' # Gblocks in the "Trimmed-Gblocks" folder. For Windows users, remember
#' # to additionally specify the Gblocks.path.
#' Filtering.align.Gblocks(input = "TempDir", LessStringent = TRUE,
#' output = "TrimmedGblocks", Type = "d", remove.empty.align = TRUE)
#'
#' # Run the function from the TempDir folder but selecting only one alignment
#' # "Mafftfftns1_Alig_co1.fas", export the Gblocks alignement using the
#' # stringent selection in a folder called "TrimmedGblocks_1file".
#' # For Windows users, remember to additionally specify the Gblocks.path.
#' Filtering.align.Gblocks(input = "TempDir", target.file = "Mafftfftns1_Alig_co1.fas",
#' LessStringent = FALSE, output = "TrimmedGblocks_1file",
#' Type = "d", remove.empty.align = TRUE)
#'
#'
#' # To remove the file created while running the example do the following:
#' # Remove the temporary folder
#' unlink("TempDir", recursive = TRUE)
#' # Remove the folder with the Gblocks outputs
#' unlink("TrimmedGblocks", recursive = TRUE)
#' unlink("TrimmedGblocks_1file", recursive = TRUE)
#' }
#'
#' @export Filtering.align.Gblocks
#'
#' @references Castresana 2000,
#' DOI: 10.1093/oxfordjournals.molbev.a026334

Filtering.align.Gblocks = function(input = NULL, target.file = NULL, LessStringent = NULL,
                                   Type = NULL, output = NULL,
                                   remove.empty.align = NULL, Gblocks.path = NULL) {


    if(is.null(target.file)){
      b = list.files(input)
      bb = b[grep(".fas", b, fixed = TRUE)] # if target.file is NULL then
      # all the file in the input folder with .fas extension are used.
    } else {
      bb = target.file
    }


    # Detect the empty alignments
    # Check for empty alignment checking the disk space for each file.
    Filespace = vector()
    i = 1
    for(i in 1:length(bb)){
      Filespace = c(Filespace, file.info(paste(input, "/", bb[i], sep = ""))[[1]])
    }
    pbfile = bb[which(Filespace == 0)]
    # If remove.empty.align = TRUE, then the empty alignments are eliminated.
    if(remove.empty.align){
      if(length(which(Filespace == 0))>0){
        bb = bb[-which(Filespace == 0)]
      }
    } else {
      # Check if some file are empty, if not the function carries on.
      if(length(pbfile) > 0){
      warning(paste(pbfile, collapse = "\n"))
      stop(paste("Some alignments might be empty; empty alignments must be remove from the folder before running the function", "\n",
                 "Look at in priority the files targeted by the warning message", "\n",  sep = ""))
      }
    }


    # Create the output folder if an output name is provided,
    # if it is null then the output will be exported in the same folder as the input folder.
    if(is.null(output)){
      output = input
    } else {
      dir.create(output, showWarnings = FALSE)
    }

    # Check the os plateform
    os <- .Platform$OS

    if(os == "windows"){
      if(missing(Gblocks.path)){
        stop("The path to the Gblocks executable must be provided in Gblocks.path")
      } else {
        src.dir = paste(getwd(), "/", input, sep = "")
        file.names = list.files(src.dir)

        oldwd = getwd()
        i = 1
        for(i in 1:length(file.names)){
          # copy the file into the directory wher Gblocks executable is found
          file.copy(from = paste(src.dir, "/", file.names[i], sep = ""),
                    to = paste(Gblocks.path, "/", file.names[i], sep = ""))
          # Change
          setwd(Gblocks.path)
          Align = seqinr::read.fasta(file.names[i], as.string = TRUE)
          MinNbSeqFlankPos = (length(labels(Align))/2) + 1  # We follow the approach of the Gblocks server
          # also used by Seaview (Gouy et al. 2010, DOI: 10.1093/molbev/msp259) to set up the
          # Minimum Number Of Sequences For A Flank Position (50% of the overall seqeunce +1 sequence).
          if(LessStringent == TRUE){
            a = paste("./Gblocks ", file.names[i], " -t=", Type, " -b2=", MinNbSeqFlankPos,
                    " -b4=5 -b5=h -e=-gbls", sep = "")
          system(a)
          out = list.files(Gblocks.path)
          file.remove(out[grep(".htm", out, fixed = TRUE)])
          out2 = list.files(Gblocks.path)
          outOri = out2[grep(".fas-gbls", out2, fixed = TRUE)]
          outRena = paste("Gblocksls_", gsub("-gbls",
                                             "", outOri, fixed = TRUE), sep = "")
          # Rename and move to output folder
          file.rename(paste(Gblocks.path, "/", outOri,
                            sep = ""), paste(oldwd, "/", output, "/", outRena, sep = ""))

          #outOri = out2[grep(".fas-gbls", out2, fixed = T)]
          #outRena = paste("Gblocksls_", gsub("-gbls", "", outOri, fixed = TRUE), sep = "")
          #file.rename(paste(Gblocks.path, "/", outOri, sep = ""), paste(output, "/", outRena,
          #                                                    sep = ""))
          } else {
            a = paste("./Gblocks ", file.names[i], " -t=", Type, " -e=-gbms", sep = "")
          system(a)
          out = list.files(Gblocks.path)
          file.remove(out[grep(".htm", out, fixed = TRUE)])
          out2 = list.files(Gblocks.path)
          outOri = out2[grep(".fas-gbms", out2, fixed = TRUE)]
          outRena = paste("Gblocksms_", gsub("-gbms", "", outOri, fixed = TRUE), sep = "")
          # Rename and move to output folder
          file.rename(paste(Gblocks.path, "/", outOri, sep = ""),
                      paste(oldwd, "/", output, "/", outRena, sep = ""))
          }
          # Change the working directory to come back to the current R working directory
          setwd(oldwd)
        }

      }
      # End the  if(os == "windows"){
    } else {

    if (LessStringent == TRUE) {
      # Gblocks runs using the less stringent selection approach.
        Gblocks.lessStr = function(x) {
            Align = seqinr::read.fasta(paste(input, "/", x, sep = ""), as.string = TRUE)
            MinNbSeqFlankPos = (length(labels(Align))/2) + 1  # We follow the approach of the Gblocks server
            # also used by Seaview (Gouy et al. 2010, DOI: 10.1093/molbev/msp259) to set up the
            # Minimum Number Of Sequences For A Flank Position (50% of the overall seqeunce +1 sequence).
            a = paste("Gblocks ", input, "/", x, " -t=", Type, " -b2=", MinNbSeqFlankPos,
                " -b4=5 -b5=h -e=-gbls", sep = "")
            system(a)
        }
        lapply(bb, Gblocks.lessStr)
        out = list.files(input)
        file.remove(paste(input, "/", out[grep(".htm", out, fixed = TRUE)], sep = ""))
        out2 = list.files(input)
        outOri = out2[grep(".fas-gbls", out2, fixed = T)]
        outRena = paste("Gblocksls_", gsub("-gbls", "", outOri, fixed = TRUE), sep = "")
        file.rename(paste(input, "/", outOri, sep = ""), paste(output, "/", outRena,
            sep = "")) # Rename and move to output folder

    } else {
        Gblocks.def = function(x) {
            # Gblocks runs using default parameter (the more stringent selection approach).
            a = paste("Gblocks ", input, "/", x, " -t=", Type, " -e=-gbms", sep = "")
            system(a)
        }
        lapply(bb, Gblocks.def)
        out = list.files(input)
        file.remove(paste(input, "/", out[grep(".htm", out, fixed = TRUE)], sep = ""))
        out2 = list.files(input)
        outOri = out2[grep(".fas-gbms", out2, fixed = TRUE)]
        outRena = paste("Gblocksms_", gsub("-gbms", "", outOri, fixed = TRUE), sep = "")
        file.rename(paste(input, "/", outOri, sep = ""), paste(output, "/", outRena,
            sep = "")) # Rename and move to output folder
    }
    }

    # To report the length of each alignment in a table, as in the Filtering.align.Trimal function.
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
                                                                    Atemp1)], sep = ""), as.string = T)
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
