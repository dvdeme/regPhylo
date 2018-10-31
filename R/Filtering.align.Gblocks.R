#' @title Remove poorly aligned nucleotides position using Gblocks software.

#' @description This function uses the software Gblocks V 0.91b (Castresana 2000)
#' to filter out poorly aligned positions and divergent regions in an alignment.
#' It allows to specify if the sequences are proteins, DNA non coding or DNA coding regions.
#' The function offers the possibility to use the default parameter (more stringent selection)
#' or the less stringent selection approach (parameter equivalent that on the Gblocks server,
#' http://molevol.cmima.csic.es/castresana/Gblocks_server.html).

#' @param input name of the folder storing the alignments in fasta format (with the
#' extension '.fas').
#' @param LessStringent if 'TRUE' allows an option for a less
#' stringent selection (equivalent of ticking the three options in Gblocks server,
#' 'Allow smaller final blocks', 'Allow gap positions within the final blocks',
#' 'Allow less strict flanking positions', the output has the extension -gbls),
#' otherwise uses the default parameters in Gblocks for a stringent selection
#' (output with the extension -gbms).
#' @param Type Type of Sequence can be Protein (p), DNA (d), or Codons (c).
#' @param output path to the folder storing the trimmed alignments.
#' The output folder is created automatically.
#' @param remove.empty.align If TRUE, the empty alignments are excluded from the computation.
#' @param Gblocks.path for Windows plateform, a character string which provides the path
#' to the Gblocks executable. For Linux the Gblocks software must be in the $PATH.
#'
#' @details The function requires, Gblocks to be installed and set up in the PATH.

#' @return The function exported the trimmed alignement directly in the same folder
#' as the input alignments. In the R environment the function returns a table with
#' the length of the different alignments for each gene region and alignment program.
#'
#' @examples # Run the function for DNA, using the less stringent
#' # selection heuristic, and removing the potential empty alignments.
#' \dontrun{
#' # To run the example we have to copy the input alignment files
#' # provided by the package to a temporary directory created into the
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
#' # Gblocks in the "Trimmed-Gblocks" folder.
#' Filtering.align.Gblocks(input = "TempDir", LessStringent = TRUE,
#' output = "TrimmedGblocks", Type = "d", remove.empty.align = TRUE)
#'
#'
#' # To clean the file created while running the example do the following:
#' # Remove the temporary folder
#' unlink("TempDir", recursive = TRUE)
#' # Remove the folder with Gblocks outputs
#' unlink("TrimmedGblocks", recursive = TRUE)
#' }
#'
#' @export Filtering.align.Gblocks
#'
#' @references Castresana 2000,
#' DOI: 10.1093/oxfordjournals.molbev.a026334

Filtering.align.Gblocks = function(input = NULL, LessStringent = NULL,
                                   Type = NULL, output = NULL,
                                   remove.empty.align = NULL, Gblocks.path = NULL) {
    b = list.files(input)
    bb = b[grep(".fas", b, fixed = TRUE)]

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
      dir.create(output)
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
}  # End of the function.
