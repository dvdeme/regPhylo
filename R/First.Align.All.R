#' @title Reverse complement and align the sequences using MAFFT and PASTA

#' @description First, this function aligns all the sequences for each distinct gene region
#' using Mafft v7.222 fftns1 (Katoh et al. 2005) with the
#' option \emph{--adjustdirection} (i.e. to automatically detect the direction of the sequences
#' and reverse complement the sequences if necessary). The function adds a '_R_' as a
#' prefix to the sequence name if it is reverse complemented. Second, the function runs Mafft fftnsi (Katoh et
#' al. 2005), and PASTA (Mirarab et al. 2015) in parallel and alphabetically orders
#' the sequences within each alignment (i.e. for each gene region).
#' \strong{Note:} PASTA is not available to run for Windows users.
#'
#' @details PASTA (as Sate I and II) may be better at handling very
#' divergent sequences and at dealing with saturated phylogenetic signal because
#' it splits the alignment into small chunks of similar sequences and performs
#' profile alignments among those smaller chunk of sequences (Roquet et al. 2013).

#' @param input path to the folder storing the files of interest. The files must have
#' the extension '.fas' and the name of the gene region must be immediately before
#' the extension separated from the rest of the name by '_' (e.g.
#' 'Alignment_co1.fas').
#'
#' @param output path to the folder storing the final alignments,
#' (the output folder is created directly by the function). The name of the
#' program appears as a prefix in the name of the alignment file using '_' as a
#' separator. The function also exports a table into the input folder with the name
#' of all the sequences per gene region that have been reverse complemented.
#'
#' @param nthread number of threads used to run the alignment software in parallel for
#' the different gene regions.
#'
#' @param methods programs used to align the sequences, the function used Mafft-FFTNS1
#' to reverse complement and align quickly the sequences, but other programs can be also
#' used simultaneously, such as  Mafft-fftnsi and PASTA, c("mafftfftnsi", "pasta").
#'
#' @param Mafft.path for the Windows platform, a character string which provides the path
#' to the mafft executable
#' (e.g. Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft")
#' (see examples below). Eventhough Mafft is not the alignment you
#' are after, mafft is necessary to reverse complement the sequences and the Mafft.path must
#' be completed in all cases. For Linux the mafft software must be in the $PATH.
#'
#' @details The 'output', 'input' and 'nthread' objects need
#' to be present in the R environment before running the function because the
#' function runs in parallel using the R package 'parallel', see example below the
#' function.

#' @details The function requires that Mafft and PASTA are installed and in the PATH.
#' \itemize{ \item For Mafft, download for the different OS and all documentation are available at
#' \url{https://mafft.cbrc.jp/alignment/software/}, (for Windows OS we have only tested the
#' "All-in-one package for Windows" version of Maftt).
#' \item For PASTA, download and instructions for installation are provided on the Github page
#' available at: \url{https://github.com/smirarab/pasta}.
#' }


#' @return The function returns the alignment in the output folder using the
#' name of the program as a prefix of the alignment (i.e. "mafftfftnsi_Alig_co1.fas").
#' The function also returns a table called "ListSeq_RevCompl_FirstAlignAll.txt"
#' into the output folder listing the sequence names that have been reverse
#' complemented by Mafft, the same table is exported into the R environment as well.
#'

#' @examples # The input, output and nthread object have to be present in
#' # the R global environment before running the function.
#' \dontrun{
#' # To demonstrate the function we provided small example alignments in the
#' # inst/extdata/FirstToAlign folder available with the package.
#' # These alignments are exported from the Seq.DF5 dataset, for
#' # the co1 and 16s regions. Notice that the 16s has only one
#' # sequence and this is explained in the warning message.
#'
#' # To run the example, copy the input alignment files
#' # provided by the package to a temporary directory in your
#' # current working directory.
#' src.dir = system.file("extdata/FirstToAlign", package = "regPhylo")
#' dir.create("TempDir.FirstToAlign")
#' # Set up the path to the TempDir folder.
#' dest.dir = paste(getwd(), "/TempDir.FirstToAlign", sep="")
#' file.names <- dir(src.dir)
#' # Copy all the files stored in regPhylo/extdata/FirstToAlign"
#' # into a temporary folder.
#' sapply(file.names, function(x) {
#' file.copy(from = paste(src.dir, x, sep = "/"),
#' to = paste(dest.dir, x, sep = "/"),
#' overwrite = FALSE) })
#'
#' # Remove the empty file created for the folder name when
#' # importing the data file
#' file.remove("TempDir.FirstToAlign/FirstAligned")
#'
#' input = "TempDir.FirstToAlign"
#' output = "TempDir.FirstToAlign/FirstAligned"
#' nthread = 2
#' methods = c("mafftfftnsi", "pasta")
#'
#' First.Align.All(input = input, output =
#' "TempDir.FirstToAlign/FirstAligned", nthread = 2,
#' methods = c("mafftfftnsi", "pasta"))
#'
#'
#' #### To run FirstAlign.All on Windows OS. ###
#'
#' input = "TempDir.FirstToAlign"
#' output = "TempDir.FirstToAlign/FirstAligned"
#' nthread = 2
#' methods = "mafftfftnsi"

#' First.Align.All(input = input, output =
#' "TempDir.FirstToAlign/FirstAligned", nthread = 2,
#' methods = "mafftfftnsi",
#' Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft")
#'
#' #### end specifications to run FirstAlign.All on Windows OS. ####
#'
#'
#' # To remove the file created while running the example do the following:
#' unlink("TempDir.FirstToAlign", recursive = TRUE)
#'
#' }

#' @references Katoh et al. 2005, DOI: 10.1093/nar/gki198
#' @references Mirarab et al. 2015 DOI: 10.1089/cmb.1014.0156
#' @references Roquet et al. 2013 DOI: 10.1111/j.1600-0587.2012.07773.x
#'
#' @export First.Align.All

First.Align.All = function(input = NULL, output = NULL, nthread = NULL, methods = NULL , Mafft.path = NULL) {
    AlignSelect = list.files(input)
    AlignSelect = AlignSelect[grep(".fas", AlignSelect)]
    AlignSelect2 = paste("revc_", AlignSelect, sep = "")
    dir.create(output)  # Create the ouput folder

    mafft = "mafft"
    os <- .Platform$OS
    if(os == "windows"){
      if(missing(Mafft.path)){
        stop("The path to the mafft executable must be provided in Mafft.path")
      }
      mafft = Mafft.path
    }

    # Run the job in parallel through a certain number of threads.
    cl <- parallel::makeCluster(nthread)  # Create the cluster.
    # Designate the functions and variables that need to be exported for the parallel version.
    parallel::clusterExport(cl, varlist = c("output", "input", "nthread", "methods"))

    # Mafft alignment using FFT-NS-1 algorithm This step also checks the necessity to
    # reverse complement the sequence. The MAFFT FFTNS1 alignment will serve as input
    # alignment for all other algorithms
    mafftfftns1.align = function(x) {
        a = paste(mafft, " --retree 1 --maxiterate 0 --adjustdirection ", input, "/",
            x, " > ", input, "/", "revc_", x, sep = "")
        system(a)
    }
    parallel::parLapply(cl, AlignSelect, mafftfftns1.align)

    # Copy the revc_ file in the output folder.
    p = list.files(input)
    # Copy in the new folder.
    file.copy(from = paste(input, "/", p[grep("revc_", p, fixed = T)], sep = ""),
        to = output)

    # Rename the revc_ files in Mafftfftns1 files.
    p1 = list.files(output)
    p2 = gsub("revc_", "Mafftfftns1_", p1, fixed = T)
    # Rename in the new folder.
    for (i in 1:length(p1)) {
        file.rename(from = paste(output, "/", p1[i], sep = ""), to = paste(output,
            "/", p2[i], sep = ""))
    }
    parallel::stopCluster(cl)  # Stop Cluster.

    cl <- parallel::makeCluster(nthread)  # Create the cluster.
    # Designate the functions and variables that need to be exported for the parallel version.
    parallel::clusterExport(cl, varlist = c("output", "input", "nthread", "methods"))

    if(length(which(methods=="mafftfftnsi"))==1){

      mafft = "mafft"
      os <- .Platform$OS
      if(os == "windows"){
        if(missing(Mafft.path)){
          stop("The path to the mafft executable must be provided in Mafft.path")
        }
        mafft = paste(Mafft.path, " ", sep = "")
      }


    # Mafft alignment using fftnsi algorithm.
    mafftfftnsi.align = function(x) {
        a = paste(mafft, " --retree 2 --maxiterate 2 ",
                  input, "/", x, " > ", output, "/",
                  "Mafftfftnsi_", x, sep = "")
        system(a)
    }
    parallel::parLapply(cl, AlignSelect2, mafftfftnsi.align)
    }

    options(warn=-1)
    if(length(which(methods=="pasta"))==1){
    # PASTA alignment (Mirarab et al. 2015) based on 'divide and conquer' approach to
    # co-infer alignment and tree, based on SATEII (Liu et al. 2011, DOI:
    # 10.1093/sysbio/syr095) and transitivity.
      if(os == "windows") {
        warning("Pasta cannot be used on a Windows plateform at the moment")
      } else {
        pasta.align = function(x) {
            a = paste("run_pasta.py -i ", input, "/", x, " -j PASTA_", " --num-cpus 1",
                      sep = "")
        system(a)
        }
        parallel::parLapply(cl, AlignSelect2, pasta.align)

        # Move the final PASTA alignments into the same folder as the other alignments
        # Define the names of the PASTA alignment files to move.
        b = list.files(input)
        b1 = paste(input, "/", b[grep(".aln", b)], sep = "")
       # Rename the files before to move them in the other folder.
        bb = gsub(".aln", ".fas", gsub("(_[0-9]*.marker001.)", "_", b1, perl = T), fixed = T)
        i = 1
        for (i in 1:length(bb)) {
          file.rename(from = b1[i], to = bb[i])
        }
        # Copy the renamed files to an appropriate output folder.
        file.copy(from = bb, to = output)
        # Remove the temporary files created by PASTA in the working directory.
        c = paste(input, "/", b[grep("PASTA_", b)], sep = "")
        file.remove(c)
      }
    }
    options(warn=0)
    parallel::stopCluster(cl)  # Stop Cluster.


    # Remove the temporary file used for the alignment
    bn = list.files(input)
    file.remove(paste(input, "/", bn[grep("revc_", bn, fixed = T)], sep = ""))

    # Remove the '_revc_' in the name of the alignment in the output folder
    b = list.files(output)
    bb = paste(output, "/", b[grep("_revc_", b, fixed = T)], sep = "")
    bb2 = gsub("_revc", "", bb, fixed = T)
    for (i in 1:length(bb)) {
        file.rename(from = bb[i], to = bb2[i])
    }

    # Automatically detect the number of alignment programs used and the number of
    # genes introduced
    a = list.files(output)
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

    # Re-order the sequences in each alignment in the same alphabetical order (using
    # the sequence name).
    x = Unigene
    listRevComp = matrix(NA, ncol = 2)[-1, ]
    i = 1
    for (i in 1:length(x)) {
        # Loop over multiple genes.
        #listgeneb = c(listgeneb, a[which(Nbgene == x[i])])
        listgeneb = a[which(Nbgene == x[i])]
        j = 1
        # Loop over multiple alignments of the same gene.
        for (j in 1:length(listgeneb)) {
          listAlig = tryCatch(seqinr::read.fasta(paste(output, "/", listgeneb[j], sep = ""),
                                                 as.string = T), error = function(e) e)
          if(inherits(listAlig, "simpleError")==FALSE) {
            SeqName = labels(listAlig)
            SeqT = vector()
            k = 1
            for (k in 1:length(listAlig)) SeqT = c(SeqT, listAlig[[k]][1])
            DFtemp = cbind(SeqName, SeqT)
            if(length(listAlig) == 1){
              DFtempord = t(as.matrix(DFtemp[order(DFtemp[, 1]), ]))
              } else {
              DFtempord = DFtemp[order(DFtemp[, 1]), ]
            }
            # Re-build and export the ordered alignments.
            Seq_Name = paste(">", paste(DFtempord[, 1], DFtempord[, 2], sep = "|_|"),
                             sep = "")
            AlignFasta = unlist(strsplit(Seq_Name, "|_|", Seq_Name, fixed = TRUE))
            write(AlignFasta, file = paste(output, "/", listgeneb[j], sep = ""))  # Write the alignment.
          } else {
            warning(paste("The alignment ", listgeneb[j], " contains no sequence, either because the original alignment contains 1 sequence only, or because PASTA couldn't find a proper tree and crashed", sep=""))
          }
        }  # End for j.

        # For each gene region, list the sequences that have been reverse complemented by
        # Mafft.

        listAli = tryCatch(seqinr::read.fasta(paste(output, "/", a[grep(x[i], a)][1], sep = ""),
            as.string = T), error = function(e) e)
        if(inherits(listAlig, "simpleError")==FALSE) { ### in case of non empty alignment
        SeqNa = labels(listAli)
        RevComp = SeqNa[grep("^_R_", SeqNa, perl = T)]  # For each gene retrieve the name of the sequences that have been reverse complemented.
        if (length(RevComp) > 0) {
            listRevComp = rbind(listRevComp, cbind(rep(Unigene[i], length(RevComp)),
                RevComp))
        } else {
            # If none of the sequence have been reverse complemented.
            listRevComp = rbind(listRevComp, cbind(rep(Unigene[i], 1), "None of the sequence has been reverse complemented"))
        }  # End if else.
        } # end for if(inherits(listAlig, "simpleError")==FALSE) {
    }  # End for i.
    # Return a table with the name of the sequences reverse complemented for each
    # gene region.
    colnames(listRevComp) = c("Gene", "SequenceName_ReverseComplemented")
    utils::write.table(listRevComp, file = paste(input, "/ListSeq_RevCompl_FirstAlignAll.txt",
        sep = ""), row.names = F, sep = "\t")
    return(listRevComp)
}  # End the function.
