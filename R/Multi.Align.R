#' @title Align sequences using multiple alignment programs

#' @description This function enables to run several alignment softwares in parallel, including Muscle v3.8.31 (Edgar 2004),
#' Mafft v7.222 fftns1, fftns2, fftnsi (Katoh et al. 2005), Prank v.140603 (Loytymoja & Goldman 2008)
#' and PASTA (Mirarab et al. 2015).
#' \strong{Note:} PASTA is not available to run for Windows users.

#' @details First, for each distinct gene region, the function aligns all the sequences
#' using Mafft v7.222 fftns1 (Katoh et al. 2005)
#' with the option \emph{--adjustdirection} to automatically detect the direction of the
#' sequences and reverse complement the sequences if necessary (the function adds
#' an '_R_' as a prefix on the sequence name if it has been reverse complemented).
#' Second, the function runs Mafft fftns2, fftnsi, Prank and PASTA, and
#' alphabetically orders the sequences within each of the gene region alignment.
#'
#' @details The function requires that Muscle, Mafft, Prank and Pasta are installed and in
#' the PATH.
#'
#' @details The 'output', 'input' and 'nthread' objects  the need
#' to be present in the R environment before running the function
#' because the function runs in parallel using the R package 'parallel' (see
#' example).

#' @param input path to the folder storing the files of interest for further alignment.
#' The files should have the extension '.fas' and the name of the gene region
#' should be immediately before the extension separated from the rest of the name
#' by '_' (e.g.  'Alignment_co1.fas').
#'
#' @param output path to the folder storing the
#' final alignments. The output folder is created directly by the function. The
#' name of the program appears as a prefix in the name of the alignment file using
#' '_' as a separator. The function also exports a table into the input folder containing
#' all the names of the sequences, per gene region, that have been reverse
#' complemented.
#'
#' @param methods programs used to align the sequences. The function used Mafft-FFTNS1
#' to reverse complement and align quickly the sequences,
#' but other programs can be used simultaneously, such as Mafft-fftns2, Mafft-fftnsi,
#' Muscle, Prank and PASTA, c("mafftfftns2", "mafftfftnsi", "muscle", "prank", "pasta").
#'
#' @param nthread number of threads used to run the alignment software in
#' parallel for the different gene regions.
#'
#' @param Mafft.path for the Windows platform, a character string which provides the path
#' to the Mafft executable
#' (e.g. Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/
#' mafft-7.409-win64-signed/mafft-win/mafft") (see examples below). Eventhough Mafft is not
#' the alignment you are after, mafft is necessary to reverse complement the sequences
#' and the Mafft.path must be completed in all cases.
#' For Linux the mafft software must be in the $PATH.

#' @param Muscle.path for the Windows platform, a character string which provides the path
#' to the Muscle executable (e.g. "C:/Users/deme/Documents/Programs/Muscle/muscle3.8.31_i86win32.exe").
#' For Linux the Muscle software must be in the $PATH.
#'
#' @param Prank.path for the Windows platform, a character string which provides the path
#' to the Prank executable
#' (e.g. "C:/Users/deme/Documents/Programs/Prank/prank.windows.140603/prank/bin/prank.exe")
#' For Linux the Prank software must be in the $PATH.
#'
#'
#' @examples # Load into the R environment the object used by the function.
#' # Here the example consists of four alignment (co1, 12srrna, cytb, rag1)
#' # for 16 species.
#' \dontrun{
#'
#' # To run the example copy the input alignment files
#' # provided by the package to a temporary directory created in the
#' # current working directory.
#' src.dir = system.file("extdata/multi.align", package = "regPhylo")
#' dir.create("TempDir.Multi.aligned")
#' # Set up the path to the TempDir folder.
#' dest.dir = paste(getwd(), "/TempDir.Multi.aligned", sep="")
#' file.names <- dir(src.dir)
#' # Copy all the files stored in regPhylo/extdata/multi.align"
#' # into a temporary folder.
#' sapply(file.names, function(x) {
#' file.copy(from = paste(src.dir, x, sep = "/"),
#' to = paste(dest.dir, x, sep = "/"),
#' overwrite = FALSE) })
#'
#' output = "multi.alignedTest"
#' nthread = 3
#' methods = c("mafftfftnsi", "pasta")
#' input = "TempDir.Multi.aligned"
#' # run the function using two additional alignment program "mafftfftnsi" and "pasta"
#' # (not availbale on Windows OS, see the example for Windows below).
#' Multi.Align(input = input, output = "multi.alignedTest", nthread = 3,
#' methods = c("mafftfftnsi", "pasta"))
#' list.files("multi.alignedTest")
#'
#' # To clean the file created while running the example do the following:
#' # Remove the temporary folder
#' unlink("TempDir.Multi.aligned", recursive = TRUE)
#' # Remove the folder with Gblocks outputs
#' unlink("multi.alignedTest", recursive = TRUE)
#'
#'
#' #### To run Multi.Align on Windows OS. ###
#'
#' output = "multi.alignedTest"
#' nthread = 3
#' methods = c("mafftfftnsi", "muscle", "prank")
#' input = "TempDir.Multi.aligned"
#'
#' # Example to run "mafftfftnsi", "muscle", and "prank" on Windows.
#' # (Mafft.path, Muscle.path, Prank.path need to be adapted to the user's
#' # computer configuration, and are provided hereas an example).
#'
#' Multi.Align(input = input, output = "multi.alignedTest", nthread = 3,
#' methods = c("mafftfftnsi", "muscle", "prank"),
#' Mafft.path = "C:/Users/deme/Documents/Programs/Mafft/mafft-7.409-win64-signed/mafft-win/mafft",
#' Muscle.path = "C:/Users/deme/Documents/Programs/Muscle/muscle3.8.31_i86win32.exe",
#' Prank.path = "C:/Users/deme/Documents/Programs/Prank/prank.windows.140603/prank/bin/prank.exe")
#'
#' #### end specifications to run Multi.Align on Windows OS. ####
#'
#'
#' # The output files are also present as external data in the regPhylo package
#' # and can be accessed by the following code:
#' # a = system.file("extdata/multi.align/multi.aligned", package = "regPhylo")
#' # list.files(a)
#' }
#'
#' @export Multi.Align
#'
#' @references Edgar 2004, DOI: 10.1093/nar/gkh340
#' @references Katoh et al. 2005, DOI: 10.1093/nar/gki198
#' @references Loytymoja & Goldman 2008, DOI: 10.1126/science.1158395
#' @references Mirarab et al. 2015, DOI: 10.1089/cmb.1014.0156


Multi.Align = function(input = NULL, output = NULL, nthread = NULL, methods = NULL, Mafft.path = NULL,
                       Muscle.path = NULL, Prank.path = NULL) {
    AlignSelect = list.files(input)
    AlignSelect = AlignSelect[grep(".fas", AlignSelect)]
    AlignSelect2 = paste("revc_", AlignSelect, sep = "")
    dir.create(output)  # Create the ouput folder.

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

    # Copy the revc_ file in the output folder
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
    if(length(which(methods=="muscle"))==1){

      muscle = "muscle"
      os <- .Platform$OS
      if(os == "windows"){
        if(missing(Muscle.path)){
          stop("The path to the muscle executable must be provided in Muscle.path")
        }
        muscle = Muscle.path
      }

    # Muscle alignment:
    muscle.align = function(x) {
        a = paste(muscle, " -in ", input, "/", x, " -out ", output, "/", "Muscle_",
            x, sep = "")
        system(a)
    }
    parallel::parLapply(cl, AlignSelect2, muscle.align)
    }

    if(length(which(methods=="mafftfftns2"))==1){

      mafft = "mafft"
      os <- .Platform$OS
      if(os == "windows"){
        if(missing(Mafft.path)){
          stop("The path to the mafft executable must be provided in Mafft.path")
        }
        mafft = Mafft.path
      }

    # Mafft alignment using FFT-NS-2 algorithm:
    mafftfftns2.align = function(x) {
        a = paste(mafft, " --retree 2 --maxiterate 0 ",
                  input, "/", x, " > ", output, "/",
                  "Mafftfftns2_", x, sep = "")
        system(a)
    }
    parallel::parLapply(cl, AlignSelect2, mafftfftns2.align)
    }

    if(length(which(methods=="mafftfftnsi"))==1){

      mafft = "mafft"
      os <- .Platform$OS
      if(os == "windows"){
        if(missing(Mafft.path)){
          stop("The path to the mafft executable must be provided in Mafft.path")
        }
        mafft = Mafft.path
      }

    # Mafft alignment using fftnsi algorithm (increase in alignment accuracy for
    # distantly related sequences)
    mafftfftnsi.align = function(x) {
        a = paste(mafft, " --retree 2 --maxiterate 2 ",
                  input, "/", x, " > ", output, "/",
                  "Mafftfftnsi_", x, sep = "")
        system(a)
    }
    parallel::parLapply(cl, AlignSelect2, mafftfftnsi.align)
    }

    if(length(which(methods=="prank"))==1){
    # PRANK alignment (Phylogeny aware alignment approach distinguishing the
    # insertion and deletion evolutionary event [especially good for distantly
    # related species]).

      prank = "prank"
      os <- .Platform$OS
      if(os == "windows"){
        if(missing(Prank.path)){
          stop("The path to the prank executable must be provided in Prank.path")
        }
        prank = Prank.path
      }


    prank.align = function(x) {
        a = paste(prank, " -d=", input, "/", x, " -o=", output, "/", "Prank_", x, sep = "")
        system(a)
    }
    parallel::parLapply(cl, AlignSelect2, prank.align)

    # Rename the Prank output alignment.
    d = list.files(output)
    d1 = paste(output, "/", d[grep(".fas.best.fas", d, fixed = T)], sep = "")
    d2 = gsub(".fas.best.fas", ".fas", d1, fixed = T)
    i = 1
    for (i in 1:length(d1)) {
        file.rename(from = d1[i], to = d2[i])
    }
    } ### end if(length(which(methods=="prank"))==1){

    if(length(which(methods=="pasta"))==1){
    # PASTA alignment (Mirarab et al. 2015) based on the 'divide and conquer'
    # approach to co-infer alignment and tree, based on SATEII (Liu et al. 2011, DOI:
    # 10.1093/sysbio/syr095) and transitivity.
    pasta.align = function(x) {
        a = paste("run_pasta.py -i ", input, "/", x, " -j PASTA_", " --num-cpus 1",
            sep = "")
        system(a)
    }
    parallel::parLapply(cl, AlignSelect2, pasta.align)
    # Move the final PASTA alignments into the same folder as the other alignments.
    # Define the names of the PASTA alignment files to move.
    b = list.files(input)
    b1 = paste(input, "/", b[grep(".aln", b)], sep = "")
    # Rename the files before moving them into the other folder
    bb = gsub(".aln", ".fas", gsub("(_[0-9]*.marker001.)", "_", b1, perl = T), fixed = T)
    i = 1
    for (i in 1:length(bb)) {
      file.rename(from = b1[i], to = bb[i])
    }
    # Copy the renamed files into the appropriate output folder.
    file.copy(from = bb, to = output)
    # Remove the temporary files created by Pasta in the working directory.
    c = paste(input, "/", b[grep("PASTA_", b)], sep = "")
    file.remove(c)
    }

    parallel::stopCluster(cl)  # Stop Cluster.


    # Remove the temporary files used for the alignment.
    bn = list.files(input)
    file.remove(paste(input, "/", bn[grep("revc_", bn, fixed = T)], sep = ""))

    # Remove the '_revc_' in the name of the alignment in the output folder.
    b = list.files(output)
    bb = paste(output, "/", b[grep("_revc_", b, fixed = T)], sep = "")
    bb2 = gsub("_revc", "", bb, fixed = T)
    for (i in 1:length(bb)) {
        file.rename(from = bb[i], to = bb2[i])
    }

    # Automatically detect the number of alignment programs used and the number of
    # genes introduced.
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

    # Order the sequences alphabetically in each alignment (using the sequence name).
    x = Unigene
    # listgeneb = vector()
    listRevComp = matrix(NA, ncol = 2)[-1, ]

    # Loop over multiple genes.
    i = 1
    for (i in 1:length(x)) {

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
          # DFtempord = DFtemp[order(DFtemp[, 1]), ]  # Re-order the sequences in alphabetical order.
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
          listRevComp = rbind(listRevComp, cbind(rep(Unigene[i], 1), NA))
        }  # End if else.
      } # end for if(inherits(listAlig, "simpleError")==FALSE) {
    }  # End for i.
    # Return a table with the name of the sequences reverse complemented for each
    # gene region.
    colnames(listRevComp) = c("Gene", "SequenceName")
    utils::write.table(listRevComp, file = paste(input, "/ListSeq_ReverseComplemented.txt",
        sep = ""), row.names = F, sep = "\t")
    return(listRevComp)
}  # End the function Multi.Align.R
