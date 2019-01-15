#' @title Estimate the best number of partitions and substitution model for each partition: wrapper for PartitionFinder2

#' @description This function is a wrapper running PartitionFinder2 (Lanfear et al. 2017)
#' for a supermatrix (or concatenated alignment) to estimate the best
#' number of partitions and their best associated substitution models.
#'
#' @details Using a fasta alignment and partition file following the
#' RAxML requirements (as exported by \code{\link{Align.Concat}} function) this function
#' prepares an alignment in phylip format ('.phy'), it builds a .cfg input file
#' describing the settings for PartitionFinder2 and calls the program. For details
#' about software installation see
#' \url{http://www.robertlanfear.com/partitionfinder/assets/Manual_v2.1.x.pdf}

#' @return The function returns:
#' \itemize{
#' \item 1) a new partition file in RAxML format according to the best
#' partition scheme proposed by PartitionFinder2 (=output file with a suffix
#' of '_PF2' and the name of the model used in PartitionFinder2),
#' \item 2) an alignment in nexus format with a new block at the end
#' describing the new partitions provided by PartitionFinder2
#' (=name of the input nexus file with
#' suffix '_PF2.nex'). This nexus file is ready to be read by BEAST2 for instance.
#' All the normal PartitionFinder2 outputs are stored in the folder called
#' 'analysis'.
#' }
#'
#' @param input the name (with the path, if necessary) of the input alignment in
#' fasta format.
#'
#' @param Partition the name (with the path, if necessary) of the file
#' providing the gene partition of the alignment file. This file follows the RAxML
#' format, see \url{http://sco.h-its.org/exelixis/resource/download/NewManual.pdf}.
#'
#' @param codon according to the file providing the gene partition in the alignment,
#' codon is the number of the row in this file corresponding to a coding gene that
#' will be split into first, second, and third codon positions, to prepare the
#' .cfg PartitionFinder2 inputfile.
#'
#' @param nexus.file the name (with the path, if
#' necessary) of the concatenated nexus alignment.
#'
#' @param Path.PartiF2 path pointing to
#' where the file 'PartitionFinder.py' is installed (e.g.
#' '/home/davidpc/Programs/PartitionFinder2/partitionfinder-2.1.1/PartitionFinder.py').
#' Note: Apply to both Linux and Windows OS.
#'
#' @param branchlengths 'linked' or 'unlinked'.
#' @param models 'all', 'allx', 'beast',
#' 'mrbayes', 'gamma', 'gammai', '<list>'
#' @param model_selection 'AIC','AICc','BIC'.
#' @param search 'all', 'greedy','rcluster', 'rclusterf', 'hcluster', 'kmeans', 'user'.
#' Type of searching algorithms among the potential substitution models.
#' @param Raxml 'TRUE', 'FALSE', PartitionFinder2 can use RAxML instead of PhyML to get a
#' quicker output, but RAxML implements only three substitution models GTR, GTR+G,
#' GTR+G+I, of the many more available in PhyML.
#' @param nthread number of threads used by PartitionFinder2 to run.
#' @param rcluster_percent See description for rcluster-max, below.
#' @param rcluster_max rcluster-max and rcluster-percent control
#' the thoroughness of the relaxed clustering algorithm together. Setting either
#' of them higher will tend to make the search more thorough and slower. Setting
#' them lower will tend to make the search quicker but less thorough.

#' @details "The rcluster algorithm works by finding the rcluster-max most similar
#' pairs of data blocks, OR the top rcluster-percent of similar datablocks, whichever
#' is smaller. It then calculates the information score (e.g. AICc) of all of these data
#' blocks and keeps the best one. Setting --rcluster-max to 1000 and --rcluster-percent
#' to 10 (i.e. the default values) is usually sufficient to ensure that PartitionFinder2 will
#' estimate a robust partitioning scheme, even on very large datasets in which
#' there may be millions of possible pairs of data blocks."
#' (\url{http://www.robertlanfear.com/partitionfinder/assets/Manual_v2.1.x.pdf})

#' @details For detailed descriptions about Partitionfinder2 options, see
#' \url{http://www.robertlanfear.com/partitionfinder/assets/Manual_v2.1.x.pdf}.

#' @details The function requires that PartitionFinder2 is installed and
#' the path of the python script "PartitionFinder.py" must be provided
#' for the Path.PartiF2 parameter.
#' Currently the function uses a temporary bash script to run Partitionfinder2
#' on Linux plateforms.
#'
#' @references Lanfear et al. 2017, DOI: 10.1093/molbev/msw260
#'
#' @examples
#' \dontrun{
#'
#' # To run the example copy the input files
#' # provided by the package to a temporary directory created into the
#' # current working directory.
#' src.dir = system.file("extdata/multi.align/ForPartiFinder2", package = "regPhylo")
#' dir.create("TempDir.ForPartiFinder2")
#' # Set up the path to the TempDir folder.
#' dest.dir = paste(getwd(), "/TempDir.ForPartiFinder2", sep="")
#' file.names <- dir(src.dir)
#' # Copy all the files stored in regPhylo/extdata/multi.align/ForPartiFinder2"
#' # into a temporary folder.
#' sapply(file.names, function(x) {
#' file.copy(from = paste(src.dir, x, sep = "/"),
#' to = paste(dest.dir, x, sep = "/"),
#' overwrite = FALSE) })
#'
#' input = "TempDir.ForPartiFinder2/Concat.fas"
#' Partition = "TempDir.ForPartiFinder2/Partitions_Concat.txt"
#' # Open the convtab.txt document to determine which genes need
#' # to be partitioned for the first, second, and third codon position.
#' read.delim("TempDir.ForPartiFinder2/convtab.txt", sep = "\t", header = TRUE)
#'
#' # Run the function using RAxML software, BIC criteria,
#' # the fast "rcluster" algorithm and the classic default parameters.
#' PartiFinder2(input = "TempDir.ForPartiFinder2/Concat.fas",
#' Partition = "TempDir.ForPartiFinder2/Partitions_Concat.txt",
#' codon = c(2:4), nexus.file = "TempDir.ForPartiFinder2/Concat.nex",
#' Path.PartiF2 = "/home/davidpc/Programs/PartitionFinder2/partitionfinder-2.1.1/PartitionFinder.py",
#' branchlengths = "linked", models = "all", model_selection = "BIC", search = "rcluster",
#' Raxml = "TRUE", nthread = 5, rcluster_percent = 10, rcluster_max = 1000)
#'
#' # Detailed results of the analysis are stored in the newly created folder "analysis".
#' # The nexus file of the concatenated alignment (supermatrix) including the best
#' # partitioning scheme provided by PartitionFinder2 is called "Concat_PF2.nex",
#' # and the partition file compatible with RAxML is
#' # called "Partitions_Concat.txt_PF2_all.txt", and both are
#' # stored in the temporary directory "TempDir.ForPartiFinder2".
#'
#' # To remove the files created while running the example do the following:
#' # Remove the folder "analysis".
#' unlink("analysis", recursive = TRUE)
#' # Remove the Temporary folder
#' unlink("TempDir.ForPartiFinder2", recursive = TRUE)
#'
#' # Remove the files created by, or for, PartitionFinder2
#' file.remove("Concat.phy")
#' file.remove("log.txt")
#' file.remove("partition_finder.cfg")
#'
#' }
#' @export PartiFinder2

PartiFinder2 = function(input = NULL, Partition = NULL, codon = NULL, nexus.file = NULL,
    Path.PartiF2 = NULL, branchlengths = NULL, models = NULL, model_selection = NULL,
    search = NULL, Raxml = NULL, nthread = NULL, rcluster_percent = NULL, rcluster_max = NULL) {

  inputtree = input
    # Convert a fasta file in phylip format (only format read by PartitionFinder2).
    Concat = ape::read.FASTA(inputtree)
    OriSeqName = labels(Concat)
    CodeSeqName = paste("Seq", seq(1, length(labels(Concat)), by = 1), sep = "")
    names(Concat) = CodeSeqName  # Change the names of the sequences to a simplified code that avoids any complication with the Phylip format due to sequence name length.
    bb = strsplit(inputtree, "/", fixed = T)
    if (length(bb[[1]]) > 1) {
        inputtree = unlist(bb)[length(unlist(bb))]
    }

    ape::write.dna(Concat, file = paste(gsub(".fas", "", inputtree, fixed = T), ".phy",
        sep = ""), format = "sequential", nbcol = -1, colsep = "")
    inputtree2 = paste(gsub(".fas", "", inputtree, fixed = T), ".phy", sep = "")


    namecfg = "partition_finder.cfg"
    # Prepare the .cfg file.
    cat(file = namecfg, "## ALIGNMENT FILE ##", "\n", paste("alignment = ", inputtree2,
        ";", sep = ""), "\n", "\n", sep = "")
    cat(file = namecfg, append = T, "## BRANCHLENGTHS: linked | unlinked ##", "\n",
        paste("branchlengths = ", branchlengths, ";", sep = ""), "\n", "\n", sep = "")
    cat(file = namecfg, append = T, "## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##",
        "\n", paste("models = ", models, ";", sep = ""), "\n", "\n", sep = "")
    cat(file = namecfg, append = T, "# MODEL SELECCTION: AIC | AICc | BIC #", "\n",
        paste("model_selection = ", model_selection, ";", sep = ""), "\n", "\n",
        sep = "")
    cat(file = namecfg, append = T, "## DATA BLOCKS: see manual for how to define ##",
        "\n", "[data_blocks]", "\n", sep = "")
    parti = readLines(Partition)
    parti = gsub("DNA, g", "G", parti, fixed = T)
    nbgene = seq(1, length(parti), by = 1)
    j = 1
    for (j in 1:length(nbgene)) {
        if (is.na(match(nbgene[j], codon)) == "TRUE") {
            cat(file = namecfg, append = T, parti[j], ";", "\n", sep = "")
        } else {
            # Delineate the starting point of the codon position.
            startcod = gsub(" ", "", strsplit(strsplit(parti[j], "=", fixed = T)[[1]][2],
                "-", fixed = T)[[1]][1], fixed = T)
            as.numeric(as.character(startcod)) + 1
            cat(file = namecfg, append = T, paste(strsplit(parti[j], " =", fixed = T)[[1]][1],
                "_pos1 = ", as.numeric(as.character(startcod)), "-", strsplit(strsplit(parti[j],
                  "=", fixed = T)[[1]][2], "-", fixed = T)[[1]][2], "\\3;", sep = ""),
                "\n", sep = "")
            cat(file = namecfg, append = T, paste(strsplit(parti[j], " =", fixed = T)[[1]][1],
                "_pos2 = ", as.numeric(as.character(startcod)) + 1, "-", strsplit(strsplit(parti[j],
                  "=", fixed = T)[[1]][2], "-", fixed = T)[[1]][2], "\\3;", sep = ""),
                "\n", sep = "")
            cat(file = namecfg, append = T, paste(strsplit(parti[j], " =", fixed = T)[[1]][1],
                "_pos3 = ", as.numeric(as.character(startcod)) + 2, "-", strsplit(strsplit(parti[j],
                  "=", fixed = T)[[1]][2], "-", fixed = T)[[1]][2], "\\3;", sep = ""),
                "\n", sep = "")
        }  # End of if else
    }  # End of for
    cat(file = namecfg, append = T, "\n", "## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##",
        "\n", "[schemes]", "\n", "search = ", search, ";", "\n", sep = "")


    if (is.null(rcluster_percent)) {
        rcluster_percent = 10 # default value
    }
    if (is.null(rcluster_max)) {
        rcluster_max = 1000 # default value
    }

    os <- .Platform$OS
    #if os== "windows"
    if(os == "windows"){
      if(missing(Path.PartiF2)){
        stop("The path to the PartitionFinder.py file must be provided in Path.PartiF2")
      }

      # Run PartitionFinder2 partitioning with RAxML
      # (quicker but has only three possible substitution models, GTR, GTR+G, GTR+G+I).
      if (Raxml == "TRUE") {
      a = paste("python ", Path.PartiF2, " ", getwd(), " -p ", nthread,
                " --raxml --rcluster-percent ", rcluster_percent, " --rcluster-max ",
                rcluster_max, sep = "")
      system(a)
      } else {
        # Run PartitionFinder2 with Phyml, a lot of more substitution models, but runs
        # slowly.
        a = paste("python ", Path.PartiF2, " ", getwd(), " -p ", nthread,
                " --rcluster-percent ", rcluster_percent, " --rcluster-max ",
                rcluster_max, sep = "")
        system(a)
      }
      # Else run on Linux
    } else {
    # Run PartitionFinder2. We create a temporary bash script to run partitionfinder2
    # from the root of the computer.  Run PartitionFinder2 partitioning with RAxML
    # (quicker but has only three possible substitution models, GTR, GTR+G, GTR+G+I).
    if (Raxml == "TRUE") {
        cat(file = "Test.sh", "#!/bin/bash", "\n", "\n", "cd ", "\n", "python ",
            Path.PartiF2, " ", getwd(), " -p ", nthread, " --raxml --rcluster-percent ",
            rcluster_percent, " --rcluster-max ", rcluster_max, "\n", sep = "")
    } else {
        # Run PartitionFinder2 with Phyml, a lot of more substitution models, but runs
        # slowly.
        cat(file = "Test.sh", "#!/bin/bash", "\n", "\n", "cd ", "\n", "python ",
            Path.PartiF2, " ", getwd(), " -p ", nthread, " --rcluster-percent ",
            rcluster_percent, " --rcluster-max ", rcluster_max, "\n", sep = "")
    }
    a = "chmod u+x Test.sh"  # Attribute the rights to execute the bash script.
    system(a)
    a = "./Test.sh"  # Run the bash script.
    system(a)
    file.remove("Test.sh")  # Remove the temporary bash script used to call PartitionFinder2.
}

    # Extract the output
    parti = readLines("analysis/best_scheme.txt")
    parti = parti[-which(parti == "")]
    # Follow the RAxML format
    parti2 = parti[c(grep("^RaxML-style", parti, perl = TRUE):grep("^MrBayes block",
        parti, perl = TRUE))]
    parti3 = parti2[-c(1, 2, length(parti2))]
    cat(file = paste(Partition, "_PF2_", models, ".txt", sep = ""), parti3, sep = "\n")

    options(warn=-1)
    # Prepare a nexus file ready for BEAST2.
    parti2_b = parti[c(grep("^Nexus formatted character sets", parti, perl = TRUE):
                         grep("^Nexus formatted character sets for IQtree",
                              parti, perl = TRUE))]
    parti3_b = parti2_b[-c(1, length(parti2_b))]
    nexusf = readLines(nexus.file)
    cat(file = paste(gsub(".nex", "", nexus.file, fixed = T), "_PF2.nex", sep = ""),
        nexusf, sep = "\n")
    cat(file = paste(gsub(".nex", "", nexus.file, fixed = T), "_PF2.nex", sep = ""),
        append = T, "\n", parti3_b, "\n", sep = "\n")
    options(warn=0)
}  # End of the function.
