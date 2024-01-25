#' @title Help to detect outlier sequences in an alignment and remove them.

#' @description This function helps to detect some potential outlier sequences
#' (i.e. mis-aligned sequences, which might be caused by different problems such
#' as gene or species annotation problems, presence of paralogous sequences, etc)
#' that should be removed from the pool of sequences (per species and gene regions)
#' before selecting the best sequence.  This function is only provided as a tool
#' to HELP detect potential outlier sequences, great care must be taken to
#' ensure sequences are not erroneously omitted as potential outliers,
#' or remain included undetected outlier sequences. We recommend checking
#' alignments by eye to ensure your confidence in the retained sequences.


#' @details This function uses the following general strategy to detect outlier
#' sequences.  First, the function detects and removes the duplicated sequences to work considering only the haplotypes. Second, a distance matrix between all pairwise sequences is built considering different strategies implemented in the dist.DNA function from the ape R package.
#' Third the function computes a summary statistic (can be the "median" or the "min" for the minimal distance) for each focal sequence with the genetic distance with all the other sequences using the Sum.Stat option. Fourth, the distribution of summary statistic is compared to an outlier threshold defined by the Outlier.Th option. Outlier.Th is a multiplier of the interquartile range distance as defined by the "range" parameter from the boxplot R function. All the sequences with a summary statistic above the threshold are considered as potential outliers and are removed from the final alignment. Then, the deletion gap only sites are removed automatically and the retained sequences can be re-aligned using a MAFFT FFT NS1 algorithm if desired.
#' In order to speed up the detection of outlier sequences, we also provide an approximate option to compute the pairwise genetic distances among sequences. The option "Approxim = TRUE" allows to perform such approximative detection by randomly selecting a number of sequences (set by the option nb.seq) several times (set by the option nb.rep) to compare the average summary statistic of each focal sequences. We can here higly parallelized this approximate search. By default, the function selects 25 random sequences 5 times, increasing the number of random sequences increases the accuracy but increase the computational time and the memory load. Based on several tests, we considered that approximate estimate are worth trying for large dataset above  5000 / 10000 sequences, because the parallel computing start to subtantially decrease the computational time for minimal differences in outlier sequences detection. 



#' @return The function return a list of 4 objects using the following headers:
#' \itemize{
#' \item 'Outlier.Sequences' which provides a vectorwith teh names of all the sequences considered as potential outlier.
#' \item 'Clean.Haplo.Alignment' is a DNAbin object of the alignment of only the haplotypes sequences after removing the outlier sequences.
#' \item 'Clean.All.Seq.Alignment' is a DNAbin object of the alignment of all the sequences after removing the outlier sequences.
#' \item 'First.Dist.Mat' is the summary statistic (either median or minimal distance, this parameter is defined by the 'Sum.Stat' option) of the genetic distance from the each focal sequences in comparison to all the other.
#' }


#' @param inputal an object of class "alignment" (seqinr R package),
#' or "DNAbin" (ape R package) or the name (including the path if necessary) of the input fasta file. The input file must be a DNA alignement file already aligned.

#' @parameter Outlier.Th is a number (by default 1.5) defining the distance to consider an outlier sequences. The Outlier.Th is a multiplier of the interquartile range from the genetic distances distribution using the Sum.stat (eg. median). This parameter is equivalent of the "range" parameter from the boxplot function from the stat R package.


#' @param output name of the two alignment exported. The name can also include the path the folder.
#' Two separate alignment are exported in a fasta format, the alignment with the ".Clean.haplo.fas" siffix includes only the haplotypes after removing outlier sequences, the alignment with the ".Clean.AllSeq.fas" suffix includes all the sequences after removing the outlier sequences.

#' @param Evol.Model evolutionary model used to compute the DNA distance among sequences, default "TN93" (Tamura and Nei 1993), This options accepts all the models used in dist.dna function from the ape R package.  

#' @param Sum.Stat type of symmary statistic ("median" by default)  to compare the DNA distance of the focal sequence with the other, by defult we use the median distance but it can also be the "min" for the minimal dNA distance between 1 sequence and its closest neighbor.

#' @param Iter.optim if TRUE then an iterative outlier detection appraoch is performed until no longer outlier sequenece are detected considering the outier threshold defined by the Outlier.Th option.

#' @param Re.Align the sequences are realigned after the indel gaps were removed using fast mafft.fftns1 algorithm, but this steps can be long especially combined the automatic iteration of the optimisation approach that can be perfomed multiple times.

#' @param Approxim if TRUE then the function perform an approximate estimate of the outlier sequence detection by using a random selection of nb.seq (25 by default) sequences nb.rep (5 by default) times for each sequence, this options is worth using for large alignement including more than 10 000 sequences.

#' @param  nb.seq is a number defining how many sequence need to be selected randomly if Approxim is TRUE, increasing the number of sequence will considerably slow dow the algorithm. default 25.

#' @param nb.rep is a number defining how many replicate of the randomw selection of sequence must be performed for each sequence of interest (default 5 times)

#' @param nbthreads a number (by default 1) defining the number of threads used to perfoimed the Approximate detection of outlier sequences 



#' @examples # Load the alignment file (class "alignment"), this is the first
#' # object of the list called Example_16S_outlier.
#' data(Example_16S_outlier)
#' Example_16S_outlier_align = Example_16S_outlier[[1]]
#' \dontrun{
#'
#' # Running the function with an exact computation approach performing all the 
#' pairwise genetic comparisons among sequences. We used "raw" model to compute 
#' the genetic distances and considered a median distance as summary statistic 
#' to define the outlier threshold. We used a classic 1.5 interquartile range 
#' distance threshold to consider the outlier sequences.
#' S16_MisAlign = Clean.Seq(inputal = Example_16S_outlier_align, Outlier.Th = 1.5, 
#' output = "Example_16S_outlier_align.rm.out.", Evol.Model = "raw", Sum.Stat = "median", 
#' Iter.optim = FALSE, Re.Align = TRUE, nbthreads = 1, Approxim = FALSE, nb.seq = 25, 
#' nb.rep = 5) 

#' S16_MisAlign[[1]]  # 36 sequences are considered as outlier sequences the model detect

#' # To remove the file created while running the example do the following:
#' file.remove(c( "Example_16S_outlier_align.rm.out..Clean.AllSeq.fas", 
#' "Example_16S_outlier_align.rm.out..Clean.haplo.fas"))
#'

#' Running the function with an approximate computation approach performing all 
#' the pairwise genetic comparisons among sequences. We used "raw" model to compute 
#' the genetic distances and considered a median distance as summary statistic to 
#' define the outlier threshold. We used a classic 1.5 interquartile range distance 
#' threshold to consider the outlier sequences. For the approximate comparison we used 
#' 25 sequences randomly repeated 5 times as suggested by default.
#' S16_MisAlign.approx = Clean.Seq(inputal = Example_16S_outlier_align, Outlier.Th = 1.5, 
#' output = "Example_16S_outlier_align.rm.out.approx", Evol.Model = "raw", Sum.Stat = "median", 
#' Iter.optim = FALSE, Re.Align = TRUE, nbthreads = 5, Approxim = TRUE, nb.seq = 25, nb.rep = 5) 

#' S16_MisAlign.approx[[1]]  # 36 sequences are considered as outlier sequences the model detect

#' setdiff(S16_MisAlign.approx[[1]], S16_MisAlign[[1]]) # excat same finding considering the exact 
#' or the approximate search of outlier sequences.


#' To remove the file created while running the example do the following:
#' file.remove(c( "Example_16S_outlier_align.rm.out.approx.Clean.AllSeq.fas", 
#' "Example_16S_outlier_align.rm.out.approx.Clean.haplo.fas"))
#'
#' }
#'

#' @export Clean.Seq



Clean.Seq = function(inputal = NULL,
                     Outlier.Th = 1.5,
                     output = NULL,
                     Evol.Model = "TN93",
                     Sum.Stat = "median",
                     Iter.optim = TRUE,
                     Re.Align = FALSE,
                     nbthreads = 1,
                     Approxim = FALSE,
                     nb.seq = 25,
                     nb.rep = 5)
{
  # if the class of the inputal object is "alignment" object recognised by ape and seqinr r packages.
  if (class(inputal) == "alignment") {
    cytb2 = inputal
    cytb2a = ape::as.DNAbin(inputal) # Convert in a DNAbin object.
  }
  # if the path to the alignment file in fasta format is provided as inputal
  if (class(inputal) == "character") {
    cytb2 = seqinr::read.alignment(inputal, format = "fasta")
    cytb2a = ape::as.DNAbin(cytb2) # Convert in a DNAbin object.
  }
  
  # if the class of the inputal object is a "DNAbin" object.
  if (class(input) == "DNAbin") {
    if(length(input[[1]]) == 1){
      input = ape::as.alignment(input) # convert to an "alignment" object.
    } else {
      input = DNAbin2alignment(input)
    }
  }
  
  
  # Remove duplicated sequences before to find the outlier sequences
  cytb2.uniq = cytb2
  cytb2.uniq$seq = cytb2$seq[which(base::duplicated(cytb2$seq) == FALSE)]
  cytb2.uniq$nam = cytb2$nam[which(base::duplicated(cytb2$seq) == FALSE)]
  cytb2.uniq$nb = length(cytb2.uniq$seq)
  cytb2Sa = ape::as.DNAbin(cytb2.uniq)
  Seq.name.origin = cytb2.uniq$nam # to keep track of the original list of sequence name, before removing
  
  
  Res.Tot = detect.out.seq(
    input = cytb2Sa,
    Evol.Model = Evol.Model,
    Sum.Stat = Sum.Stat,
    Outlier.Th = Outlier.Th,
    Approxim = Approxim,
    nb.seq = nb.seq,
    nb.rep = nb.rep,
    nbthreads = nbthreads
  )
  
  if (length(Res.Tot[[1]]) > 0) {
    # remove the previously outlier sequences identified
    pos.2.rm = match(Res.Tot[[1]], cytb2.uniq$nam)
    cytb2.uniq$nam = cytb2.uniq$nam[-pos.2.rm]
    cytb2.uniq$seq = cytb2.uniq$seq[-pos.2.rm]
    cytb2.uniq$nb = length(cytb2.uniq$seq)
    # Remove the gap only site
    cytb2.uniq = rm.del.gap(input = cytb2.uniq)
    
    if (Re.Align == TRUE) {
      #Export a temporary file
      if (class(cytb2.uniq) == "alignment") {
        cytb2Sa = ape::as.DNAbin(cytb2.uniq[[1]])  # convert to DNAbin object and export as fasta format.
      } else {
        cytb2Sa = cytb2.uniq[[1]]
      }
      ape::write.dna(
        cytb2Sa,
        file = "Temp.realig.fas",
        format = "fasta",
        nbcol = -1,
        colsep = ""
      )
      
      mafft = "mafft"
      os <- .Platform$OS
      if (os == "windows") {
        if (missing(Mafft.path)) {
          stop("The path to the mafft executable must be provided in Mafft.path")
        }
        mafft = Mafft.path
      }
      a = paste(
        mafft,
        " --retree 1 --maxiterate 0 --thread ",
        nbthreads,
        " Temp.realig.fas > ",
        "Temp.mafft.realign.fas",
        sep = ""
      )
      system(a)
      
      # import the new alignement
      cytb2.uniq = seqinr::read.alignment("Temp.mafft.realign.fas", format = "fasta")
      cytb2Sa = ape::as.DNAbin(cytb2.uniq) # Convert in a DNAbin object.
      
    }
  }
  
  
  if (class(cytb2.uniq) == "DNAbin") {
    cytb2Sa = cytb2.uniq
    if(length(cytb2.uniq[[1]]) == 1){
      cytb2.uniq = ape::as.alignment(cytb2.uniq) # convert to an "alignment" object.
    } else {
      cytb2.uniq = DNAbin2alignment(cytb2.uniq)
    }
  } else {
    cytb2Sa = ape::as.DNAbin(cytb2.uniq)
  }
  
  
  if (Iter.optim == TRUE) {
    if (length(Res.Tot[[1]]) < 1) {
      print(
        "No outlier sequence detected, either Outlier.Th value (default 1.5) could be decrease to increase the chance to get an outlier sequence, and /or the Evol.Model should be changed"
      )
    } else {
      repeat {
        Res.optim = detect.out.seq(
          input = cytb2Sa,
          Evol.Model = Evol.Model,
          Sum.Stat = Sum.Stat,
          Outlier.Th = Outlier.Th,
          Approxim = Approxim,
          nb.seq = nb.seq,
          nb.rep = nb.rep,
          nbthreads = nbthreads
        )
        
        if (length(Res.optim[[1]]) > 0) {
          # remove the previously outlier sequences identified
          pos.2.rm = match(Res.optim[[1]], cytb2.uniq$nam)
          cytb2.uniq$nam = cytb2.uniq$nam[-pos.2.rm]
          cytb2.uniq$seq = cytb2.uniq$seq[-pos.2.rm]
          cytb2.uniq$nb = length(cytb2.uniq$seq)
          # Remove the gap only site
          cytb2.uniq = rm.del.gap(input = cytb2.uniq)
          
          if (Re.Align == TRUE) {
            #Export a temporary file
            if (class(cytb2.uniq[[1]]) == "alignment") {
              cytb2Sa = ape::as.DNAbin(cytb2.uniq[[1]])  # convert to DNAbin object and export as fasta format.
            } else {
              cytb2Sa = cytb2.uniq[[1]]
            }
            ape::write.dna(
              cytb2Sa,
              file = "Temp.realig.fas",
              format = "fasta",
              nbcol = -1,
              colsep = ""
            )
            
            a = paste(
              mafft,
              " --retree 1 --maxiterate 0 --thread ",
              nbthreads,
              " Temp.realig.fas > ",
              "Temp.mafft.realign.fas",
              sep = ""
            )
            system(a)
            
            # import the new alignement
            cytb2.uniq = seqinr::read.alignment("Temp.mafft.realign.fas", format = "fasta")
            cytb2Sa = ape::as.DNAbin(cytb2.uniq) # Convert in a DNAbin object.
            
          }
        }
        
        if (length(Res.optim[[1]]) == 0) {
          break
        }
      }
    }
  }
  
  # Remove the outlier sequence from the global alignment including the duplicated sequences.
  
  Outlier.Seqs = setdiff(Seq.name.origin, cytb2.uniq$nam)
  if (length(Outlier.Seqs) == 0) {
    cytb2.toExport = cytb2
    Outlier.Seqs = NA
  } else {
    cytb2.toExport = cytb2
    cytb2.toExport$seq = cytb2.toExport$seq[-match(Outlier.Seqs, cytb2.toExport$nam)]
    cytb2.toExport$nam = cytb2.toExport$nam[-match(Outlier.Seqs, cytb2.toExport$nam)]
    cytb2.toExport$nb = length(cytb2.toExport$seq)
  }
  cytb2.toExport = ape::as.DNAbin(cytb2.toExport)
  
  # export the final alignment
  if (!is.null(output)) {
    ape::write.dna(
      cytb2Sa,
      file = paste(output, ".Clean.haplo.fas", sep = ""),
      format = "fasta",
      nbcol = -1,
      colsep = ""
    )
    
    ape::write.dna(
      cytb2.toExport,
      file = paste(output, ".Clean.AllSeq.fas", sep = ""),
      format = "fasta",
      nbcol = -1,
      colsep = ""
    )
  }
  
  if (Re.Align == TRUE) {
    # remove the temporary file
    unlink(c("Temp.realig.fas", "Temp.mafft.realign.fas"))
  }
  
  return(
    list(
      Outlier.Sequences = Outlier.Seqs,
      Clean.Haplo.Alignment = cytb2Sa,
      Clean.All.Seq.Alignment = cytb2.toExport,
      First.Dist.Mat = Res.Tot[[2]]
    )
  )
}

#' For debug
#' rm(cytb2Sa, NewAlign.Noindels, NewAlign.Noindels2, Res1, cytb2.uniq, Res.optim, Seq.name.origin, cytb2, cytb2a, cytb2.toExportOutlier.Seqs)


#' @title This function detects outlier sequences and is used internally by the Clean.Seq function.

#' @description This function helps to detect sequences that are abnormally distant from all other sequences in 
#' the alignment. See the Clean.Seq function for a more detailed description.

#' @param input is a DNAbin alignment object.

#' @param Evol.Model evolutionary model used to compute the DNA distance among sequences, 
#' default "TN93" (Tamura and Nei 1993), it allows all the models used in dist.dna function from the ape R package. 

#' @param Sum.Stat type of symmary statistic to compare the DNA distance of the focal sequence with the other, 
#' by defult we use the median distance but it can also be the "min" for the minimal dNA distance between 1 sequence and its closest neighbor.

#' @param Outlier.Th is a number (by default 1.5) defining the distance to consider an outlier sequences. 
#' The Outlier.Th is a multiplier of the interquartile range from the genetic distances distribution using the 
#' Sum.stat (eg. median). This parameter is equivalent of the "range" parameter from the boxplot function from the stat R package.

#' @param Approxim if TRUE (defaults FALSE) the functions uses an approximation approches to estimate the outlier sequences. 
#' If sequence is compare to a random selection of other sequences (defined by the parameter nb.seq) x times (the number 
#' of repetition of randoms election is performed by the parameter  nb.rep). This function can speed up a lot the discovery 
#' of potential oulier sequences for lareg alignment (more than 5000 - 10000 sequences) becasue the function can run 
#' in parallel. However we need to keep in mind that this option provide only an approximation of the outlier 
#' sequences when Approxim  = TRUE.

#' @param nb.seq a number (by default 25) defining the number of sequences than are randomly choosen too detect an 
#' oultier sequences. This option works only if Approxim = TRUE. To increase the precision of the outlier detection 
#' increasing this number help, but it will be slower to run!


#' @param nb.rep a number (by default 5) precising how many time the random selection of sequences will be performed. 
#' This option works only if Approxim = TRUE. 

#' @param nbthreads a number (by default 1) defining the number of threads used to performed the Approximate 
#' detection of outlier sequences 

#' @export detect.out.seq

detect.out.seq = function(input = NULL,
                          Evol.Model = Evol.Model,
                          Sum.Stat = "median",
                          Outlier.Th = 1.5,
                          Approxim = FALSE,
                          nb.seq = 25,
                          nb.rep = 5,
                          nbthreads = 1) {
  if (Approxim == TRUE) {
    # cpu0 = Sys.time()
    Res1 = do.call(rbind,
                   parallel::mclapply(1:dim(input)[1], mc.cores = nbthreads, function(x) {
                     nb.pot = seq(1, dim(input)[1])[-x]
                     repl = replicate(nb.rep, sample(nb.pot, nb.seq))
                     
                     res.ind = unlist(lapply(1:dim(repl)[2], function(i) {
                       input.repl = rbind(input[x, ], input[repl[, i], ])
                       
                       matdis = ape::dist.dna(input.repl,
                                              model = Evol.Model,
                                              pairwise.deletion = TRUE)
                       matdis.mat = as.matrix(matdis)
                       
                       diag(matdis.mat) = NA # remove the diagonal
                       if (Sum.Stat == "median") {
                         # Compute the mediane distance for each sequence in comparison to all the other sequences
                         medist = apply(matdis.mat, 2, median, na.rm = T)
                       }
                       if (Sum.Stat == "min") {
                         # Compute the mediane distance for each sequence in comparison to all the other sequences
                         medist = apply(matdis.mat, 2, min, na.rm = T)
                       }
                       medist[1]
                     }))
                     c(labels(input)[x], mean(res.ind))
                   }))
    # cpu1 = Sys.time()
    # cpu1 - cpu0
    
    Res1 = as.data.frame(Res1)
    Res1[, 2] = as.numeric(Res1[, 2])
    medist.DF = Res1[order(Res1[, 2], decreasing = T), ]
    
  } else {
    matdis = ape::dist.dna(input, model = Evol.Model, pairwise.deletion = TRUE)
    matdis.mat = as.matrix(matdis)
    diag(matdis.mat) = NA # remove the diagonal
    if (Sum.Stat == "median") {
      # Compute the mediane distance for each sequence in comparison to all the other sequences
      medist = apply(matdis.mat, 2, median, na.rm = T)
    }
    if (Sum.Stat == "min") {
      # Compute the mediane distance for each sequence in comparison to all the other sequences
      medist = apply(matdis.mat, 2, min, na.rm = T)
    }
    medist.DF = data.frame(names(sort(medist, decreasing = T)), sort(medist, decreasing = T))
  }
  
  
  # define the threshold distance
  TH.out = quantile(medist.DF[, 2], 0.75) + Outlier.Th * IQR(medist.DF[, 2])
  
  # identify the outlier sequence.
  Res2 = medist.DF[which(medist.DF[, 2] > TH.out), 1]
  
  
  return(list(Outlier.seq = Res2, All.Distance.DF = medist.DF))
}













#' @title Function to remove sequences (or loci) with rare insertion, or containing too many ambiguous nucleotides, 
#' or sequences that are too short.

#' @description This function helps to detect and remove sequences and loci with rare insertions, or containing too many ambiguous nucleotides, or sequences that are too short. 

#' @param input an object of class "alignment" (seqinr R package),
#' or "DNAbin" (ape R package) or the name (including the path if necessary) of the input fasta file.

#' @param Remove.del.gaps if TRUE (by default) it removes the deletion gap only site.

#' @param Remove.Seq.Indels if TRUE (default FALSE) removes the sequence with a rare insertion creating 
#' an insertion gap for most of the other sequences.

#' @param Nb.Seq.Indels is an integer precising the maximal number of sequences that could include 
#' a rare insertion blocks creating an insertion gap for most of the other sequences (valide only if Remove.Seq.Indels = TRUE).

#' @param Min.Perc.Seq.With.Info minimal percentage (on the scale between 0 and 1) of sequence 
#' with a documented nucleotide for a loci preceding or following the indel blocks, when considering an
#' insertion as a potential sequencing error and causing and indel blocks for the other sequence.
#' This value by default is set up to 0.1, to avoid the few longer sequences to be considered as 
#' problematic sequences causing indels.

#' @param Remove.Indels.Only if TRUE (default FALSE) the indels blocks detected are removed form the 
#' alignment but the sequence is maintained in the alignment. This option works only when the option Remove.Seq.Indels is TRUE.

#' @param Remove.Short.Seq if TRUE (default FALSE) it removes sequences that are too short below a threshold in percentage 
#' of the of the number of missing nucleotides of the longuest sequence (default 15% of missing nucleotides). 
#' So sequences with more than 50% of missing nucleotides are removed.

#' @param Max.Per.missing a number precising the maximal percentage (by default 50%) of missing nucleotides allow in a sequence 
#' to be retained. If the sequence contains more missing data than the Max.Per.missing threshold, then the sequences is 
#' removed from the alignment. This functions works only if Remove.Short.Seq is set to TRUE.

#' @param Remove.Seq.TooManyAmbig if TRUE (default FALSE) then the function removes sequences with too many ambiguous 
#' nucleotides (N,R, Y, S, W, K, M, B, D, H, V and N nucleotide code) from the IUPAC nucleotide code.

#'  @param Percent.Ambig a number precising the maximal percentage (by default 50%) of ambiguous nucleotides per 
#'  sequence (by default 30%). This option works only if Percent.Ambig is set to TRUE.

#'  @param Remove.Loc.Low.Freq if TRUE (default FALSE), this function remove the loci of the alignment with less than 
#'  X percent of sequences with some informtaion, the percentage of information is provide by the option 
#'  Minimal.Locus.Freq (by default 30%). This option is used typically to trim the beginning and the end of the alignment 
#'  to avoid few extra long sequences.


#'  @param Minimal.Locus.Freq the percentage of minimal sequences with some nucleotides used to retained a locus in 
#'  the alignment. By default at least the locus should have 30% of sequences with some nucleotides to be retained. 
#'  This option works only if Remove.Loc.Low.Freq is TRUE.

#' @export rm.del.gap

rm.del.gap = function(input = NULL,
                      Remove.del.gaps = TRUE,
                      Remove.Seq.Indels = FALSE,
                      Nb.Seq.Indels = 2,
                      Min.Perc.Seq.With.Info = 0.1,
                      Remove.Indels.Only = FALSE,
                      Remove.Short.Seq = FALSE,
                      Max.Per.missing = 50,
                      Remove.Seq.TooManyAmbig = FALSE,
                      Percent.Ambig = 30,
                      Remove.Loc.Low.Freq = FALSE,
                      Minimal.Locus.Freq = 30) {
  # if the class of the inputal object is a "DNAbin" object.
  if (class(input) == "DNAbin") {
    if(length(input[[1]]) == 1){
    input = ape::as.alignment(input) # convert to an "alignment" object.
    cytb2Smat = seqinr::as.matrix.alignment(input)  # convert in a matrix
  } else {
    input = DNAbin2alignment(input)
    cytb2Smat = alignment2matrix(input)
  }
  } else {
    cytb2Smat = seqinr::as.matrix.alignment(input)  # convert in a matrix
  }
  
  
  Seq.names.original = row.names(cytb2Smat)
  
  if (Remove.del.gaps == "TRUE") {
    # remove the deletion gap site only.
    nbindels = apply(as.matrix(seq(1, dim(cytb2Smat)[2], 1)), 1, function(x) {
      sum(stringr::str_count(cytb2Smat[, x], "-"))
    })
    
    a = which(nbindels == dim(cytb2Smat)[1])
    
    if (length(a) > 0) {
      cytb2Smat = cytb2Smat[,-a]
    }
  }
  
  
  if (Remove.Seq.Indels == "TRUE") {
    nbindels = apply(as.matrix(seq(1, dim(cytb2Smat)[2], 1)), 1, function(x) {
      sum(stringr::str_count(cytb2Smat[, x], "-"))
    })
    a = which(nbindels >= (dim(cytb2Smat)[1] - Nb.Seq.Indels))
    
    
    if (length(a) > 0) {
      aa = seqle(a)
      Starting = aa$values
      Ending = aa$values + (aa$length - 1)
      
      # For each indels blocks
      Seq.to.remove.Indel = do.call(rbind, lapply(1:length(Starting), function(x) {
        if (Starting[x] == 1 | Ending[x] == dim(cytb2Smat)[2]) {
          resa = 0
        } else {
          Test.preceding = sum(stringr::str_count(cytb2Smat[, (Starting[x] - 1)], "-")) # test if the preceding loci of the starting position of the insertion include enougth informed loci.
          Test.Ending = sum(stringr::str_count(cytb2Smat[, (Ending[x] + 1)], "-")) # test if the following position of the ending position of the insertion include enougth informed loci.
          nb.seq = dim(cytb2Smat)[1]
          if (Test.preceding < nb.seq * (1-Min.Perc.Seq.With.Info) & Test.Ending < nb.seq * (1-Min.Perc.Seq.With.Info)) {
            resa = 1
          } else {
            resa = 0
          }
        }
        
        # detect the sequence causing the true indels
        if (resa == 1) {
          seq.name.indel = which(!cytb2Smat[, Starting[x]] == "-")
          data.frame(
            Seq.name = seq.name.indel,
            Starting.Indel = rep(Starting[x], length(seq.name.indel)),
            Ending.Indel = rep(Ending[x], length(seq.name.indel))
          )
        }
      }))
      
      
      if (dim(Seq.to.remove.Indel)[1] > 0) {
        if (Remove.Indels.Only == TRUE) {
          bb = unique(Seq.to.remove.Indel[, c(2, 3)])
          Pos.ToRemove = unlist(lapply(1:dim(bb)[1], function(x) {
            seq(bb[x, 1], bb[x, 2])
          }))
          Seq.to.remove.indels = row.names(cytb2Smat)[unique(Seq.to.remove.Indel[, 1])]
          cytb2Smat = cytb2Smat[,-Pos.ToRemove]
          
        } else {
          Seq.to.remove.indels = row.names(cytb2Smat)[unique(Seq.to.remove.Indel[, 1])]
          cytb2Smat = cytb2Smat[-unique(Seq.to.remove.Indel[, 1]), ]
          
        }
        
        # Remove the deletion gap site only, create by the suppression of the sequence with insertion.
        nbindels = apply(as.matrix(seq(1, dim(cytb2Smat)[2], 1)), 1, function(x) {
          sum(stringr::str_count(cytb2Smat[, x], "-"))
        })
        
        aa = which(nbindels == dim(cytb2Smat)[1])
        
        if (length(aa) > 0) {
          cytb2Smat = cytb2Smat[,-aa]
        }
        
        
      } else {
        Seq.to.remove.indels = NA
        
      }
      
      
    } else {
      Seq.to.remove.indels = NA
    }
    
  } else {
    Seq.to.remove.indels = NA
  }
  
  
  
  # Remove sequence that are too short
  if (Remove.Short.Seq == "TRUE") {
    nbmissing = apply(as.matrix(seq(1, dim(cytb2Smat)[1], 1)), 1, function(x) {
      sum(stringr::str_count(cytb2Smat[x,], "-"))
    })
    
    Seq.To.remove = which(((nbmissing / dim(cytb2Smat)[2]) * 100) >= Max.Per.missing)
    
    if (length(Seq.To.remove) > 0) {
      Seq.To.remove.short = row.names(cytb2Smat)[Seq.To.remove]
      cytb2Smat = cytb2Smat[-Seq.To.remove,]
    } else {
      Seq.To.remove.short = NA
    }
    
  } else {
    Seq.To.remove.short = NA
  }
  
  
  # Remove sequence with too many ambiguities (N,R, Y, S, W, K, M, B, D, H, V and N nucleotide code) 
  # from the IUPAC nucleotide code.
  if (Remove.Seq.TooManyAmbig == "TRUE") {
    nb.loc = dim(cytb2Smat)[2]
    
    nbAmbig = unlist(lapply(1:dim(cytb2Smat)[1], function(x) {
      nb.loc.noindel = nb.loc - sum(stringr::str_count(cytb2Smat[x,], "-"))
      ((nb.loc.noindel - (
        sum(stringr::str_count(cytb2Smat[x,], "a")) + sum(stringr::str_count(cytb2Smat[x,], "c")) + sum(stringr::str_count(cytb2Smat[x,], "t")) + sum(stringr::str_count(cytb2Smat[x,], "g")) +  sum(stringr::str_count(cytb2Smat[x,], "u"))
      )) / nb.loc.noindel) * 100
    }))
    
    Seq.To.remove.Too.Ambig = which(nbAmbig > Percent.Ambig)
    
    if (length(Seq.To.remove.Too.Ambig) > 0) {
      Seq.To.remove.amb = row.names(cytb2Smat)[Seq.To.remove.Too.Ambig]
      cytb2Smat = cytb2Smat[-Seq.To.remove.Too.Ambig,]
    } else {
      Seq.To.remove.amb = NA
    }
    
  } else {
    Seq.To.remove.amb = NA
  }
  
  
  if (Remove.Loc.Low.Freq == TRUE) {
    # remove the locus in low frequency 
    nbindels = apply(as.matrix(seq(1, dim(cytb2Smat)[2], 1)), 1, function(x) {
      sum(stringr::str_count(cytb2Smat[, x], "-"))
    })
    
    a = which(nbindels >= round(dim(cytb2Smat)[1] / 100 * Minimal.Locus.Freq, digits = 0))
    
    if (length(a) > 0) {
      cytb2Smat = cytb2Smat[,-a]
    }
    
  }
  
  
  
  # convert the matrix in string for all the sequences
  NewAlign.Noindels2 = vector()
  i = 1
  for (i in 1:dim(cytb2Smat)[1]) {
    NewAlign.Noindels2 = c(NewAlign.Noindels2, seqinr::c2s(cytb2Smat[i,]))
  }
  output = list()
  output$nb = dim(cytb2Smat)[1]
  output$nam = row.names(cytb2Smat)
  output$seq = NewAlign.Noindels2
  output$com = NA
  class(output) = class(input)
  #input$seq = NewAlign.Noindels2
  
  # Convert to DNAbin object and export as fasta format.
  cytb2Sa = ape::as.DNAbin(output)
  
  Seq.names.final = row.names(cytb2Smat)
  
  Seq.removed = setdiff(Seq.names.original, Seq.names.final)
  
  Seq.removed.DF = data.frame(Sequence.names.removed = Seq.removed,
                              Removal.type = rep(NA, length(Seq.removed)))
  
  if (Remove.Seq.Indels == "TRUE" &
      is.na(Seq.to.remove.indels[1]) == FALSE &
      Remove.Indels.Only == FALSE) {
    if (length(Seq.to.remove.indels) > 0) {
      Seq.removed.DF[na.omit(match(Seq.to.remove.indels, Seq.removed.DF[, 1])), 2] = "Indels"
    }
  }
  
  if (Remove.Seq.Indels == "TRUE" &
      is.na(Seq.to.remove.indels[1]) == FALSE &
      Remove.Indels.Only == TRUE) {
    if (length(Seq.to.remove.indels) > 0) {
      Seq.removed.DF.add = data.frame(
        Sequence.names.removed = Seq.to.remove.indels,
        Removal.type =  rep("Indels", length(Seq.to.remove.indels))
      )
      
      Seq.removed.DF = data.frame(rbind(Seq.removed.DF, Seq.removed.DF.add))
      
    }
  }
  
  
  
  if (Remove.Short.Seq == "TRUE" &
      is.na(Seq.To.remove.short[1]) == FALSE) {
    if (length(Seq.To.remove.short) > 0) {
      Seq.removed.DF[na.omit(match(Seq.To.remove.short, Seq.removed.DF[, 1])), 2] = "Too.short"
    }
  }
  
  if (Remove.Seq.TooManyAmbig == "TRUE" &
      is.na(Seq.To.remove.amb[1]) == FALSE) {
    if (length(Seq.To.remove.amb) > 0) {
      Seq.removed.DF[na.omit(match(Seq.To.remove.amb, Seq.removed.DF[, 1])), 2] = "Too.Ambiguous"
    }
  }
  
  return(list(Alignment = cytb2Sa, Sequences.Removed = Seq.removed.DF))
}


# for Debug
# rm(Remove.Seq.TooManyAmbig, Remove.Short.Seq, Remove.Seq.Indels, Seq.removed, Seq.removed.DF, output, NewAlign.Noindels2, Seq.names.final, Seq.removed, Seq.To.remove.short, Seq.to.remove.indels, Seq.To.remove.amb, Seq.To.remove, nbAmbig,nb.loc, Seq.To.remove.Too.Ambig, nbmissing, cytb2Smat, bb, Pos.ToRemove, resa, Ending, Starting, seq.name.indel, Test.Ending, Test.preceding, Seq.to.remove.Indel, aa, nbindels)



#' @title Detect automatically the reading frame and sequences with the presence of stop codon 
#' in a coding DNA sequence alignment.
#' 
#' @description This function automatically detects the reading frame of a DNA sequences coding for a protein, 
#' and then identifies the sequences with the presence of stop codons. 

#' @param input a "DNAbin" or an "alignment" object containing the multiple DNA sequences already aligned.
#' 
#' @param genet.code a number precising the genetic code of the DNA sequence (by default 5 for invertebrate mitochondrial). See the help page of the function
#' \emph{translate} from the \emph{seqinr} R package available at https://cran.r-project.org/web/packages/seqinr/seqinr.pdf 
#' eg: 1 is standard, 2 is vertebrate.mitochondrial, 3 is yeast.mitochondrial, 4 protozoan.mitochondrial+mycoplasm, 
#' 5 invertebrate.mitochondrial etc...

#' @return The functions returns a data.frame with the sequence names including stop codons and the number of stops codon detected.
#' 
#' @export Detect.stop.codon

Detect.stop.codon = function(input = NULL, genet.code = 5){
  
  if (class(input) == "DNAbin") {
    if(length(input[[1]]) == 1){
      input = ape::as.alignment(input) # convert to an "alignment" object.
    } else {
      input = DNAbin2alignment(input)
    }
  }
  
  # define the codonstart
  frame.vec = c(0,1,2)
  ReadingFrame = do.call(cbind, lapply(1:3, function(i){
    unlist(lapply(1:length(input$seq), function(x){
      aa.trans = paste(seqinr::translate(seq = seqinr::s2c(input$seq[[x]]), frame = frame.vec[i], numcode = genet.code), collapse = "")
      stringr::str_count(aa.trans, stringr::fixed("*"))
    }))
  }))
  # interesting the function trans from the ape R package provide slightly different results that are not in agreement with seaview software for instance, while the translate function from the seqinr R package provide the appropriate behavior!
  
  pos.codonstart = which.min(apply(ReadingFrame, 2,sum))
  Seq.With.stop.codon = which(ReadingFrame[,pos.codonstart] > 0)
  Nb.stop.codon = ReadingFrame[which(ReadingFrame[,pos.codonstart] > 0), pos.codonstart]
  data.frame(Seq.with.stop.codon = input$nam[Seq.With.stop.codon], Nb.stop.codon = Nb.stop.codon)
}




#' @title Detect suites of consecutive numbers in a vector
#' @description a function to detect suites of consecutive numbers in a vector from 
#' https://stackoverflow.com/questions/8400901/group-integer-vector-into-consecutive-runs/8402950#8402950
#' This function is used internally by the \emph{rm.del.gap} function.
#' @param x the vector of interest with the suite of numbers.

#' @export seqle


seqle <- function(x,incr=1) { 
  if(!is.numeric(x)) x <- as.numeric(x) 
  n <- length(x)  
  y <- x[-1L] != x[-n] + incr 
  i <- c(which(y|is.na(y)),n) 
  list(lengths = diff(c(0L,i)),
       values = x[head(c(0L,i)+1L,-1L)]) 
} 






#' @title Convert a DNAbin object (ape format) to an alignment object (seqinr format)
#' 
#' @description This function converts a DNAbin object (ape format) loaded with the ape::read.FASTA 
#' to an alignment object (seqinr format). It avoids the conversion from character 
#' to numeric in the DNAbin object.
#' @param input a DNAbin object
#'
#' @export DNAbin2alignment

DNAbin2alignment = function(input){
  nb = length(input)
  seq = as.character(input)
  nam = names(input)
  com = NA
  out.al = list(nb = nb, seq = seq, nam = nam, com = com)
  class(out.al) = "alignment"
  return(out.al)
}

#' @title Convert an alignment object (seqinr format) to a matrix

#' @description This function converts an alignment object (seqinr format) to a matrix, when the initial 
#' alignment was converted from a DNAbin object loaded with the ape::read.FASTA function.
#' It solves some recent probelm I encounter with the as.matrix.alignment function 
#' from the seqinr package
#' @param input an alignment object (seqinr format)
#' 
#' 

alignment2matrix = function(input){
  out.mat = do.call(rbind, lapply(input$seq, function(x){
    x
  }))
  row.names(out.mat) = input$nam
  return(out.mat)
}


