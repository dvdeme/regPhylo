#' @title Include CladeAge blocks into an xml file formatted by BEAUTi
#' to calibrate the tree in an absolute time frame using BEAST2

#' @description This function edits an .xml file following CladeAge requirements to
#' time calibrate the tree, based on the age estimate of the first fossil occurrence for
#' the monophyletic clades retained for the calibration.
#' The CladeAge approach allows the user to objectively determine the shape of the prior distribution
#' of the calibration points used to calibrate the tree in an absolute time frame.
#' @details CladeAge is based on the oldest fossil occurrence of the
#' clades used for the calibration, the net diversification rate, diversification
#' turnover and the fossil sampling rates of the group under investigation.  For a
#' tutorial about CladeAge analysis see the original paper Matschiner et al. 2017
#' and the Tutorial "A Rough guide to CladeAge"
#' (available at: https://www.beast2.org/tutorials/).
#'
#' @return The output .xml file is ready for analysis in BEAST2.

#' @param xml.input an xml input file including topological constraints set-up with
#' BEAUTi and/or edited using the function MultiTopoCont.EditXML4BEAST2.
#' @param CalPointTable a data.frame with 6 columns: the first column reports the name
#' of the clade used for CladeAge, the second reports the taxonomic reference for
#' the clade, the third column reports the support for the clade, the fourth
#' column reports the justification of the first fossil occurrence for the clade,
#' the fifth and sixth columns report the minimum and maximum age of the fossil,
#' respectively. If only one age estimate is available for the fossil then the
#' same age must be reported in both columns. (Only the columns 1, 5 and 6 are
#' used by the function).
#' @param MinDivRate minimum net diversification rate of the clade.
#' @param MaxDivRate maximum net diversification rate of the clade.
#' @param MinTurnoverRate minimum turnover diversification rate of the clade.
#' @param MaxTurnoverRate maximum turnover diversification rate of the clade.
#' @param MinSamplingRate minimum sampling rate of fossils for the clade.
#' @param MaxSamplingRate maximum sampling rate of fossils for the clade.
#' @param inputTaxono taxonomic table with species names in the first column
#' followed by different levels of the classification (e.g.  genus, family,
#' superfamily, order, class, Division....). If taxa without DNA sequences are
#' also present in the analysis then the 'output.new.TaxoTable' table exported
#' by the \code{\link{MultiTopoConst.EditXML4BEAST2}} function must be provided instead.
#' @param input.tree the same rooted tree used by the
#' \code{\link{MultiTopoConst.EditXML4BEAST2}} function.
#' If taxa without DNA sequences are present in the analysis then the 'output.new.tree' tree
#' exported by the \code{\link{MultiTopoConst.EditXML4BEAST2}} function must be provided
#' instead.
#' @param xmltreename name of the xml tree in the input xml file.
#' @param output name of the xml file output ready for BEAST2 analysis
#' (with the path, if necessary).
#' @param Partitions If "TRUE" then the alignment included in the xml files contains
#' more than 1 partition.


#' @examples # To run the example, copy the input files
#' # provided by the package to a temporary directory created in the
#' # current working directory.
#' \dontrun{
#' src.dir = system.file("extdata/TopoConstraints", package = "regPhylo")
#' dir.create("TempDir.CladeAge")
#' # Set up the path of the TempDir folder.
#' dest.dir = paste(getwd(), "/TempDir.CladeAge", sep="")
#' file.names <- dir(src.dir)
#' # Copy all the files stored in regPhylo/extdata/TopoConstraints"
#' # into a temporary folder.
#' sapply(file.names, function(x) {
#' file.copy(from = paste(src.dir, x, sep = "/"),
#' to = paste(dest.dir, x, sep = "/"),
#' overwrite = FALSE) })
#'
#' # We include 4 calibration constraints based on 4 clades (Elopomorpha,
#' # Anguilliformes, Stomiati, Perciformes)
#' # Import the table (i.e. "CalPointTable") listing the 4 clades constrained
#' # with the occurrence of the first fossil.
#' CalibrationTable4clades = read.delim("TempDir.CladeAge/Calib_CA_Fossil_4cl.csv",
#' sep="\t", header = TRUE)
#'
#'
#' #### Example restricted to taxa with at least a DNA sequence in the supermatrix.
#'
#' # Load the classification table (the same table used for the
#' # ConstraintTaxo2newick function), one of two ways:
#' # either through the .Rdata;
#' data(TopoConstraints) # the second object of the list is the classification table
#' dim(TopoConstraints[[2]]) # 16 by 23.
#' # Or the classification table that has been loaded into the temporary directory,
#' # and can be loaded into the R environment by doing the following:
#' ClassifDF = read.csv("TempDir.CladeAge/Classif16sp.csv", header = TRUE)
#' dim(ClassifDF) # 16 by 23
#'
#' # Load the re-rooted tree (the same tree as for the
#' # ConstraintTaxo2newick function) in R (note that a rooted tree is available in the
#' # package and has been loaded in the temporary directory).
#' require(ape)
#' TreeRooted = read.nexus("TempDir.CladeAge/RAxML_bipartitions.Concat_7GTR_Allconst_autoMRE_ReRooted")
#'
#' # All diversification/turnover/sampling rates are from Matschiner et al. 2017.
#' CladeAgeCalib.xml(xml.input = "TempDir.CladeAge/SimpleXml_Wcont.xml", input.tree = TreeRooted,
#' output="TempDir.CladeAge/SimpleXml_ReadyForBEAST.xml", CalPointTable = CalibrationTable4clades,
#' MinDivRate = 0.041, MaxDivRate = 0.081, MinTurnoverRate = 0.0011,
#' MaxTurnoverRate = 0.37, MinSamplingRate = 0.0066, MaxSamplingRate = 0.01806,
#' xmltreename = "Subset1", inputTaxono = ClassifDF, Partitions = "TRUE")
#'
#'
#'
#' #### Example including taxa without DNA sequences in the supermatrix.
#'
#' # Load the new classification table including the two additional taxa without DNA
#' # exported by the function MultiTopoConst.EditXML4BEAST2.
#' NewClassifDF = read.delim("TempDir.CladeAge/Classif18sp_2NoDNA.csv", sep = "\t", header = TRUE)
#'
#' # Load the new rooted "RAxML" tree including the two additional taxa and also the
#' # bootstrap values for each node exported by the function MultiTopoConst.EditXML4BEAST2
#' require(ape)
#' NewTree = read.tree("TempDir.CladeAge/BackboneTreeAll_2spNoDNA.txt")
#'
#' # We then load the calibration table.
#' CalibrationTable4clades = read.delim("TempDir.CladeAge/Calib_CA_Fossil_4cl.csv",
#' sep="\t", header = TRUE)
#'
#' # Run the function with all other setting and options unchanged.
#' CladeAgeCalib.xml(xml.input = "TempDir.CladeAge/SimpleXml_2SpNoDNA_Wcont.xml",
#' input.tree = NewTree,
#' output = "TempDir.CladeAge/SimpleXml_2SpNoDNA_ReadyForBEAST.xml",
#' CalPointTable = CalibrationTable4clades,
#' MinDivRate = 0.041, MaxDivRate = 0.081, MinTurnoverRate = 0.0011,
#' MaxTurnoverRate = 0.37, MinSamplingRate = 0.0066, MaxSamplingRate = 0.01806,
#' xmltreename = "Subset1", inputTaxono = NewClassifDF, Partitions = "TRUE")
#'
#' # To remove the files created while running the example do the following:
#' unlink("TempDir.CladeAge", recursive = TRUE)
#'
#' }
#'
#' @references Matschiner et al. 2017, DOI:10.1093/sysbio/syw076
#'
#' @export CladeAgeCalib.xml


CladeAgeCalib.xml = function(xml.input = NULL, input.tree = NULL, output = NULL,
    CalPointTable = NULL, MinDivRate = NULL, MaxDivRate = NULL, MinTurnoverRate = NULL,
    MaxTurnoverRate = NULL, MinSamplingRate = NULL, MaxSamplingRate = NULL, xmltreename = NULL,
    inputTaxono = NULL, Partitions = NULL) {

    # Load the xml.input file.
    xml.input = readLines(xml.input)

    # If multiple partitions, then set the xmltreename up using the following
    # approach.
    if (Partitions == "TRUE") {
        ggg = xml.input[grep("<tree id=\"Tree.t:", xml.input, perl = T)]
        xmltreename = substr(ggg, (regexpr("Tree.t:", ggg, fixed = T) + 7), (regexpr("\" name=",
            ggg, fixed = T) - 1))
    }

    # Setup the initial popsize of the tree to 100 (this helps the starting tree and
    # avoid Likelihood = Infinity).
    ggg1 = xml.input[grep("name=\"popSize\">", xml.input, perl = T)]
    xml.input[grep("name=\"popSize\">", xml.input, perl = T)] = paste(substr(ggg1,
        1, regexpr(">1.0<", ggg1, fixed = T)), "100</parameter>", sep = "")

    # Conversion table between the clade name already present in the xml.file and the
    # starting and ending row of the xml file to ease the extraction of the clades of
    # interest later on.
    StartClade = grep("<distribution id=\"Node", xml.input)
    e = grep("</taxonset>", xml.input)
    EndClade = c(StartClade[-1] - 1, e[length(e)] + 1)
    CladeNumber = vector()
    i = 1
    for (i in 1:length(StartClade)) {
        CladeNumber = c(CladeNumber, substring(xml.input[StartClade[i]], regexpr("id=\"Node",
            xml.input[StartClade[i]]) + 4, regexpr(".prior\" spec=", xml.input[StartClade[i]],
            fixed = T) - 1))
    }
    # Store information about the clade of the xml file, name of the clade, number of
    # the clade, and row starting the taxonset and row ending the taxonset.
    DFCladexml = cbind(CladeNumber, gsub("Node", "", CladeNumber), StartClade, EndClade)
    colnames(DFCladexml) = c("CladeName", "CladeNumber", "StartClade", "EndClade")


    # Start to feed a temporary file with the CLADEAGE xml block.
    Temp2 = "TempCladeAge.xml"
    cat(file = Temp2, "         <distribution id=\"fossilCalibrations\" spec=\"util.CompoundDistribution\">",
        "\n", sep = "")

    # Estimate the number of 'RealParameter' that are already implemented in the xml
    # file.
    RealParameter = xml.input[grep("RealParameter.", xml.input)]
    CountRealParameter = vector()
    k = 1
    for (k in 1:length(RealParameter)) {
        t1 = 0
        t2 = 0
        if (length(grep(" estimate=\"", RealParameter[k], fixed = T)) == 1) {
            t1 = 1
        }
        if (length(grep(" lower=\"", RealParameter[k], fixed = T)) == 1) {
            t2 = 2
        }
        ttt = t1 + t2
        if (ttt == 1) {
            CountRealParameter = c(CountRealParameter, substring(RealParameter[k],
                regexpr("RealParameter.", RealParameter[k]) + 14, regexpr(" estimate=\"",
                  RealParameter[k]) - 2))
        }
        if (ttt == 2) {
            CountRealParameter = c(CountRealParameter, substring(RealParameter[k],
                regexpr("RealParameter.", RealParameter[k]) + 14, regexpr(" lower=\"",
                  RealParameter[k]) - 2))
        }
        if (ttt == 3) {
            CountRealParameter = c(CountRealParameter, substring(RealParameter[k],
                regexpr("RealParameter.", RealParameter[k]) + 14, regexpr(" estimate=\"",
                  RealParameter[k]) - 2))
        }
    }
    MaxCountRealParameter = max(as.numeric(as.character(CountRealParameter)))

    # A matrix to prepare the number of the realParameter in the xml file.
    RealParameterMat = matrix(MaxCountRealParameter + seq(1, 4 * dim(CalPointTable)[1],
        by = 1), ncol = 4, nrow = dim(CalPointTable)[1], byrow = "TRUE")

    # Remove the species names from the taxonomic table.
    inputTaxob = inputTaxono[, -1]

    # Implement the CLADEAGE constraints in the xml file.
    rowsToRemoved = vector()
    logFossilPriorName = vector()

    i = 1
    for (i in 1:dim(CalPointTable)[1]) {
        # loop over each clade
        target = as.character(CalPointTable[i, 1])

        aa = vector()
        j = 1
        for (j in 1:dim(inputTaxob)[2]) {
            # Loop over all the hierarchical levels of the Linnaean classification. aa=c(aa,
            aa = c(aa, which(inputTaxob[, j] == target))
        }  # End for j.


        if (i == 1) {
            idref = "id="
        } else {
            idref = "idref="
        }

        Tipnames = vector()
        if (length(aa) > 1) {
            # First, set-up the clades with at least two taxa.
            Tipnames = c(Tipnames, as.character(inputTaxono[aa, 1]))
            # Test if the tips of the clade defined in CalPointTable defined a monophyletic
            # group in our phylogenetic tree.
            DescenTip = stats::na.omit(input.tree$tip.label[phytools::getDescendants(input.tree, ape::getMRCA(input.tree,
                Tipnames))])
            if (setequal(DescenTip, Tipnames))
                {
                  # If the group is also monophyletic then implement the CLADEAGE xml code,
                  # enforcing the monophyly of the group, because monophyly of the group is
                  # supported by previous studies.  If this condition is FALSE, that means, our
                  # tree is discordant with the clade delineated in the CalPointTable and our
                  # taxonomic table which is based on DeepFinV4 (Betancur et al. 2017, DOI:
                  # 10.1186/s12862-017-0958-3, and Roberts et al. 2017 available at :
                  # https://www.tepapa.govt.nz/sites/default/files/checklist_of_the_fishes_of_new_zealand_v_1_0_july_2017.pdf).
                  # In this case, we cannot use the constraint, because the nature of the clade has
                  # changed, and potentially the occurrence of the first fossil as well.

                  # Test the presence of the clade in the existing xml file.
                  NumberofTheClade = ape::getMRCA(input.tree, Tipnames) - length(input.tree$tip.label)
                  if (is.na(match(NumberofTheClade, DFCladexml[, 2])) == "FALSE") {
                    # Test if the clade is already constrained to be monophyletic in the xml file.
                    ww = DFCladexml[match(NumberofTheClade, DFCladexml[, 2]), ]
                    rowsToRemoved = c(rowsToRemoved, paste(as.character(ww[3]), ":",
                      as.character(ww[4]), sep = ""))
                    # Prepare the CLADEAGE xml block.
                    logFossilPriorName = c(logFossilPriorName, paste("\"", as.character(paste(target,
                      ww[1], sep = "_")), ".fossilprior\"", sep = ""))
                    cat(file = Temp2, append = T, "             <distribution id=\"",
                      as.character(paste(target, ww[1], sep = "_")), ".fossilprior\" spec=\"beast.math.distributions.FossilPrior\" monophyletic=\"true\" tree=\"@Tree.t:",
                      xmltreename, "\">", "\n", sep = "")
                    cat(file = Temp2, append = T, "                 <fossilDistr id=\"FossilCalibration.",
                      i - 1, "\" spec=\"beast.math.distributions.FossilCalibration\">",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                      RealParameterMat[i, 1], "\" name=\"minOccuranceAge\">", as.character(CalPointTable[i,
                        5]), "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                      RealParameterMat[i, 2], "\" name=\"maxOccuranceAge\">", as.character(CalPointTable[i,
                        6]), "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"minDivRate\" name=\"minDivRate\">", MinDivRate, "</parameter>",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"maxDivRate\" name=\"maxDivRate\">", MaxDivRate, "</parameter>",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"minTurnoverRate\" name=\"minTurnoverRate\">", MinTurnoverRate,
                      "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"maxTurnoverRate\" name=\"maxTurnoverRate\">", MaxTurnoverRate,
                      "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"minSamplingRate\" name=\"minSamplingRate\">", MinSamplingRate,
                      "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"maxSamplingRate\" name=\"maxSamplingRate\">", MaxSamplingRate,
                      "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                      RealParameterMat[i, 3], "\" name=\"minSamplingGap\">0.0</parameter>",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                      RealParameterMat[i, 4], "\" name=\"maxSamplingGap\">0.0</parameter>",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                 </fossilDistr>",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                 <taxonset id=\"",
                      as.character(paste(target, ww[1], sep = "_")), "\" spec=\"TaxonSet\">",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, xml.input[(as.numeric(as.character(ww[3])) +
                      2):as.numeric(as.character(ww[4]))], sep = "\n")

                  } else {
                    # Monophyletic clade retrieves in the RAxML tree but not constrained (meaning
                    # that the bootstrap is below the threshold specified by the function
                    # MultiTopoCont.EditXML4BEAST2.R). In this case, we enforce also the monophyly of
                    # this group because there is a good support for this group being monophyletic in
                    # previous studies.

                    # Prepare the CLADEAGE xml block.
                    logFossilPriorName = c(logFossilPriorName, paste("\"", as.character(target),
                      ".fossilprior\"", sep = ""))
                    cat(file = Temp2, append = T, "             <distribution id=\"",
                      as.character(target), ".fossilprior\" spec=\"beast.math.distributions.FossilPrior\" monophyletic=\"true\" tree=\"@Tree.t:",
                      xmltreename, "\">", "\n", sep = "")
                    cat(file = Temp2, append = T, "                 <fossilDistr id=\"FossilCalibration.",
                      i - 1, "\" spec=\"beast.math.distributions.FossilCalibration\">",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                      RealParameterMat[i, 1], "\" name=\"minOccuranceAge\">", as.character(CalPointTable[i,
                        5]), "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                      RealParameterMat[i, 2], "\" name=\"maxOccuranceAge\">", as.character(CalPointTable[i,
                        6]), "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"minDivRate\" name=\"minDivRate\">", MinDivRate, "</parameter>",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"maxDivRate\" name=\"maxDivRate\">", MaxDivRate, "</parameter>",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"minTurnoverRate\" name=\"minTurnoverRate\">", MinTurnoverRate,
                      "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"maxTurnoverRate\" name=\"maxTurnoverRate\">", MaxTurnoverRate,
                      "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"minSamplingRate\" name=\"minSamplingRate\">", MinSamplingRate,
                      "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter ",
                      idref, "\"maxSamplingRate\" name=\"maxSamplingRate\">", MaxSamplingRate,
                      "</parameter>", "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                      RealParameterMat[i, 3], "\" name=\"minSamplingGap\">0.0</parameter>",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                      RealParameterMat[i, 4], "\" name=\"maxSamplingGap\">0.0</parameter>",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                 </fossilDistr>",
                      "\n", sep = "")
                    cat(file = Temp2, append = T, "                 <taxonset id=\"",
                      as.character(target), "\" spec=\"TaxonSet\">", "\n", sep = "")
                    # Prepare the taxonset block.
                    k = 1
                    for (k in 1:length(Tipnames)) {
                      cat(file = Temp2, append = T, "                    <taxon idref=\"",
                        Tipnames[k], "\" spec=\"Taxon\"/>", "\n", sep = "")
                    }  # End for k
                    cat(file = Temp2, append = T, "                </taxonset>",
                      "\n", "            </distribution>", "\n", sep = "")
                  }  # End else if(is.na(match(NumberofTheClade, DFCladexml[,2]))=='FALSE'){

                }  # End if(setequal(DescenTip, Tipnames)){.

        } else {
            # End of if(length(aa)>0){

            if (length(aa) == 1)
                {
                  # The clade is represented by only one taxon, but according to Matschiner et al.
                  # 2017 this taxon will help the calibration of the more inclusive nodes. A single
                  # taxa clade is never implemented as a topological constraint in the previous
                  # step, so we set up everything from scratch.
                  Tipnames = c(Tipnames, as.character(inputTaxono[aa, 1]))

                  # Prepare the CLADEAGE xml block.
                  logFossilPriorName = c(logFossilPriorName, paste("\"", as.character(paste(target,
                    Tipnames, sep = "_")), ".fossilprior\"", sep = ""))
                  cat(file = Temp2, append = T, "             <distribution id=\"",
                    as.character(paste(target, Tipnames, sep = "_")), ".fossilprior\" spec=\"beast.math.distributions.FossilPrior\" monophyletic=\"true\" tree=\"@Tree.t:",
                    xmltreename, "\">", "\n", sep = "")
                  cat(file = Temp2, append = T, "                 <fossilDistr id=\"FossilCalibration.",
                    i - 1, "\" spec=\"beast.math.distributions.FossilCalibration\">",
                    "\n", sep = "")
                  cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                    RealParameterMat[i, 1], "\" name=\"minOccuranceAge\">", as.character(CalPointTable[i,
                      5]), "</parameter>", "\n", sep = "")
                  cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                    RealParameterMat[i, 2], "\" name=\"maxOccuranceAge\">", as.character(CalPointTable[i,
                      6]), "</parameter>", "\n", sep = "")
                  cat(file = Temp2, append = T, "                     <parameter ",
                    idref, "\"minDivRate\" name=\"minDivRate\">", MinDivRate, "</parameter>",
                    "\n", sep = "")
                  cat(file = Temp2, append = T, "                     <parameter ",
                    idref, "\"maxDivRate\" name=\"maxDivRate\">", MaxDivRate, "</parameter>",
                    "\n", sep = "")
                  cat(file = Temp2, append = T, "                     <parameter ",
                    idref, "\"minTurnoverRate\" name=\"minTurnoverRate\">", MinTurnoverRate,
                    "</parameter>", "\n", sep = "")
                  cat(file = Temp2, append = T, "                     <parameter ",
                    idref, "\"maxTurnoverRate\" name=\"maxTurnoverRate\">", MaxTurnoverRate,
                    "</parameter>", "\n", sep = "")
                  cat(file = Temp2, append = T, "                     <parameter ",
                    idref, "\"minSamplingRate\" name=\"minSamplingRate\">", MinSamplingRate,
                    "</parameter>", "\n", sep = "")
                  cat(file = Temp2, append = T, "                     <parameter ",
                    idref, "\"maxSamplingRate\" name=\"maxSamplingRate\">", MaxSamplingRate,
                    "</parameter>", "\n", sep = "")
                  cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                    RealParameterMat[i, 3], "\" name=\"minSamplingGap\">0.0</parameter>",
                    "\n", sep = "")
                  cat(file = Temp2, append = T, "                     <parameter id=\"RealParameter.",
                    RealParameterMat[i, 4], "\" name=\"maxSamplingGap\">0.0</parameter>",
                    "\n", sep = "")
                  cat(file = Temp2, append = T, "                 </fossilDistr>",
                    "\n", sep = "")
                  cat(file = Temp2, append = T, "                 <taxonset id=\"",
                    as.character(paste(target, Tipnames, sep = "_")), "\" spec=\"TaxonSet\">",
                    "\n", sep = "")
                  # Prepare the taxonset block.
                  cat(file = Temp2, append = T, "                    <taxon idref=\"",
                    Tipnames, "\" spec=\"Taxon\"/>", "\n", sep = "")
                  cat(file = Temp2, append = T, "                </taxonset>", "\n",
                    "            </distribution>", "\n", sep = "")
                }  # End if(length(aa)==1) becasue if(length(aa)==0) there is no taxa, so we do nothing.
        }  # End else if(length(aa)>0){

    }  # End for i

    ## Remove from the original xml file the clades that are already constrained and
    ## also include a CLADEAGE calibration prior.
    xml.new = cbind(seq(1, length(xml.input), by = 1), xml.input)
    tocheck = vector()
    i = 1
    for (i in 1:length(rowsToRemoved)) {
        To.rm = seq(as.numeric(as.character(strsplit(rowsToRemoved[i], ":")[[1]][1])),
            as.numeric(as.character(strsplit(rowsToRemoved[i], ":")[[1]][2])), by = 1)
        tocheck = c(tocheck, seq(as.numeric(as.character(strsplit(rowsToRemoved[i],
            ":")[[1]][1])), as.numeric(as.character(strsplit(rowsToRemoved[i], ":")[[1]][2])),
            by = 1))
        xml.new = xml.new[-match(To.rm, xml.new[, 1]), ]
    }
    xml.new = xml.new[-grep("<distribution id=\"fossilCalibrations\"", xml.new[,
        2], fixed = T), ]
    xml.new2 = xml.new[, -1]



    ## In the modified xml file include the CLADEAGE blocks.
    CladeAgeXmlBlock = readLines(Temp2)

    # length(CladeAgeXmlBlock)

    cat(file = output, xml.new2[c(1:(grep("<operator id=", xml.new2)[1] - 3))], sep = "\n",
        "\n")
    cat(file = output, append = T, CladeAgeXmlBlock, sep = "\n")
    cat(file = output, append = T, "          </distribution>", "\n", sep = "")
    cat(file = output, append = T, "    </distribution>", "\n", "\n", sep = "")
    logidrefpos = grep("<log idref=", xml.new2)
    cat(file = output, append = T, xml.new2[c(grep("<operator id=", xml.new2)[1]:logidrefpos[(length(logidrefpos) -
        3)])], sep = "\n")
    cat(file = output, append = T, paste("        <log idref=", logFossilPriorName,
        "/>", sep = ""), sep = "\n")
    cat(file = output, append = T, xml.new2[c((logidrefpos[(length(logidrefpos) -
        3)] + 1):length(xml.new2))], sep = "\n")

    # Remove the temporary files.
    file.remove(Temp2)
}  # End of the function.
