#' @title Implements multiple monophyletic constraints in the xml file
#' for BEAST2 using a rooted RAxML tree with bootstrap supports

#' @description This function implements multiple monophyletic constraints in the xml file for
#' BEAST2 using a rooted RAxML tree with bootstrap support. Taxa without DNA can be included
#' in the tree if required. These taxa can be constrained to move freely within a clade,
#'  but the constraints must be specified in the TaxaNoDNA data.frame,
#'  and a taxonomic table must be provided (similar table to the one used for the
#'  function \code{\link{ConstraintTaxo2newick}}).
#'
#' @details If the inclusion of a taxon without DNA applies to a clade including subclades
#' that were constrained in the previous maximum likelihood analysis (RAxML tree),
#' then those constraints are removed to allow the taxon without DNA to move
#' freely across the entire diversity of the clade (for instance, if a taxon
#' without DNA is added into a family X, and the family X contains two genera
#' already constrained in the Maximum likelihood analysis, then these two
#' constraints of monophyly will be removed, to allow the taxon without DNA to
#' also enter within those two genera).

#' @param inputtree a rooted phylogenetic tree of class 'phylo' (usually a RAxML tree
#' already setup with monophyletic constraints) and with bootstrap values reported
#' for each node (we recommand using Dendroscope 3.5.9 to root the tree and to ensure bootstrap
#' values are associated with  the appropriate nodes).
#'
#' @param output name of xml file output ready for BEAST2 analysis
#' (with the path, if necessary).
#'
#' @param bootstrapTH a bootstrap threshold (between [0,100]) used to determine
#' the node that must be constrained (all nodes with a bootstrap equal or
#' above the threshold value, will be constrained to form monophyletic groups).
#'
#' @param xmltreename name of the xml tree in the input xml file.
#'
#' @param input.xml name of the input xml file (with the path, if necessary)
#' edited by beauti without any topological constraints.
#'
#' @param Partitions if 'Yes' indicate the presence of multiple nucleotide partitions in
#' the alignment.
#'
#' @param TaxaNoDNA a three column dataframe (i.e. class "data.frame"), with the species name in
#' the first column (Species names, without space, '?', '-' or '.'), a second
#' column with the hierarchical level where the constraint applies ('Family',
#' 'Genus'...), and a third column with the name of the hierarchical constraint
#' (e.g.  'Aplodactylidae', 'Aplodactylus'). When specifying
#' the hierarchical level and the name of the constraints, ensure these names refer
#' to the levels and names refer to the levels and names already present
#' in the 'TaxoTable' (see below). If a dataframe is provided in the TaxaNoDNA option,
#' then the next three arguments must be specified.
#'
#' @param TaxoTable a classification table of class 'data.frame'. The first
#' column is the species names (or the name used as tip.label by the phylogenetic
#' tree), then the following columns are the different hierarchical levels of the
#' Linnaean classification.
#'
#' @param output.new.TaxoTable name of the new taxonomic table
#' (with the path, if necessary) including taxa without DNA.
#'
#' @param output.new.tree name of the phylogenetic tree
#' (with the path, if necessary) that is going to be
#' exported, including taxa without DNA, in a Newick format.

#' @return An .xml file including all the blocks related to the constrained taxa.
#'
#' @examples
#' \dontrun{
#' # To run the example copy the input files
#' # provided by the package to a temporary directory created in the
#' # current working directory.
#'
#' src.dir = system.file("extdata/TopoConstraints", package = "regPhylo")
#' dir.create("TempDir.TC")
#' # Set up the path of the TempDir folder.
#' dest.dir = paste(getwd(), "/TempDir.TC", sep="")
#' file.names <- dir(src.dir)
#' # Copy all the files stored in regPhylo/extdata/TopoConstraints"
#' # into a temporary folder.
#' sapply(file.names, function(x) {
#' file.copy(from = paste(src.dir, x, sep = "/"),
#' to = paste(dest.dir, x, sep = "/"),
#' overwrite = FALSE) })
#'
#' # All the inputs, including the re-rroted tree are directly available in the
#' # "TempDir.TC" directory. Step one and two described above offer
#' # a detailed description about how a rooted phylogenetic tree with bootstrap support
#' # can be built from R using RAxML and Dendroscope (these steps can be skipped).
#'
#'
#' # First, we build a tree with RAxML using a guiding tree to constraints
#' # the topology considering all constraints.
#'
#' # To run this example RAxML has to be installed and in the $PATH
#' # (for Windows user see the note below).
#'
#' nthread = 5
#' input = "TempDir.TC/Concat.fas"
#' output = "Concat_7GTR_Allconst_autoMRE"
#' ConstraintTree = "TempDir.TC/BackboneTreeAll.txt"
#' PartitionFile = "TempDir.TC/Partitions_Concat.txt_PF2_all.txt"
#'
#' ## in our case the RAxML executable is called "raxmlHPC-PTHREADS-AVX"
#' a=paste("raxmlHPC-PTHREADS-AVX -g ", ConstraintTree,
#' " -f a -x 22345 -p 12345 -# autoMRE -q ", PartitionFile,
#' " -m GTRCAT -T ", nthread, " -s ", input, " -n ",
#' output, sep="")
#' system(a)
#'
#'
#' # Note for Windows user: RAxML executable for Windows can be found at:
#' # https://github.com/stamatak/standard-RAxML/tree/master/WindowsExecutables_v8.2.10
#' # and once the appropriate executbale is downloaded and extracted, RAxML can be run
#' # from that folder.
#' # All input files must be present in the same folder.
#' # RAxML can be run through R using similar type of command that
#' # the one described above for Linux.
#'
#'
#' # Second, we use Dendroscope v.3.5.9 to root the tree to ensure
#' # that each node has the appropriate bootstrap values after re-rooting the tree.
#' # Open the RAxML tree including
#' # the bootstrap supports
#' # (‘RAxML_bipartitions.Concat_7GTR_Allconst_autoMRE_ReRooted’)
#' # in Dendroscope 3.5.9 (available at: http://dendroscope.org/),
#' # specifying that the internal nodes should be interpreted as edge labels.
#' # Then re-rooted the tree using the MRCA of "Halosauropsis_macrochir",
#' # Aldrovandia_affinis", "Bassanago_bulbiceps", and "Conger_verreauxi",
#' # all belonging to the "outgroup" of Elopomorpha.
#'
#' # Load the re-rooted tree into R (a rooted tree is available in the
#' # package and has been loaded into the temporary directory)
#' require(ape)
#' TreeRooted = read.nexus("TempDir.TC/RAxML_bipartitions.Concat_7GTR_Allconst_autoMRE_ReRooted")
#'
#'
#' #### Example including taxa with at least one DNA sequence in the supermatrix.
#'
#' # The input.xml file, called SimpleXml.xml, was prepareb in Beauti. It is available
#' # in the folder "TempDir.TC/xmlfiles". This xml file includes
#' # the basic template of the xml file (including the supermatrix,
#' # the partitions, the substitution and clock models, the priors,
#' # the MCMC parameters). Although the hard constraints are not
#' # included in this input.xml file, they will be incorporated by running the
#' # MultiTopoConst.EditXML4BEAST2 function below.
#'
#' MultiTopoConst.EditXML4BEAST2(inputtree = TreeRooted,
#' output = "TempDir.TC/SimpleXml_Wcont.xml",
#' bootstrapTH = 100, xmltreename = "Subset1",
#' input.xml = "TempDir.TC/SimpleXml.xml",
#' Partitions = "TRUE")
#'
#'
#' #### Example including taxa without DNA in the supermatrix.
#'
#' # Load a table with the new taxa without DNA and the topological constraints.
#' # This table will be used to fill the option "TaxaNoDNA".
#' TaxaNoDNA = cbind(c("Titi_titi", "Toto_toto"), c("Family", "Genus"), c("Congridae", "Scorpaena"))
#' colnames(TaxaNoDNA) = c("SpeciesName", "hier.level", "ConstraintName")
#' TaxaNoDNA = as.data.frame(TaxaNoDNA)
#'
#' # Load the classification table (the same table as used for the
#' # ConstraintTaxo2newick function), there are two way to do this:
#' # either through the .Rdata
#' data(TopoConstraints) # the second object of the list is the classification table
#' dim(TopoConstraints[[2]]) # 16 by 23.
#' # or the classification table that has been loaded into the temporary directory,
#' # and can be loaded into the R environment doing the following
#' ClassifDF = read.csv("TempDir.TC/Classif16sp.csv", header = TRUE)
#' dim(ClassifDF) # 16 by 23
#'
#'
#'
#' # Run the function, but this time the input.xml file called "SimpleXml_2SpNoDNA.xml"
#' # generated by Beauti must be selected. This xml file template was built using
#' # an alignment (called "Concat_2spNoDNA_PF2.nex") including two empty sequences
#' # (i.e. "Titi_titi, "Toto_toto") built by the Align.Concat function.
#' # This alignment has been completed by the nexus block related to the best
#' # partition scheme defined by partitionFinder2 (here we used the same
#' # partition scheme as defined for the 16 taxa with at least one DNA sequence).
#'
#' MultiTopoConst.EditXML4BEAST2(inputtree = TreeRooted,
#' output = "TempDir.TC/SimpleXml_2SpNoDNA_Wcont.xml",
#' bootstrapTH = 100, xmltreename = "Subset1",
#' input.xml = "TempDir.TC/SimpleXml_2SpNoDNA.xml",
#' Partitions = "TRUE", TaxaNoDNA = TaxaNoDNA,
#' TaxoTable = ClassifDF,
#' output.new.TaxoTable = "TempDir.TC/Classif18sp_2NoDNA.csv",
#' output.new.tree = "TempDir.TC/BackboneTreeAll_2spNoDNA.txt")
#'
#' # Plot the newly constrained tree including the 2 taxa without DNA
#' NewTree = read.tree("TempDir.TC/BackboneTreeAll_2spNoDNA.txt")
#' plot(NewTree)
#'
#' # To see the new classification table including the 2 taxa without DNA.
#' NewClassifDF = read.delim("TempDir.TC/Classif18sp_2NoDNA.csv",
#' sep = "\t", header = TRUE)
#' dim(NewClassifDF) # now there are 18 rows and still 23 columns.
#'
#' # To remove the files created while running the example, do the following:
#' unlink("TempDir.TC", recursive = TRUE)
#' a = list.files()
#' a1 = a[grep("Concat_7GTR_Allconst_autoMRE", a)]
#' file.remove(a1)
#'
#' }
#'
#'
#' @export MultiTopoConst.EditXML4BEAST2


MultiTopoConst.EditXML4BEAST2 = function(inputtree = NULL, output = NULL, bootstrapTH = NULL,
    xmltreename = NULL, input.xml = NULL, Partitions = NULL, TaxaNoDNA = NULL, TaxoTable = NULL,
    output.new.TaxoTable = NULL, output.new.tree = NULL) {

    if (is.data.frame(TaxaNoDNA))
        {
            # If a dataframe with taxa without DNA is provided.

            # Add each new taxa without DNA into the RAxML tree.
            TaxoTableNew = TaxoTable

            # Re-order the TaxaNoDNA dataframe according to the node number where the
            # constraint applies. Higher number first and then decreasing in number.
            NodeNb = vector()
            i = 1
            for (i in 1:dim(TaxaNoDNA)[1]) {
                n = match(TaxaNoDNA[i, 2], colnames(TaxoTable))
                if (length(n) == 0) {
                  stop("Hierarchical level not found in the TaxoTable")
                }
                spn = as.character(TaxoTable[which(TaxoTable[, n] == as.character(TaxaNoDNA[i,
                  3])), 1])
                if (length(spn) == 0) {
                  stop("Name of the hierarchical level not found in the table. The name of the clade constraint must be already present in the TaxoTable")
                }
                Nodenb = ape::getMRCA(inputtree, spn)
                if (is.null(Nodenb)) {
                  Nodenb = "NA"
                }
                NodeNb = c(NodeNb, Nodenb)
            }
            TaxaNoDNA = cbind(TaxaNoDNA, NodeNb)
            options(warn = -1)
            TaxaNoDNA = TaxaNoDNA[order(as.numeric(as.character(TaxaNoDNA[, 4])),
                decreasing = F), ]  # Re-order the TaxaNoDNA dataframe in decreasing order of the node number where the constraints apply.
            options(warn = 0)

            # Add the taxa without DNA into the tree within the appropriate clade defined by
            # the constraints.
            i = 1
            for (i in 1:dim(TaxaNoDNA)[1]) {
                n = match(TaxaNoDNA[i, 2], colnames(TaxoTable))
                spn = as.character(TaxoTable[which(TaxoTable[, n] == as.character(TaxaNoDNA[i,
                  3])), 1])
                SpChoice = sample(spn, 1)  # choose a species at random within the clade of interest as an anchor to create a new tip.
                options(warn = -1)  # Disable warning for bind.tip function
                inputtree = phytools::bind.tip(inputtree, tip.label = as.character(TaxaNoDNA[i,
                  1]), where = which(inputtree$tip.label == SpChoice))  # Bind a new tip for the taxa without DNA at random close to a terminal tip belonging to a selected group present in the object 'spn'.
                options(warn = 0)
                inputtree = ape::compute.brlen(inputtree)  # Compute the branch lengths of the new tree (in our case branch length doesn't matter, just the topology and bootstrap values for the nodes are important).

                nnb = ape::getMRCA(inputtree, c(spn, as.character(TaxaNoDNA[i, 1])))  # Extract the node number of the MRCA of the new clade including the taxa without DNA.
                inputtree$node.label[nnb - length(inputtree$tip.label)] = "100"  # Impose a bootstrap value above the threshold to constrain this node in order to allow the bayesian search to randomly position the tip of the taxon without DNA within this clade.
                desc = phytools::getDescendants(inputtree, nnb)  # Extract the descendant node and tips.
                ndesc = desc[which(desc > length(inputtree$tip.label))]  # Select the nodes only.
                inputtree$node.label[ndesc - length(inputtree$tip.label)] = "NA"  # All the nodes within the group have an assigned bootstrap value of NA in order to remove the constraints within the clade allowing the taxa to move freely across all the diversity of the clade (so it removes the potential nested constraints applied at lower hierarchical levels or supported by the Maximum likelihood analysis with RAxML).

                # Feed the taxonomic table the taxa without DNA.
                a = rep("", dim(TaxoTable)[2])
                a[1] = as.character(TaxaNoDNA[i, 1])
                a[2:n] = as.character(as.vector(as.matrix(TaxoTable[which(TaxoTable[,
                  1] == SpChoice), c(2:n)])))
                if (!colnames(TaxoTable)[n] == "Genus") {
                  # If the constraints apply to hierarchical level other than the genus, then we
                  # feed the new taxoTable extracting the genus information from the species name
                  # of the taxa.
                  a[which(colnames(TaxoTable) == "Genus")] = strsplit(as.character(TaxaNoDNA[i,
                    1]), "_")[[1]][1]
                }
                a = t(as.matrix(a))
                colnames(a) = colnames(TaxoTableNew)
                TaxoTableNew = rbind(TaxoTableNew, as.matrix(a))  # Include the taxonomic information of the taxa without DNA in the new taxonomic table.
            }
            utils::write.table(TaxoTableNew, file = output.new.TaxoTable, sep = "\t", row.names = F,
                quote = F)
            ape::write.tree(inputtree, file = output.new.tree)
        }  # End is.table(TaxaNoDNA)){



    ## Load the xml file.
    input.xml1 = readLines(input.xml)

    ## Extract the clades with a bootstrap support above the threshold that must be
    ## constrained.

    ## Make node labels.
    NodeName = ape::makeNodeLabel(inputtree, method = "number")
    # Make a table of correspondence between the node labels and the bootstrap
    # values.
    DFtab = as.matrix(cbind(NodeName$node.label, inputtree$node.label, nodeNb = seq(1,
        length(inputtree$node.label), by = 1)))
    # Selection of all the nodes with bootstrap support >= to the bootstrap threshold
    # (bootstrapTH).
    options(warn = -1)
    NodeBoot100 = DFtab[which(as.numeric(as.character(DFtab[, 2])) >= as.numeric(as.character(bootstrapTH))),
        ]
    options(warn = 0)

    # Create a list of tip names for every node in the tree.
    Tipnamelist = caper::clade.members.list(NodeName, tip.labels = TRUE)
    SpNamAlP = vector()  # A vector to store the tip names is already present in the previous constraints.
    Temp1 = "xmlConstraintOnly.Temp1.xml"

    # Create a xml block with all the taxonset.
    cat(file = Temp1, "     <taxonset id=\"Set_All\" spec=\"TaxonSet\">", "\n", sep = "")
    i = 1
    for (i in 1:length(inputtree$tip.label)) {
        cat(file = Temp1, append = T, "         <taxon id=\"", inputtree$tip.label[i],
            "\" spec=\"Taxon\"/>", "\n", sep = "")
    }
    cat(file = Temp1, append = T, "     </taxonset>", "\n", sep = "")

    # If the alignment has multiple partitions then set up the xmltree name using the
    # following approach.
    if (Partitions == "TRUE") {
        ggg = input.xml1[grep("<tree id=\"Tree.t:", input.xml1, perl = T)]
        xmltreename = substr(ggg, (regexpr("Tree.t:", ggg, fixed = T) + 7), (regexpr("\" name=",
            ggg, fixed = T) - 1))
    }


    # Create a second block with the specific monophyletic groups.
    Temp1b = "xmlConstraintOnly.Temp1b.xml"
    i = 1
    for (i in 1:dim(NodeBoot100)[1]) {
        cat(file = Temp1b, append = T, "            <distribution id=\"", NodeBoot100[i,
            1], ".prior\" spec=\"beast.math.distributions.MRCAPrior\" monophyletic=\"true\" tree=\"@Tree.t:",
            xmltreename, "\">", "\n", "                <taxonset id=\"", NodeBoot100[i,
                1], "\" spec=\"TaxonSet\">", "\n", sep = "")
        splistClade = Tipnamelist[[as.numeric(as.character(NodeBoot100[i, 3]))]]
        j = 1
        for (j in 1:length(splistClade)) {
            cat(file = Temp1b, append = T, "                    <taxon idref=\"",
                splistClade[j], "\" spec=\"Taxon\"/>", "\n", sep = "")
        }  # End for j
        cat(file = Temp1b, append = T, "                </taxonset>", "\n", "            </distribution>",
            "\n", sep = "")
        SpNamAlP = c(SpNamAlP, splistClade)
        SpNamAlP = unique(SpNamAlP)
    }  # End for i

    # Insert the xml file storing the topological constraints previously built into a
    # new xml file.
    ConstraintXml = readLines(Temp1b)
    ConstraintXml_1 = readLines(Temp1)
    cat(file = output, input.xml1[1:(grep("<run id=", input.xml1)[1] - 1)], sep = "\n")
    cat(file = output, append = T, ConstraintXml_1, sep = "\n", "\n")
    cat(file = output, append = T, input.xml1[(grep("<run id=", input.xml1)[1]):(grep("<operator id=",
        input.xml1)[1] - 3)], sep = "\n")
    cat(file = output, append = T, ConstraintXml, sep = "\n")
    cat(file = output, append = T, input.xml1[(grep("<operator id=", input.xml1)[1] -
        2):length(input.xml1)], sep = "\n")
    file.remove(c(Temp1, Temp1b))  # Remove the temporary files including the topological constraints only.
}  # End of the function.
