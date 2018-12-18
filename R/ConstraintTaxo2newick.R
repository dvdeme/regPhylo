#' @title Build a multifurcating topological constraining tree for RAxML
#'
#' @description The function builds a multifurcating phylogenetic tree
#' from a classification table and a table of phylogenetic constraints ready to be
#' used by RAxML as constraint tree (option -g in RAxML) to guide the
#' reconstruction of the molecular phylogenetic tree.

#' @param inputTaxo a classification table, the first column is the species names (or
#' the name used as tip.label by the phylogenetic tree), then the following
#' columns are the different hierarchical levels of the Linnean classification.
#' @param inputConst is a two column table, the first column refers to the hierarchical
#' level of the topological constraints (e.g.= 'Family', or 'Order', or
#' 'Subdivision'...), the name of the hierarchical level have to be the same as
#' the headers of the classification table, the second column refers to the name
#' of the taxa with a constraint of monophyly (e.g. 'Aplodactylidae',
#' Aulopiformes', 'Percomorphaceae'....).
#' @param outputNewick name of the output multifurcating newick tree that will
#' be exported in a .txt file, (can also include the path to the folder).

#' @return This function exports into the R environment a list of 2 objects; the first
#' object is the taxonomic table modified to include the constraints, and the
#' second object is the multifurcating tree converted in a 'phylo' object.
#' The function also exports a newick tree in a txt document that can be used to constrain
#' the topology in RAxML.
#'
#' @details Warnings: branch lengths of the multifurcating tree are misleading, only the
#' topology matters.
#'
#' @examples # Load a tables listing the topological constraints (first object of the list)
#' # and the classification table (second object of the list).
#' \dontrun{
#' data(TopoConstraints)
#' # The table storing the constraints include 22 topological constraints overall
#' # including constraints at the Family, Order, Series, Subdivision, Division,
#' # Subsection, Subcohort, Cohort, Supercohort, Infraclass, Subclass.
#' #
#' # The Classification table include 16 species from the New Zealand marine
#' # ray-finned fish species list.
#'
#' # Create a Temporary folder to store the outputs of the function.
#' dir.create("TempDir.TopoConstraints")
#' # Run the function considering all the constraints
#' BackBoneTreeAll = ConstraintTaxo2newick(inputTaxo = TopoConstraints[[2]],
#' inputConst = TopoConstraints[[1]], outputNewick = "TempDir.TopoConstraints/BackboneTreeAll")
#'
#' # plot the constraining tree (the branch length does not matter, only the topology matters).
#' plot(BackBoneTreeAll[[2]], cex=0.8)
#'
#' # Use only the constraint at the Family level
#' FamilyConst=TopoConstraints[[1]][TopoConstraints[[1]][,1]=="Family",]
#'
#' # Run the function considering only the constraints at the family level.
#' BackBoneTreeFamily = ConstraintTaxo2newick(inputTaxo = TopoConstraints[[2]],
#' inputConst = FamilyConst, outputNewick = "TempDir.TopoConstraints/BackboneTreeFamily")
#'
#' # plot the constraining tree (the branch length does not matter,
#' # only the topology matters), notice that only constrained taxa
#' # are present on the guiding tree, the unconstrained taxa will
#' # be positioned on the tree based on their molecular affinities.
#' plot(BackBoneTreeFamily[[2]], cex=0.8)
#'
#' # Use only the constraint at the Family and Series levels.
#' FamilySeriesConst=TopoConstraints[[1]][c(which(TopoConstraints[[1]][,1] == "Family"),
#' which(TopoConstraints[[1]][,1] == "Series")),]
#'
#' # Run the function considering only the constraints at the family and order levels.
#' BackBoneTreeFamilySeries = ConstraintTaxo2newick(inputTaxo = TopoConstraints[[2]],
#' inputConst = FamilySeriesConst, outputNewick = "TempDir.TopoConstraints/BackboneTreeFamilySeries")
#'
#' # plot the constraining tree (the branch length does not matter,
#' # only the topology matters). notice that only constrained taxa
#' # are present on the guiding tree, the unconstrained taxa will
#' # be positioned on the tree based on their molecular affinities.
#' plot(BackBoneTreeFamilySeries[[2]], cex=0.8)
#'
#' # To clean the files created while running the example do the following:
#' unlink("TempDir.TopoConstraints", recursive = TRUE)
#' }

#' @export ConstraintTaxo2newick

ConstraintTaxo2newick = function(inputTaxo = NULL, inputConst = NULL, outputNewick = NULL) {

    # Evaluate the number of hierarchical levels with topological constraints.
    NbHier = length(unique(inputConst[, 1]))
    Hier = as.character(unique(inputConst[, 1]))

    # Extract the column of the taxonomic table targeting the constraints.
    taxo1 = as.matrix(inputTaxo[, match(Hier, colnames(inputTaxo))])
    row.names(taxo1) = inputTaxo[, 1]
    colnames(taxo1) = colnames(inputTaxo)[match(Hier, colnames(inputTaxo))]

    # Check if some species are unconstrained by the topological constraints, in
    # which case unconstrained species are removed from the classification table
    # (they are not needed by RAxML and free to move everywhere according to their
    # molecular affinities).
    Speciestotal = vector()
    i = 1
    for (i in 1:NbHier) {
        ## Check within each hierarchical level.
        p = inputConst[which(inputConst[, 1] == Hier[i]), 2]
        j = 1
        for (j in 1:length(p)) {
            Speciestotal = c(Speciestotal, names(which(taxo1[, match(Hier[i], colnames(taxo1))] ==
                p[j])))
        }
    }
    UnconstraintSp = setdiff(row.names(taxo1), unique(Speciestotal))
    if (length(UnconstraintSp) > 0) {
        taxo1 = taxo1[-match(UnconstraintSp, row.names(taxo1)), ]  # Remove the unconstrained species from the taxonomic table.
    }

    if(class(taxo1)=="character"){
      taxo1 = as.matrix(taxo1)
      colnames(taxo1) = colnames(inputTaxo)[match(Hier, colnames(inputTaxo))]
    }

    # Replace missing hierarchical levels with NA.
    taxo1[which(taxo1 == "")] = NA

    i = 1
    for (i in 1:NbHier) {
        # Within each hierarchical level.
        a = setdiff(unique(taxo1[, i]), inputConst[which(inputConst[, 1] == Hier[i]),
            2])  # Identify the taxa that are not monophyletic within a hierarchical level.

        if (i > 1) {
            # Check if we are dealing with the first hierarchical level.

            if (length(a) > 0)
                {
                  # Check if some of the clades within a hierarchical level are not constrained.
                  j = 1
                  for (j in 1:length(a)) {
                    if (is.na(a[j])) {
                      # In case a = NA.
                      taxo1[as.vector(which(is.na(taxo1[, i]) == TRUE)), i] = paste(Hier[i],
                        "_", as.vector(taxo1[as.vector(which(is.na(taxo1[, i]) ==
                          TRUE)), i - 1]), sep = "")
                    } else {
                      taxo1[which(taxo1[, i] == a[j]), i] = paste(Hier[i], "_", as.vector(taxo1[which(taxo1[,
                        i] == a[j]), i - 1]), sep = "")
                    }  ## End else if(a[j]==NA){.
                  }  # End for j.
                }  # End if(length(a)>0){
        } else {
            # In case i=1 (first hierarchical level), need to select the row.names of the
            # classification table (=tip.label of the phylogenetic tree).
            if (length(a) > 0)
                {
                  j = 1
                  for (j in 1:length(a)) {
                    if (is.na(a[j])) {
                      # In case a = NA.
                      taxo1[as.vector(which(is.na(taxo1[, i]) == TRUE)), i] = paste(Hier[i],
                        "_", row.names(taxo1)[as.vector(which(is.na(taxo1[, i]) ==
                          TRUE))], sep = "")
                    } else {
                      taxo1[which(taxo1[, i] == a[j]), i] = paste(Hier[i], "_", row.names(taxo1)[which(taxo1[,
                        i] == a[j])], sep = "")
                    }  ## End else if(a[j]==NA){
                  }  # End for j
                }  # End if(length(a)>0){
        }  # End if(i>1){ else
    }  # End for i
    taxo1 = as.data.frame(taxo1)  # Convert into a dataframe to allow columns to be converted as factors.


    # Taxo2newick function is a modified function of the function taxo2phylog from
    # ade4 R package used by ConstraintTaxo2newick to generate the newick tree.
    # taxo= is a taxonomic classification table with the species names as row names
    # and each hierarchical level in a separate column. The classification level has
    # to be factors.

    taxo2newick = function(taxo, abbrev = TRUE) {
      nc <- ncol(taxo)  ## NB. of hierarchical levels
      for (k in 1:nc) {
        w <- as.character(k)
        w <- paste("l", w, sep = "")
        w1 <- levels(taxo[, k])
        if (abbrev)
          w1 <- abbreviate(w1)
        levels(taxo[, k]) <- paste(w, w1, sep = "")
      }
      leaves.names <- row.names(taxo)
      res <- paste(";", sep = "")
      x <- taxo[, nc]
      xred <- as.character(levels(x))
      w <- "("
      w = paste(w, paste(xred, collapse = ","), sep = "")
      res <- paste(w, ")", res, sep = "")
      res <- sub(",\\)", "\\)", res)
      for (j in nc:1) {
        x <- taxo[, j]
        if (j > 1)
          y <- taxo[, j - 1] else y <- as.factor(leaves.names)
          for (k in 1:nlevels(x)) {
            w <- "("
            old <- as.character(levels(x)[k])
            yred <- unique(y[x == levels(x)[k]])
            yred <- as.character(yred)
            w = paste(w, paste(yred, collapse = ","), sep = "")
            w <- paste(w, ")", sep = "")
            w <- sub(",\\)", "\\)", w)
            res <- gsub(old, w, res)
          }
      }
      return(res)
    }

    # Convert a taxonomic table in newick using taxo2newick function.
    treeA = taxo2newick(taxo1, abbrev = TRUE)
    treeA1 = phytools::read.newick(text = treeA)  # Read a newick tree and convert it into a phylo object.
    treeA2 = ape::collapse.singles(treeA1)  # Collapse single node.
    ape::write.tree(treeA2, file = paste(outputNewick, ".txt", sep = ""))  # Export the tree in newick format.
    treeA3 = ape::read.tree(paste(outputNewick, ".txt", sep = ""))
    # Export a list of two objects
    # 1)= the modified classification table taking into account the topological constraints,
    # 2)= a final newick tree.
    return(list(TaxonomicTable = taxo1, NewickConstraintTree = treeA3))
}
