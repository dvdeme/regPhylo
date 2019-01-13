#' @title Estimate the species overlap between each per of selected gene region
#'
#' @description This function estimates the number of species overlapping between each pair of
#' selected gene regions.
#' @return The function returns a first matrix providing the raw
#' number of species overlapping between pairs of gene regions, and a second
#' matrix that provides the proportion of species overlap between pairs of gene
#' regions, the proportion is computed as a percentage of the richest gene region
#' (i.e. the gene region that is represented in the greatest number of species).


#' @param input the species-by-gene matrix which is provided as the first object in the
# list returned by the function SpeciesGeneMat.Bl.R.
#' @param gene.Sel a vector of a selected gene region
#' (gene names have to be consistent with the header of the
#' table with the suffix '_CleanDataset.txt' exported by the function
#' SpeciesGeneMat_Bl.R).  Example (Cytochrome c oxydase subunit 1, should be
#' written 'co1' and not 'COI' or 'COX1').

#' @examples # Load a species-by-gene matrix for instance exported by SpeciesGeneMat.Bl function.
#' data(Seq.DF4) ## the first object of the list is the Species-by-gene matrix
#'
#' # Run the function.
#' Matrix.Overlap(input = Seq.DF4[[1]], gene.Sel = c("co1", "16srrna"))
#'
#'
#' @export Matrix.Overlap


Matrix.Overlap = function(input = NULL, gene.Sel = NULL) {
    input = input[, match(gene.Sel, colnames(input))]
    nbi = seq(1, dim(input)[2], 1)
    nbj = seq(1, dim(input)[2], 1)
    ResMat = matrix(NA, ncol = length(nbj), nrow = length(nbi))
    ResMatp = matrix(NA, ncol = length(nbj), nrow = length(nbi))
    i = 1
    for (i in 1:length(nbi)) {
        spnami = row.names(input)[which(as.numeric(as.character(input[, nbi[i]])) >
            0)]
        j = 1
        for (j in 1:length(nbj)) {
            spnamj = row.names(input)[which(as.numeric(as.character(input[, nbj[j]])) >
                0)]
            ResMat[nbi[i], nbj[j]] = length(intersect(spnami, spnamj))
            ResMatp[nbi[i], nbj[j]] = length(intersect(spnami, spnamj))/max(c(length(spnami),
                length(spnamj))) * 100
        }
    }
    colnames(ResMat) = colnames(input)
    row.names(ResMat) = colnames(input)
    colnames(ResMatp) = colnames(input)
    row.names(ResMatp) = colnames(input)
    return(list(NumberOfSpecies = ResMat, PercentageOfSpecies = ResMatp))
}
