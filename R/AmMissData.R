#' @title Estimate the percentage of missing gene regions within a species-by-gene matrix
#'
#' @description This function extracts the percentage of missing data in the species-by-gene
#' matrix according to a list of selected gene regions.
#'
#' @param input the species-by-gene matrix which is provided as the first object in the
#' list returned by the function SpeciesGeneMat.Bl.R.
#' @param gene.list a vector of the selected gene regions (gene names
#'  have to be consistent with the header of the table with the suffix
#'  '_CleanDataset.txt' exported by the function SpeciesGeneMat_Bl.R).
#'  For example: Cytochrome c oxidase subunit 1, should be written 'co1'
#'  and not 'COI' or 'COX1'.
#'
#' @return The function returns the percentage of missing gene regions.
#'
#' @examples # Load a species-by-gene matrix, such as that exported by SpeciesGeneMat.Bl function.
#' data(Seq.DF4) ## the first object of the list is the species-by-gene matrix
#'
#' # Run the function
#' AmMissData(input = Seq.DF4[[1]], gene.list = c("co1", "16srrna"))
#'
#' @export AmMissData
#'
AmMissData = function(input = NULL, gene.list = NULL) {
    DF = as.matrix(input[, stats::na.omit(match(gene.list, names(as.data.frame(input))))])
    mode(DF) = "numeric"
    DFpa = DF
    DFpa[which(DFpa[] > 1)] <- 1  # Turn abundance into presence/absence data.
    PourMissing = 100 - (sum(DFpa) * 100)/(dim(DFpa)[1] * dim(DFpa)[2])
    return(PourMissing)
}
