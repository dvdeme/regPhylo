#' @title Estimate the minimal number of gene regions to maximise species coverage

#' @description This function returns additional information including the number of genes
#' required to have 100% species coverage, the list of those genes, and a species
#' accumulation curve from a gene perspective giving an idea of the contribution
#' of each gene to maximise the species coverage.

#' @details For a pre-selected number of
#' genes, the function returns a list of the gene names maximizing the species
#' coverage and the species coverage. For a preselected set of genes the function
#' provides the species coverage, the list of species potentially missing, and the
#' minimum list of genes enabling us to have 100% coverage.

#' @param input the species-by-gene matrix which is provided as the first object in the
# list returned by the function SpeciesGeneMat.Bl.R.
#' @param NBGene either a number of genes that we would like to include optimizing
#' the species coverage, or a vector of pre-defined gene regions
#' or NULL when the function is turned off.

#' @examples # Load a Species-by-gene matrix exported from the SpeciesGeneMat.Bl function
#' data(Seq.DF4) ## the first object of the list is the Species-by-gene matrix
#'
#' # Run the function without a pre-selection of gene regions.
#' SelGene.MaxSpCov(input = Seq.DF4[[1]])
#' # Run the function with a pre-selection of gene regions.
#' SelGene.MaxSpCov(input = Seq.DF4[[1]], NBGene = c("co1", "12srrna"))
#' SpbyGeneMat=rbind(as.matrix(Seq.DF4[[1]]), c("Titi_titi",0, 2, 1), c("Toto_toto", 0, 0, 4))
#' row.names(SpbyGeneMat)=SpbyGeneMat[,1]
#' SpbyGeneMat=as.data.frame(SpbyGeneMat)
#' SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = c("co1", "12srrna"))
#' SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = c("co1", "16srrna"))
#' SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = c("co1"))
#' SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = c("12srrna"))
#' SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = c("16srrna"))
#' SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = 1)
#' SelGene.MaxSpCov(input = SpbyGeneMat, NBGene = 2)


#' @export SelGene.MaxSpCov


SelGene.MaxSpCov = function(input = NULL, NBGene = NULL) {
  if(dim(input)[1]==1){ # In case there is only one species.)
    return(list("There is only one species in the species-by-gene matrix, all the gene regions are present",
                  input))
  } else {

    if (dim(input)[2] == 2) {
        # In the case where there is only a single gene.
        return(paste("A single gene covers 100% of the species", colnames(input)[2],
            sep = ""))
    } else {
        # In the case where there is more than one gene.

        # Species coverage accumulation curve for the genes
        transfoNumeric = function(x) as.numeric(as.character(x))  # Convert the species-by-gene matrix into 'numeric' class.
        Sp.DNAMatA <- apply(input[, -1], 2, transfoNumeric)
        Sp.DNAMatPA = Sp.DNAMatA
        Sp.DNAMatPA[Sp.DNAMatPA > 1] <- 1  # Transform the abundance to Presence/Absence data.

        # In case there is only one gene selected.
        if(is.character(NBGene) && length(NBGene) == 1 | is.numeric(NBGene) && NBGene == 1){
          # if the NBGene is a character
          if (is.character(NBGene)) {
            Sp.DNAMatPA = Sp.DNAMatPA[, sort(match(NBGene, colnames(Sp.DNAMatPA)))]
            Sp.DNAMatPA = as.matrix(Sp.DNAMatPA)
            row.names(Sp.DNAMatPA) = row.names(input)
            colnames(Sp.DNAMatPA) = NBGene
            CovPer = round((colSums(Sp.DNAMatPA)/dim(Sp.DNAMatPA)[1])*100, 2)
            return(list(paste("The selected gene region ", NBGene, " covers ", CovPer,
                              " % of the species", sep = ""), Sp.DNAMatPA))
          }

          # if the NBGene is a number
          if(is.numeric(NBGene)){
            sumcol = colSums(Sp.DNAMatPA, na.rm = TRUE)
            BestGene = colnames(Sp.DNAMatPA)[which(sumcol == max(sumcol))]
            res = as.matrix(Sp.DNAMatPA[, which(sumcol == max(sumcol))])
            colnames(res) = BestGene
            row.names(res) = row.names(input)
            return(list(paste("The ", BestGene, " equally maximize(s) the species coverage.",
                              sep = ""),
                        paste("Selecting the gene region maximizing the species coverage covers ",
                              round((max(sumcol)/dim(Sp.DNAMatPA)[1])*100, 2), " % of species",
                              sep = ""), res))
          }
        } # end if(is.character(NBGene) && length(NBGene)==1 | is.numeric(NBGene) && NBGene==1){

        #if(dim(input)[1]==1){ # In case there is only one species.
        #  Sp.DNAMatPA=t(as.data.frame(Sp.DNAMatPA))
        #  row.names(Sp.DNAMatPA)=row.names(input)
        #}

        if (is.character(NBGene)) {
          Sp.DNAMatPA = Sp.DNAMatPA[, sort(match(NBGene, colnames(Sp.DNAMatPA)))]
          row.names(Sp.DNAMatPA) = row.names(input)
        }

        accuCurve = vector()
        accuCurve = c(accuCurve, sum(Sp.DNAMatPA[, 1]))
        i = 1
        for (i in 2:dim(Sp.DNAMatPA)[2]) {
            a = rowSums(Sp.DNAMatPA[, c(1:i)])  # Select the columns with the DNA fragment you are interested in.
            a[which(a > 1)] <- 1
            accuCurve = c(accuCurve, sum(a))
        }


        # Select the minimum number of genes to improve the species coverage.
        BringNsp = vector()
        BringNsp = c(BringNsp, sum(Sp.DNAMatPA[, 1]))
        for (i in 2:length(accuCurve)) {
            BringNsp = c(BringNsp, accuCurve[i] - accuCurve[i - 1])
        }
        info1 = length(which(BringNsp > 0))
        BringNsp = cbind(BringNsp, genename = colnames(Sp.DNAMatPA))

        # Determine the minimum number of genes required
        BringNsp2 = BringNsp[which(as.numeric(as.character(BringNsp[, 1])) > 0),
            ]

        ## If a single gene gets 100% species coverage the class(BringNsp2) will be a
        ## 'character', otherwise it will be a 'matrix' if more than one gene is required
        ## to get 100% species coverage.
        if (class(BringNsp2) == "matrix") {
            Sp.DNAMatPA2 = Sp.DNAMatPA[, match(BringNsp2[c(1:dim(BringNsp2)[1]),
                2], colnames(Sp.DNAMatPA))]

            # Determine the best arrangement of genes that minimise the number of genes but
            # maximising the species coverage.
            NBGeneOpti = vector()
            GeneOpti = list(2:dim(Sp.DNAMatPA2)[2])
            i = 1
            for (i in 2:dim(Sp.DNAMatPA2)[2]) {
                Sp.DNAMatPA2r = cbind(Sp.DNAMatPA2[, -i], Sp.DNAMatPA2[, i])
                colnames(Sp.DNAMatPA2r) = c(colnames(Sp.DNAMatPA2)[-i], colnames(Sp.DNAMatPA2)[i])

                accuCurve2 = vector()
                accuCurve2 = c(accuCurve2, sum(Sp.DNAMatPA2r[, 1]))
                j = 1
                for (j in 2:dim(Sp.DNAMatPA2r)[2]) {
                  a = rowSums(Sp.DNAMatPA2r[, c(1:j)])  # Select the columns with the DNA fragment you are interested in.
                  a[which(a > 1)] <- 1
                  accuCurve2 = c(accuCurve2, sum(a))
                }  ## End for j.

                BringNsp2 = vector()
                BringNsp2 = c(BringNsp2, sum(Sp.DNAMatPA2r[, 1]))
                for (z in 2:length(accuCurve2)) {
                  BringNsp2 = c(BringNsp2, accuCurve2[z] - accuCurve2[z - 1])
                }  ## End for z.
                NBGeneOpti = c(NBGeneOpti, length(which(BringNsp2 > 0)))
                GeneOpti[[i]] = cbind(BringNsp2, genename = colnames(Sp.DNAMatPA2r))
            }  ## End for i.
            GeneOpti = GeneOpti[-1]  # Removed the first element of the list not needed.

            BestComb = which(NBGeneOpti == min(NBGeneOpti))  # Determine the best combination of genes.
            if (length(BestComb > 1))
                {
                  BestComb = BestComb[1]
                }  # When multiple combinations lead to the same results, we arbitrarily select the first combination.
            GeneOpti[[BestComb]]

            ## The best arrangement of the species-by-gene matrix that minimizes the number of
            ## genes, but maximize the species coverage.
            Sp.DNAMatPA3 = cbind(Sp.DNAMatPA[, match(GeneOpti[[BestComb]][, 2], colnames(Sp.DNAMatPA))],
                Sp.DNAMatPA[, -match(GeneOpti[[BestComb]][, 2], colnames(Sp.DNAMatPA))])
            colnames(Sp.DNAMatPA3) = c(GeneOpti[[BestComb]][, 2], colnames(Sp.DNAMatPA)[-match(GeneOpti[[BestComb]][,
                2], colnames(Sp.DNAMatPA))])

            accuCurve = vector()
            accuCurve = c(accuCurve, sum(Sp.DNAMatPA3[, 1]))
            i = 1
            for (i in 2:dim(Sp.DNAMatPA3)[2]) {
                a = rowSums(Sp.DNAMatPA3[, c(1:i)])  # Select the columns with the DNA fragment you are interested in.
                a[which(a > 1)] <- 1
                accuCurve = c(accuCurve, sum(a))
            }

            # Select the minimum number of genes to improve the species coverage.
            BringNsp = vector()
            BringNsp = c(BringNsp, sum(Sp.DNAMatPA3[, 1]))
            for (i in 2:length(accuCurve)) {
                BringNsp = c(BringNsp, accuCurve[i] - accuCurve[i - 1])
            }
            info1 = length(which(BringNsp > 0))
            BringNsp = cbind(BringNsp, genename = colnames(Sp.DNAMatPA3))


            # Accumulation curve after ordering the contribution of each gene to maximize the
            # species coverage.
            BringNspMat = cbind(GeneOrder = seq(1, dim(Sp.DNAMatPA)[2]), BringNsp)
            BringNspMatOrda = BringNspMat[order(as.numeric(as.character(BringNspMat[,
                2])), decreasing = T), ]

            BringNspMatOrd = apply(BringNspMatOrda[, -3], 2, transfoNumeric)
            info2 = BringNspMatOrda[which(BringNspMatOrd[, 2] > 0), 3]

            NewAccuCurve = vector()
            NewAccuCurve = c(NewAccuCurve, sum(Sp.DNAMatPA[, 1]))
            for (i in 2:dim(BringNspMatOrd)[1]) {
                NewAccuCurve = c(NewAccuCurve, sum(BringNspMatOrd[c(1:i), 2]))
            }

            # A function that optimizes the selection of genes maximising the species
            # coverage for a given number of genes.
            if (is.numeric(NBGene)) {
                info3 = BringNspMatOrda[c(1:NBGene), 3]
                info4 = NewAccuCurve[NBGene]
                info5 = (NewAccuCurve[NBGene]/max(NewAccuCurve)) * 100
                graphics::plot(NewAccuCurve, main = "Species coverage accumulation curve from a gene perspective",
                  xlab = "Gene index", ylab = "Number of species", type = "l", col = "red",
                  lwd = 2)
                graphics::abline(h = max(NewAccuCurve), lty = 2, lwd = 1.5, col = "black")
                graphics::abline(v = NBGene, lty = 2, lwd = 1.5, col = "blue")
                graphics::abline(v = info1, lty = 2, lwd = 1, col = "grey")
                graphics::abline(h = info4, lty = 2, lwd = 1.5, col = "green")
                graphics::legend("bottomright", c("Species accumulation curve", "Total number of species",
                  paste("Number of genes selected=", NBGene, sep = ""), paste("Minimum number of genes for 100% species coverage=",
                    info1, sep = ""), paste("Species coverage for selected genes=",
                    info4, " (", round(info5, 2), "%)", sep = "")), seg.len = 2,
                  lwd = c(2, 1.5, 1.5, 1, 1.5), lty = c(1, 2, 2, 2, 2), col = c("red",
                    "black", "blue", "grey", "green"), box.lty = 1, box.lwd = 1,
                  box.col = "black", inset = 0.05, cex = 0.75)
                return(list(Minimum_Number_of_Gene_With_Full_Species_Coverage = info1,
                  List_Gene_Name_Full_Species_Coverage = info2, List_Gene_Name_Selected_Genes = info3,
                  Species_Coverage_Selected_Genes = info4, Species_Coverage_Selected_Genes_percentage = info5))
            }

            # A function that determines the species coverage for a given set of genes.
            if (is.character(NBGene)) {
              graphics::plot(NewAccuCurve, main = "Species coverage accumulation curve from a gene perspective",
                  xlab = "Gene index", ylab = "Number of species", ylim = c(0, dim(input)[1]),
                  type = "l", col = "red", lwd = 2)
              graphics::abline(h = max(NewAccuCurve), lty = 2, lwd = 1.5, col = "red")
              graphics::abline(h = dim(input)[1], lty = 2, lwd = 1.5, col = "blue")
              graphics::abline(v = info1, lty = 2, lwd = 1, col = "grey")
              graphics::legend("bottom", c("Species accumulation curve", "Total number of species",
                  paste("Species coverage for selected genes=", max(NewAccuCurve),
                    " (", round(max(NewAccuCurve)/dim(input)[1] * 100, 2), "%)",
                    sep = ""), paste("Minimum nb. of selected genes maximising the species coverage: ",
                    info1, sep = "")), seg.len = 2, lwd = c(2, 1.5, 1.5, 1), lty = c(1,
                  2, 2, 2), col = c("red", "blue", "red", "grey"), box.lty = 1, box.lwd = 1,
                  box.col = "black", inset = 0.05, cex = 0.75)
                return(list(SpeciesCoverage = max(NewAccuCurve)/dim(input)[1] * 100,
                  GeneContributionCoverage = BringNspMatOrda[, c(2, 3)], MissingSpecies = as.vector(input[which(rowSums(Sp.DNAMatPA) ==
                    0), 1]), MinimumGeneMaximisingSpeciesCoverageForGeneSelection = BringNspMatOrda[c(1:info1),
                    3]))
            }

            # A function that determines the minimum number of genes offering 100% species
            # coverage.
            if (is.null(NBGene))
                {
                  # Drawing a species accumulation curve from a gene perspective.
              graphics::plot(NewAccuCurve, main = "Species coverage accumulation curve from a gene perspective",
                    xlab = "Gene index", ylab = "Number of species", type = "l",
                    col = "red", lwd = 2)
              graphics::abline(h = max(NewAccuCurve), lty = 2, lwd = 1.5, col = "blue")
              graphics::abline(v = info1, lty = 2, lwd = 1.5, col = "green")
              graphics::legend("bottomright", c("Species accumulation curve", "Total number of species",
                    paste("Minimum number of genes for 100% species coverage=", info1,
                      sep = "")), seg.len = 2, lwd = c(2, 1.5, 1.5), lty = c(1, 2,
                    2), col = c("red", "blue", "green"), box.lty = 1, box.lwd = 1,
                    box.col = "black", inset = 0.05, cex = 0.75)
                  return(list(Minimum_Number_of_Gene_With_Full_Species_Coverage = info1,
                    List_Gene_Name_Full_Species_Coverage = info2))
                }  # End if(is.null(NBGene))

        } else {
            # If a single gene gets 100% species coverage.
            Sp.DNAMatPA2 = Sp.DNAMatPA[, match(BringNsp2[2], colnames(Sp.DNAMatPA))]

            if (is.null(NBGene)) {
                return(c("A single gene covers 100% of species", BringNsp2[[2]]))
            }

            if (is.numeric(NBGene)) {
                return(c("A single gene covers 100% of species", BringNsp2[[2]]))
            }

            if (is.character(NBGene)) {
                return(list(c(paste("A single gene in the selection covers ", round((as.numeric(as.character(BringNsp2[[1]]))/dim(input)[1]) *
                  100, 2), "% of species, the potential that other genes do not provide genetic information for any additional species",
                  sep = ""), BringNsp2[[2]]), Matrix_PresenceAbsence_Of_Gene_Per_Species = Sp.DNAMatPA))
            }
        }  # end of if(class(BringNsp2)=='matrix'){ } else {}
    }  # end if(dim(input)[2]==2){
  } # end if(dim(input)[1]==1){
}  # end of the function
