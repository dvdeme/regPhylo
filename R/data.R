#' Sequences and metadata extracted from the Barcode of Life Database (BOLD), GenBank (NCBI) and the repository of the New Zealand National
#' Fish Collection for the species "Diastobranchus capensis".
#'
#'
#' This example data provides a list of tables storing sequences and metadata for the species "Diastobranchus capensis".
#' Each table corresponds to the output table of the following functions: \code{\link{GetSeq_BOLD}}, \code{\link{GetSeqInfo_NCBI_taxid}},
#' \code{\link{Congr.NCBI.BOLD.perReposit}}, \code{\link{GeoCoord.WGS84}}, \code{\link{GeoCodeName}}.
#' @format A list of six tables:
#' \itemize{
#' \item 1) This table contains the raw output of the function \code{\link{GetSeq_BOLD}} from the regPhylo package (i.e. sequences and metadata extracted from the Barcode of Life Database (BOLD)
#' for the species "Diastobranchus capensis" the 9/10/2018).
#' The table "Seq.BOLD" is a dataframe with 8 rows and 82 columns.
#' \item 2) This table contains the raw output of the function \code{\link{GetSeqInfo_NCBI_taxid}} from the regPhylo package (i.e. sequences and metadata extracted from GenBank (NCBI)
#' for the species "Diastobranchus capensis" the 9/10/2018).
#' The table "Seq.NCBI" is a dataframe with 4 rows and 23 columns.
#' \item 3) This table contains 1 DNA sequence and the associated metadata present in the repository of
#' the New Zealand National Fish Collection of the National Museum of New Zealand, Te Papa Tongarewa.
#' The table "Seq.Perso" is a dataframe with 1 row and 24 columns
#' \item 4) This table contains the raw output of the function \code{\link{Congr.NCBI.BOLD.perReposit}} including data
#' coming from NCBI, BOLD, and a personal repository, for the species "Diastobranchus capensis".
#' The table "Seq.DF" is a data frame with 11 rows and 25 columns.
#' \item 5) This table contains the raw output of the function \code{\link{GeoCoord.WGS84}} which splits
#' the Lat_Lon field into two distinct fields "Latitude" and "Longitude" and converts all the geographic
#' coordinates into WGS84 decimal degrees using +- North, and +- East (i.e. -45.10, -53.23).
#' The table "Seq.DF1" is a dataframe with 11 rows and 26 columns
#' \item 6) This table contains the raw output of the function \code{\link{GeoCodeName}} which extracts the geographic
#' coordinates from nominatim openstreetmap API for all the sequences with only place name information.
#' The table provide two additional columns: the place name used for the query; and information about the precision of the geographic coordinates.
#' The table "Seq.DF3" is a data frame with 11 rows and 28 columns.
#' }
"Seq.Diastocapen"


#' Output of the function SpeciesGeneMat.Bl and Select.DNA for the species "Diastobranchus capensis"
#'
#' A data file for use in function examples. It contains a list of 6 objects:
#' the first five elements are the output of the \code{\link{SpeciesGeneMat.Bl}} function,
#' the sixth object is the table output of the \code{\link{Select.DNA}} function.
#' \itemize{
#' \item The first table provides the species-by-gene matrix.
#' \item The second table provides a summary table from a gene region perspective
#' (i.e. number of species and sequences for each gene region).
#' \item The third table provides a summary table from the
#' species perspective (i.e. number of gene regions and sequences for each species).
#' \item The fourth object provides an empty
#' string listing the species lost when removing microsatellites, unassigned DNA,
#' or blacklisted sequences.
#' \item The fifth element is a table including all the sequences and the metadata exported by the function
#' \code{\link{SpeciesGeneMat.Bl}} after removing microsatellites, unassigned DNA and potentially blacklisted sequences.
#' This table also includes two additional columns, one reporting the entry order of the sequences,
#' and the other reporting the name of the gene region used for the species-by-gene matrix.
#' \item The last object is a table exported by the function \code{\link{Select.DNA}}.
#' }
#'
#' @format A list of five data frames and one string:
#' \itemize{
#' \item 1) The first table provides the species-by-gene matrix
#' a data frame with 1 row and 4 columns.
#' \item 2) The second table provides a summary table from a gene region perspective (i.e. number of species and sequences for each gene region).
#' a data frame with 3 rows and 3 columns
#' \item 3) The third table provides a summary tables from the species perspective (i.e. number of gene regions and sequences for each species).
#' a data frame with 1 row and 3 columns
#' \item 4) The fourth object provides an empty string listing the species lost when removing microsatellites, unassigned DNA, or blacklisted sequences
#' an empty string.
#' \item 5) The fifth object is the fourth table (i.e. "CleanDataTable"), it provides all the sequences and the metadata exported by the function \code{\link{SpeciesGeneMat.Bl}} after removing microsatellites, unassigned DNA and potentially blacklisted sequences.
#' a data frame with 1 rows and 30 columns
#' \item 6) This data frame including all the sequences and metadata for selected gene regions.
#' a data frame with 11 rows and 29 columns (the columns with the entry order of the sequences has been removed).
#' }
"Seq.DF4"


#' Example of an alignment including 319 sequences of 16S mitochondrial DNA, and a table of potential outlier sequences detected by the Detect.Outlier.Seq function.
#'
#'
#' The aim of this example is to illustrate the use of the \code{\link{Detect.Outlier.Seq}} and the \code{\link{Rm.OutSeq.Gap}} functions. The .Rdata object is a list of two objects:
#' the first is an "alignment" of 319 sequences of 16S mitochondrial DNA; the second object "Table.Outlier.Seq.16S" is a list of primary and secondary
#' potential outlier sequences detected using the \code{\link{Detect.Outlier.Seq}} function, based on the the following settings:
#' Table.Outlier.Seq.16S = Detect.Outlier.Seq(inputal = Example_16S_outlier, Strat.DistMat = "Comb", Dist.Th = 0.6, output = "Example_16S_outliers_2.txt", Second.Outlier = "Yes", Bitsc.Th = 0.8)
#' @format A list of two objects 1) an object of class "alignment" including a list of 319 DNA sequences of 16S mitochondrial DNA called "16S_example_outliers.fas",
#' 2) A table with 85 rows and 16 columns called "Table.Outlier.Seq.16S".
"Example_16S_outlier"


#' A list of two tables: a table of the topological constraints to build a multifurcating tree, and a classification table, for 16 species of NZ fish.
#'
#'
#' This list contains two tables for use in the function examples. The first table
#' contains information regarding topological constraints: the first column refers to the hierarchical
#' level of the topological constraints (e.g.= 'Family', or 'Order', or 'Subdivision'...),
#' and the second column refers to the name of the hierarchy (i.e. Aplodactylidae, Anguilliformes).
#' The second table contains classification information: the first column provides the species names
#' (or the name used as tip.label by the phylogenetic tree), then the following columns are the different
#' hierarchical levels of the Linnean classification.
#' @format A list of two tables: the first table has 22 rows and 2 columns; the second table has 16 rows and 23 columns.
"TopoConstraints"






















