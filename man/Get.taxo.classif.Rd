% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get.taxo.classif.R
\name{Get.taxo.classif}
\alias{Get.taxo.classif}
\title{Extract the species list, classification and check for the presence of synonyms in the 
NCBI taxonomy, itis, gbif, bold and taxref databases for a taxa.}
\usage{
Get.taxo.classif(
  input.id = NULL,
  db = NULL,
  downto = "species",
  api_key = NULL,
  local.db = FALSE,
  input.ID.Rank = "Genus"
)
}
\arguments{
\item{input.id}{a vector with the taxid of the taxa of interest (extract all the descendant
species, from the ncbi, bold, itis, gbif and taxref databases).}

\item{db}{a vector with the names of the database to query such as "ncbi", "itis", "gbif", "bold",
or "taxref" (TAXREFV1.5, from the rtaxref R package).}

\item{downto}{a hierachical level to look for descent, by default it is "species". see function
downstream from the taxize R package.}

\item{api_key}{a string with the api key of the database.}

\item{local.db}{[currently not working but work in progress] if TRUE (by default FALSE), it uses
local imported database to query the classification and perfomed request, can be much faster for
important request. This is working with "ncbi", "gbif" (not for the search of synonyms). if FALSE
the function use the API to query the online server.}

\item{input.ID.Rank}{the Rank of the ID requested eg. "Genus" or "Species"}
}
\value{
The function returns a table with the following fields
c("species.valid", "species.syn", "genus", "family", "order", "class", "txid.species.valid",
"txid.species.syn", "txid.genus", "txid.family", "txid.order", "txid.class", "database",
"extraction.date").
}
\description{
This functions extracts the species list, the classification and the check for the
presence of synonyms species names in the ncbi taxonomy, itis, gbif, bold, taxref (MNHN NPN)
databases from a taxid of a taxon from a given hierarchy within the Linnean classification, from
genus or higher (eg. family , order..), obtain by the get_ids("plecoptera", db = "ncbi") function
for instance.
}
\details{
This function works only if the api_key of the user is provided in the option
\emph{api_key}. This function used the \emph{rtaxref} R package to acces the TAXREF database.
}
