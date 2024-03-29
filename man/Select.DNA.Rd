% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Select.DNA.R
\name{Select.DNA}
\alias{Select.DNA}
\title{Extract all the sequences and metadata from the data table for a list of selected gene regions}
\usage{
Select.DNA(
  input = NULL,
  gene.list = NULL,
  output = NULL,
  timeout = 10,
  Seqinr = FALSE
)
}
\arguments{
\item{input}{The table created using the function \code{\link{SpeciesGeneMat.Bl}}
with the suffix '_CleanDataset.txt'.}

\item{gene.list}{a vector of the selected gene regions (gene names have
to be consistent with the header of the table with the suffix '_SpDNA_Mat.txt'
exported by the function SpeciesGeneMat_Bl.R). For example (Cytochrome c oxidase subunit 1, should be
written 'co1' and not 'COI' or 'COX1').}

\item{output}{the name of the output table exported into the working directory.}

\item{timeout}{the timeout in seconds for socketConnection used in, the
choosebank function of the seqinr R package. The default is 10 seconds.
It might be necessary to increase the timeout (i.e. the time to get an answer from
the server) if the function cannot retrieve any DNA sequence for
a certain species when DNA sequences are known to be available in GenBank.
Alternatively, if the server connection is quick the timeout can be decreased to 5 seconds
(the default in choosebank) to speed-up the function.}

\item{Seqinr}{by default FALSE, so the rentrez R package is used to query Genbank, 
if TRUE, the old version using seqinr R package is used.}
}
\description{
This function selects the sequences and associated information for a list
of selected gene regions for all the species present in the data table. For
full mitochondrial genomes, or all sequences with a sequence length > 5000 bp
listed as 'TooLong', the function extracts the sequences regions of interest (from GenBank)
 using the query function of the seqinr R package.
}
\examples{
# Load the data table exported by the SpeciesGeneMat.Bl function
data(Seq.DF4) # the table of sequences and metadata is called "CleanDataTable"
# and is the fifth object of the list "Seq.DF4".

# Run the function
Seq.DF5=Select.DNA(input = Seq.DF4$CleanDataTable, gene.list = c("co1", "16srrna"),
output = "Seq.DF4.dataTable")

}
