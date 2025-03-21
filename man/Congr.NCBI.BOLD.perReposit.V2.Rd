% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Congr.NCBI.BOLD.perReposit.R
\name{Congr.NCBI.BOLD.perReposit.V2}
\alias{Congr.NCBI.BOLD.perReposit.V2}
\title{Merge outputs from GenBank, BOLD (local queries) and a personal repository
(if available) into a single table and remove potential duplicate
sequences.}
\usage{
Congr.NCBI.BOLD.perReposit.V2(
  input.NCBI = NULL,
  input.BOLD = NULL,
  output = NULL,
  input.perReposit = NULL,
  perReposit = NULL
)
}
\arguments{
\item{input.NCBI}{table coming from the function \code{\link{GetSeqInfo_NCBI_taxid}}.}

\item{input.BOLD}{table coming from the function \code{\link{GetSeq_BOLD}}.}

\item{output}{name of the output table in txt format.}

\item{input.perReposit}{a personal repository in a table format loaded
in the R environment with the following fields: 'TaxaName', 'AccessBold',
'AccessNCBI', 'Sequence', 'SeqLength', 'Definition', 'OrganismClassif',
'Source', 'Title', 'Authors', 'Journal', 'Pubmed', 'Year', 'Organism',
'Organelle', 'Mol_type', 'Db_xref', 'Product', 'Genes', 'Location',
'isolation_source', 'Lat_lon', 'Collection_date', 'Date_Extract'.  The fields
must follow the same definition as above, except that 'AccessBold',
'AccessNCBI' should remain empty, and the internal unique id of the sequence
must be reported in the field 'Db_xref'. Any other information in this field
must be provided after the unique sequence id and separated by '; '.}

\item{perReposit}{if a name is provided then compare the NCBI and BOLD output to the
information contained in the personal repository, in case some sequences in the
personal repository have been already deposited in NCBI or BOLD.}
}
\value{
The information contained in the output file includes:
\itemize{
\item 'TaxaName' (Binomial species name included in the second column of the NCBI
and BOLD species list),
\item 'AccessBold' (concatenation of the BOLD processid and the marker_code = sequenceID in BOLDV5
web portal),
\item 'AccessNCBI' (NCBI sequence accession number),
\item 'Sequence' (The nucleotide sequence itself, or 'Too_Long' if more than 5000bp),
\item 'SeqLength' (Sequence length),
\item 'Definition' (NCBI definition field, or for BOLD a combination of
'TaxaName', the string 'BOLD' and the BOLD 'processid'),
\item 'OrganismClassif' (Classification of the organism),
\item 'Source' (NCBI source or for BOLD combination of 'TaxaName', 'Voucher_type', 'processid', 
'marker_code'),
\item 'Title' (Title provided by NCBI),
\item 'Authors' (Authors provided by NCBI, or for BOLD 'collectors'),
\item 'Journal' (Journal provided by NCBI),
\item 'Pubmed' (Pubmed references provided by NCBI),
\item 'Year' (Year provided by NCBI, or for BOLD year of the 'collection_date_start'),
\item 'Organism' (Organism provided by NCBI or for BOLD 'TaxaName'),
\item 'Organelle' (Organelle provided by NCBI),
\item 'Mol_type' (Mol_type provided by NCBI or 'genomic DNA' for BOLD),
\item 'Db_xref' (Db_xref for NCBI, or for BOLD combination of 'processid',
 'specimenid', sampleid', 'recordID', 'museumid', collection_code', 'inst', 'sovereign_inst' separated by '; '),
\item 'Product' (Product provided by NCBI),
\item 'Genes' (Gene provided by NCBI or, for BOLD 'marker_code'),
\item 'Location' (Field Country in NCBI, or for BOLD combination of 'country.ocean', 'province.state',
'region', 'sector', 'site'),
\item 'isolation_source' (isolation_source provided by NCBI),
\item 'Lat_lon' (Lat_lon provided by NCBI or for 'coord' for BOLD),
\item 'Collection_date' (Collection_date provided by NCBI, or for BOLD 'collection_date_start'),
\item 'Date_Extract' (Date of the extraction of the data),
\item 'OriginDatabase' (Name of the Original database).
}
}
\description{
This function assembles a large common dataframe for all the data retrieved
from GenBank (through the NCBI) and BOLD (from a local database file 
that can be downloaded from https://boldsystems.org/data/data-packages/) ), and 
also from a personal repository, if available.
The function checks for the presence of duplicated sequences and
selects the most relevant information (i.e. selects the longest sequence, and
maximise the metadata information across the different sources for the
location, geographic coordinates, and collection date).
}
