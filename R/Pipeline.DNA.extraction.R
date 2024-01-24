#' @title A pipeline to automate the retrieval and organisation of DNA sequences and metadata
#' from NCBI and BOLD databases from a list fo taxa of interest.


#' @description This function allows to extract DNA sequences and metadata from NCBI and BOLD,
#' dereplicates the data from both databases, improve the georeferencing of DNA sequences
#' using openstreetmap API, built a species by gene matrix, select the DNA region with the best
#' taxonomic coverage, export DNA sequences and the associated metadata for selected Gene region,
#' select the best DNA sequences per taxa and perform first alignment allowing automatic detection
#' of sequences entered in a wrong direction.



#' @details this function will require a specific NCBI personal API key in order to
#' get a faster access to the NCBI server.
# NCBI personal API key: 4XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# I give my API key to access NCBI server allow a faster processing of the data.
# set_entrez_key("4XXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
# Sys.getenv("ENTREZ_KEY")


#' @details This function build a pipeline of the following function "GetSeqNumber_NCBI_taxid" (to detect the number of sequence per taxon and identify the model species with too many
#' sequence to download), "GetSeqInfo_NCBI_taxid" (get the DNA and the metadata from NCBI),
#' "GetSeq_BOLD" (get the DNA and the metadata from BOLD), "Congr.NCBI.BOLD.perReposit" (dereplicate
#' and improve metadata from NCBI and BOLD), "GeoCoord.WGS84" (homogenize the geographic coordinates),
#' "GeoCodeName" (improve the georeferencing of the data with just locality names using the
#' openstreetmap API), "SpeciesGeneMat.Bl" (build the species-by-gene matrix),
#' "SelGene.MaxSpCov" (identify the DNA region maximising the species coverage), "Select.DNA" (export
#' the metadata of the selected DNA region, including the extraction of the DNA region within larger
#' gene region such as the COI within the mitogenome), "SelBestSeq" (select the DNA sequences based on
#' several criteria including some geographic proximity to some geographic region of interest, and
#' export the alignments in fasta format), "First.Align.All" (detection of sequences that are entered
#' in the wrong direction and first alignments using Mafft and pasta software).
#' When a taxa is detected with too many sequences set by the max.seq option (10 000sequences by
#' default), these taxa are momentarily excluded and the species-by-gene matrix is performed without
#' them to identify the DNA region with the highest taxonomic coverage, then once those region
#' identified, the DNA sequences and metadata of those region are extracted for the taxa initially
#' recognized with too many sequences. This function allows to run different steps, see the definition of the "Steps" parameter.



#' @param input.ncbi a two column table with the unique NCBI taxonomic ID (i.e. taxid) in
#' the first column (incluidng the taxonomic ID of synonyms species names)
#' and the accepted binomial species name in the second column
#' which will be retained as 'TaxaName' in the output file.

#' @param input.bold a two column table with the first column containing all
#' the binomial species names including synonyms and the second column contains
#' only the names that will be reported in the output table as 'species_name'.

#' @param gene can be a particular gene region or 'ALL' for all entries
#' of that species. Please note that this will only retrieve genes
#' with the same name in Genbank. For example, "COX1" will not return the
#' same entries as "COI". Default "ALL"

#' @param Path.output the path to the output folder within which all the intermediate and final
#' results files will be stored.

#' @param Nb.DNA.marker can be a numeric ,eg. 4, to select the first gene with the highest species
#' coverage, it can be a gene name ,eg. "co1", or it an be "Min.All.Sp.cov" in order to include the
#' minimum gene region to get 100 percent species coverage.

#' @param Remove.Large.Nuc.Fragment if TRUE (by default) long nuclear DNA assembly fragment such,
#' whole nuclear genome assembly, chromosomes assembly are removed from the NCBI ("GenBank") 
#' output data file of the GetSeqInfo_NCBI_taxid function. So far all fragment longer than 50000 pb
#' are removed.

#' @param Steps is a vector of numeric value from 1 to 3, with 1 extract the DNA sequences and
#' associated metadata from NCBI and BOLD, dereplicate the data, improve the metadata location, and
#' build the Species by Gene matrix, step 2 selects the gene of interest and exports the metadata
#' and the alignments in fasta format, step 3 performs a first alignment approach of the selected gene
#' of interest to detect sequence entered in wrong direction and performed automatic reverse
#' complement. The steps are performed sequentially so the steps 2 cannot be performed alone
#' without performing the steps 1.

#' @param nthread, if steps 3 is present in the option Step (eg. Steps = c(1,2,3)) then we can choose
#' how many thread we would like to use to perform the first alignment steps in case there are multiple
#'  markers to align.

#' @param nthread.per.align if steps 3 is present in the option Step (eg. Steps = c(1,2,3)) then we
#' can choose how many threads we would like to use to perform the first alignment step. This option
#' is interesting when few alignments have many sequence > 1000 DNA sequences.

#' @param methods if steps 3 is present in the option Step (eg. Steps = c(1,2,3)) then we can choose
#' the method used for the alignment approach, it can be "mafft_auto", "mafftfftnsi", "mafftfftns2",
#' and "pasta". See the function "First.Align" from the regPhylo R package. It can also be a
#' combination of several methods eg. c("mafftfftnsi", "mafftfftns2", "pasta").

#' @return this function returns a set of intermediate and final files and folders within the folder
#' precised in the Path.output option. The structure of the folders is organised as follow,
#' The function return a list of 4 objects using the following headers:
#' \itemize{
#' \item the "Raw.Data" folder containing the raw data files of NCBI, "Seq.NCBI.txt", and BOLD
#' ,Seq.BOLD.txt", DNA and metadata extraction, the combine file of DNA sequences from NCBI an BOLD
#' after dereplication ,"Seq_NCBI_BOLD.txt", the intermediate files with all the improve
#' georeferencing "Seq_NCBI_BOLD_Geo1.txt", "Seq_NCBI_BOLD_Geo2.txt", the files of the
#' "species-"by-gene" matrix "Matrix.Sp.DNA.NaiveSpDNA_Mat.txt" including files with all the
#' metadata "Matrix.Sp.DNA.Naive_CleanDataset.txt", and summary files of the species-by-gene matrix
#' from a DNA perspective "Matrix.Sp.DNA.NaiveSummary_DNApers.txt" and from a species perspective
#' "Matrix.Sp.DNA.NaiveSummary_SPECIESpers.txt".
#' \item the "Selected.DNA" folder containing the non-aligned fasta files of the selected DNA regions
#' and the metadata of all the selected DNA sequences "Align_InfoTabSelSeq.txt", and the folder
#' "FirstAlign"if the step 3 of the pipeline was performed. This folder contains the alignments of
#' the selected DNA region from the different alignment software specify in the option "methods".
#' \item the folder "Selected.META" contains a file of the Metadata of the sequences of the
#' selected DNA regions "MetaData.Select.DNA.Select.DNA.txt".
#' }


#' @export Pipeline.DNA.extraction



Pipeline.DNA.extraction = function(input.ncbi = NULL,
                                   input.bold = NULL,
                                   gene = "ALL",
                                   max.seq = 10000,
                                   Path.output = NULL,
                                   Nb.DNA.marker = 1,
                                   Remove.Large.Nuc.Fragment = TRUE,
                                   Steps = c(1, 2),
                                   nthread = 1,
                                   methods = "mafftfftns2",
                                   nthread.per.align = NULL) {
  # 1_Get the number of sequence per taxon for NCBI
  NCBI.nb.seq = as.data.frame(GetSeqNumber_NCBI_taxid(splist = input.ncbi, gene = gene))
  NCBI.nb.seq$NbSequences = as.numeric(NCBI.nb.seq$NbSequences)
  
  
  # 2_Exclude the taxon with too many sequences, with a threshold define by max.seq
  a = which(NCBI.nb.seq[, "NbSequences"] > max.seq)
  if (length(a) > 0) {
    input.ncbi.short = input.ncbi[match(NCBI.nb.seq[-a, 1], input.ncbi[, 1]), ]
  } else {
    input.ncbi.short = input.ncbi
  }
  
  dir.create(paste(Path.output, "/Raw.Data", sep = ""))
  
  Path.output.raw = paste(Path.output, "/Raw.Data/", sep = "")
  
  
  # 3_Extract the DNA sequence and Metadata from GenbBank.
  NCBI.stats = GetSeqInfo_NCBI_taxid(
    splist = input.ncbi.short,
    gene = gene,
    filename = paste(Path.output.raw, "Seq.NCBI.txt", sep = "")
  )
  
  # Load the NCBI sequence and metadata in R
  NCBI.data = read.delim(
    paste(Path.output.raw, "Seq.NCBI.txt", sep = ""),
    sep = "\t",
    header = T
  )
  
  # Remove DNA entry related to very large nuclear fragment assembly 
  # such as the whole genome assembly, chromosomes assembly etc... 
  if(Remove.Large.Nuc.Fragment == TRUE){
    # First remove all the accession number corresponding to Whole genome assembly projet
    a = which(nchar(NCBI.data[,"AccessNb"]) > 8)
    if(length(a)> 0){
      NCBI.data =  NCBI.data[-a,]
    }
    # Convert the sequence length into numeric values
    NCBI.data$SeqLength = as.numeric(as.character(NCBI.data$SeqLength))
    
    # detect sequence possibly longer than  50 000 pb.
    b = which(NCBI.data[,"SeqLength"] > 50000)
    if(length(b) > 0){
      NCBI.data =  NCBI.data[-b,]
    }
  }
  
  
  
  
  # 4_Extract DNA sequences and Metadata from BOLD.
  
  BOLD.data = GetSeq_BOLD(
    splist = input.bold,
    filename = paste(Path.output.raw, "Seq.BOLD.txt", sep = ""),
    bold.id = FALSE
  )
  BOLD.stats = as.data.frame(BOLD.data[[1]])
  
  # 5_Merge and dereplication of GenBank and BOLD Data.
  AllSeqDF = Congr.NCBI.BOLD.perReposit(
    input.NCBI = NCBI.data,
    input.BOLD = BOLD.data[[2]],
    output = paste(Path.output.raw, "Seq_NCBI_BOLD.txt", sep = ""),
    input.perReposit = NULL
  )
  
  
  # 6_Homogeneize the Geographic information related to the sequences. (not very important at that stage).
  AllSeqDF1 = GeoCoord.WGS84(
    input = AllSeqDF,
    output = paste(Path.output.raw, "Seq_NCBI_BOLD_Geo1.txt", sep = "")
  )
  
  # 7_Improve the geographic information (might not be very reliable, but you do not care really).
  AllSeqDF2 = GeoCodeName(
    input = AllSeqDF1,
    output = paste(Path.output.raw, "Seq_NCBI_BOLD_Geo2.txt", sep = ""),
    AutoCorrNZ = FALSE
  )
  
  
  # 8_Build the Species by gene matrix (That's very Important !!).
  Sp.DNA.Mat.naive = SpeciesGeneMat.Bl(
    input = AllSeqDF2,
    output = paste(Path.output.raw, "Matrix.Sp.DNA.Naive", sep = "")
  )
  
  
  
  # 9_If some species were excluded from the NCBI search because there were too many sequences, we can extract the sequence from the selected gene regions of interest only for those taxa.
  
  a = which(NCBI.nb.seq[, "NbSequences"] > max.seq)
  if (length(a) > 0) {
    # Select the Taxa with too many sequences.
    input.ncbi.short2 = input.ncbi[match(NCBI.nb.seq[a, 1], input.ncbi[, 1]), ]
    
    # 9.1_Select the gene of interest based on the initial search, to get those gene regions for the
    # taxa with too many sequences.
    # Minimum gene region to get 100% species coverage
    Min.Gene.Max.cov = SelGene.MaxSpCov(input = Sp.DNA.Mat.naive$Species.Gene_matrix)
    if (class(Nb.DNA.marker) == "numeric") {
      GeneSelection = Sp.DNA.Mat.naive$Summary_DNA[1:Nb.DNA.marker, 1]
    } else {
      if (Nb.DNA.marker == "Min.All.Sp.cov") {
        if (class(Min.Gene.Max.cov) == "character") {
          GeneSelection = Min.Gene.Max.cov[[2]]
        } else {
          GeneSelection = Min.Gene.Max.cov$List_Gene_Name_Full_Species_Coverage
        }
        
      } else {
        GeneSelection = Nb.DNA.marker
      }
    }
    
    
    # Homogenize the gene names according to NCBI requirements for the request.
    gene.listB = toupper(GeneSelection)
    gene.listB = gsub("CO1", "COI", gene.listB, fixed = TRUE)
    gene.listB = gsub("SRRNA", "S ribosomal RNA", gene.listB, fixed = TRUE)
    
    ## check the presence of the CO1 gene in the selection
    co1.pres = match("co1", GeneSelection)
    
    if (!is.na(co1.pres)) {
      ## include additional gene nomenclature for the CO1.
      gene.listB = c(gene.listB, "CO1", "COX1", "COXI")
    }
    
    
    # check if all the alternative name are associated to DNA sequences.
    NCBI.nb.seq.V2 = do.call(rbind, lapply(1:length(gene.listB), function(x) {
      df = as.data.frame(GetSeqNumber_NCBI_taxid(splist = input.ncbi.short2, gene = gene.listB[x]))
      c(gene.listB[x], sum(as.numeric(df$NbSequences), na.rm = T))
    }))
    
    #Remove all the gene name without DNAsequence retrieved
    gene.listB.short = NCBI.nb.seq.V2[which(as.numeric(NCBI.nb.seq.V2[, 2]) > 0), 1]
    
    # Replace by tolower and remove space
    gene.listB.export = tolower(gsub(" ", "", gsub(
      "S ribosomal RNA", "SRRNA", gene.listB.short
    )))
    
    # Export the new DNA sequences.
    NCBI.data2 = do.call(rbind, lapply(1:length(GeneSelection), function(x) {
      NCBI.stats2 = GetSeqInfo_NCBI_taxid(
        splist = input.ncbi.short2[2:4, ],
        gene = gene.listB[x],
        filename = paste(
          Path.output.raw,
          "Seq.NCBI.",
          gene.listB.export[x],
          ".txt",
          sep = ""
        )
      )
      read.delim(
        paste(
          Path.output.raw,
          "Seq.NCBI.",
          gene.listB.export[x],
          ".txt",
          sep = ""
        ),
        sep = "\t",
        header = T
      )
    }))
    
    # Remove the potential duplicated sequences
    To.remove = which(duplicated(NCBI.data2$AccessNb) == TRUE)
    if (length(To.remove) > 0) {
      NCBI.data2 = NCBI.data2[-To.remove, ]
    }
    
    
    # Merge and dereplication of GenBank and BOLD Data.
    AllSeqDF.V2 = Congr.NCBI.BOLD.perReposit(
      input.NCBI = NCBI.data2,
      input.BOLD = NULL,
      output = paste(Path.output.raw, "Seq_NCBI.Congr.txt", sep = ""),
      input.perReposit = NULL
    )
    
    # Homogeneize the Geographic information related to the sequences. (not very important at that stage).
    AllSeqDF1.V2 = GeoCoord.WGS84(
      input = AllSeqDF.V2,
      output = paste(Path.output.raw, "Seq_NCBI.Congr_Geo1.txt", sep = "")
    )
    
    # Improve the geographic information (might not be very reliable, but you do not care really).
    AllSeqDF2.V2 = GeoCodeName(
      input = AllSeqDF1.V2,
      output = paste(Path.output.raw, "Seq_NCBI.Congr_Geo2.txt", sep = ""),
      AutoCorrNZ = FALSE
    )
    
    
    AllSeqDF3 = as.data.frame(rbind(AllSeqDF2, AllSeqDF2.V2))
    # 8_Build the Species by gene matrix (That's very Important !!).
    Sp.DNA.Mat.naive = SpeciesGeneMat.Bl(
      input = AllSeqDF3,
      output = paste(Path.output.raw, "Matrix.Sp.DNA.Naive", sep = "")
    )
    
  }
  
  if (length(Steps) == 1) {
    results = Sp.DNA.Mat.naive
    
    return(results)
  }
  
  
  if (match(2, Steps) == 2) {
    ## PERFORM THE STEPS 2
    
    # 9.1_Select the gene of interest based on the initial search, to get those gene regions for the
    # taxa with too many sequences.
    # Minimum gene region to get 100% species coverage
    Min.Gene.Max.cov = SelGene.MaxSpCov(input = Sp.DNA.Mat.naive$Species.Gene_matrix)
    if (class(Nb.DNA.marker) == "numeric") {
      GeneSelection = Sp.DNA.Mat.naive$Summary_DNA[1:Nb.DNA.marker, 1]
    } else {
      if (Nb.DNA.marker == "Min.All.Sp.cov") {
        if (class(Min.Gene.Max.cov) == "character") {
          GeneSelection = Min.Gene.Max.cov[[2]]
        } else {
          GeneSelection = Min.Gene.Max.cov$List_Gene_Name_Full_Species_Coverage
        }
        
      } else {
        GeneSelection = Nb.DNA.marker
      }
    }
    
    # 11_Export all DNA sequences and metadata for the selected gene regions
    # Extract all the sequences (including from long annotated DNA sequences) and
    # associated information.
    PoolSeq = read.delim(
      paste(
        Path.output.raw,
        "Matrix.Sp.DNA.Naive_CleanDataset.txt",
        sep = ""
      ),
      sep = "\t",
      h = T,
      check.names = F
    )
    #PoolSeq = read.delim("../Output/Extract_DNA/Sp.DNA.Mat.Naive__CleanDataset.txt", sep="\t", h=T)
    
    ## Create a folder to export the metadata
    dir.create(paste(Path.output, "/Selected.META", sep = ""))
    Path.output.Select.Meta = paste(Path.output, "/Selected.META", sep = "")
    
    # Export the metadata
    PoolAllSeq.Select.Genes = Select.DNA(
      input = PoolSeq,
      gene.list = GeneSelection,
      output = paste(Path.output.Select.Meta, "/MetaData.Select.DNA", sep = ""),
      timeout = 15,
      Seqinr = FALSE
    )
    
    
    # 12_Extraction of all the sequences per species for the selected gene regions
    
    ## Create a folder to export the selected DNA alignments
    dir.create(paste(Path.output, "/Selected.DNA", sep = ""))
    Path.output.Select.DNA = paste(Path.output, "/Selected.DNA/Align", sep = "")
    
    
    Export.AllSeq.5genes = SelBestSeq(
      input = PoolAllSeq.Select.Genes,
      output = Path.output.Select.DNA,
      RefPoint = cbind(0, 0),
      Alignment = T,
      MaxSeq = "ALL",
      gene.list = GeneSelection,
      SeqChoice = "Median"
    )
    
    if (length(Steps) == 2) {
      results = c(
        list(DNA.Seq.And.Metadata.Selected.DNA = PoolAllSeq.Select.Genes),
        list(GeneSelection = GeneSelection),
        Sp.DNA.Mat.naive
      )
      
      return(results)
    }
    
  }
  
  if (match(3, Steps) == 3) {
    ## PERFORM THE STEPS 3
    
    dir.create(paste(Path.output, "/Selected.DNA/FirstAlign", sep = ""))
    
    # 13_First alignment of all sequences and detection of sequences that should be reverse complemented.
    
    # Some variable must be defined in the R Global environment if the alignment software (MAfft and Pasta) run the multiple DNA markers in parallel. (each alignment is performed by a single thread).
    input = paste(Path.output, "/Selected.DNA", sep = "") #"../Output/Extract_DNA/AllSeqNaive_5genes"
    output = paste(Path.output, "/Selected.DNA/FirstAlign", sep = "")
    nthread = nthread
    methods = methods
    
    Res.align = First.Align.All(
      input = input,
      output = output,
      nthread = nthread,
      methods = methods,
      nthread.per.align = nthread.per.align
    )
    
    results = c(
      list(Alignment.log.file = Res.align),
      list(DNA.Seq.And.Metadata.Selected.DNA = PoolAllSeq.Select.Genes),
      list(GeneSelection = GeneSelection),
      Sp.DNA.Mat.naive
    )
    return(results)
  }
  
  
} # End of the function




#' @title Get the number of DNA sequences available for taxa of interest in NCBI.

#' @description This function provides the number of DNA sequences available for the taxa of interest
#' in th NCBI database, this function is used internally in the Pipeline.DNA.extraction function

#' @param splist a two column table with the unique NCBI taxonomic ID (i.e. taxid) in
#' the first column and the accepted binomial species name in the second column
#' which will be retained as 'TaxaName' in the output file.

#' @param gene can be a particular gene region or 'ALL' for all entries
#' of that species. Please note that this will only retrieve genes
#' with the same name in Genbank. For example, "COX1" will not return the
#' same entries as "COI". Default "ALL"

#' @param chunk_size number of records to be downloaded at a time from Genbank.
#' If connection drops try decreasing this value to 20 or smaller. Default 50.


GetSeqNumber_NCBI_taxid = function(splist = NULL,
                                   gene = "ALL",
                                   chunk_size = 50) {
  # Extract the character of string from the right.
  substrRight <- function(x, n) {
    substr(x, nchar(x) - n + 1, nchar(x))
  }
  
  print(paste("Starting processing:", date()))
  # Empty file to start sequence deposition.
  #write(paste("Starting processing:", date(), filename))
  
  # Summary table providing the number of sequences for each taxa.
  TabSpGenSum = matrix("NA", ncol = 2)[-1,]
  colnames(TabSpGenSum) = c("TaxaName", "NbSequences")
  
  
  # Loop over species.
  for (k in 1:dim(splist)[1]) {
    if (gene == "ALL") {
      # If gene='ALL' all the DNA sequences of that species are recorded.
      oo = rentrez::entrez_search(
        db = "nucleotide",
        term = paste("txid", splist[k, 1], " [Organism:exp]", sep = ""),
        use_history = TRUE
      )
      
    } else {
      # If the name of the gene is provided, the sequences of this gene of interest are recorded.
      oo = rentrez::entrez_search(
        db = "nucleotide",
        term = paste(
          "txid",
          splist[k, 1],
          " [Organism:exp]",
          " AND ",
          gene ,
          "[All Fields] ",
          sep = ""
        ),
        use_history = TRUE
      )
    }
    # Filling up the summary table
    imax = as.numeric(oo$count)
    TabSpGenSum = rbind(TabSpGenSum, c(as.character(splist[k, 1]), imax))
    print(paste(as.character(splist[k, 2]), ": ", imax, "seq."))
  }
  return(TabSpGenSum)
}
