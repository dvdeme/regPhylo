#' @title Extract all the sequences and metadata from the data table for a list of selected gene regions

#' @description This function selects the sequences and associated information for a list
#' of selected gene regions for all the species present in the data table. For
#' full mitochondrial genomes, or all sequences with a sequence length > 5000 bp
#' listed as 'TooLong', the function extracts the sequences regions of interest (from GenBank)
#'  using the query function of the seqinr R package.
#'

#' @param input The table created using the function \code{\link{SpeciesGeneMat.Bl}}
#' with the suffix '_CleanDataset.txt'.
#'
#' @param gene.list a vector of the selected gene regions (gene names have
#' to be consistent with the header of the table with the suffix '_SpDNA_Mat.txt'
#' exported by the function SpeciesGeneMat_Bl.R). For example (Cytochrome c oxidase subunit 1, should be
#' written 'co1' and not 'COI' or 'COX1').
#' @param output the name of the output table exported into the working directory.

#' @param timeout the timeout in seconds for socketConnection used in, the
#' choosebank function of the seqinr R package. The default is 10 seconds.
#' It might be necessary to increase the timeout (i.e. the time to get an answer from
#' the server) if the function cannot retrieve any DNA sequence for
#' a certain species when DNA sequences are known to be available in GenBank.
#' Alternatively, if the server connection is quick the timeout can be decreased to 5 seconds
#'(the default in choosebank) to speed-up the function.
#'
#'@param Seqinr by default FALSE, so the rentrez R package is used to query Genbank, 
#'if TRUE, the old version using seqinr R package is used.

#' @examples # Load the data table exported by the SpeciesGeneMat.Bl function
#' data(Seq.DF4) # the table of sequences and metadata is called "CleanDataTable"
#' # and is the fifth object of the list "Seq.DF4".
#'
#' # Run the function
#' Seq.DF5=Select.DNA(input = Seq.DF4$CleanDataTable, gene.list = c("co1", "16srrna"),
#' output = "Seq.DF4.dataTable")
#'

#' @export Select.DNA

Select.DNA = function(input = NULL, gene.list = NULL, output = NULL, timeout = 10, Seqinr = FALSE) {

    # Check the dimensions (number of columns) of the input table, if the input table
    # also includes the sequence order in the 29th column, this last column is
    # removed.
    if (dim(input)[2] == 30)
        input = input[, -29]

    # Select all the rows with all the sequences of the gene.list in the table.
    r.pos = vector()
    i = 1
    for (i in 1:length(gene.list)) {
        if(gene.list[i] == "mll") {
            # In the case that the 'mll' gene is in the gene.list, ensure that only this one
            # is called and not 'mll2' or 'mll4'
            r.pos = c(r.pos, which(input[, 29] == gene.list[i]))
        } else {
            if(gene.list[i] == "h3") {
                # In the case that the 'h3' gene is in the gene.list, ensure that only this one
                # is called and not gnrh3 or sh3px3.
                r.pos = c(r.pos, which(input[, 29] == gene.list[i]))
            } else {
                r.pos = c(r.pos, grep(gene.list[i], as.character(input[, 29]), fixed = T))
            }  # End else if(gene.list[i]=='h3')
        }  # End else gene.list[i]=='mll')
    }  # End for i
    r.pos.uni = unique(r.pos)

    # Suitable selection
    ShortTab = input[r.pos.uni, ]

    # Extract the accession numbers of the sequences called 'TooLong' in the field
    # 'Sequence'.
    Ori.multSeq.po = grep("Too_Long", ShortTab[, 4])  # Identify the position of the sequence with multiple genes associated in the ShortTab.

    if(Seqinr == TRUE){
    
    if (length(Ori.multSeq.po) > 0)
        {
            # If at least one sequence need to be retrieved from NCBI.

            # Homogenize the gene names according to NCBI requirements for the request.
            gene.listB = toupper(gene.list)
            gene.listB = gsub("CO1", "COI", gene.listB, fixed = TRUE)
            gene.listB = gsub("SRRNA", "S ribosomal RNA", gene.listB, fixed = TRUE)

            # Connect to the Genbank database.
            seqinr::choosebank("genbank", timeout = timeout)
            DFadd = matrix(NA, ncol = dim(ShortTab)[2])[-1, ]
            i = 1
            for (i in 1:length(gene.list)) {
                # For multiple genes of interest.
                j = 1
                Addtab = as.matrix(ShortTab[Ori.multSeq.po, ])
                for (j in 1:length(Ori.multSeq.po)) {
                  # For all the multiple sequences.
                  seqinr::autosocket()
                  ee = tryCatch(seqinr::query("ee", paste("AC=", ShortTab[Ori.multSeq.po[j], 3],
                    " and k=", gene.listB[i], sep = "")), error = function(e) e)
                  if (is.null(ee$nelem) == "TRUE")
                    ese = 0 else if (ee$nelem == 0)
                    ese = 0 else ese = 1

                  if (ese == 0) {
                    # If the first trial failed, it can be due to a wrong CO1 syntax.
                    if (gene.listB[i] == "COI") {
                      ee = tryCatch(seqinr::query("ee", paste("AC=", ShortTab[Ori.multSeq.po[j],
                        3], " and k=", "COX1", sep = "")), error = function(e) e)  # 2nd CO1 syntax option.
                      if (is.null(ee$nelem) == "TRUE")
                        ese = 0 else if (ee$nelem == 0)
                        ese = 0 else ese = 1
                      if (ese == 0) {
                        ee = tryCatch(seqinr::query("ee", paste("AC=", ShortTab[Ori.multSeq.po[j],
                          3], " and k=", "CO1", sep = "")), error = function(e) e)  # 3rd CO1 syntax option.
                        if (is.null(ee$nelem) == "TRUE")
                          ese = 0 else if (ee$nelem == 0)
                          ese = 0 else ese = 1
                      }
                      if (ese == 0) {
                        ee = tryCatch(seqinr::query("ee", paste("AC=", ShortTab[Ori.multSeq.po[j],
                          3], " and k=", "COXI", sep = "")), error = function(e) e)  # 4th CO1 syntax option.
                        if (is.null(ee$nelem) == "TRUE")
                          ese = 0 else if (ee$nelem == 0)
                          ese = 0 else ese = 1
                      }
                      if (ese == 0) {
                        ee = tryCatch(seqinr::query("ee", paste("AC=", ShortTab[Ori.multSeq.po[j],
                          3], " and k=", "cytochrome c oxidase subunit I", sep = "")),
                          error = function(e) e)  # 5th CO1 syntax option.
                        if (is.null(ee$nelem) == "TRUE")
                          ese = 0 else if (ee$nelem == 0)
                          ese = 0 else ese = 1
                      }
                      if (ese == 0) {
                        ee = tryCatch(seqinr::query("ee", paste("AC=", ShortTab[Ori.multSeq.po[j],
                          3], " and k=", "cytochrome oxidase subunit I", sep = "")),
                          error = function(e) e)  # 6th CO1 syntax option.
                        if (is.null(ee$nelem) == "TRUE")
                          ese = 0 else if (ee$nelem == 0)
                          ese = 0 else ese = 1
                      }


                      if (ese == 1) {
                        Addtab[j, 4] = paste(seqinr::getSequence(ee$req)[[1]], collapse = "")
                        Addtab[j, 3] = ee$req[[1]][1]  # Accession number.
                        Addtab[j, 5] = seqinr::getLength(ee$req)  # Sequence length.
                        Addtab[j, 18] = gene.listB[i]
                        Addtab[j, 19] = gene.listB[i]
                        Addtab[j, 29] = gene.list[i]
                      } else {
                        # End if(ese==1) within COI==TRUE
                        Addtab[j, 4] = NA
                        Addtab[j, 3] = NA
                        Addtab[j, 5] = NA
                      }  # End when ese==0 when CO1==TRUE
                    } else {
                      #
                      Addtab[j, 4] = NA
                      Addtab[j, 3] = NA
                      Addtab[j, 5] = NA
                    }  # End if(gene.listB[i]=='CO1'){ else
                  } else {
                    # End if ese==0
                    Addtab[j, 4] = paste(seqinr::getSequence(ee$req)[[1]], collapse = "")
                    Addtab[j, 3] = ee$req[[1]][1]  # Accession number.
                    Addtab[j, 5] = seqinr::getLength(ee$req)  # Sequence length.
                    Addtab[j, 18] = gene.listB[i]
                    Addtab[j, 19] = gene.listB[i]
                    Addtab[j, 29] = gene.list[i]
                  }  # End if(ese==1)
                }  # End for j
                Addtab2 = Addtab[which(is.na(as.numeric(as.character(Addtab[, 5]))) ==
                  "FALSE"), ]
                DFadd = rbind(DFadd, Addtab2)
            }  # End for i

            ShortTab2 = ShortTab[-Ori.multSeq.po, ]
            ShortTab = rbind(ShortTab2, DFadd)
            CheckOK = tryCatch(seqinr::closebank(), error = function(e) e)
            if(inherits(CheckOK, "error")){
              stop(paste("The function could not retrieve the desired sequence(s), because the default ", "\n",
                         " time-lag (timeout = 10 sec.) to connect to the seqinr server is too short:", "\n",
                         " increasing the timeout to 15 or 20 seconds should solve the problem.", sep = ""))
            }

        }  # End if(length(Ori.multSeq.po)>0){
    } else {
      if (length(Ori.multSeq.po) > 0){
        # If at least one sequence need to be retrieved from NCBI.
        
        # Homogenize the gene names according to NCBI requirements for the request.
        gene.listB = toupper(gene.list)
        gene.listB = gsub("CO1", "COI", gene.listB, fixed = TRUE)
        gene.listB = gsub("SRRNA", "S ribosomal RNA", gene.listB, fixed = TRUE)
        
        
        # Connect to the Genbank database.
        DFadd = matrix(NA, ncol = dim(ShortTab)[2])[-1, ]
        i = 1
        for (i in 1:length(gene.list)) {
          # For multiple genes of interest.
          j = 1
          Addtab = as.matrix(ShortTab[Ori.multSeq.po, ])
          for (j in 1:length(Ori.multSeq.po)) {
            # For all the multiple sequences.
            
            # new way to get the sequence 
            oo = rentrez::entrez_search(db="nucleotide", 
                                        term =paste(ShortTab[Ori.multSeq.po[j], 3], "[ACCESSION]", sep = ""), 
                                        use_history= TRUE)
            
            #oo = rentrez::entrez_search(db="nucleotide", 
            #term =paste("ON838225", "[ACCESSION]", sep = ""), 
            #use_history= TRUE)
            
            
            
            # Fetch the sequences for records in a chunk as a Genbank XML
            oo.seqs = rentrez::entrez_fetch(db='nucleotide',
                                            web_history = oo$web_history, 
                                            rettype = 'gb', retmode = 'xml',
                                            parsed = TRUE)
            oo.seqs.xml = XML::xmlToList(oo.seqs)
            
            posF = which(names(oo.seqs.xml[[1]]) == "GBSeq_feature-table")
            posSeq = which(names(oo.seqs.xml[[1]]) == "GBSeq_sequence")
            # Build an annotaed table of the gene sequence 
            Annotated.DF = unique(as.data.frame(do.call(rbind, 
                                                        lapply(1:length(oo.seqs.xml[[1]][[posF]]), function(x)
                                                        {
                                                          
                                                          feat1 = oo.seqs.xml[[1]][[posF]][[x]]
                                                          a = feat1$GBFeature_intervals$GBInterval$GBInterval_from # starting Position
                                                          b = feat1$GBFeature_intervals$GBInterval$GBInterval_to # ending position
                                                          
                                                          c = unlist(lapply(1:length(feat1$GBFeature_quals), function(i){
                                                            feat1$GBFeature_quals[[i]][[2]]}))
                                                          
                                                          if(length(c) == 1){ # x = 75
                                                            c1 = c(c, NA)
                                                          }
                                                          
                                                          if(length(c) == 3){ # x = 3
                                                            c1 = c[c(1,2)]
                                                          }
                                                          
                                                          if(length(c) == 5){ # x = 1
                                                            c1 = c[c(1,4)]
                                                          }
                                                          
                                                          if(length(c) == 6){ # x = 72
                                                            c1 = c[c(1,4)]
                                                          }
                                                          
                                                          if(length(c) == 8){ #x = 76
                                                            c1 = c[c(1,6)]
                                                          }
                                                          
                                                          
                                                          if(length(c) == 2 | length(c) == 4 | length(c) == 7 | length(c) > 8){
                                                            c1 = c(NA, NA)
                                                          }
                                                          
                                                          c(a,b,c1)
                                                        }))))
            
            Annotated.DF = Annotated.DF[-1,]
            
            
            if (gene.listB[i] == "COI"){
              
              select = which(Annotated.DF[,3] == gene.listB[i] | 
                               Annotated.DF[,4] == gene.listB[i]) 
              
              if(length(select) == 0){
                select = which(Annotated.DF[,3] == "CO1" | Annotated.DF[,4] == "CO1") 
              }
              
              if(length(select) == 0){
                select = which(Annotated.DF[,3] == "COX1" | Annotated.DF[,4] == "COX1") 
              }
              
              if(length(select) == 0){
                select = which(Annotated.DF[,3] == "COXI" | Annotated.DF[,4] == "COXI") 
              }
              
              if(length(select) == 0){
                select = which(Annotated.DF[,3] == "cytochrome c oxidase subunit I" |
                                 Annotated.DF[,4] == "cytochrome c oxidase subunit I") 
              }
              
              if(length(select) == 0){
                select = which(Annotated.DF[,3] == "cytochrome oxidase subunit I" | 
                                 Annotated.DF[,4] == "cytochrome oxidase subunit I") 
              } 
              
              ### for other gene than COI
            } else {
              
              select = which(Annotated.DF[,3] == gene.listB[i] | 
                               Annotated.DF[,4] == gene.listB[i]) 
              
            }
            
            if(length(select) > 0){
              Boundary.From = as.numeric(Annotated.DF[select[1], 1])
              Boundary.To = as.numeric(Annotated.DF[select[1], 2])
              if(Boundary.From < Boundary.To){
                dna1 = substring(oo.seqs.xml[[1]][[posSeq]], Boundary.From, Boundary.To)
              } else {
                # in case the start and end of the sequence are inverted 
                dna1 = substring(oo.seqs.xml[[1]][[posSeq]], Boundary.To, Boundary.From)
              }
              
              Addtab[j, 4] = dna1
              Addtab[j, 3] = ShortTab[Ori.multSeq.po[j], 3]  # Accession number.
              Addtab[j, 5] = nchar(dna1)  # Sequence length.
              Addtab[j, 18] = gene.listB[i]
              Addtab[j, 19] = gene.listB[i]
              Addtab[j, 29] = gene.list[i]
              
            } else {
              Addtab[j, 4] = NA
              Addtab[j, 3] = NA
              Addtab[j, 5] = NA
              
            }
          } # end fo j
          Addtab2 = Addtab[which(is.na(as.numeric(as.character(Addtab[, 5]))) ==
                                   "FALSE"), ]
          DFadd = rbind(DFadd, Addtab2)
        } # end for i
        
        ShortTab2 = ShortTab[-Ori.multSeq.po, ]
        ShortTab = rbind(ShortTab2, DFadd)
        
      }
      
    }

    utils::write.table(ShortTab, file = paste(output, ".Select.DNA.txt", sep = ""), sep = "\t",
        row.names = FALSE, quote = FALSE)
    return(ShortTab)

}  # End function.
