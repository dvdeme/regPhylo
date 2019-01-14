#' @title Extract all the seqeunces and metadata from the data table for a list of selected gene regions

#' @description This function selects the sequences and associated information for a list
#' of selected gene regions for all the species present in the data table. For
#' full mitochondrial genomes, or all sequences with a sequence length > 5000 bp
#' listed as 'TooLong', the function extracts the sequences regions of interest
#'  using the query function of the seqinr R package.
#'

#' @param input The clean table coming from the function SpeciesGeneMat_Bl.R. which is
#' the output with the suffix '_CleanDataset.txt'.
#' @param gene.list a vector of the selected gene regions (gene names have
#' to be consistent with the header of the table with the suffix '_SpDNA_Mat.txt'
#' exported by the function SpeciesGeneMat_Bl.R). For example (Cytochrome c oxydase subunit 1, should be
#' written 'co1' and not 'COI' or 'COX1').
#' @param output the name of the output table exported in the working directory.

#' @examples # Load the data table exported by the SpeciesGeneMat.Bl function
#' data(Seq.DF4) # the clean table of sequences and metadata is called "CleanDataTable"
#' # and is the fifth object of the list "Seq.DF4".
#'
#' # Run the function
#' Seq.DF5=Select.DNA(input = Seq.DF4$CleanDataTable, gene.list = c("co1", "16srrna"),
#' output = "Seq.DF4.dataTable")

#' @export Select.DNA

Select.DNA = function(input = NULL, gene.list = NULL, output = NULL) {

    # Check the dimensions (number of columns) of the input table, if the input table
    # also includes the sequence order in the 29th column, this last column is
    # removed.
    if (dim(input)[2] == 30)
        input = input[, -29]

    # Select all the rows with all the sequences of the gene.list in the table.
    r.pos = vector()
    i = 1
    for (i in 1:length(gene.list)) {
        if (gene.list[i] == "mll") {
            # In the case that the 'mll' gene is in the gene.list, ensure that only this one
            # is called and not 'mll2' or 'mll4'
            r.pos = c(r.pos, which(input[, 29] == gene.list[i]))
        } else {
            if (gene.list[i] == "h3") {
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


    if (length(Ori.multSeq.po) > 0)
        {
            # If at least one sequence need to be retrieved from NCBI.

            # Homogenize the gene names according to NCBI requirements for the request.
            gene.listB = toupper(gene.list)
            gene.listB = gsub("CO1", "COI", gene.listB, fixed = TRUE)
            gene.listB = gsub("SRRNA", "S ribosomal RNA", gene.listB, fixed = TRUE)

            # Connect to the Genbank database.
            seqinr::choosebank("genbank")
            DFadd = matrix(NA, ncol = dim(ShortTab)[2])[-1, ]
            i = 1
            for (i in 1:length(gene.list)) {
                # For multiple genes of interest.
                j = 1
                Addtab = as.matrix(ShortTab[Ori.multSeq.po, ])
                for (j in 1:length(Ori.multSeq.po)) {
                  # For all the multiple sequences.
                  seqinr::autosocket()
                  ee = tryCatch(seqinr::query(paste("AC=", ShortTab[Ori.multSeq.po[j], 3],
                    " et k=", gene.listB[i], sep = "")), error = function(e) e)
                  if (is.null(ee$nelem) == "TRUE")
                    ese = 0 else if (ee$nelem == 0)
                    ese = 0 else ese = 1

                  if (ese == 0) {
                    # If the first trial failed, it can be due to a wrong CO1 syntax.
                    if (gene.listB[i] == "COI") {
                      ee = tryCatch(seqinr::query(paste("AC=", ShortTab[Ori.multSeq.po[j],
                        3], " et k=", "COX1", sep = "")), error = function(e) e)  # 2nd CO1 syntax option.
                      if (is.null(ee$nelem) == "TRUE")
                        ese = 0 else if (ee$nelem == 0)
                        ese = 0 else ese = 1
                      if (ese == 0) {
                        ee = tryCatch(seqinr::query(paste("AC=", ShortTab[Ori.multSeq.po[j],
                          3], " et k=", "CO1", sep = "")), error = function(e) e)  # 3rd CO1 syntax option.
                        if (is.null(ee$nelem) == "TRUE")
                          ese = 0 else if (ee$nelem == 0)
                          ese = 0 else ese = 1
                      }
                      if (ese == 0) {
                        ee = tryCatch(seqinr::query(paste("AC=", ShortTab[Ori.multSeq.po[j],
                          3], " et k=", "COXI", sep = "")), error = function(e) e)  # 4th CO1 syntax option.
                        if (is.null(ee$nelem) == "TRUE")
                          ese = 0 else if (ee$nelem == 0)
                          ese = 0 else ese = 1
                      }
                      if (ese == 0) {
                        ee = tryCatch(seqinr::query(paste("AC=", ShortTab[Ori.multSeq.po[j],
                          3], " et k=", "cytochrome c oxidase subunit I", sep = "")),
                          error = function(e) e)  # 5th CO1 syntax option.
                        if (is.null(ee$nelem) == "TRUE")
                          ese = 0 else if (ee$nelem == 0)
                          ese = 0 else ese = 1
                      }
                      if (ese == 0) {
                        ee = tryCatch(seqinr::query(paste("AC=", ShortTab[Ori.multSeq.po[j],
                          3], " et k=", "cytochrome oxidase subunit I", sep = "")),
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

        }  # End if(length(Ori.multSeq.po)>0){

    utils::write.table(ShortTab, file = paste(output, ".Select.DNA.txt", sep = ""), sep = "\t",
        row.names = FALSE, quote = FALSE)
    return(ShortTab)
}  # End function.
