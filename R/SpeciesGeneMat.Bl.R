#' @title Clean gene nomenclature and build a species-by-gene matrix

#' @description This function cleans the gene names with a special focus on 37 DNA
#' markers commonly used in if fish molecular phylogenies
#' and then build a species-by-gene matrix, and export a summary tables from a DNA
#' and species perspectives reporting the species or the gene coverage respectively.
#' The function also exports a clean table of the dataset
#' reporting the clean name for the gene, and also the original gene name, in
#' addition the 'blacklisted' sequences, the microsatellites, and unassigned DNA
#' are also removed from this table.
#'
#' @details The 37 gene regions targeted in priority are: 12S, 16S, COI,
#' CytB, ENC1, Glyt (=gtdc2), myh6, plagl2, Ptr (=ptchd1), rag1, SH3PX3 (=snx33),
#' sreb2 (=gpr85), tbr1 (=tbr1b), zic1, rhodopsin (=rh), tmo4c4, bmp4, rag2, irbp,
#' egr1, mll, ficd, hoxc6a, Kiaa1239, panx2, ripk4, sidkey, vcpip, RNF213,
#' SLC10A3, ube3A-like, ZNF503, ZNF536, EHHADH, GCS1, GPR61, ube3A, that are
#' commonly used in fish molecular phylogenies (Li et al. 2007; Santini et al.
#' 2009, 2013; Near et al. 2012, 2013; Betancur et al. 2013, 2015; Rabosky et al.
#' 2013, 2018).

#' @return The function computes and exports into the working directory
#' the following object: 1) a species-by-gene matrix (with the number of sequences
#' for each species and each DNA type also called gene or gene region for
#' simplification). The order of the genes and species in the matrix is as follow:
#' it starts from the gene with the maximum species coverage and finish with the
#' gene with the minimum species coverage (genes in columns); for the species (in
#' rows), it starts with the species with the minimal gene coverage and finishes
#' with the species with the maximum gene coverage.  2) a table with the number of
#' sequences for each gene and number of species for each gene (gene perspective).
#' 3) a table with the number of sequences for each species and number of genes
#' for each species (species perspective).  4) a clean table of the dataset
#' reporting the clean name for the gene, and also the original gene name, in
#' addition the 'blacklisted' sequences, the microsatellites, and unassigned DNA
#' are also removed from this table.  The function also returns a list of objects
#' in the R environment, the table 1, 2, and 3 and the last object is a list of
#' species that have been removed after cleaning the dataset.

#' @param input a table which is an output of the GeoCodeName() function
#' @param output a list of 4 tables exported in the working directory 1) the species-by-gene
#' matrix (file name suffix = '_Mat.txt'), 2) the table from the gene perspective
#' (file name suffix = '_DNApers.txt'), 3) the table from the species perspective
#' (file name suffix = '_SPECIESpers.txt', 4) a table with clean gene names or
# discarding the blacklisted sequences and any microsatellites and unassigned DNA
# (file name suffix = 'CleanDataset.txt').
#' @param NCBI.Trash a vector of the NCBI accession numbers that are
#' blacklisted and need to be removed from the dataset.
#' @param BOLD.Trash a vector of the BOLD 'sequenceID' (or 'recordID') that are
#' blacklisted and need to be removed from the dataset.

#' @examples # Load a table exported by GeoCodeName()
#' data(Seq.Diastocapen)
#' Seq.DF3 = Seq.Diastocapen$Seq.DF3
#'
#' # Run the function
#' Seq.DF4=SpeciesGeneMat.Bl(input=Seq.DF3, output="Seq.DF4.")


#' @export SpeciesGeneMat.Bl


SpeciesGeneMat.Bl = function(input = NULL, output = NULL, NCBI.Trash = NULL, BOLD.Trash = NULL) {

    # Species list before cleaning the sequences.
    list1 = unique(as.character(input[, 1]))

    ## Removing the external blacklisted sequences for NCBI
    if (is.null(NCBI.Trash) == FALSE) {
        if (length(match(intersect(unique(NCBI.Trash), as.character(input[, 3])),
            as.character(input[, 3]))) > 0) {
            input = input[-match(intersect(unique(NCBI.Trash), as.character(input[,
                3])), as.character(input[, 3])), ]
        }
    }
    ## BOLD
    if (is.null(BOLD.Trash) == FALSE) {
        if (length(match(intersect(unique(BOLD.Trash), as.character(input[, 2])),
            as.character(input[, 2]))) > 0) {
            input = input[-match(intersect(unique(BOLD.Trash), as.character(input[,
                2])), as.character(input[, 2])), ]
        }
    }

    if (length(grep("ProductGeneClean", colnames(input))) == 0) {
        # In case the cleaning of gene region names has to be done.

        ## Correction for mll gene.  Distinguishing the mll1 and mll2.
        mll1 = c("AY363630", "AY363653", "AF137234", "AY363666", "AY363658", "AF137241",
            "AY363657", "AY363659", "AF137244", "AY363638", "AY363656", "EF589659")
        Genes = as.character(input[, 19])
        Genes[match(mll1, input[, 3])] <- "mll1"

        mll4 = c("HM071049", "HM071050")
        Genes[match(mll4, input[, 3])] <- "mll4"
        # MLL2 called in Dettai and Lecointre 2005, has been changed to MLL4, see
        # Lautredou et al. 2013, and other copies of the MLL gene has been called MLL2 in
        # 2013 by Lautredou et al. 2013. The huge majority of these sequences belong to
        # MLL4.

        # Add the sequence order.
        input = cbind(input[, c(1:18)], Genes, input[, c(20:dim(input)[2])], SeqOrder = seq(1,
            dim(input)[1], by = 1))

        # Remove rows without DNA sequences attached (sometimes happens in the BOLD
        # database)
        if (length(grep("^$", input[, 4], perl = T)) > 0) {
            input = input[-grep("^$", input[, 4], perl = T), ]
        }

        # Find all the sequences with no (i.e. 'NA') DNA annotation in the 'Product' and
        # 'gene' fields.
        input2 = input[intersect(which(is.na(input[, 18]) == "TRUE"), which(is.na(input[,
            19]) == "TRUE")), ]
        # Remove all the microsatellites.
        if (length(grep("[Mm]icrosatellite", as.character(input2[, 6]), perl = TRUE)) >
            0) {
            input3 = input2[-grep("[Mm]icrosatellite", as.character(input2[, 6]),
                perl = TRUE), ]
        } else {
            input3 = input2
        }

        # Remove the 'unassigned DNA', in Mol_Type field
        if (length(grep("unassigned DNA", as.character(input3[, 16]), fixed = TRUE)) >
            0) {
            input3 = input3[-grep("unassigned DNA", as.character(input3[, 16]), fixed = TRUE),
                ]
        } else {
            input3 = input3
        }

        # When the DNA annotation is only present in the field called 'Definition',
        # extract the name of the DNA fragment contained the field 'Definition'

        # Clean some inconsistent wording in the field 'Definition'
        HHH = gsub("D-loop", "D-loop", input3[, 6], perl = TRUE)
        HHH = gsub("control[ -]region", "D-loop", HHH, perl = TRUE)
        HHH = gsub("mitochondrial genes for tRNA-Pro", "tRNA-Pro", HHH, fixed = TRUE)
        HHH = gsub("[Cc]ytochrome [c]?[ ]?ox[yi]dase subunit I", "co1", HHH, perl = TRUE)
        HHH = gsub("[ ]?ribosomal", "rRNA", HHH, perl = TRUE)
        HHH = gsub("[Aa][Tt][Pp]ase subunit 8-like gene", "atp8", HHH, perl = TRUE)
        HHH = gsub("[Aa][Tt][Pp]ase subunit 6-like gene", "atp6", HHH, perl = TRUE)
        HHH = gsub("ptr[ -]like", "ptr-like", HHH, perl = TRUE)
        HHH = gsub("[Ii]nternal [Tt]ranscribed [Ss]pacer[ ]?", "its", HHH, perl = TRUE)
        HHH = gsub("[Cc][iy]tochrome b", "cytb", HHH, perl = TRUE)
        HHH = gsub("histone[ -]?III", "histone 3", HHH, perl = TRUE)
        HHH = gsub("ultra conserved element locus", "ucel", HHH, perl = TRUE)
        HHH = gsub("[Bb]eta[ -]actin", "beta-actin", HHH, perl = TRUE)
        HHH = gsub("[Mm]exico1", "", HHH, perl = TRUE)
        HHH = tolower(HHH)

        # Molecule names to check (usually the name included in the field 'Definition'
        # when nothing is provided in the fields 'Product' or 'Gene').  Ucel is for ultra
        # conserved element locus.
        Tocheck = c("ucel", "d-loop", "12srrna", "18srrna", "16srrna", "5srrna",
            "5.8srrna", "atp8", "atp6", "its1", "its2", "histone 3", "ptr", "cytb",
            "co1", "tRNA-pro", "tRNA-thr", "tRNA-phe", "beta-actin", "creatine kinase",
            "non-LTR retrotransposon")

        DNAfrag = matrix(NA, ncol = 1)[-1, ]
        i = 1
        for (i in 1:length(HHH)) {
            DNAfragi = vector()
            j = 1
            for (j in 1:length(Tocheck)) {
                if (length(grep(Tocheck[j], HHH[i])) > 0) {
                  DNAfragi = paste(DNAfragi, Tocheck[j], sep = "; ")
                }
            }  # End for j
            if (length(DNAfragi) > 0) {
                DNAfrag[i] = DNAfragi
            } else {
                DNAfrag[i] = NA
            }
        }  # End for i

        # Copy and paste the field 'Genes' into the object 'product' when the field
        # 'Product' is NA.
        product = as.character(input[, 18])
        product[which(is.na(input[, 18]) == "TRUE")] = as.character(input[which(is.na(input[,
            18]) == "TRUE"), 19])


        options(warn = -1)  # Remove the warnings that appear if DNAfrag has a length greater to one, we can use length because we report the data in a matrix with a single column only.
        if (length(DNAfrag) == 1 & is.na(DNAfrag) == TRUE) {
            product = product
        } else {
            # Removed the '; ' in front of each row.
            DNAfrag = gsub("^; ", "", DNAfrag, perl = TRUE)
            input4 = cbind(input3, DNAfrag)
            # Copy and paste the name of the DNA fragment contained in the field 'Definition'
            # when the fields 'Genes' and 'Product' are empty ('NA')
            product[as.numeric(as.character(input4[, 29]))] = as.character(input4[,
                30])
        }
        options(warn = 0)  # Remove 'disable' warnings.


        # Fill the overall table.
        inputb = cbind(input, product_All = product)

        # Remove the microsatellites and the unassigned DNA or other ambiguous or very
        # rare DNA fragments in the overall table
        if (length(which(is.na(product) == "TRUE")) > 0) {
            inputCl = inputb[-which(is.na(product) == "TRUE"), ]
        } else {
            inputCl = inputb
        }

        # Clean the field called 'Genes' abbreviation terminology
        GeneClean = gsub("COI-5P", "co1", inputCl[, 19], fixed = TRUE)
        GeneClean = gsub("COI-3P", "co1", GeneClean, fixed = TRUE)
        GeneClean = gsub("III", "3", GeneClean, fixed = TRUE)
        GeneClean = gsub("II", "2", GeneClean, fixed = TRUE)
        GeneClean = gsub("IV", "4", GeneClean, fixed = TRUE)
        GeneClean = gsub("COI", "co1", GeneClean, fixed = TRUE)
        GeneClean = tolower(GeneClean)  # Convert all the gene abbreviation in lowercase.
        GeneClean = gsub("[ ]", "", GeneClean, perl = TRUE)  # Remove any space.
        GeneClean = gsub("cox", "co", GeneClean, fixed = TRUE)
        GeneClean = gsub("coi", "co1", GeneClean, fixed = TRUE)
        GeneClean = gsub("tmo[-]?4c4", "tmo4c4", GeneClean, perl = TRUE)
        GeneClean = gsub("rag[ -]?", "rag", GeneClean, perl = TRUE)
        GeneClean = gsub("cyt+[ ]?b", "cytb", GeneClean, perl = TRUE)
        GeneClean = gsub("^rh[1]?[a]?", "rhod", GeneClean, perl = TRUE)
        GeneClean = gsub("rhodo[d]?", "rhod", GeneClean, perl = TRUE)
        GeneClean = gsub("rhod2", "rh2", GeneClean, perl = TRUE)
        GeneClean = gsub("gnrhod3-3", "gnrh3-3", GeneClean, perl = TRUE)
        GeneClean = gsub("gnrhod", "gnrh", GeneClean, perl = TRUE)
        GeneClean = gsub("sbgnrhod", "sbgnrh", GeneClean, perl = TRUE)
        GeneClean = gsub("atpase", "atp", GeneClean, perl = TRUE)
        GeneClean = gsub("snx33", "sh3px3", GeneClean, perl = TRUE)
        GeneClean = gsub("sreb$", "sreb2", GeneClean, perl = TRUE)
        GeneClean = gsub("ptchd1", "ptr", GeneClean, perl = TRUE)
        GeneClean = gsub("gpr85", "sreb2", GeneClean, perl = TRUE)
        GeneClean = gsub("tbr1b", "tbr1", GeneClean, perl = TRUE)
        GeneClean = gsub("gtdc2", "glyt", GeneClean, perl = TRUE)
        GeneClean = gsub("s7ribosomalprotein", "s7", GeneClean, perl = TRUE)
        GeneClean = gsub("mot4c4", "tmo4c4", GeneClean, perl = TRUE)
        GeneClean = gsub("co1 (co1)", "co1", GeneClean, fixed = TRUE)
        GeneClean = gsub("nadh", "nd", GeneClean, perl = TRUE)
        GeneClean = gsub("ets-", "ets", GeneClean, fixed = TRUE)
        GeneClean = gsub("egr-", "egr", GeneClean, fixed = TRUE)
        GeneClean = gsub("bmp[-]?4", "bmp4", GeneClean, perl = TRUE)
        GeneClean = gsub("^at8$", "atp8", GeneClean, perl = TRUE)
        GeneClean = gsub("^at6$", "atp6", GeneClean, perl = TRUE)
        GeneClean = gsub(";at8", "atp8", GeneClean, fixed = TRUE)
        GeneClean = gsub(";at6", "atp6", GeneClean, fixed = TRUE)
        GeneClean = gsub("rps7", "s7", GeneClean, perl = TRUE)
        GeneClean = gsub("mll$", "mll4", GeneClean, perl = TRUE)
        GeneClean = gsub("mll41", "mll1", GeneClean, perl = TRUE)
        GeneClean = gsub("mll42", "mll2", GeneClean, perl = TRUE)
        GeneClean = gsub("irbp", "irbp2", GeneClean, perl = TRUE)  # The second copy is present in all fishes whereas the first copy is present in a very limited number of species, in the absence of more precise information, we assumed that irbp belonged to the second copy (see Dettai & Lecointre 2008 MPE and Chen et al. 2008 Gene).

        ## Clean the name of the product_All (include the fields 'Product' and
        ## 'Definition').
        ProductCl = gsub("III", "3", inputCl[, 30], fixed = TRUE)
        ProductCl = gsub("II", "2", ProductCl, fixed = TRUE)
        ProductCl = gsub("IV", "4", ProductCl, fixed = TRUE)
        ProductCl = gsub(" idase ", "", ProductCl, fixed = TRUE)
        ProductCl = gsub("[ ]?subunit I", "1", ProductCl, fixed = TRUE)
        ProductCl = gsub("[ ]?[Cc]ytochrome[ ]?[c]?[ ]?ox[yi]dase subunit[ ]?", "co",
            ProductCl, perl = TRUE)
        ProductCl = gsub("[ ]?[Cc]ytochrome[ ]?[c]?[ ]?ox[yi]dase 1", "co1", ProductCl,
            perl = TRUE)
        ProductCl = gsub("COI", "co1", ProductCl, perl = TRUE)
        ProductCl = tolower(ProductCl)  # convert all letter in lowercase.

        # Exclude all the rows with multiple genes in the DNA fragment from the automatic
        # cleaning.
        All.pos = seq(1, dim(inputCl)[1], by = 1)
        W.pos = All.pos[-grep(";", GeneClean)]

        if (length(W.pos) > 0)
            {
                # Automatic cleaning for the single gene per row.
                i = 1
                for (i in 1:length(W.pos)) {
                  if (is.na(GeneClean[W.pos[i]]) == "FALSE") {
                    if (is.na(inputCl[W.pos[i], 30]) == "FALSE") {
                      ProductCl[W.pos[i]] = gsub(ProductCl[W.pos[i]], GeneClean[W.pos[i]],
                        ProductCl[W.pos[i]], fixed = TRUE)
                    }
                  }
                }
            }  # end if(length(W.pos)>0){

        # Additional manual cleaning due to inconsistent wording of the gene.  Clean the
        # synonymies of the gene names in these fields.
        ProductClean = gsub("^[ ]?", "", ProductCl, perl = TRUE)
        ProductClean = gsub("coi", "co1", ProductClean, perl = TRUE)
        ProductClean = gsub("co1-5p", "co1", ProductClean, perl = TRUE)
        ProductClean = gsub("glucocortico1d", "glucocorticoid", ProductClean, perl = TRUE)
        ProductClean = gsub("[ ]?small subunit ribosomal rna", "rrna", ProductClean,
            perl = TRUE)
        ProductClean = gsub("[ ]?large subunit ribosomal rna", "rrna", ProductClean,
            perl = TRUE)
        ProductClean = gsub("[ ]?ribosomal rna subunit", "rrna", ProductClean, perl = TRUE)
        ProductClean = gsub("[ ]?ribosomal rna", "rrna", ProductClean, perl = TRUE)
        ProductClean = gsub("[ ]?rrna", "rrna", ProductClean, perl = TRUE)
        ProductClean = gsub("[ ]?small subunitrrna", "rrna", ProductClean, perl = TRUE)
        ProductClean = gsub("tmo+[-]?4c4", "tmo4c4", ProductClean, perl = TRUE)
        ProductClean = gsub("cyt[h]?oc[h]?rome[ ]?b", "cytb", ProductClean, perl = TRUE)
        ProductClean = gsub("nadh dehydrogenase subunit[ ]?", "nd", ProductClean,
            perl = TRUE)
        ProductClean = gsub("nadh dehydrogenase[ ]?", "nd", ProductClean, perl = TRUE)
        ProductClean = gsub("nadhdehydrogenasesubunit", "nd", ProductClean, perl = TRUE)
        ProductClean = gsub("nadh[ ]?", "nd", ProductClean, perl = TRUE)
        ProductClean = gsub("recombination[- ]?activating protein[ ]?", "rag", ProductClean,
            perl = TRUE)
        ProductClean = gsub("recombination[- ]?activation protein subunit[ ]?", "rag",
            ProductClean, perl = TRUE)
        ProductClean = gsub("recombinase activating gene[ ]?", "rag", ProductClean,
            perl = TRUE)
        ProductClean = gsub("recombinase activating protein[ ]?", "rag", ProductClean,
            perl = TRUE)
        ProductClean = gsub("zinc finger protein of the cerebellum 1", "zic1", ProductClean,
            perl = TRUE)
        ProductClean = gsub("zic family member 1", "zic1", ProductClean, perl = TRUE)
        ProductClean = gsub("zic family member protein 1", "zic1", ProductClean,
            perl = TRUE)
        ProductClean = gsub("internal transcribed spacer[ ]?", "its", ProductClean,
            perl = TRUE)
        ProductClean = gsub("r[h]?od[ -]?opsin", "rhod", ProductClean, perl = TRUE)
        ProductClean = gsub("rh1[ -]?opsin", "rhod", ProductClean, perl = TRUE)
        ProductClean = gsub("rh1[-]?a", "rhod", ProductClean, perl = TRUE)
        ProductClean = gsub("rh2a opsin", "rh2a", ProductClean, perl = TRUE)
        ProductClean = gsub("rhodo[d]?", "rhod", ProductClean, perl = TRUE)
        ProductClean = gsub("atp synthase [f]?[0o]?[ ]?subunit[ ]?", "atp", ProductClean,
            perl = TRUE)
        ProductClean = gsub("atp synthetase [f]?[0o]?[ ]?subunit[ ]?", "atp", ProductClean,
            perl = TRUE)
        ProductClean = gsub("atpase subunit[ ]?", "atp", ProductClean, perl = TRUE)
        ProductClean = gsub("atpase", "atp", ProductClean, perl = TRUE)
        ProductClean = gsub("atps[ ]?", "atp", ProductClean, perl = TRUE)
        ProductClean = gsub("atp synthase", "atp", ProductClean, perl = TRUE)
        ProductClean = gsub("atp[ ]?", "atp", ProductClean, perl = TRUE)
        ProductClean = gsub("adenosine triphos[p]?hatase subunit[ ]?", "atp", ProductClean,
            perl = TRUE)
        ProductClean = gsub("^at8$", "atp8", ProductClean, perl = TRUE)
        ProductClean = gsub("^at6$", "atp6", ProductClean, perl = TRUE)
        ProductClean = gsub("rag-", "rag", ProductClean, perl = TRUE)
        ProductClean = gsub("[c]?[a]?[r]?[d]?[i]?[a]?[c]?[ ]?[m]?[u]?[s]?[c]?[l]?[e]?[ ]?myosin heavy chain 6 alpha",
            "myh6", ProductClean, perl = TRUE)
        ProductClean = gsub("^pleiomorphic adenoma protein-like 2$", "plagl2", ProductClean,
            perl = TRUE)
        ProductClean = gsub("^pleiomorphic adenoma protein-like 2 protein$", "plagl2",
            ProductClean, perl = TRUE)
        ProductClean = gsub("^interphotoreceptor retinoid-binding protein$", "irbp2",
            ProductClean, perl = TRUE)  # The second copy is present in all fishes whereas the first copy is present in a very limited number of species, in the absence of more precise information, we assumed that irbp belonged to the second copy (see Dettai & Lecointre 2008 MPE and Chen et al. 2008 Gene).
        ProductClean = gsub("^mixed[ -]?lineage leukemia-like protein$", "mll4",
            ProductClean, perl = TRUE)
        ProductClean = gsub("mixed[ -]?lineage leukemia-like protein[ ]?", "mll",
            ProductClean, perl = TRUE)
        ProductClean = gsub("sh3 and px[3]? domain-containing 3-like protein", "sh3px3",
            ProductClean, perl = TRUE)
        ProductClean = gsub("sh3 and px domain-containing 3 protein", "sh3px3", ProductClean,
            perl = TRUE)
        ProductClean = gsub("ectodermal[ -]?neural cortex 1-like protein", "enc1",
            ProductClean, perl = TRUE)
        ProductClean = gsub("early growth response[ ]?", "egr", ProductClean, perl = TRUE)
        ProductClean = gsub("egr-", "egr", ProductClean, fixed = TRUE)
        ProductClean = gsub("t-box brain 1b protein", "tbr1", ProductClean, perl = TRUE)
        ProductClean = gsub("t[-]?box brain[ ]?[1]?$", "tbr1", ProductClean, perl = TRUE)
        ProductClean = gsub("receptor-interacting serine-threonine kinase[ ]?[4]?",
            "ripk4", ProductClean, perl = TRUE)
        ProductClean = gsub("glycosyltransferase", "glyt", ProductClean, perl = TRUE)
        ProductClean = gsub("si:dkey-174m14.3", "sidkey", ProductClean, perl = TRUE)
        ProductClean = gsub("brain super conserved receptor[ ],[2]?", "sreb2", ProductClean,
            perl = TRUE)
        ProductClean = gsub("myoglobin", "mb", ProductClean, perl = TRUE)
        ProductClean = gsub("leucine-rich repeat and wd repeat-containing", "kiaa1239",
            ProductClean, perl = TRUE)
        ProductClean = gsub("loc562320", "kiaa1239", ProductClean, perl = TRUE)
        ProductClean = gsub("ubiquitin protein ligase e3a-like protein", "ube3a",
            ProductClean, perl = TRUE)
        ProductClean = gsub("^ubiquitin protein ligase e3a$", "ube3a", ProductClean,
            perl = TRUE)
        ProductClean = gsub("^ring finger protein 213$", "rnf213", ProductClean,
            perl = TRUE)
        ProductClean = gsub("valosin-containing protein p97/p47 complete", "vcpip",
            ProductClean, perl = TRUE)
        ProductClean = gsub("zinc finger protein 503", "znf503", ProductClean, perl = TRUE)
        ProductClean = gsub("peroxisomal enoyl-coa", "ehhadh", ProductClean, perl = TRUE)
        ProductClean = gsub("g protein-coupled receptor 61", "gpr61", ProductClean,
            perl = TRUE)
        ProductClean = gsub("fic domain protein", "ficd", ProductClean, perl = TRUE)
        ProductClean = gsub("homeo box c6a", "hoxc6a", ProductClean, perl = TRUE)
        ProductClean = gsub("interleukin[- ]?1[ ]?beta", "il1-beta", ProductClean,
            perl = TRUE)
        ProductClean = gsub("interleukin[ -]?8", "il8", ProductClean, perl = TRUE)
        ProductClean = gsub("interleukin[ -]?10", "il10", ProductClean, perl = TRUE)
        ProductClean = gsub("ubiquitin protein ligase e3a-like protein", "ube3a",
            ProductClean, perl = TRUE)
        ProductClean = gsub("zinc finger protein 536", "znf536", ProductClean, perl = TRUE)
        ProductClean = gsub("control[ -]region", "d-loop", ProductClean, perl = TRUE)
        ProductClean = gsub("creatine[ -]kinase", "creatine kinase", ProductClean,
            perl = TRUE)
        ProductClean = gsub("prolactin receptor", "prolactin", ProductClean, perl = TRUE)
        ProductClean = gsub("s7ribosomalprotein", "s7", ProductClean, perl = TRUE)
        ProductClean = gsub("histone[ ]?3", "histone h3", ProductClean, perl = TRUE)
        ProductClean = gsub("ef1alpha", "ef1a", ProductClean, perl = TRUE)
        ProductClean = gsub("elongation factor 1-alpha", "ef1a", ProductClean, perl = TRUE)
        ProductClean = gsub("co1 (co1)", "co1", ProductClean, fixed = TRUE)
        ProductClean = gsub("cytochrome ox[yi]dase i", "co1", ProductClean, perl = TRUE)
        ProductClean = gsub("cytochrome ox[yi]dase[ ]?", "co", ProductClean, perl = TRUE)
        ProductClean = gsub("nwd2", "kiaa1239", ProductClean, perl = TRUE)
        ProductClean = gsub("bmp[-]?4", "bmp4", ProductClean, perl = TRUE)
        ProductClean = gsub("bone morphogenetic protein 4", "bmp4", ProductClean,
            fixed = TRUE)
        ProductClean = gsub("s7ribosomalprotein", "s7", ProductClean, perl = TRUE)
        ProductClean = gsub("ribosomal protein s7", "s7", ProductClean, perl = TRUE)
        ProductClean = gsub("mot4c4", "tmo4c4", ProductClean, perl = TRUE)
        ProductClean = gsub(",", ".", ProductClean, fixed = TRUE)
        ProductClean = gsub("ets[- ]2 oncogene", "ets2", ProductClean, perl = TRUE)
        ProductClean = gsub("ets[ -]?2", "ets2", ProductClean, perl = TRUE)
        ProductClean = gsub("rps7", "s7", ProductClean, perl = TRUE)
        ProductClean = gsub("histone h3", "h3", ProductClean, perl = TRUE)

        # Include ProductClean, GeneClean as new columns in the table.
        inputCl1 = cbind(inputCl, ProductClean, GeneClean)
        ProductGene_Clean = paste(inputCl1[, 31], inputCl1[, 32], sep = ";")  # Merge ProductClean and GeneClean in a single vector.

        # Remove the duplicates among ProductClean and GeneClean
        GeneNameHomog = matrix(NA, ncol = 1)[-1, ]
        i = 1
        a = gsub("; ", ";", as.character(ProductGene_Clean), fixed = TRUE)
        for (i in 1:dim(inputCl1)[1]) {
            GeneNameHomog[i] = paste(unique(unlist(strsplit(as.character(a[i]), ";"))),
                collapse = "; ")
        }
        GeneNameHomog = gsub("; NA", "", GeneNameHomog, fixed = TRUE)  # Removed the '; NA' pattern.
        GeneNameHomog = gsub("^rrna$", "rrna_unidentified", GeneNameHomog, perl = T)
        GeneNameHomog = gsub("^irbp22$", "irbp2", GeneNameHomog, perl = T)  # Correct some over-corrected gene names.
        GeneNameHomog = gsub("^irbp21$", "irbp1", GeneNameHomog, perl = T)  # Correct some over-corrected gene names, this one is the irbp1 present in minority of fish species see Dettai & Lecointre 2008.

        # Create the unique list of DNA fragments.
        DNAFragUniq = unique(unlist(strsplit(GeneNameHomog, "; ")))

        # A simplified and clean overall data table
        inputCl1 = cbind(inputCl1[, c(1:29)], ProductGeneClean = GeneNameHomog)


        # Check for the presence of NA in the field 'Sequence', if some exist, use seqinr
        # to extract them (although this should not have happened).
        tt = intersect(which(is.na(inputCl1[, 4])), which(inputCl1[, 28] == "NCBI"))
        if (length(tt) > 0)
            {
                seqinr::choosebank("genbank")
                SeqTocomplete = as.character(inputCl1[, 4])
                i = 1
                for (i in 1:length(tt)) {
                  seqinr::autosocket()
                  ee = tryCatch(seqinr::query(paste("AC=", as.character(inputCl1[tt[i], 3]),
                    sep = "")), error = function(e) e)
                  SeqTocomplete[tt[i]] = paste(seqinr::getSequence(ee$req)[[1]], collapse = "")
                }  # End for i
                seqinr::closebank()
                inputCl1 = cbind(inputCl1[, c(1:3)], Sequence = SeqTocomplete, inputCl1[,
                  c(5:30)])
            }  # End if(length(tt)>0).

    } else {
        # In case the cleaning of the gene region has already been done.
        inputCl1 = input
        DNAFragUniq = unique(unlist(strsplit(as.character(inputCl1[, grep("ProductGeneClean",
            colnames(inputCl1))]), "; ")))
    }


    # Report the name of potential species lost when removing the missing sequences.
    list2 = unique(as.character(inputCl1[, 1]))
    MissingSpecies = setdiff(list1, list2)

    ## Build the species-by-gene (DNA fragment= Product+Definition+Gene) matrix.
    sp = unique(as.character(inputCl1[, 1]))  # Species list.

    # Define the unique set of DNA fragment= unique homogenisation of
    # Product+Definition+Genes.
    SpeciesDNAFragMat = matrix(NA, ncol = length(DNAFragUniq), nrow = length(sp))
    i = 1
    for (i in 1:length(sp)) {
        DFtemp = inputCl1[which(inputCl1[, 1] == sp[i]), ]
        j = 1
        for (j in 1:length(DNAFragUniq)) {
            a = length(grep(DNAFragUniq[j], DFtemp[, grep("ProductGeneClean", colnames(DFtemp))],
                fixed = TRUE))  ## NB of sequences of the gene j of the species i.
            SpeciesDNAFragMat[i, j] = a
        }
    }
    colnames(SpeciesDNAFragMat) = DNAFragUniq
    rownames(SpeciesDNAFragMat) = sp


    # Transform the species-by-gene matrix into a presence/absence matrix for each
    # species.
    SpeciesDNAFragMatPA = SpeciesDNAFragMat
    SpeciesDNAFragMatPA[which(SpeciesDNAFragMatPA[, ] > 1)] <- 1

    if (dim(SpeciesDNAFragMatPA)[2] > 1) {
        # If there is more than one gene in the species-by-gene matrix order the original
        # species-by-gene matrix.
        SpeciesDNAFragMatTemp = rbind(SpeciesDNAFragMat, colSums(SpeciesDNAFragMatPA))
        SpeciesDNAFragMatTemp = cbind(SpeciesDNAFragMatTemp, c(rowSums(SpeciesDNAFragMatPA),
            1))
        rownames(SpeciesDNAFragMatTemp) = c(sp, "Cov.Gene")
        colnames(SpeciesDNAFragMatTemp) = c(DNAFragUniq, "Cov.Sp")
        # Order the matrix from the species with the least number of genes to the species with the maximum number of genes.
        SpeciesDNAFragMatTempOrd = SpeciesDNAFragMatTemp[order(SpeciesDNAFragMatTemp[,
            dim(SpeciesDNAFragMatTemp)[2]]), ][, -dim(SpeciesDNAFragMatTemp)[2]]
        tSpeciesDNAFragMatTempOrd = t(SpeciesDNAFragMatTempOrd)
        # Order the matrix from the gene with the maximum species coverage to the gene with the minimum species coverage.
        SpeciesDNAFragMat = t(tSpeciesDNAFragMatTempOrd[order(tSpeciesDNAFragMatTempOrd[,
            grep("Cov.Gene", colnames(tSpeciesDNAFragMatTempOrd))], decreasing = TRUE),
            ][, -grep("Cov.Gene", colnames(tSpeciesDNAFragMatTempOrd))])
        if(dim(SpeciesDNAFragMat)[1]==1){
          row.names(SpeciesDNAFragMat)=sp

          # All the sequences per DNA fragment, from the gene with maximum of sequences to gene with minimum of sequence per species.
          aa = as.data.frame(as.matrix(sort(colSums(SpeciesDNAFragMat), decreasing = T)))
          colnames(aa)=sp

          # Order the matrix from the genes with the maximum number of species coverage to
          # the minimum species coverage.
          bb = as.data.frame(as.matrix(sort(colSums(SpeciesDNAFragMatPA), decreasing = T)))
          colnames(bb)=sp
        } else {
          # All the sequences per DNA fragment.
          aa = as.data.frame(as.matrix(sort(colSums(SpeciesDNAFragMat), decreasing = T)))

          # Order the matrix from the genes with the maximum number of species coverage to
          # the minimum species coverage.
          bb = as.data.frame(as.matrix(sort(colSums(SpeciesDNAFragMatPA), decreasing = T)))
        }

        # DNA fragment perspective.
        DNASeq_Summary = merge(cbind(Name_DNA = row.names(bb), NB_Species = bb[,
            1]), cbind(Name_DNA = row.names(aa), NB_Seq = aa[, 1]), by.x = 1, by.y = 1,
            all.x = TRUE)
        DNASeq_SummarySort = DNASeq_Summary[order(as.numeric(as.character(DNASeq_Summary[,
            2])), decreasing = T), ]

        # Species perspective.
        aa2 = as.data.frame(as.matrix(sort(rowSums(SpeciesDNAFragMat), decreasing = T)))
        bb2 = as.data.frame(as.matrix(sort(rowSums(SpeciesDNAFragMatPA), decreasing = T)))

        SpeciesSeq_Summary = merge(cbind(Name_Species = row.names(bb2), NB_TypeDNA = bb2[,
            1]), cbind(Name_Species = row.names(aa2), NB_Seq = aa2[, 1]), by.x = 1,
            by.y = 1, all.x = TRUE)
        SpeciesSeq_SummarySort = SpeciesSeq_Summary[order(as.numeric(as.character(SpeciesSeq_Summary[,
            2])), decreasing = T), ]

        if(length(which(colnames(SpeciesDNAFragMat)=="NA"))==0){
          # Export the species-by-gene (DNA fragments) matrix.
          SpeciesDNAFragMat = cbind(Species_Name = row.names(SpeciesDNAFragMat), SpeciesDNAFragMat)
        } else {
          # Export the species-by-gene (DNA fragments) matrix, and removing a potential marker called NA.
          SpeciesDNAFragMat = cbind(Species_Name = row.names(SpeciesDNAFragMat),
                                    SpeciesDNAFragMat[,-which(colnames(SpeciesDNAFragMat)=="NA")])
        }

        # Export the species-by-gene (DNA fragments) matrix.
        # SpeciesDNAFragMat = cbind(Species_Name = row.names(SpeciesDNAFragMat), SpeciesDNAFragMat[,
        #    -dim(SpeciesDNAFragMat)[2]])  ## Remove the NA.

    } else {
        # If there is only one gene.
        SpeciesDNAFragMat = as.matrix(SpeciesDNAFragMat[order(SpeciesDNAFragMat[,
            dim(SpeciesDNAFragMat)[2]]), ])
        colnames(SpeciesDNAFragMat) = DNAFragUniq[1]

        # All the sequences per DNA fragment.
        aa = as.data.frame(as.matrix(sort(colSums(SpeciesDNAFragMat), decreasing = T)))
        bb = as.data.frame(as.matrix(sort(colSums(SpeciesDNAFragMatPA), decreasing = T)))

        # DNA (gene) perspective.
        DNASeq_SummarySort = cbind(Name_DNA = row.names(bb), NB_Species = bb[, 1],
            NB_Seq = colSums(SpeciesDNAFragMat))

        # Species perspective.
        aa2 = as.data.frame(as.matrix(sort(rowSums(SpeciesDNAFragMat), decreasing = T)))
        bb2 = as.data.frame(as.matrix(sort(rowSums(SpeciesDNAFragMatPA), decreasing = T)))

        SpeciesSeq_Summary = merge(cbind(Name_Species = row.names(bb2), NB_TypeDNA = bb2[,
            1]), cbind(Name_Species = row.names(aa2), NB_Seq = aa2[, 1]), by.x = 1,
            by.y = 1, all.x = TRUE)
        SpeciesSeq_SummarySort = SpeciesSeq_Summary[order(as.numeric(as.character(SpeciesSeq_Summary[,
            2])), decreasing = T), ]

        # Export the species-by-gene (DNA fragments) matrix.
        SpeciesDNAFragMat = cbind(Species_Name = row.names(SpeciesDNAFragMat), SpeciesDNAFragMat)
    }
    SpeciesDNAFragMat=as.data.frame(SpeciesDNAFragMat)
    # The species-by-gene matrix
    utils::write.table(SpeciesDNAFragMat, file = paste(output, "SpDNA_Mat.txt", sep = ""),
        sep = "\t", row.names = FALSE)
    # Summary DNA perspective.
    utils::write.table(DNASeq_SummarySort, file = paste(output, "Summary_DNApers.txt", sep = ""),
        sep = "\t", row.names = FALSE)
    # Summary species perspective
    utils::write.table(SpeciesSeq_SummarySort, file = paste(output, "Summary_SPECIESpers.txt",
        sep = ""), sep = "\t", row.names = FALSE)
    # Export a clean table.
    utils::write.table(inputCl1, file = paste(output, "_CleanDataset.txt", sep = ""), sep = "\t",
        row.names = FALSE)

    return(list(Specie.Gene_matrix = SpeciesDNAFragMat, Summary_DNA = DNASeq_SummarySort,
        Summary_Species = SpeciesSeq_SummarySort, MissingSpecies_WithoutSequences = MissingSpecies))
}  # End of the function.


# Example Sp.DNAMat=SpeciesGeneMat.Bl(input=AllSeqDF3, output='SpAll.DNA.Mat_',
# NCBI.Trash=NCBI.SeqTrash, BOLD.Trash=BOLD.SeqTrash)


# Help for Debug.  rm(SpeciesDNAFragMat, DNASeq_SummarySort,
# SpeciesSeq_SummarySort, aa, bb, aa2, bb2, a, SpeciesDNAFragMatPA,
# SpeciesDNAFragMat, sp, ProductClean, inputCl1, inputCl2, inputCl3, input,
# input2, input3, input4, ProductCl, GeneClean, DNAfrag, product, All.pos, W.pos,
# GeneNameHomog, a, tSpeciesDNAFragMatTempOrd, tSpeciesDNAFragMatTemp,
# tSpeciesDNAFragMatTemp, SpeciesDNAFragMatTemp, DNAFragUniq, inputb,
# NCBI.SeqTrash, RAG1.Trash, zic1.Trash, myh6.Trash, plagl2.Trash, rhod.Trash)

