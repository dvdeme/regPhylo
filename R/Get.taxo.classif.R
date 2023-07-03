#' @title Extract a taxonomic information and check for synonyms in the NCBI, itis, gbif, bold and
#' taxref (MNHN) databases

#' @description This functions extracts the species list, the classification and the synonyms 
#' species name of the accepted species names from the ncbi taxonomy, itis, gbif, bold, taxref (MNHN
#' INPN) databases from a taxid of a taxon from a hierarchical order from species level to higher
#' taxonomic levels (eg. "Genus", "Family" , "Order"..), obtain by the \emph{get_ids("plecoptera", db =
#' "ncbi")} function for instance.

#' @details This function allows to check for the presence of the synonyms found in a database, their
#' are also cross checked in the other databases. In this case, the function makes the distinction between 
#' the original search ("Ori") in each database, form the secondary search ("Sec") when the cross checking of synonyms
#' is performed, this information is available through the column \emph{search.type}.
#' This function allow to loop the search over multiple species or genera or family at the 
#' same time over the different databases.

#' @param Taxa.list the taxid of the taxa of interest, (extract all the descendant species, from the
#' ncbi, itis, gbif, bold and taxref database.

#' @param dbs a vector with the database to query either "ncbi", "itis", "gbif", "bold", or "taxref" (TAXREFV1.5,
#' from the rtaxref R package).

#' @param downto hierarchical level to look for descent, by default it is species. see function
#' downstream from the taxize R package.

#' @param local.db THIS OPTION IS CURRENTLY INVALID, WORK IN PROGRESS, if TRUE (by default FALSE),it use local imported database to query the
#' classification and performed request, can be much faster for important request. 
#' This is working with "ncbi","gbif" (not for the search of synonyms). 
#' If FALSE the function use the API to query the online server.

#' @param ID.Rank the name of the Rank of the Taxa.list, it can be either "Species", or any other
#' levels (by default, the function considers the "Genus"), however in the results only the Species,
#' Genus, Family, Order and Class hierarchical levels are reported.

#' @param input.ID.Rank the Rank of the ID requested

#' @param Save.Intermediate.ID.rank the path of the folder to store the results of the different
#' Taxa.list. 


#' @return The function returns a table with the following fields
#' c("species.valid", "species.syn", "genus", "family", "order", "class", "txid.species.valid",
#' "txid.species.syn", "txid.genus", "txid.family", "txid.order", "txid.class", "database", "search.type",
#' "extraction.date").


#' @example
#' # List the database
#' \dontrun{
#' dbs = c("ncbi", "bold", "itis", "gbif", "taxref")
#' # List the API keys (needed for NCBI).
#' api_keys = c("4XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX", "NULL", "NULL", "NULL", "NULL")

#' # List of the different genera to search
#' Taxa.list = c("Metanoea", "Anomalopterygella", "Cryptothrix", "Monocentra", "Leptodrusus")

#' # Run the function
#' Included.In.Drusus = Get.taxo.classif.Multi(Taxa.list = Taxa.list, dbs = dbs, api_keys = api_keys, ncbi.api.key = "472778dfff3834888186b7b0335e635b0e08", ID.Rank = "Genus", Save.Intermediate.ID.rank = "Results/Genus")
#' 
#' }

#' @export Get.taxo.classif.Multi


Get.taxo.classif.Multi = function(Taxa.list = NULL, dbs = NULL, api_keys = NULL, ncbi.api.key = NULL,  
                                  ID.Rank = "Genus", downto = "species", Save.Intermediate.ID.rank = NULL){

Overall.Results = do.call(rbind, lapply(1:length(Taxa.list), function(i){
# get the id
input.ids = unlist(lapply(1:length(dbs), function(x){
get.db.ids(taxon.name = Taxa.list[i], db = dbs[x])
}))


### Get the species list, the classification, and the synonyms
res.DF = do.call(rbind, lapply(1:length(dbs), function(x){
#x= 1
#for (x in 1:length(dbs)){
if(is.na(input.ids[x]) == TRUE){
c(rep("NA", 12),dbs[x], NA)
} else {
Get.taxo.classif(input.id = input.ids[x], db = dbs[x], downto = downto, api_key = api_keys[x], input.ID.Rank = ID.Rank)
}
}))



### Remove the taxa included as NA (absent in the database).
ToremoveNA = which(res.DF[,1]== "NA")
if(length(ToremoveNA) > 0){
which(is.na(res.DF[,1]))
res.DF = res.DF[-ToremoveNA,]
}
#dbs = unique(res.DF$database)


if(dim(res.DF)[1] > 0){
if(length(dbs) > 1){
### Determine the synonyms names that are not found in the other DB.
  ### THE PROBLEM LAy HERE, TO GLOBAL as a method to compare the synonyms.
Syno.Sp = setdiff(res.DF[,"species.syn"], res.DF[,"species.valid"])

if(length(Syno.Sp) > 0){ 

Toremove = which(is.na(Syno.Sp))
if(length(Toremove) > 0){
Syno.Sp = Syno.Sp[-Toremove]
}

#### Check if among the different synonyms recovered in the different DB, in that case we look for in the 
# other database with the different synonym name.
Res.DF3 = do.call(rbind, lapply(1:length(Syno.Sp), function(j){
#j=1
#for(j in 98:length(Syno.Sp)){
a = res.DF[which(res.DF[,"species.syn"] == Syno.Sp[j]), c("species.valid", "species.syn", "database")]

## Check in which Database the synonyms species should be look for
BD.2.look = setdiff(unique(res.DF[,"database"]), a[,"database"])

if(length(BD.2.look) > 0){ ### in case the synonyms are already reported in all the db but with distinct valid species names

# get taxid 
input.ids2 = unlist(lapply(1:length(BD.2.look), function(x){
get.db.ids(taxon.name = a[1, "species.syn"], db = BD.2.look[x], input.ID.Rank = "Species")
}))

### Check the DB to keep for additional search
RemainDB = BD.2.look[-which(is.na(input.ids2))]
remainId = input.ids2[-which(is.na(input.ids2))]

if(length(RemainDB) > 0){
### check the presence of the ncbi DB in that case the proper api_key need to to be setup.
ncbi.pres = which(RemainDB == "ncbi")
api_keys2 = rep("NULL", length(RemainDB))
if(length(ncbi.pres) > 0){
api_keys2[ncbi.pres] = ncbi.api.key
}

### Get the species list, the classification, and the synonymes
res.DF2 = do.call(rbind, lapply(1:length(RemainDB), function(x){

if(RemainDB[x] == "taxref"){
Get.taxo.classif(input.id = as.numeric(as.character(remainId[x])), db = RemainDB[x], downto = "species", 
                 api_key = api_keys2[x], input.ID.Rank = "Species")
} else {
Get.taxo.classif(input.id = remainId[x], db = RemainDB[x], downto = "species", api_key = api_keys2[x], input.ID.Rank = "Species")
}
}))

res.DF2
}
}

}))



### The database name follow by ".Ori" indicates that this is the original 
# data extraction from the database, to distinguished with synonyms species 
# considered in other database that could be secondary extracted 
# in the following steps.
res.DF$database = paste(res.DF$database, ".Ori", sep="")


### Combine the results of the first search based on the genus of interest and the 
# results of the second search based on the results from the synonymised genera from the different database.
if(is.null(Res.DF3)){
  res.tot = res.DF
} else {
  ### Add the suffix "Sec" to the name of teh database to distinguish the secondary
  # search done due to the presence of additional synonymes species in other database.
  Res.DF3$database = paste(Res.DF3$database, ".Sec", sep="")
  
  ### For some reason sometimes the taxref return the taxid of the genus and not from the species because
  # that species is not present in the taxref (but the genus is present), 
  # and the then the function extract the information at the genus level while we are looking for some species.
  # So we remove those non wanted genera.
  if(ID.Rank == "Genus" ) {
    Toremove = which(is.na(Res.DF3$txid.genus))
    if(length(Toremove) >0){
      Res.DF3 = Res.DF3[-Toremove,]
    }
  }
  
  # Finally put the two data.frame together
  res.tot = data.frame(rbind(res.DF, Res.DF3))
}

} else {
res.tot = res.DF
}

} else {
res.tot = res.DF
}

## Remove Potential replicates 
Duplic.2Remove = which(duplicated(res.tot[,-14]) == TRUE)
if(length(Duplic.2Remove) > 0){
res.tot = res.tot[-Duplic.2Remove,]
}

### Separate the Original and the secondary search from the database search
res.tot$search.type  = unlist(lapply(strsplit(res.tot$database, ".", fixed = T), function(k){k[2]}))
res.tot$database  = unlist(lapply(strsplit(res.tot$database, ".", fixed = T), function(k){k[1]}))
# re-organise the columns
res.tot = res.tot[,c(1:13, 15, 14)]

if(is.null(Save.Intermediate.ID.rank) == FALSE){
saveRDS(res.tot, file = paste(Save.Intermediate.ID.rank, "/", Taxa.list[i],".", paste(dbs, collapse="."), ".RDS", sep = ""))
}

} else {
warning(paste("None of the dbs could find ", Taxa.list[i], sep=""))
}


res.tot
}))

return(Overall.Results)
}


# For Debug 
# rm(res.tot, res.DF, Res.DF3,res.DF2, a, b, c, remainId, remainDB, Syno.Sp, BD.2.look, input.ids2, input.ids)



##########################################################"
#' @title From a taxonomic name get the id of taxon from different open-acces databases, ncbid, gbif, bold, itis, and taxref

#' @description This function from a taxonomic name (eg. "Plecopetra") gets the id of a taxon from
#' different open acess data base, ncbid, gbif, bold, itis, and taxref (this function is used
#' internally by Get.taxo.classif.Multi

#' @param taxon.name = the name of the taxon of interest eg. "Drusus"

#' @param db =  the names of the taxonomic data base, among "itis", "ncbi", "tropicos", "gbif", "eol", "nbn", "pow", "bold", and taxref (via rtaxref Rpackage to get acess to the taxref from the French INPN database)).

#' @param ncbi.api.key a string with the ncbi api key (personal to the user, it should be created via NCBI website)

#' @param input.ID.Rank the Rank of the ID requested eg. "Genus" or "Species"

#' @return the function returns a vector of the taxonomic id of the corresponding databases.

#' @export get.db.ids


get.db.ids = function(taxon.name = NULL, db = NULL, ncbi.api.key = NULL, input.ID.Rank = "Genus"){
if(db == "bold"){
res.id = tryCatch(as.character(taxize::bold_search(taxon.name)[1,1]), error=function(e) "error")
if(res.id == "error" | res.id == taxon.name){
res.id = NA
}

} else {
if(db == "taxref"){### get the taxid of the TAXREF data base
  if(input.ID.Rank == "Species"){
    res.id = tryCatch(as.data.frame(rtaxref::rt_taxa_search(sciname = taxon.name)), error=function(e) "error")
    if(dim(res.id)[2] == 1){
      class(res.id) = "character"
    }
  } else {
res.id = tryCatch(as.data.frame(rtaxref::rt_taxa_fuzzymatch(taxon.name)), error=function(e) "error")
}
if(class(res.id) == "character"){
res.id = NA
} else {
res.id = res.id[1,1]
}

} else {
### For ncbi, gbif, itis

if(db == "itis"){
res.id = tryCatch(as.character(taxize::get_tsn(taxon.name, accepted = TRUE, rows = 1)[[1]][1]), error=function(e) "error")

} else {
res.id = tryCatch(as.character(taxize::get_ids(taxon.name, db = db, rows = 1)[[1]][1]), error=function(e) "error")

#if(res.id == "error"){
#res.id = NA
#}

}
}
}


if(is.na(res.id) == TRUE){
warning(paste(taxon.name, " cannot be found in ", db, sep = ""))
}

res.id
}




#################################################################################################

#' @title Extract the species list, classification and check for the presence of synonyms in the 
#' NCBI taxonomy, itis, gbif, bold and taxref databases for a taxa.

#' @description This functions extracts the species list, the classification and the check for the
#' presence of synonyms species names in the ncbi taxonomy, itis, gbif, bold, taxref (MNHN NPN)
#' databases from a taxid of a taxon from a given hierarchy within the Linnean classification, from
#' genus or higher (eg. family , order..), obtain by the get_ids("plecoptera", db = "ncbi") function
#' for instance.

#' @details This function works only if the api_key of the user is provided in the option
#' \emph{api_key}. This function used the \emph{rtaxref} R package to acces the TAXREF database.

#' @param input.id a vector with the taxid of the taxa of interest (extract all the descendant
#' species, from the ncbi, bold, itis, gbif and taxref databases).

#' @param db a vector with the names of the database to query such as "ncbi", "itis", "gbif", "bold",
#' or "taxref" (TAXREFV1.5, from the rtaxref R package).

#' @param downto a hierachical level to look for descent, by default it is "species". see function
#' downstream from the taxize R package.

#' @param local.db [currently not working but work in progress] if TRUE (by default FALSE), it uses
#' local imported database to query the classification and perfomed request, can be much faster for
#' important request. This is working with "ncbi", "gbif" (not for the search of synonyms). if FALSE
#' the function use the API to query the online server.

#' @param api_key a string with the api key of the database. 

#' @param input.ID.Rank the Rank of the ID requested eg. "Genus" or "Species"

#' @return The function returns a table with the following fields
#' c("species.valid", "species.syn", "genus", "family", "order", "class", "txid.species.valid",
#' "txid.species.syn", "txid.genus", "txid.family", "txid.order", "txid.class", "database",
#' "extraction.date").

#' @export Get.taxo.classif


Get.taxo.classif = function(input.id = NULL, db = NULL, downto = "species", api_key = NULL, 
                            local.db = FALSE, input.ID.Rank = "Genus"){

### Test if the input.id correspond to a species level or higher taxonomic rank.
if(input.ID.Rank == "Species"){
if(db == "taxref"){

Child.ID.DF = as.data.frame(do.call(rbind, lapply(1:length(input.id), function(x){
rtaxref::rt_taxa_id(input.id[x])
})))

} else {

### if other db accessible from taxize are accessible at the species level
if(local.db == TRUE){
sp.Info = taxizedb::taxid2name(input.id, db = db)

} else {

if(db == "itis"){
# for "itis" when the taxa is invalid and recognised as a synonyms we also provide the tsn of the valid taxa.
sp.Info = lapply(1:length(input.id), function(k){
z = taxize::itis_getrecord(input.id[k])
# check if the taxa is  avalid taxa or a synonyms in itis DB.
if(class(z$acceptedNameList$acceptedNames) == "data.frame"){
data.frame(id = z$taxRank$tsn, name = z$scientificName$combinedName, rank = gsub("[ ]*", "", z$taxRank$rankName, perl = T), acceptedName = z$acceptedNameList$acceptedNames$acceptedName, acceptedTsn =  z$acceptedNameList$acceptedNames$acceptedTsn)
} else {
data.frame(id = z$taxRank$tsn, name = z$scientificName$combinedName, rank = gsub("[ ]*", "", z$taxRank$rankName, perl = T), acceptedName = NA, acceptedTsn = NA)
}
})


} else {
sp.Info = taxize::id2name(input.id, db = db)
}

}

if(dim(sp.Info[[1]])[2] == 5){
child.ID = data.frame(do.call(rbind, lapply(sp.Info, function(x){
c(x$id, x$name, x$rank, x$acceptedName, x$acceptedTsn)
})))
names(child.ID) = c("tsn", "taxonname", "rankname", "acceptedName", "acceptedTsn")
} else {
child.ID = data.frame(do.call(rbind, lapply(sp.Info, function(x){
c(x$id, x$name, x$rank)
})))
}



if(db == "gbif") names(child.ID) = c("key", "name", "rank")
if(db == "bold") names(child.ID) = c("id", "name", "rank")
if(db == "ncbi") names(child.ID) = c("childtaxa_id", "childtaxa_name", "rank")

child.ID = list(child.ID)
}

#### If the taxonomic rank of interest is higher than the species level.
} else {

### Get all the descendents down to the species level
if(local.db == TRUE){
child.ID = taxizedb::downstream(input.id, db = db, downto = downto)
} else {
if(db == "taxref"){
Child.ID.DF = tryCatch(as.data.frame(rtaxref::rt_taxa_children(as.numeric(as.character(input.id)))), error=function(e) "error")

### If the taxonomic names correspond to a synonyms of another validated taxa name.
if(!class(Child.ID.DF) == "data.frame"){
taxa.info = as.data.frame(rtaxref::rt_taxa_id(as.numeric(as.character(input.id))))
Child.ID.DF = as.data.frame(rtaxref::rt_taxa_children(as.numeric(as.character(taxa.info$referenceId))))
}

} else {
child.ID = taxize::downstream(input.id, db = db, downto = downto)
}
}
}



#### GBIF data base.
if(db == "gbif"){
if(dim(child.ID[[1]])[1] == 0){

if(local.db == TRUE){
sp.Info = taxizedb::taxid2name(input.id, db = db)
} else {
sp.Info = taxize::id2name(input.id, db = db)
}

child.ID.DF = data.frame(sp.Info[[1]]$id, sp.Info[[1]]$name, sp.Info[[1]]$rank)

} else {

child.ID.DF = data.frame(child.ID[[1]]$key, child.ID[[1]]$name, child.ID[[1]]$rank)

}

### Get the Classification
if(local.db == TRUE){
Classif.tax = taxizedb::classification(child.ID.DF[,1], db = db)
} else {
Classif.tax = taxize::classification(child.ID.DF[,1], db = db)
}

Classif.Summary = Get.classif(input = Classif.tax, Tax.list = child.ID.DF[,c(2,1)])
Classif.Summary = data.frame(Classif.Summary, database = rep(db, dim(Classif.Summary)[1]), extraction.date = rep(Sys.time(), dim(Classif.Summary)[1]))
Classif.Summary = data.frame(species.valid = Classif.Summary[,1],  species.syn = Classif.Summary[,1], Classif.Summary[,c(2:5)], txid.species.valid = Classif.Summary[,6], txid.species.syn = Classif.Summary[,6], Classif.Summary[, c(7:12)])
Classif.Summary$extraction.date = as.character(Classif.Summary$extraction.date)

### Get the synonyms
Complement.Syno = do.call(rbind, lapply(1:length(Classif.Summary[,1]), function(x){
#x = 1
#for(x in 1:length(Classif.Summary[,1])){
ab = as.data.frame(rgbif::name_lookup(Classif.Summary[x,1], rank = "species", status = "SYNONYM", limit = 1000, verbose = FALSE)$data)

if(dim(ab)[1] == 0){
Classif.Summary[x,]

} else {

abb = tryCatch(ab[,c("key", "canonicalName")], error=function(e) "error")
if(class(abb) == "character"){
abb = tryCatch(ab[,c("key", "species")], error=function(e) "error")
}


DF = do.call(rbind, lapply(1:dim(abb)[1], function(i){
Classif.Summary[x,]
}))

DF$species.syn = abb[,2]
DF$txid.species.syn = abb[,1]
rbind(Classif.Summary[x,], DF)


}
}))
Classif.Summary = unique(Complement.Syno)
}


#### BOLD database
if(db == "bold"){
if(dim(child.ID[[1]])[1] == 0){ 
sp.Info = id2name(input.id, db = db)
child.ID.DF = data.frame(sp.Info[[1]]$id, sp.Info[[1]]$name, sp.Info[[1]]$rank)
} else {
child.ID.DF = data.frame(child.ID[[1]]$id, child.ID[[1]]$name, child.ID[[1]]$rank)
}

Classif.tax = taxize::classification(child.ID.DF[,1], db = db)
Classif.Summary = Get.classif(input = Classif.tax, Tax.list = child.ID.DF[,c(2,1)])
Classif.Summary = data.frame(Classif.Summary, database = rep(db, dim(Classif.Summary)[1]), extraction.date = rep(Sys.time(), dim(Classif.Summary)[1]))
Classif.Summary = data.frame(species.valid = Classif.Summary[,1],  species.syn = Classif.Summary[,1], Classif.Summary[,c(2:5)], txid.species.valid = Classif.Summary[,6], txid.species.syn = Classif.Summary[,6], Classif.Summary[, c(7:12)])
Classif.Summary$extraction.date = as.character(Classif.Summary$extraction.date)
### BOLD does not keep track of the synonyms so we cannot check for those !
warning("BOLD does not provide any information about synonymes, so \"species.syn\" is identical to \"species.valid\"")
}


### ITIS database
if(db == "itis"){

if(is.null(child.ID[[1]]$tsn[1])){
child.ID[[1]] = data.frame(tsn = "No data")
}

if(child.ID[[1]]$tsn[1] == "No data"){ ### if the input.id refer to a species already
sp.Info = taxize::id2name(input.id, db = db)
child.ID.DF = data.frame(sp.Info[[1]]$id, sp.Info[[1]]$name, sp.Info[[1]]$rank)
names(child.ID.DF) = c("tsn", "taxonname", "rankname")

} else {

if(input.ID.Rank == "Species"){
child.ID.DF = data.frame(child.ID[[1]]$tsn, child.ID[[1]]$taxonname, child.ID[[1]]$rankname, child.ID[[1]]$acceptedName, child.ID[[1]]$acceptedTsn)
names(child.ID.DF) = c("tsn", "taxonname", "rankname", "acceptedName", "acceptedTsn")
} else {
child.ID.DF = data.frame(child.ID[[1]]$tsn, child.ID[[1]]$taxonname, child.ID[[1]]$rankname)
names(child.ID.DF) = c("tsn", "taxonname", "rankname")
}
}


### Get the Classification
if(input.ID.Rank == "Species"){
if(is.na(child.ID.DF$acceptedTsn)){
Classif.tax = taxize::classification(child.ID.DF[,1], db = db)
} else {
Classif.tax = taxize::classification(child.ID.DF$acceptedTsn, db = db)
}
} else {
Classif.tax = taxize::classification(child.ID.DF[,1], db = db)
}

Classif.Summary = Get.classif(input = Classif.tax, Tax.list = child.ID.DF[,c(2,1)])
Classif.Summary = data.frame(Classif.Summary, database = rep(db, dim(Classif.Summary)[1]), extraction.date = as.character(rep(Sys.time(), dim(Classif.Summary)[1])))
Classif.Summary = data.frame(species.valid = Classif.Summary[,1],  species.syn = child.ID.DF[,"taxonname"], Classif.Summary[,c(2:5)], txid.species.valid = Classif.Summary[,6], txid.species.syn = child.ID.DF[,"tsn"], Classif.Summary[, c(7:12)])

### Get the synonyms from the "itis" database.
Synonymes.itis = taxize::synonyms(Classif.Summary$species.valid, db="itis")
x=1
Classif.Summary = do.call(rbind, lapply(1:length(Synonymes.itis), function(x){
if(length(Synonymes.itis[[x]]) > 0){
res = do.call(rbind, lapply(1:dim(Synonymes.itis[[x]])[1], function(i){
a = unlist(c(Classif.Summary[x,1], Synonymes.itis[[x]][i,4], Classif.Summary[x,c(3:7)], Synonymes.itis[[x]][i,5], Classif.Summary[x, c(9:13)], as.character(Classif.Summary[x,14])))
names(a) = names(Classif.Summary)
rbind(Classif.Summary[x,], a)
}))
} else {
res = Classif.Summary[x,]
}
res
}))
Classif.Summary = unique(Classif.Summary) ### remove potential duplicates 

}



##### NCBI
if(db == "ncbi"){
if(is.null(child.ID[[1]]$childtaxa_id)){ ### if the input.id refer to a species already

if(local.db == TRUE){
sp.Info = taxizedb::taxid2name(input.id, db = db)
} else {
sp.Info = taxize::id2name(input.id, db = db)
}

child.ID.DF = data.frame(sp.Info[[1]]$id, sp.Info[[1]]$name, sp.Info[[1]]$rank)
} else {
child.ID.DF = data.frame(child.ID[[1]]$childtaxa_id, child.ID[[1]]$childtaxa_name, child.ID[[1]]$rank)
}

if(local.db == TRUE){
#### NEED TO WORK ON THE NEW VERSION OF THE get_classif.ncbi FUNCTION TO ALLOW REQUEST FROM THE LOCAL DOWNLOAD OF THE NCBI TAXONOMY DATABASE.
} else {
Classif.Summary = get_classif.ncbi(input = child.ID.DF[,1], api_key = api_key)
}
}



#### TAXREF
if(db == "taxref"){

### get the info about the classification 
Classif.summary = Child.ID.DF[,c("scientificName", "scientificName", "genusName", "familyName", "orderName", "className", "referenceId", "id")]


Valid.species = as.data.frame(do.call(rbind, lapply(1:length(Classif.summary$referenceId), function(x){
rtaxref::rt_taxa_id(Classif.summary$referenceId[x])
})))
Valid.species = Valid.species[,c("scientificName", "referenceId")]

hier.levels = c("Genre", "Famille", "Ordre", "Classe")

Valid.species.Taxonomy = as.data.frame(do.call(rbind, lapply(1:length(Valid.species$referenceId), function(x){
b = tryCatch(as.data.frame(rt_taxa_classification(as.numeric(as.character(Valid.species$referenceId[x])))), error=function(e) "error")
if(class(b) == "character"){
c= rep(NA, 4)
} else {
c = c(b[match(hier.levels, b[, "rankName"]), "referenceId"])
}
unlist(c(Valid.species[x,], c, "taxref", as.character(rep(Sys.time()))))
})))
names(Valid.species.Taxonomy) = c("species.valid", "referenceId", "txid.genus", "txid.family", "txid.order", "txid.class", "database", "extraction.date")


Classif.Summary = do.call(rbind, lapply(1:dim(Classif.summary)[1], function(x){
a = tryCatch(as.data.frame(rtaxref::rt_taxa_synonyms(Classif.summary[x,c("referenceId")])), error=function(e) "error")

if(class(a) == "character"){
res1 = data.frame(Classif.summary[x,], Valid.species.Taxonomy[x,c("txid.genus", "txid.family", "txid.order", "txid.class", "database", "extraction.date")])  ## there is no synonymes so the full names is provided
} else {
res1 = a[,c("scientificName", "genusName", "familyName", "orderName", "className", "id")] ## there is no synonymes so the full names is
res2 = data.frame(scientificName = rep(Valid.species[x, "scientificName"], dim(res1)[1]),  
                  res1[,-dim(res1)[2]], 
                  referenceId = rep(Valid.species[x, "referenceId"], dim(res1)[1]), 
                  id = res1[, dim(res1)[2]], 
                  txid.genus = rep(Valid.species.Taxonomy[x, "txid.genus"], dim(res1)[1]), 
                  txid.family = rep(Valid.species.Taxonomy[x, "txid.family"], dim(res1)[1]), 
                  txid.order = rep(Valid.species.Taxonomy[x, "txid.order"], dim(res1)[1]), 
                  txid.class = rep(Valid.species.Taxonomy[x, "txid.class"], dim(res1)[1]), 
                  database = rep(Valid.species.Taxonomy[x, "database"],dim(res1)[1]), 
                  extraction.date = rep(Valid.species.Taxonomy[x, "extraction.date"], dim(res1)[1]))
res0 = data.frame(Classif.summary[x,], Valid.species.Taxonomy[x,c("txid.genus", "txid.family", "txid.order", "txid.class", "database", "extraction.date")])  ## there is no synonymes so the full names is provided
res2 = data.frame(rbind(res0, res2))
res1 = res2
}
res1
}))

# Remove potential duplicate (without accounting the date and timing of extraction).
ToRemove = which(duplicated(Classif.Summary[,-14]))
if(length(ToRemove) > 0){
Classif.Summary = Classif.Summary[-ToRemove,]
}

}

names(Classif.Summary) = c("species.valid", "species.syn" , "genus", "family", "order", "class", "txid.species.valid", "txid.species.syn", "txid.genus", "txid.family", "txid.order", "txid.class", "database", "extraction.date")


return(Classif.Summary)
}


## For debug 
# rm(Classif.Summary, a,Synonymes.itis, res, child.ID.DF)



#################################################################################################

#' @title Get the classification table and synonyms from NCBI taxonomy database

#' @description This function gets the classification and the synonyms taxonomic name of the taxa of
#' interest from the NCBI taxonomy database using the rentrez R package. This function is used
#' internally by \emph{Get.taxo.classif} and \emph{Get.taxo.classif.Multi} functions.

#' @details This function works only if the api_key of the user is provided in the option
#' \emph{api_key}.

#' @param input a vector with the ncbi taxonomic ID of the taxa of interest.

#' @param api_key a vector with the personal api key to get a better access to the NCBI API.

#' @return The function returns a table with the following fields
#' c("species.valid", "species.syn", "genus", "family", "order", "class", "txid.species.valid",
#' "txid.species.syn", "txid.genus", "txid.family", "txid.order", "txid.class", "database",
#' "extraction.date").

#' @export get_classif.ncbi

get_classif.ncbi = function(input = NULL, api_key = NULL){
#lapply(1:length(input), function(x){
# Extract the info from ncbi, based on the taxid.
tax_rec <- rentrez::entrez_fetch(db="taxonomy", id=input, rettype="xml", parsed=TRUE, api_key = api_key)
# convert the xml classification file into a list
xml_data <- XML::xmlToList(tax_rec)

info = data.frame(do.call(rbind, lapply(1:length(xml_data), function(x){
#x= 1
#for(x in 1:length(xml_data)){

lineageEx = which(names(xml_data[[x]]) == "LineageEx")

if(length(which(names(xml_data[[x]]) == "OtherNames"))> 0){ 
# Aut.name = xml_data[[x]][[which(names(xml_data[[x]]) == "OtherNames")]]$Name$DispName
if(length(which(names(xml_data[[x]]$OtherNames) == "Synonym"))> 0){
nb.syn = which(names(xml_data[[x]]$OtherNames) == "Synonym")
Synonyms = unlist(lapply(1:length(nb.syn), function(i){
xml_data[[x]]$OtherNames[[nb.syn[i]]]
}))
} else {
Synonyms = xml_data[[x]]$ScientificName
}
} else {
Synonyms = xml_data[[x]]$ScientificName
}

### get the classification information
Classif =do.call(rbind, lapply(1:length(xml_data[[x]][[lineageEx]]), function(i){
unlist(xml_data[[x]][[lineageEx]][[i]])
}))

Hier.2.keep = c("class", "order", "family", "genus", "species")
aa = c(xml_data[[x]]$ScientificName, xml_data[[x]]$ScientificName, Classif[match(rev(Hier.2.keep)[-1], Classif[,3]),2], xml_data[[x]]$TaxId, xml_data[[x]]$TaxId, Classif[match(rev(Hier.2.keep)[-1], Classif[,3]),1], "ncbi", as.character(Sys.time()))

if(length(Synonyms) > 0){
aa = do.call(rbind, lapply(1:length(Synonyms), function(j){
c(aa[1], Synonyms[j], aa[3:length(aa)])
}))
}
aa
})))

names(info) = c("species.valid", "species.syn", "genus", "family", "order", "class", "txid.species.valid", "txid.species.syn", "txid.genus", "txid.family", "txid.order", "txid.class", "database", "extraction.date")
return(info)
}

### For debug !
# rm(info, aa, Classif, Synonyms, xml_dat, xml, nb.syn, lineageEx, tax_rec)



#################################################################################################
#' @title Arrange the classification information extracted from the itis, gbif and bold databases. 

#' @description This function arranges the classification from the ncbi or itis database using the taxize R package. It keeps the class, order, family, genus and species names with their associated txid. This function work internally for the \emph{Get.taxo.classif} function.

#' @param input is the output of the /emph{taxize::classification} function, which is 
#' a list of data.frames with the taxonomic classifiation of the supplied taxa.

#' @param Tax.list is a data frame with two columns the first include the accepted taxonomic name 
#' and the second the synonyms taxonomic names.

#' @export Get.classif

Get.classif = function(input = NULL, Tax.list = NULL){

Hier.2.keep = c("class", "order", "family", "genus", "species")
res = data.frame(do.call(rbind, lapply(1:length(input), function(x){
a = input[[x]]
if(dim(a)[1] > 0){
b= a[match(Hier.2.keep, a[,2]), ]
c(rev(b[,1]), rev(b[,3]))
} else {
c(Tax.list[x,1], rep(NA, 4), Tax.list[x,2], rep(NA, 4))
}
})))
names(res) = c(rev(Hier.2.keep), paste("txid.", rev(Hier.2.keep), sep = ""))
return(res)
}



##############################################################################################
###### Utils functions used by the rtaxref R package, present in the rtaxref.utils.R document.




#### Function used to extract the classification from the TAXREF via the API.
#### this function is derived from the rtaxref R package https://github.com/Rekyt/rtaxref

#' @title Extract the classification of taxa from the TAXREF via the API.

#' @description This function extracts the classification of a taxa from the TAXREF via the API.
#' t=This function is derived from the rtaxref R package https://github.com/Rekyt/rtaxref

#' @param id is a numeric value of the TAXREF id of the taxa.

#' @export rt_taxa_classification

rt_taxa_classification = function(id){

  check_required_arg(id, "retrieve the classification of a taxon based on its id")

  stopifnot("'id' must be a numeric" = is.numeric(id))

  api_query = rt_GET("taxa/", id, "/classification")

  parse_taxa(api_query)

}




# Other functions below are directly extracted from the https://github.com/Rekyt/rtaxref/blob/main/R/utils.R


rt_base_url = function() {
    "https://taxref.mnhn.fr/"
}

#' rtaxref User Agent
#'
#' @noRd
rt_ua <- function() {
  paste0("http://github.com/Rekyt/rtaxref R package rtaxref/",
         utils::packageVersion("rtaxref"))
}

#' @importFrom httr GET status_code add_headers
rt_GET <- function(..., query = NULL) {
  if(!is.null(query)) {
    query <- rt_flatten_query(query)
  }

  GET(rt_base_url(), config = add_headers("user-agent" = rt_ua()),
      path = paste0("api/", ...), query = query)
}

#' Flatten query
#'
#' Flatten a query, meaning that when an argument of the query is an array
#' the function gives back a list with as many elements in the array named the
#' same.
#'
#' @param query {`list(1)`}\cr{}
#'              a list that represents the query argument from `rt_GET()`
#'
#' @noRd
rt_flatten_query <- function(query) {
  flat_query <- lapply(names(query), function(el) {
    trans_list <- as.list(query[[el]])
    names(trans_list) <- rep(el, length(query[[el]]))
    return(trans_list)
  })

  unlist(flat_query, recursive = FALSE)
}



#' Simple Check
#'
#' Check that argument is not null nor equal to ""
#'
#' @param arg actual variable\cr{}
#'            give the variable to be used and checked
#' @param stop_message {`character(1)`}\cr{}
#'                     The `stop()` message to display explaining why this
#'                     argument is needed
#' @noRd
check_required_arg = function(arg, stop_message) {
  if (is.null(arg) | (length(arg) != 0 && arg == "")) {
    stop("'", substitute(arg), "' argument is needed to ", stop_message,
         call. = FALSE)
  }
}

check_arg_in_list = function(arg, list, with_null = TRUE) {
  if (
    ((isFALSE(arg %in% list) | is.null(arg)) & isFALSE(with_null)) |
    (isTRUE(with_null) & !is.null(arg) & isFALSE(arg %in% list))
  ) {
    stop(
      "'", substitute(arg), "' argument should be in '",
      paste(list, collapse = ", "), "'", ifelse(with_null, " or NULL", ""),
      call. = FALSE
    )
  }
}

#' @importFrom httr content http_error http_status status_code
parse_taxa = function(api_query, cut_names = TRUE) {

  reason = http_status(api_query)$reason

  if (status_code(api_query) == 400 & reason == "Bad Request") {

    stop("The query is invalid. Please try another query.", call. = FALSE)

  } else if (status_code(api_query) == 404 & reason == "Not Found"){

    stop("The query returned no results. Please try another query",
         call. = FALSE)

  } else if (http_error(api_query)) {

    stop("TAXREF is down. Please try again later.", call. = FALSE)

  }

  raw_response = content(api_query, type = "application/json",
                         encoding = "UTF-8", simplifyDataFrame = TRUE,
                         flatten = TRUE)

  if (!("_embedded" %in% names(raw_response))) {
    # If there is a single answer

    not_links = setdiff(names(raw_response), "_links")

    # If there is a single response
    response = as.data.frame(lapply(raw_response[not_links], function(x) {
      ifelse(is.null(x), NA, x)
    }))

  } else {
    # If there are several answers
    response = raw_response[["_embedded"]][[1]]

    # Tidy names
    all_names = colnames(response)

    name_categories = grepl(".", all_names, fixed = TRUE)

    # if there is composed names
    if (any(name_categories)) {
      name_cat = unique(lapply(all_names[name_categories], function(x) {
        strsplit(x, ".", fixed = TRUE)[[1]][1]
      }))

      # If taxon put it in front
      if ("taxon" %in% name_cat) {
        other_names = setdiff(seq_along(all_names), grep("taxon", all_names))

        response = response[, c(grep("taxon", all_names), other_names)]
      }

      if (cut_names) {
        # Delete category from original column names
        colnames(response) = gsub(paste(name_cat, "\\.", sep = "",
                                        collapse = "|"), "", colnames(response))
      }
    }
  }

  if (identical(dim(response), c(0L, 0L))) {
    stop("The query returned no results. Please try another query",
         call. = FALSE)
  }

  tibble::as_tibble(response, .name_repair = "universal")
}


