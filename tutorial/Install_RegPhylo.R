remove.packages("regPhylo") # for updating regphylo with new github commits need to remove then reinstall.


install.packages("seqinr")
install.packages("ape")
install.packages("RJSONIO")
install.packages("stringr")
install.packages("fields")
install.packages('devtools')
install.packages('geomedb')
install.packages('phytools')
install.packages('caper')
install.packages('bold')
install.packages('httr')
install.packages('rentrez')
install.packages('curl')
install.packages('xml2')
install.packages('crul')

#Might help if having problems with some of the above
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")





# The "geomedb" requires the latest version available on Github, to download it, 

install.packages("devtools")
library(devtools)
install_github("biocodellc/fimsR-access")


# Install the package from GitHub
#Published version
####install_github("dvdeme/regPhylo")

# Install fork of original with some updated tools
install_github("ignacio3437/regPhylo")

#Install Blast+, Gblocks, trimAl and PartitionFinder2, mumsa (https://msa.sbc.su.se/cgi-bin/msa.cgi) if needed on command line


library(regPhylo)
library(bold)
library(httr)
library(rentrez)
library(seqinr)
library(ape)
library(RJSONIO)
library(stringr)
library(fields)
library(parallel)
library(caper)
library(phytools)
library(geomedb)


#help(package=regPhylo)

