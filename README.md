# regPhylo

The regPhylo R package provides functions to assist in building phylogenetic trees appropriate for use in ecological studies, particularly for a regional species pool.

[Paper](https://doi.org/10.1111/geb.12986)


# Installation

regPhylo requires that the following package must be installed and loaded in the R environment:
*bold, seqinr, ape, geomedb, worrms, RJSONIO, stringr, fields, parallel, caper, phytools*
```{r, eval=FALSE}
install.packages(c("bold", "seqinr", "ape", "geomedb", "worrms", 
                   "RJSONIO", "stringr", "fields", "parallel", 
                   "caper", "phytools"))
library(bold)
library(seqinr)
library(ape)
library(geomedb)
library(worrms)
library(RJSONIO)
library(stringr)
library(fields)
library(parallel)
library(caper)
library(phytools)
```

To install the *regPhylo* R package do the following:
```{r, eval=FALSE}
install.packages("devtools")
library(devtools)

# Install the package from GitHub
install_github("dvdeme/regPhylo")
```
Load the *regPhylo* package.
```{r, eval=TRUE}
library(regPhylo)
```
&nbsp;

To see the list and short description of all the functions and data available in Rdata format.
```{r, eval=TRUE}
help(package=regPhylo)
```
&nbsp;

# Tutorial using a case study.
[A tutorial is available](https://doi.org/10.5061/dryad.2v4h7bv) to help the user to navigate through the different functions. The tutorial is based on an example aiming to build a posterior distribution of time-calibrated trees for the New Zealand marine ray-finned fishes.
