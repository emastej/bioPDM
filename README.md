bioPDM
======

### Introduction

`bioPDM`  is a downstream omics analysis method that uses mediation to explore the causal pathway
between an exposure and an outcome through high dimensional, multivariate
omic data. The goal of bioPDM is to provide a method to perform mediation analysis
on omic data where $p >> n$. The overall objective is to find groups of biological features,
or principal directions of mediation (PDM) that mediate the effect of an exposure on a phenotype
of interest. 

bioPDM assumes the user has an exposure $X$, a phenotype $Y$, and multivariate omic data 
$M$ measured on the same subjects. bioPDM cannot accomodate missing data at this time. Please use your
favorite imputation method before using bioPDM. 

Please perform these functions prior to analyzing your data with bioPDM: 

1. Imputation of missing values
2. Transformation, normalization, and batch correction


Install via Github:

    if (!require("devtools")) install.packages("devtools")
    devtools::install_github("emastej/bioPDM")

