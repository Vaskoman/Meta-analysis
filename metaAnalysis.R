# Meta-analysis for top DO and GO-BP categories
# Datasets and genes that belong to these categories
# Novel biologically significant genes

## FUNCTIONS ####

## setwd, load libraries, & load data ####
setwd("/home/vassil/Documents/Bcells/BioInfo/Meta-analysis")
source ("http://www.bioconductor.org/biocLite.R")
biocLite("clusterProfiler") # also choose to update all packages
library(clusterProfiler)



## DO categories with the most datasets ####
# extract datasets belonging to each category
# extract gene names belonging to each category
# extract genes from master tables per category
# go over to select novel relevant significant genes (just up but not down)


## GO BP categories with the most datasets ####
# extract datasets belonging to each category
# extract gene names belonging to each category
# extract genes from master tables per category
# go over to select novel relevant significant genes (just up but not down)
