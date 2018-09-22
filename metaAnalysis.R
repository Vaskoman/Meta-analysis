# Meta-analysis for top DO and GO-BP categories
# Datasets and genes that belong to these categories
# Novel biologically significant genes

## FUNCTIONS ####

## setwd, load libraries, & load data ####
setwd("/home/vassil/Documents/Bcells/Meta-analysis")
source ("http://www.bioconductor.org/biocLite.R")
biocLite("clusterProfiler") # also choose to update all packages
library(clusterProfiler)


## DO categories with the most datasets ####
load ("ckAllDO")
dat <- as.data.frame(ckAll) # transform to dataframe
top20 <- table (dat$ID) # get the top 20 enriched DO categories
top20order <- top20[order(top20, decreasing = TRUE)] # sort by number of appearances
top20catDO <- names(top20order)[1:20] # select only the top 20 DO categories
dat20 <- dat[dat$ID %in% top20catDO,] # subset enriched DO for the top 20 only

# extract datasets belonging to each category:
cat20 <- rep(top20catDO, top20order[1:20])
a <- dat20[match(cat20, dat20$ID),]
dat20[order(dat20$ID),]

df <- data.frame(name=letters[1:4], value=c(rep(TRUE, 2), rep(FALSE, 2)))
target <- c("b", "c", "a", "d")
df[match(target, df$name),]
# extract gene names belonging to each category
# extract genes from master tables per category
# go over to select novel relevant significant genes (just up but not down)


## GO BP categories with the most datasets ####
# extract datasets belonging to each category
# extract gene names belonging to each category
# extract genes from master tables per category
# go over to select novel relevant significant genes (just up but not down)
