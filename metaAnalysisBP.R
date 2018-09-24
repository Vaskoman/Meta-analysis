# Meta-analysis for top GO-BP categories
# Datasets and genes that belong to these categories
# Novel biologically significant genes

## FUNCTIONS ####
extractDOgenes <- function (DO, table=dat20) {
  # returns a list of gene IDs corresponding to a
  # selected DO category from a selected enrichment table
  print (paste0 ("Extracting gene IDs for ", DO, " category 
                 from ", table, " enrichment table."))
  genes <- table[table$ID == DO, "geneID"]
  genesList <- strsplit(genes, "/", fixed = TRUE)
  genes <- unlist(genesList)
  genes <- unique(genes)
  print("Extraction complete!")
  return(genes)
}


## setwd, load libraries, & load data ####
setwd("/home/vassil/Documents/Bcells/Meta-analysis")
# source ("http://www.bioconductor.org/biocLite.R")
# biocLite("clusterProfiler") # also choose to update all packages
library(clusterProfiler)


## BP categories with the most datasets ####
load ("ckAllBP")
dat <- as.data.frame(ckAll) # transform to dataframe
categories <- table (dat$ID) # get the top 20 enriched DO categories
categories.order <- categories[order(categories, decreasing = TRUE)] # sort by number of appearances

dat.order <- data.frame(matrix(ncol=length(names(dat)), nrow=0))
colnames(dat.order) <- names(dat)

# extract datasets belonging to each category:
for (i in names(categories.order)){
  a <- dat[dat$ID==i,]
  dat.order <- rbind(dat.order,a)
}

dat20 <- dat.order[dat.order$ID %in% names(categories.order)[1:20],]
save(dat20, file = "dat20BP")
write.table(dat20, file="BP.top20.txt", sep="\t", row.names=FALSE)

rm (list = c("ckAll", "a", "dat", "dat.order","categories",
             "categories.order","i"))

