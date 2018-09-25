# Meta-analysis for top DO
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


## DO categories with the most datasets ####
load ("ckAllDO")
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
save(dat20, file = "dat20DO")
write.table(dat20, file="DO.top20.txt", sep="\t", row.names=FALSE)

rm (list = c("ckAll", "a", "dat", "dat.order","categories",
             "categories.order","i"))

## Extract DO genes for each category (unique) ####

# Define a vector with all categories of interest
DOs <- c("DOID:1936", "DOID:552", "DOID:2320", "DOID:824",
         "DOID:1575", "DOID:854", "DOID:16", "DOID:0050338",
         "DOID:4001")
names(DOs) <- c("ATHEROSCLEROSIS", "PNEUMONIA", "OBSTRUCTIVE.LUNG.DISEASE",
                "PERIODONTITIS", "RHEUMATIC.DISEASE", "COLLAGEN.DISEASE",
                "INTEGUMENTARY.SYSTEM.DISEASE",
                "PRIMARY.BACTERIAL.INFECTIOUS.DISEASE", "OVARIAN.CARCINOMA")
DOsNames <- c("ATHEROSCLEROSIS", "PNEUMONIA", "OBSTRUCTIVE.LUNG.DISEASE",
                "PERIODONTITIS", "RHEUMATIC.DISEASE", "COLLAGEN.DISEASE",
                "INTEGUMENTARY.SYSTEM.DISEASE",
                "PRIMARY.BACTERIAL.INFECTIOUS.DISEASE", "OVARIAN.CARCINOMA")

## Extract genes from master tables for each DO ####
for (i in 1:length(DOs)){
  genes <- extractDOgenes (DOs[i])
  assign(paste0(DOsNames[i],"genes"), genes)
  rm (genes)
}

## Create DO directories and save genes IDs ####
for (i in DOsNames){
  if (!file.exists(i)) {
    dir.create(i)
  }
  filename <- paste0(i, "genes")
  savefile <- get(filename)
  dirname <- paste0 (getwd(), "/", i, "/", filename)
  save (list=savefile, file = dirname)
}

## Load master tables ####

load ("masterTableUp")
load ("masterTableDown")

## Get tables per DO category and save ####

for (i in DOsNames) {
  geneNames <- get(paste0(i, "genes"))
  tabUp <- masterTableUp[masterTableUp$Gene.ID %in% geneNames,]
  tabDown <- masterTableDown[masterTableDown$Gene.ID %in% geneNames,]
  assign(paste0(i,"tabUp"), tabUp)
  assign(paste0(i,"tabDown"), tabDown)
  filename <- paste0(i,"tabUp")
  save (list=filename, file = paste0(i, "/",i,"tabUp"))
  write.table(get(paste0(i,"tabUp")),
              file = paste0(i, "/",i,"tabUp.txt"),
              sep="\t", row.names = FALSE)
  filenameDown <- paste0(i,"tabDown")
  save (list=filenameDown, file = paste0(i, "/",i,"tabDown"))
  write.table(get(paste0(i,"tabDown")),
              file = paste0(i, "/",i,"tabDown.txt"),
              sep="\t", row.names = FALSE)
}

rm (list = c("tabUp", "tabDown", "filenameDown",
             "filename","geneNames"))
rm (list = c("dirname","i"))

#########################################################################

rm (list = ls())
load("/home/vassil/Documents/Bcells/Meta-analysis/ATHEROSCLEROSIS/ATHEROSCLEROSIStabDown")
load("/home/vassil/Documents/Bcells/Meta-analysis/ATHEROSCLEROSIS/ATHEROSCLEROSIStabUp")

tabUp <- ATHEROSCLEROSIStabUp
tabDown <- ATHEROSCLEROSIStabDown


# Split per cell type:
# Get non-significant genes:
# Check for regulation in opposite table
# Table with newly significant genes



genesNOsig <- tabUpo[tabUpo$P.Value>0.05, "Gene.ID"]
genesNOsig.u <- unique(genesNOsig)
tabUpo.non.sig <- tabUpo[tabUpo$Gene.ID %in% genes2,]

genes2 <- genesNOsig.u[!(genesNOsig.u %in% ATHEROSCLEROSIStabDown$Gene.ID)]
Epith <- grep ("^E.", tabUpo.non.sig$dataset)

tabUpo.non.sig.E <- tabUpo.non.sig[Epith,]

# Get a list of genes that appear as statistically non-significant (>0.05)
# Separate them based on dataset cell type
# Remove genes that appear as regulated in the opposite direction in
# same-cell-type datasets
# BP association of these genes and GO diagram (BPs and MFs)
# Select novel genes (not regulated by VD) and look up their disease function
