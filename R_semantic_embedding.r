
#libraries
library(GO.db)
library(topGO)
library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(GOSemSim)

file <- "yeast_uetz"

db <- org.Sc.sgd.db
mapping <- "org.Sc.sgd.db"
ID <- "ENSEMBL"

##load all community gene lists
setwd(sprintf("/home/david/Documents/ghsom/"))

#background gene list
backgroundFilename <- sprintf("%s_all_genes.txt", file)
allGenes <- scan(backgroundFilename, character())

scGO <- godata(mapping, ont="BP", keytype=ID)

semanticDistances <- mgeneSim(allGenes, semData=scGO, measure="Resnik", combine="BMA", verbose=FALSE)

semanticDistances

geneSim(allGenes, semData=scGO, measure="Wang", combine="BMA", verbose=FALSE)


