
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

write.csv(semanticDistances, file = sprintf("%s_resnik_similarity.csv", file))

wangSemanticDistances <- mgeneSim(allGenes, semData=scGO, measure="Wang", combine="BMA", verbose=FALSE)

wangSemanticDistances

write.csv(wangSemanticDistances, file = sprintf("%s_wang_similarity.csv", file))


