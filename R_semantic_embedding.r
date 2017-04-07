
#libraries
library(GO.db)
library(topGO)
library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(GOSemSim)

file <- "yeast_union"

ont="BP"
db <- org.Sc.sgd.db
mapping <- "org.Sc.sgd.db"
ID <- "ENSEMBL"

##load all community gene lists
setwd(sprintf("/home/david/Documents/ghsom/"))

#background gene list
backgroundFilename <- sprintf("%s_all_genes.txt", file)
allGenes <- scan(backgroundFilename, character())

scGO <- godata(mapping, ont=ont, keytype=ID)

semanticDistances <- mgeneSim(allGenes, semData=scGO, measure="Resnik", combine="BMA", verbose=FALSE)

ncol(semanticDistances)

semanticDistances <- 1 - semanticDistances

head(semanticDistances)

semanticDistances2 <- sapply(allGenes, function(i) {
    sapply(allGenes, function(j) {
        geneSim(i, j, semData=scGO, measure="Wang", combine="BMA")["geneSim"]
    })
})

length(allGenes)

nrow(semanticDistances2)

semanticDistances2[is.na(semanticDistances2)] <- 0

write.table(semanticDistances, sep=",", file = sprintf("%s_wang_similarity.csv", file), row.names=FALSE, col.names=FALSE)

library(GOSim)

setOntology(ont, loadIC=FALSE)
setEvidenceLevel(evidences="all",organism=org.Sc.sgdORGANISM, gomap=org.Sc.sgdGO)

head(allGenes)

allGeneSims <- 1 - getGeneSim(allGenes, similarity="funSimMax", 
                          similarityTerm="relevance", normalization = TRUE)

allGeneSims[is.na(allGeneSims)] <- 1

head(allGeneSims)

write.table(allGeneSims, sep=",", file = sprintf("%s_rel_similarity_GOSim.csv", file), row.names=TRUE, col.names=FALSE)
