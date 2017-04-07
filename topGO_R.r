
library(topGO)
library(org.Sc.sgd.db)

t <- org.Sc.sgd.db

columns(t)

g1 <- sample(keys(t), 10)
g2 <- sample(keys(t), 10)

g1

g2

library(GOSemSim)
scGO <- godata("org.Sc.sgd.db", ont="BP", keytype="ENSEMBL")

library(topGO)

allGenes <- keys(t)

geneList <- factor(as.integer(allGenes %in% g1))
names(geneList) <- keys(t)

GOdata <- new("topGOdata", description=sprintf("topGO object"),
              ontology = "BP", allGenes = geneList,
              annotationFun = annFUN.org, mapping = "org.Sc.sgd.db", 
              ID = "ENSEMBL", nodeSize = 10)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

sampson <- load("sampson.RData")

sampson

Sampson

data(hansell)

install.packages(hansell)


