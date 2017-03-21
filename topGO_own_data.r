
library(GO.db)
library(topGO)
library(org.Sc.sgd.db)

columns(GO.db)

keys <- sample(keys(GO.db), 5)

keys

?sample

select(GO.db, keys=keys, keytype="GOID", columns=c("TERM", "DEFINITION"))

setwd('/home/david/Documents/ghsom')
allGenes <- read.table("Y2H_union.txt")$V1
allGenes <- unique(allGenes)
length(allGenes)

setwd('/home/david/Documents/ghsom')
myInterestingGenes <- read.table("community_one.txt")$V1
myInterestingGenes

geneList <- factor(as.integer(allGenes %in% myInterestingGenes))
# names(geneList) <- ensembl2GO[allGenes]
# names(geneList) <- ensembl2ORF[allGenes]
names(geneList) <- allGenes
geneList

GOdata <- new("topGOdata", description="Practice GOData object UNION", ontology = "BP", allGenes = geneList,
              annotationFun = annFUN.org, mapping = "org.Sc.sgd.db", ID = "ENSEMBL", nodeSize = 10)
GOdata

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher

resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
resultFisher.elim

allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   elimFisher = resultFisher.elim,
                   orderBy = "classicFisher", ranksOf = "elimFisher", topNodes = 10)
allRes

showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')

graph(GOdata)

setwd('/home/david/Documents/ghsom')
myInterestingGenes2 <- read.table("community_two.txt")$V1
myInterestingGenes2

geneList2 <- factor(as.integer(allGenes %in% myInterestingGenes2))
# names(geneList) <- ensembl2GO[allGenes]
# names(geneList) <- ensembl2ORF[allGenes]
names(geneList2) <- allGenes
geneList2

GOdata2 <- new("topGOdata", description="Practice GOData object UNION", ontology = "BP", allGenes = geneList2,
              annotationFun = annFUN.org, mapping = "org.Sc.sgd.db", ID = "ENSEMBL", nodeSize = 10)
GOdata2

resultFisher2 <- runTest(GOdata2, algorithm = "classic", statistic = "fisher")
resultFisher2

resultFisher2.elim <- runTest(GOdata2, algorithm = "elim", statistic = "fisher")
resultFisher2.elim

allRes2 <- GenTable(GOdata2, classicFisher = resultFisher,
                   elimFisher = resultFisher.elim,
                   orderBy = "classicFisher", ranksOf = "elimFisher", topNodes = 10)
allRes2

showSigOfNodes(GOdata2, score(resultFisher2), firstSigNodes = 5, useInfo ='all')

printGraph(GOdata, resultFisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

printGraph(GOdata2, resultFisher2, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

g1 <- names(score(resultFisher)[score(resultFisher) < 0.05])
g2 <- names(score(resultFisher2)[score(resultFisher2) < 0.05])

source("https://bioconductor.org/biocLite.R")
biocLite("GOSemSim")

library(GOSemSim)
scGO <- godata('org.Sc.sgd.db', ont="BP")

semSim <- mgoSim(g1, g2, semData=scGO, measure="Wang", combine=NULL)
mean(semSim)

semSim


