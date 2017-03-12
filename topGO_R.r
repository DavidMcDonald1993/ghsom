
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("ALL")
biocLite("topGO")
biocLite("hgu95av2.db")
biocLite("Rgraphviz")
biocLite("colMap")

library(topGO)
library(ALL)
data(ALL)
data(geneList)

affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)

sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)

sampleGOdata

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultFisher

resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS

resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
resultKS.elim

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
allRes

pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)

sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
cbind(termStat(sampleGOdata, sel.go),
      elim = pValue.elim[sel.go],
      classic = pValue.classic[sel.go])

showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo ='all')


