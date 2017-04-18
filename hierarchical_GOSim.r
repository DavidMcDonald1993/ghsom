
library(GO.db)
library(topGO)
library(GOSim)
library(org.Sc.sgd.db)
library(igraph)

file <- "yeast_uetz"

ont <- "BP"
p <- 0.5
init <- 1

db <- org.Sc.sgd.db
mapping <- "org.Sc.sgd.db"
ID <- "ENSEMBL"

##load all community gene lists
setwd(sprintf("/home/david/Documents/ghsom/%s_hierarchy_communities_%s_%s", file, p, init))

setOntology(ont, loadIC=TRUE)
setEvidenceLevel(evidences="all",organism=org.Sc.sgdORGANISM, gomap=org.Sc.sgdGO)

generateMap <- function(filename){
    map <- as.matrix(read.csv(filename, sep=",", header = F))
    communities <- map[,1]
    map <- map[,2:ncol(map)]
    rownames(map) <- communities
    colnames(map) <- communities
    return (map)
}

#background gene list
backgroundFilename <- "all_genes.txt"
allGenes <- scan(backgroundFilename, character())

#shortest path files
shortestPathFiles  <- list.files(pattern="*shortest_path*")

#shortest paths list
shortestPaths <- sapply(shortestPathFiles, generateMap)
names(shortestPaths) <- sapply(names(shortestPaths), function(name) strsplit(name, "_")[[1]][[1]])

#communitiy assignemtns
assignments <- as.matrix(read.csv("assignment_matrix.csv", sep=",", header=F))
rownames(assignments) <- allGenes
colnames <- sapply(1:ncol(assignments), function(i) as.character(i-1))
colnames(assignments) <- colnames

getDepth <- function(com) {
    return(which(apply(assignments, 2, function(i) any(i == com))))
}

getGenes <- function(com){
    return(names(which(assignments[,getDepth(com)] == com)))
}

getSubCommunities <- function(com){
    return(try(as.character(unique(assignments[getGenes(com), getDepth(com) + 1]))))
}

getSuperCommunity <- function(com){
    return(try(as.character(unique(assignments[getGenes(com), getDepth(com) - 1]))))
}

getShortestPath <- function(com){
    return (try(shortestPaths[[com]]))
}

allGenesInDB <- keys(db)
allGenes <- allGenes[allGenes %in% allGenesInDB]
enrichmentResults <- sapply(1:max(assignments), function(i) {

    genesOfInterest <- getGenes(i)
    genesOfInterest <- genesOfInterest[genesOfInterest %in% allGenesInDB]
    GOenrichment(genesOfInterest, allGenesInDB, cutoff=0.05, method="weight01")
}
)

rownames(enrichmentResults) <- c("terms","p-values","genes")
colnames(enrichmentResults) <- 2:max(assignments)

communitySimilarity <- function(community) {
    termSims <- getTermSim(termlist = names(community), method = "Lin", verbose = F)
    if (length(termSims) > 1) {
        return(mean(termSims[upper.tri(termSims)]))
    } else {
        return (NaN)
    }
}

communitySimilarity(enrichmentResults[["p-values", "27"]])

layerSimilarity <- function(layer) {
    pvalueList <- enrichmentResults["p-values", unique(assignments[,layer][assignments[,layer] != -1]) - 1]
    communitiesSimilarity <- sapply(pvalueList, communitySimilarity)
    communitiesSimilarity <- communitiesSimilarity[!is.na(communitiesSimilarity)]
    return(mean(communitiesSimilarity))
}

layerMeanSimilarities <- sapply(colnames, layerSimilarity)

layerMeanSimilarities

geneCommunities <- sapply(1:max(assignments), function (i) getGenes(i)[getGenes(i) %in% allGenesInDB])

getSubCommunities(6)

as.list(org.Sc.sgdPATH[geneCommunities[[6]]])

as.list(org.Sc.sgdPATH2ORF[["00330"]]) %in% allGenes

length(allGenes)

geneCommunities[[1]]

allGenes <- allGenes[allGenes%in% allGenesInDB]

GOenrichment(allGenes, allGenesInDB, cutoff = 0.01, method = "weight01")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("ReactomePA")

library(ReactomePA)
data(geneList)
de <- names(geneList)[abs(geneList) > 1.5]
head(de)

x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

entrezCommunities <- sapply(1:max(assignments), function(i){
    orfs <- getGenes(i)
    orfs <- orfs[orfs%in%allGenesInDB]
    return(as.character(org.Sc.sgdENTREZID[orfs]))
})

entrezCommunities[[1]]

columns(db)

pathwayEnrichments <- sapply(entrezCommunities[2:length(entrezCommunities)],
                             function(i) enrichPathway(gene=i, organism = "yeast", universe = entrezCommunities[[1]], 
                                                                            pvalueCutoff = 0.05, readable = T))

library(ReactomePA)
getDb("human")

orfs <- getGenes(2)
orfs <- orfs[orfs%in%allGenesInDB]
entrez <- as.character(org.Sc.sgdENTREZID[orfs])

orfs

entrez

allGenes <- allGenes[allGenes%in%allGenesInDB]
entrezAll <- as.character(org.Sc.sgdENTREZID[allGenes])
entrezAllDB <- as.character(org.Sc.sgdENTREZID[allGenesInDB])

x <- enrichPathway(gene = entrezAll, organism = "yeast", universe = entrezAllDB)

head(as.data.frame(x))

print(x)

unionAllGenes <- scan(character(), file="../yeast_union_all_genes.txt")

unionAllGenes <- unionAllGenes[unionAllGenes%in%allGenesInDB]

x <- enrichPathway(gene = as.character(org.Sc.sgdENTREZID[unionAllGenes]), organism = "yeast", universe = entrezAllDB)

print (x)

source("http://bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
install.packages("WGCNA") 

library(WGCNA)

allowWGCNAThreads()

library(org.Sc.sgd.db)
uniProt <- org.Sc.sgdUNIPROT

uniProt[["P04076"]]


