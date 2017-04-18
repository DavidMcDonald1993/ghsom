
library(GO.db)
library(topGO)
library(GOSim)
library(org.Sc.sgd.db)
library(igraph)

file <- "yeast_uetz"

ont <- "BP"
p <- 0.1
init <- 1
eta <- 0.0001

##load all community gene lists
setwd(sprintf("/home/david/Documents/ghsom/%s_communities_%s_%s", file, p, init, eta))

setOntology(ont, loadIC=TRUE)
setEvidenceLevel(evidences="all",organism=org.Sc.sgdORGANISM, gomap=org.Sc.sgdGO)
# calcICs()

db <- org.Sc.sgd.db
mapping <- "org.Sc.sgd.db"
ID <- "ENSEMBL"

#background gene list
backgroundFilename <- "all_genes.txt"
allGenes <- scan(backgroundFilename, character())

#load communities from file
g <- list()
numCom <- 0
filename <- sprintf("community_%s.txt", numCom)
while (file.exists(filename)) {
    numCom <- numCom + 1
    g[[numCom]] <- scan(filename, character())
    filename <- sprintf("community_%s.txt", numCom)
}

#distances between neurons
shortest.path <- read.csv("shortest_path.csv", sep=",", header=FALSE)

length(g)

calcICs()

allGeneNames <- scan(character(), file="../yeast_uetz_all_genes.txt")
allGenes <- allGeneNames[as.integer(allGenes)]
g <- sapply(g, function(i) allGeneNames[as.integer(i)])

enrichments <- sapply(g, function(i) GOenrichment(i, allGenes, cutoff=0.05, method="weight01"))

enrichedGOTerms <- function(genes, allGenes, cutoff, correction, ont, mapping, ID, algorithm){
    interestingGenes <- factor(as.integer(allGenes %in% genes))
    names(interestingGenes) <- allGenes
    
    GOdata <- new("topGOdata", description=sprintf("topGO object"),
              ontology = ont, allGenes = interestingGenes,
              annotationFun = annFUN.org, mapping = mapping, 
              ID = ID, nodeSize = 1)
    
    result <- runTest(GOdata, algorithm = algorithm, statistic = "fisher")
    if (correction){
        GOs <- score(result)[which(p.adjust(score(result), method="BH") <= cutoff)]
    } else {
        GOs <- score(result)[score(result) <= cutoff]
    }
    
    plot <- showSigOfNodes(GOdata, score(result), firstSigNodes = 10, useInfo ='all', swPlot = FALSE)
    
    return(list(GOdata, GOs, plot))
}

enrichedGOs  <- sapply(g, enrichedGOTerms, allGenes=allGenes, 
                      cutoff=0.05, correction=FALSE, ont=ont, mapping=mapping, ID=ID, algorithm="weight01")

p.values <- enrichments[2,]

lengths(p.values)

lengths(g)

shortest.path <- shortest.path[lengths(p.values) > 0, lengths(p.values) > 0]

g <- g[sapply(p.values, function(i) length(i) > 0)]
p.values <- p.values[sapply(p.values, function(i) length(i) > 0)]

minimumSubsumers <- function(gos) sapply(gos, function(i) sapply(gos, function (j) getMinimumSubsumer(i, j)))
allMinimumSubsumers <- sapply(p.values, function(i) 
    sort(summary(as.factor(minimumSubsumers(names(i)))), decreasing = TRUE))

allMinimumSubsumers

minGO <- sapply(p.values, function(i) names(i)[which.min(i)])
maxGO <- sapply(p.values, function(i) names(i)[which.max(i)])

select(GO.db, keys=minGO, columns=c("TERM","DEFINITION"))

select(GO.db, keys=maxGO, columns=c("TERM","DEFINITION"))

head(shortest.path)

fall <- function(i) !all(is.na(i))
fany <- function(i) !any(is.na(i))
    
l <- length(g)   

geneSims <- sapply(1:l, function(i) {
    sapply(i:l, function(j){
        if (i == j){
            return(1)
        } else {
            #gene sims as dataframe
            t <- getGeneSim(g[[i]], g[[j]], similarity="max", similarityTerm="Resnik", normalization=TRUE)
            ##remove na columns and rows
            t <- t[apply(t, 1, fall), apply(t, 2, fall)]
            ##remove any remaining rows with nan
            t <- t[apply(t, 1, fany),]
            ##BMA
            return((sum(apply(t, 1, max)) + sum(apply(t, 2, max))) / (nrow(t) + ncol(t)))
        }
       
    })
})

head(geneSims)

geneSimsDF <- sapply(1:length(g), function(i)
    sapply(1:length(g), function(j) {
        
        if (i <= j){
            x <- i
            y <- j - i + 1  
        } else {
            x <- j
            y <- i - j + 1
        }
        
        return(geneSims[[x]] [[y]]) }))

names <- character()
for (i in 1:length(g)){
    names <- c(names, sprintf("Com %s", i))
}

rownames(geneSimsDF) <- names
colnames(geneSimsDF) <- names
geneSimsDF <- round(geneSimsDF, 3)
geneSimsDF

library(gridExtra)
grid.table(geneSimsDF)

rownames(shortest.path) <- names
colnames(shortest.path) <- names
head(shortest.path)

shortest.path

library(gridExtra)
grid.table(shortest.path)

library(GOSemSim)

scGO <- godata(OrgDb = mapping, keytype = ID, ont = ont)

wangClusters <- mclusterSim(clusters = g, semData = scGO, measure = "Wang", combine = "BMA")

rownames(wangClusters) <- names
colnames(wangClusters) <- names

wangClusters

library(gridExtra)
grid.table(wangClusters)

resnikClusters <- mclusterSim(clusters = g, semData = scGO, measure = "Resnik", combine = "BMA")

rownames(resnikClusters) <- names
colnames(resnikClusters) <- names

library(gridExtra)
grid.table(resnikClusters)

distances <- numeric(length = (numCom * (numCom - 1)) / 2)
semSims <- numeric(length = (numCom * (numCom - 1)) / 2)

completed <- 0

for (c1 in 1:length(g)) {
    
    for (c2 in c1:length(g)) {
        
        if (c1 == c2) next   
        
        completed <- completed + 1  
        semSims[completed] <- wangClusters[c1, c2]
            
        distances[completed] <- shortest.path[c1, c2]
        
        print(sprintf("Completed: %s", completed))
    }
}

plot(distances, semSims, xlab="Distance on Map", ylab="Wang semantic similarity")


