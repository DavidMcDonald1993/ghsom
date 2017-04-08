
library(GO.db)
library(topGO)
library(GOSim)
library(org.Sc.sgd.db)
library(igraph)

file <- "yeast_union"

ont <- "MF"
p <- 0.5
init <- 1


setOntology(ont, loadIC=FALSE)
setEvidenceLevel(evidences="all",organism=org.Sc.sgdORGANISM, gomap=org.Sc.sgdGO)

db <- org.Sc.sgd.db
mapping <- "org.Sc.sgd.db"
ID <- "ENSEMBL"

##load all community gene lists
setwd(sprintf("/home/david/Documents/ghsom/%s_communities_%s_%s", file, p, init))

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

numCom

allGeneNames <- scan(character(), file="../yeast_uetz_all_genes.txt")
allGenes <- allGeneNames[as.integer(allGenes)]
g <- sapply(g, function(i) allGeneNames[as.integer(i)])

enrichments <- sapply(g, function(i) GOenrichment(i, allGenes, cutoff=0.05, method="elim"))

p.values <- enrichments[2,]

G <- getGOGraph(names(p.values[[1]]))
G2  <-  igraph.from.graphNEL(G)
plot(G2)

allEnrichedGenes <- enrichments[3,]

getEnrichedGenes <- function(enrichedGenes){
    v <- character()
    for (i in enrichedGenes){
        v <- c(v, i[1])
    }
    return(v)
}

clusters <- sapply(allEnrichedGenes, function(i) getEnrichedGenes(i))

clusterSim <- mclusterSim(clusters = clusters, semData = scGO, measure = "Wang", combine = "BMA")

head(clusterSim)
head(shortest.path)

lengths(p.values)

shortest.path <- shortest.path[lengths(p.values) > 0, lengths(p.values) > 0]

p.values <- p.values[sapply(p.values, function(i) length(i) > 0)]

lengths(p.values)

minGO <- sapply(p.values, function(i) names(i)[which.min(i)])
maxGO <- sapply(p.values, function(i) names(i)[which.max(i)])

select(GO.db, keys=minGO, columns=c("TERM","DEFINITION"))

select(GO.db, keys=maxGO, columns=c("TERM","DEFINITION"))

head(shortest.path)

colnames(shortest.path) <- maxGO
rownames(shortest.path) <- maxGO

library(GOSemSim)

scGO <- godata(ont = ont, OrgDb = mapping, keytype = ID)

semSimGO <- mgoSim(maxGO, maxGO, semData=scGO, measure="Wang", combine=NULL)

head(semSimGO)
head(shortest.path)

termSim <- getTermSim(maxGO, method = "Resnik")

head(termSim)
head(shortest.path)

#gene sims as dataframe
t <- getGeneSim(g[[1]], g[[2]], similarity="max", similarityTerm="Resnik", normalization=TRUE)

t

fall <- function(i) !all(is.na(i))
fany <- function(i) !any(is.na(i))
##remove na columns and rows
s <- t[apply(t, 1, fall), apply(t, 2, fall)]
##remove any remaining rows with nan
s <- s[apply(s, 1, fany),]

s

##BMA
((sum(apply(s, 1, max)) + sum(apply(s, 2, max))) / (nrow(s) + ncol(s)))

fall <- function(i) !all(is.na(i))
fany <- function(i) !any(is.na(i))
    
l = 2    

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

head(shortest.path)

geneSimsDF <- sapply(1:length(geneSims), function(i) {
    sapply(1:length(geneSims), function(j) {
        if (j >= i){
            mean(geneSims[[i]][[j]])
        } else{
            mean(geneSims[[j]][[i]])
        }
        
    })
})

i <- 1
j <- 3
getGeneSim(g[[i]], g[[j]], similarity="max", similarityTerm="relevance")

maxGoSimilarity <- function(gos1, gos2) {
    max(sapply(gos1, function(i) {
        sapply(gos2, function(j) {
            getTermSim(c(i, j), method="relevance")[1, 2]
        })
    }))
}

clusterGoSimilarities <- sapply(1:length(p.values), function(i) {
    sapply(1:length(p.values), function(j){
        maxGoSimilarity(names(p.values[[i]]), names(p.values[[j]]))
    })
})

head(clusterGoSimilarities)

head(shortest.path)

distances <- numeric(length = (numCom * (numCom - 1)) / 2)
semSims <- numeric(length = (numCom * (numCom - 1)) / 2)

completed <- 0

for (c1 in 1:length(g)) {
    
    for (c2 in c1:length(g)) {
        
        if (c1 == c2) next   
        
        completed <- completed + 1  
        semSims[completed] <- clusterGoSimilarities[c1, c2]
            
        distances[completed] <- shortest.path[c1, c2]
        
        print(sprintf("Completed: %s", completed))
    }
}

plot(distances, semSims, xlab="Distance on Map", ylab="GO similairties")


