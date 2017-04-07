
library(topGO)
library(GOSim)
library(org.Sc.sgd.db)

file <- "yeast_uetz"

ont <- "BP"
p <- 0.8
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

length(g)

enrichments <- sapply(g, function(i) GOenrichment(i, allGenes, cutoff=0.01, method="elim"))

p.values <- enrichments[2,]

names((p.values[[1]]))

lengths(p.values)

clusterSimilarity <- sapply(p.values, function(i) mean(getTermSim(names(i), method="relevance")))

clusterSimilarity

geneSims12 <- getGeneSim(g[[1]], g[[2]], similarity="max", similarityTerm="Resnik", normalization=TRUE)

##remove na rows and columns
f <- function(i) !all(is.na(i))
geneSims12 <- geneSims12[apply(geneSims12,1,f), apply(geneSims12, 2, f)]

fany <- function(i) !any(is.na(i))
    
m <- matrix(c(1,NaN,3,4,5,6), nrow=3)

df <- data.frame(m)
df <- df[apply(df, 1, fany),]
df
apply(df, 1, fany)

fall <- function(i) !all(is.na(i))
fany <- function(i) !any(is.na(i))
    
geneSims <- sapply(1:length(g), function(i) {
    sapply(i:length(g), function(j){
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
            return(sum(apply(t, 1, max) + sum(apply(t, 2, max))) / (nrow(t) + ncol(t)))
        }
       
    })
})

head(geneSims)

apply(shortest.path, 1, max)

(sum(apply(shortest.path, 1, max)) + sum(apply(shortest.path, 2, max))) / 12

typeof(shortest.path)

head(shortest.path)

geneSims <- sapply(geneSims, function(i){
    sapply(i, function(j){
        j[is.na(j)] <- 0
    })
})

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


