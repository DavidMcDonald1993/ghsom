
source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")

#libraries
library(GO.db)
library(topGO)
library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(GOSemSim)

db <- org.Sc.sgd.db
mapping <- "org.Sc.sgd.db"
ID <- "ENSEMBL"
# db <- org.Hs.eg.db
# mapping <- "org.Hs.eg.db"
# ID <- "ENTREZ"

#background gene list
setwd('/home/david/Documents/ghsom')
allGenes <- scan("Uetz_screen.txt", character())
# allGenes <- scan("Y2H_union.txt", character())
# allGenes <- scan("HI-II-14.txt", character())
allGenes <- unique(allGenes) 
length(allGenes)

allGenes <- select(org.Sc.sgd.db, keys=keys(org.Sc.sgd.db), columns="ORF")$ORF
length(allGenes)

##load all community gene lists
setwd("/home/david/Documents/ghsom/uetz_communities_04")
# setwd("/home/david/Documents/ghsom/union_communities_08")
# setwd("/home/david/Documents/ghsom/hi_communities_08")

g <- list()
numCom <- 0
filename <- sprintf("community_%s.txt", numCom)
while (file.exists(filename)) {
    numCom <- numCom + 1
    g[[numCom]] <- scan(filename, character())
    filename <- sprintf("community_%s.txt", numCom)
}
numCom

#distances between neurons
shortest.path <- read.csv("shortest_path.csv", sep=",", header=FALSE)

#factor of interesting genes
geneList <- factor(as.integer(allGenes %in% g[[3]]))
names(geneList) <- allGenes

#construct topGO object
GOdata <- new("topGOdata", description=sprintf("topGO object for community 1"),
              ontology = "BP", allGenes = geneList,
              annotationFun = annFUN.org, mapping = mapping, 
              ID = ID, nodeSize = 10)

#fishers exact test classic
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

gos <- score(resultFisher)[score(resultFisher) < 0.05]

#tabulate results
allRes <- GenTable(GOdata, classicFisher = resultFisher,
              orderBy = "classicFisher", topNodes = 500)

gos

select(GO.db, keys=names(gos), columns=c("TERM", "DEFINITION"))

generateTopGOData <- function(g){
    
}

filename <- sprintf("Uetz-04.rda")

if (file.exists(filename)){
    
    print(sprintf("loading: %s", filename))
    load(filename)
    
} else {
    
    cutOff <- 0.05

    geneLists <- vector("list", numCom) 
    GOdataObjects <- vector("list", numCom) 
    resultFishers <- vector("list", numCom) 
    results <- vector("list", numCom) 
    gos <- vector("list", numCom) 

    #perform enrichment analyses
    for (c in 1:numCom){

        #factor of interesting genes
        geneList <- factor(as.integer(allGenes %in% g[[c]]))
        names(geneList) <- allGenes
        geneLists[[c]] <- geneList

        #construct topGO object
        GOdata <- new("topGOdata", description=sprintf("topGO object for community %s", c),
                      ontology = "BP", allGenes = geneList,
                      annotationFun = annFUN.org, mapping = mapping, 
                      ID = ID, nodeSize = 10)
        GOdataObjects[[c]] <- GOdata

        #fishers exact test classic
        resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        resultFishers[[c]] <- resultFisher

        #tabulate results
        allRes <- GenTable(GOdata, classicFisher = resultFisher,
                      elimFisher = resultFisher.elim,
                      orderBy = "classicFisher", topNodes = 500)
        results[[c]] <- allRes

        #go terms < cut off 
        gos[[c]] <- score(resultFisher)[score(resultFisher) < cutOff]

        print(sprintf("community %s complete", c))
    }
    
    print(sprintf("Saving data: %s", filename))
    save(geneLists, GOdataObjects, resultFishers, results, gos, file=filename)
}

##SEMATIC SIMILARITY
#construct gosemsim object
semsimfile <- sprintf("Uetz-semsimfile.rda")
if (file.exists(semsimfile)){
    load(semsimfile)
} else {
    hsGO <- godata(mapping, ont="BP", keytype="ENSEMBL")
    save(hsGO, file=semsimfile)
    print(sprintf("saved semsimfile: %s", semsimfile))
}


information_content <- function(term){
    return (goSim(term, term, semData=hsGO, measure="Resnik"))
}

most_representative_term_ic <- function(namedTerms){
    ics <- sapply(names(namedTerms), information_content)
    names(ics) <- names(namedTerms)
    return(names(sort(ics, decreasing=TRUE)[1]))
}

most_representative_term_ancestor <- function(namedTerms){
    
    counts <- numeric(length(namedTerms))
    names(counts) <- names(terms)

    for (term in terms) {
        ancestors <- as.list(GOBPANCESTOR[names(term)])
        for (ancestor in ancestors[[names(term)]]) {
            if (ancestor %in% names(counts)) {
                counts[ancestor] <- counts[ancestor] + 1
            }
        }

    }
#     return (sort(names(counts), decreasing=TRUE)[1])
    return (sort(counts / sum(counts), decreasing=TRUE))
}

representativeTerms <- sapply(gos, most_representative_term_ic)

representativeTerms

select(GO.db, keys=representativeTerms, columns=c("TERM", "DEFINITION"))

shortest.path

representativeTerms[7]

representativeTerms[8]

goSim(representativeTerms[7], representativeTerms[8], semData=hsGO, measure="Wang")

sims <- mgoSim(representativeTerms, representativeTerms, semData=hsGO, measure="Resnik", combine=NULL)

sims <- mclusterSim(g, semData=hsGO, measure="Rel", combine="BMA")

head(sims)

head(shortest.path)

simfile <- sprintf("HI-II-14-sims.rda")
if (file.exists(simfile)){
    load(simfile)
} else {
    sims <- mclusterSim(g, semData=hsGO, measure="Wang", combine="BMA")
    save(sims, file=simfile)
    print (sprintf("saved sim file: %s", simfile))
}


weighted_similarity <- function(namedTerms1, namedTerms2) {
    

    s <- mgoSim(names(namedTerms1), 
           names(namedTerms2), semData=hsGO, measure="Wang", combine=NULL)

    t <- matrix(numeric(), 
                nrow=length(namedTerms1), 
                ncol=length(namedTerms2))
    
    for (t1 in 1:length(namedTerms1)){
        w1 <- namedTerms1[t1]
        term1 <- names(w1)
        for (t2 in 1:length(namedTerms2)){
            w2 <- namedTerms2[t2]
            term2 <- names(w2)
            t[[t1, t2]] <- s[term1, term2] * w1 * w2
        }
    }
         
    return(max(t))
    
}


t <- matrix(numeric(), nrow=numCom, ncol=numCom)
for (t1 in 1:numCom) {
    namedTerms1 <- normalisedRepresentativeTerms[[t1]]
    for (t2 in 1:numCom) {
        namedTerms2 <- normalisedRepresentativeTerms[[t2]]
        t[[t1, t2]] <- weighted_similarity(namedTerms1, namedTerms2)
    }
}
# rownames(t) <- representativeTerms
# colnames(t) <- representativeTerms
head(t)

semSimTable <- mgoSim(representativeTerms, representativeTerms, semData=hsGO, measure="Wang", combine=NULL)

t <- matrix(numeric(), nrow=numCom, ncol=numCom)
for (t1 in 1:numCom) {
    term1 <- representativeTerms[t1]
    for (t2 in 1:numCom) {
        term2 <- representativeTerms[t2]
        t[[t1, t2]] <- semSimTable[term1, term2]
    }
}
rownames(t) <- representativeTerms
colnames(t) <- representativeTerms
head(t)

head(t * 100)

head(shortest.path)

distances <- numeric(length = (numCom * (numCom - 1)) / 2)
semSims <- numeric(length = (numCom * (numCom - 1)) / 2)

completed <- 0

for (c1 in 1:numCom) {
    
    t1 <- representativeTerms[c1]
#     gs1 <- g[[c1]]
#     if (length(gos[[c1]]) == 0) next
    
    for (c2 in c1:numCom) {
        
        if (c1 == c2) next
            
            t2 <- representativeTerms[c2]
            
#         if (length(gos[[c2]]) == 0) next
            
#         gs2 <- g[[c2]]    
        
        completed <- completed + 1  
        
        #compute semantic similarity of two protein clusters
#         semSims[completed] <- clusterSim(gs1, gs2, semData=scGO, measure="Wang", combine="BMA")
#         semSims[completed] <- mgoSim(gos[[c1]], gos[[c2]], semData=scGO, measure="Wang", combine="BMA")
#         semSims[completed] <- semSimTable[t1, t2]
        semSims[completed] <- t[c1, c2]
            
        distances[completed] <- shortest.path[c1, c2]
        
        print(sprintf("Completed: %s", completed))
    }
}
# distances <- distances[distances > 0]
# semSims <- semSims[semSims > 0]

plot(distances, semSims, xlab="Distance on Map", ylab="Semantic Similarity")

cor(distances, semSims)

source("https://bioconductor.org/biocLite.R")
biocLite("GOSim")

sim

library(GOSim)
setOntology("BP")
gomap <- get("gomap",env=GOSimEnv)
allgenes = sample(names(gomap), 1000) # suppose these are all genes
genesOfInterest = sample(allgenes, 20) # suppose these are all genes of interest
sim = getGeneSim(genesOfInterest,verbose=FALSE) # and these are their similarities
hc = hclust(as.dist(1-sim), method="ward") # use them to perform a clustering
plot(hc) # plot the cluster tree
cl = cutree(hc, k=3) # take 3 clusters
if(require(cluster)){
    ev = evaluateClustering(cl, sim) # evaluate the clustering
    print(ev$clusterstats) # print out some statistics
    plot(ev$clustersil,main="") # plot the cluster silhouettes
}
# investigate cluster 1 further
if(require(topGO))
    GOenrichment(genesOfInterest[cl == 1], allgenes, cutoff=0.05) # print out what cluster 1 is about

sim = getGeneSim(g[[1]], g[[2]], verbose=FALSE)

GOenrichment(g[[1]], allGenes, cutoff=0.05) # print out what cluster 1 is about


