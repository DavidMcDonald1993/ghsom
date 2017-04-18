
#libraries
library(GO.db)
library(topGO)
library(org.Sc.sgd.db)
library(GOSemSim)
library(gridExtra)

file <- "yeast_uetz"

ont <- "BP"
p <- 0.5
init <- 1

db <- org.Sc.sgd.db
mapping <- "org.Sc.sgd.db"
ID <- "ENSEMBL"

##load all community gene lists
setwd(sprintf("/home/david/Documents/ghsom/%s_hierarchy_communities_%s_%s", file, p, init))

scGO <- godata(OrgDb = mapping, keytype = ID, ont = ont)

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

communitySimilarity <- function(community) {
    geneSims <- mgeneSim(genes = getGenes(as.character(community)), 
                         semData = scGO, measure = "Wang", combine = "BMA", verbose=F)
    if (length(geneSims) > 1) {
        return(mean(geneSims[upper.tri(geneSims)]))
    } else {
        return (NaN)
    }
}

layerSimilarity <- function(layer) {
    communitiesSimilarity <- sapply(unique(assignments[,layer][assignments[,layer] != -1]), communitySimilarity)
    communitiesSimilarity <- communitiesSimilarity[!is.na(communitiesSimilarity)]
    return(mean(communitiesSimilarity))
}

layerMeanSimilarities <- sapply(colnames, layerSimilarity)

layerMeanSimilarities

plot(colnames, layerMeanSimilarities, xlab="Layer", ylab="Mean Similarity", type = "l")

mapSimilarity <- function(mapID) {
    map <- getShortestPath(as.character(mapID))
    if (is.null(map)) return (NaN)
    communities <- sapply(rownames(map), function(rowname) getGenes(rowname))
    communitySimilarities <- mclusterSim(clusters = communities, semData = scGO, measure = "Wang", combine = "BMA")
    return(mean(communitySimilarities[upper.tri(communitySimilarities)]))
}

similarityOfAllMapsOnLayer <- function(layer) {
    mapSimilarities <- sapply(unique(assignments[,layer][assignments[,layer] != -1]), mapSimilarity)
    mapSimilarities <- mapSimilarities[!is.na(mapSimilarities)]
    return(mean(mapSimilarities))
}

layerMeanMapSimilarities <- sapply(colnames[1:length(colnames)-1], similarityOfAllMapsOnLayer)

layerMeanMapSimilarities

getDistance <- function(c1, c2) {
    
    m1 <- getSuperCommunity(c1)
    m2 <- getSuperCommunity(c2)
    
    if (m1 == m2) return(getShortestPath(m1)[c1, c2])
        
    d1 <- getDepth(c1)
    d2 <- getDepth(c2)
        
    if (d1 == d2) return (2 + getDistance(m1, m2))
        
    if (d1 > d2) return (1 + getDistance(m1, c2))
        
    if (d2 > d1) return (1 + getDistance(c1, m2))
    
}

distances <- sapply(as.character(2:max(assignments)), function(i)
    sapply(as.character(2:max(assignments)), function(j) {
        return(getDistance(i, j))
    }))

head(distances)

similarities <- sapply(2:max(assignments), function(i)
    sapply(2:max(assignments), function(j) {
        return(clusterSim(cluster1 = geneCommunities[[i]], 
                          cluster2 = geneCommunities[[j]], semData = scGO, measure = "Wang", combine = "BMA"))
    }))

head(similarities)


