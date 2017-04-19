
library(ReactomePA)
library(org.Sc.sgd.db)
library(clusterProfiler)

file <- "yeast_reactome"

# ont <- "BP"
p <- 0.2

db <- org.Sc.sgd.db
# mapping <- "org.Sc.sgd.db"
# ID <- "ENSEMBL"

##load all community gene lists
setwd(sprintf("/home/david/Documents/ghsom/%s_hierarchy_communities_%s", file, p))

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

##convert all genes to ENTREZID
conversion <- select(db, allGenes, "ENTREZID", "UNIPROT")
conversion <- subset(conversion, !duplicated(conversion$UNIPROT))
allGenes <- conversion$ENTREZID

#shortest path files
shortestPathFiles  <- list.files(pattern="*shortest_path*")

#shortest paths list
shortestPaths <- sapply(shortestPathFiles, generateMap)
names(shortestPaths) <- sapply(names(shortestPaths), function(name) strsplit(name, "_")[[1]][[1]])

#communitiy assignemtns
assignments <- as.matrix(read.csv("assignment_matrix.csv", sep=",", header=F))
rownames(assignments) <- allGenes
colnames <- sapply(1:ncol(assignments), function(i) as.character(i - 1))
colnames(assignments) <- colnames
    
#filter out genes with no ENTREZID
assignments <- assignments[!is.na(rownames(assignments)),]
    
#all ORF identifers in org.Sc.sgd.db converted to EntrezID
allGenesInDB <- select(db, keys(db), "ENTREZID", "ORF")$ENTREZID
allGenesInDB <- allGenesInDB[!is.na(allGenesInDB)]
    
#communities detected
communities <- unique(as.character(assignments))
communities <- communities[communities != ""]
communities <- sort(communities)

communities

getDepth <- function(com) {
    return(which(apply(assignments, 2, function(i) any(i == com))))
}

getGenes <- function(com){
    depth <- getDepth(com)
    return(names(which(assignments[, depth] == com)))
}

getSubCommunities <- function(com){
    depth <- getDepth(com)
    genesInCommunity <- subset(assignments, assignments[,depth] == com)
    return(as.character(unique(genesInCommunity[,depth + 1])))
}

getSuperCommunity <- function(com){
    depth <- getDepth(com)
    genesInCommunity <- subset(assignments, assignments[,depth] == com)
    return(as.character(unique(genesInCommunity[,depth - 1])))
}

getShortestPath <- function(com){
    return (try(shortestPaths[[com]]))
}
                       
getNeighbours <- function(com){
    
    superCommunity <- getSuperCommunity(com)
    superCommunityMap <- getShortestPath(superCommunity)
    v <- superCommunityMap[com, ] == 1
    return (names(v[v]))
    
}

getShortestPath("1")

genesInCommunities <- sapply(communities, function(i) getGenes(i))

lengths(genesInCommunities)

enrichmentResults <- sapply(genesInCommunities, 
                            function (i) enrichPathway(gene = i, universe = allGenesInDB, organism = "yeast"))
names(enrichmentResults) <- communities

getSubCommunities("1")

x <- enrichmentResults[["1-5"]]

nrow(as.data.frame(x))

barplot(x, showCategory=10, title = "Top 10 Enriched Pathways")

dotplot(x, showCategory=15)

enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)

numbersOfEnrichedPathways <- sapply(enrichmentResults, function(i) nrow(as.data.frame(i)))
enrichedClusters <- genesInCommunities[numbersOfEnrichedPathways > 0]

res <- compareCluster(genesInCommunities[numbersOfEnrichedPathways > 0], 
                      fun="enrichPathway", universe = allGenesInDB, organism = "yeast")

png(filename=sprintf("cluster_pathway_enrichment_overall_%s.png", p), width = 2000)
plot(res)
dev.off()

com  <-  "4"

depth <- getDepth(com)
genesInCommunity <- subset(assignments, assignments[,depth] == com)
length(unique(genesInCommunity[,depth + 1]))

is.null(getSubCommunities("4"))

plotPathwayEnrichments <- function(community){
    
    subCommunities <- getSubCommunities(community)
    
    if (length(subCommunities) > 0) {

        communitiesOfInterest <- c(community, subCommunities)

        res <- compareCluster(genesInCommunities[communitiesOfInterest], 
                      fun="enrichPathway", universe = allGenesInDB, organism = "yeast")
        
        png(filename=sprintf("cluster_pathway_enrichment_%s.png", community), width = 2000)
        print(plot(res))
        dev.off()
        
    }

}

plotPathwayEnrichments("1")

viewPathway(pathName = "Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)", 
            organism = "yeast", readable = F)


