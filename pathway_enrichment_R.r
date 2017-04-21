
library(ReactomePA)
library(org.Sc.sgd.db)
library(clusterProfiler)

file <- "yeast_reactome"

# ont <- "BP"
e_sg <- 0.3
e_en <- 0.2

db <- org.Sc.sgd.db
# mapping <- "org.Sc.sgd.db"
# ID <- "ENSEMBL"

##load all community gene lists
setwd(sprintf("/home/david/Documents/ghsom/%s_hierarchy_communities_%s_%s", file, e_sg, e_en))

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
assignments <- as.matrix(read.csv("assignment_matrix.csv", sep=",", header=F, colClasses="character"))
assignments[assignments == ""] <- NA 
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
    if (depth < ncol(genesInCommunity)){
        return(as.character(unique(genesInCommunity[,depth + 1])))
    } else {
        return (NULL)
    }
    
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

genesInCommunities <- sapply(communities, function(i) getGenes(i))

communities

length(communities)

lengths(genesInCommunities)

enrichmentResults <- sapply(genesInCommunities, 
                            function (i) enrichPathway(gene = i, universe = allGenesInDB, organism = "yeast"))
names(enrichmentResults) <- communities

x <- enrichmentResults[["01-05"]]

nrow(as.data.frame(x))

barplot(x, showCategory=10, title = "Top Enriched Pathways")

dotplot(x, showCategory=15)

enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)

numbersOfEnrichedPathways <- sapply(enrichmentResults, function(i) nrow(as.data.frame(i)))
enrichedCommunities <- genesInCommunities[numbersOfEnrichedPathways > 0]

res <- compareCluster(enrichedCommunities, 
                      fun="enrichPathway", universe = allGenesInDB, organism = "yeast")

png(filename=sprintf("community_pathway_enrichment_all_communities.png"), width=1500)
plot(res)
dev.off()

plotPathwayEnrichments <- function(community){
    
    subCommunities <- getSubCommunities(community)
    
    if (!is.null(subCommunities) && !any(is.na(subCommunities) > 0)) {

        communitiesOfInterest <- c(community, subCommunities)
        print(communitiesOfInterest)
        genesOfInterest <- enrichedCommunities[communitiesOfInterest]
        genesOfInterest <- genesOfInterest[!is.na(names(genesOfInterest))]
        print (genesOfInterest)
        
        if (length(genesOfInterest) > 2) {
            res <- compareCluster(genesOfInterest,
            fun="enrichPathway", universe = allGenesInDB, organism = "yeast")

            png(filename=sprintf("community_pathway_enrichment_%s.png", community), width=500 + length(genesOfInterest) * 150)
            print(plot(res))
            dev.off()
        } 
        
    }

}

sapply(communities, plotPathwayEnrichments)

viewPathway(pathName = "Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)", 
            organism = "yeast", readable = F)


