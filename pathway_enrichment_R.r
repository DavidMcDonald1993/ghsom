
library(ReactomePA)
library(org.Sc.sgd.db)
library(clusterProfiler)
library(GOSim)
library(topGO)

file <- "yeast_reactome"

ont <- "BP"
e_sg <- 0.8
e_en <- 0.3

db <- org.Sc.sgd.db
mapping <- "org.Sc.sgd.db"
ID <- "ENTREZID"

##load all community gene lists
setwd(sprintf("/home/david/Documents/ghsom/hierarchical_exploration_10000/%s_hierarchy_communities_%s_%s", file, e_sg, e_en))
# setwd(sprintf("/home/david/Desktop/%s_hierarchy_communities_%s", file, e_sg))

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
assignments <- as.matrix(read.csv("assignment_matrix.csv", sep=",", header=F, row.names=1, colClasses="character"))
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

communities

length(communities)

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
                       
getAllSubCommunities <- function(com){
    
    subCommunities <- getSubCommunities(com)
    if (NA %in% subCommunities){
        return(NULL)
    }
    q <- as.list(subCommunities)
    allSubCommunities <- subCommunities
    
    while (length(q) > 0){
        com <- q[[1]]
        q <- q[-1]
        subCommunities <- getSubCommunities(com)
        if (!NA %in% subCommunities){
            q <- append(q, subCommunities)
            allSubCommunities <- append(allSubCommunities, subCommunities)
        }
    }
    
    return(allSubCommunities)
    
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

allGenes <- allGenes[!is.na(allGenes)]

enrichmentResultsFile <- "enrichmentResults.rda"
if (!file.exists(enrichmentResultsFile)) {
       enrichmentResults <- sapply(genesInCommunities, 
                            function (genes) enrichPathway(gene = genes, organism = "yeast", minGSSize = 5,  
                                                           pAdjustMethod = "none"))
       names(enrichmentResults) <- communities
       save(enrichmentResults, file=enrichmentResultsFile)
        print("saving")
} else {
       load(enrichmentResultsFile)    
        print("loaded")
}

df <- as.data.frame(enrichmentResults[["01"]])

head(df)

pathways <- df$Description

genes <- df$geneID

genes <- sapply(genes, function(i) strsplit(i, "/"))

genes <- sapply(genes, function(i) select(db, i, "UNIPROT", "ENTREZID")$UNIPROT)
names(genes) <- pathways    

genes

cat(sapply(pathways, toString), file="pathways.csv", sep="\n")

cat(sapply(genes, toString), file="pathway_genes.csv", sep="\n")

write.csv(genes, file="pathway_genes.csv", sep=",")

# nrow(as.data.frame(x))

# barplot(x, showCategory=10, title = "Top Enriched Pathways")

# dotplot(x, showCategory=15)

# enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)

numbersOfEnrichedPathways <- sapply(enrichmentResults, function(i) nrow(as.data.frame(i)))
enrichedCommunities <- genesInCommunities[numbersOfEnrichedPathways > 0 & lengths(genesInCommunities) > 3]

numbersOfEnrichedPathways

data.frame(numbersOfEnrichedPathways, lengths(genesInCommunities))

res <- compareCluster(enrichedCommunities, 
                      fun="enrichPathway", universe = allGenesInDB, organism = "yeast", minGSSize = 5,  
                                                           pAdjustMethod = "none")

png(filename=sprintf("community_pathway_enrichment_all_communities.png"), width=1500)
plot(res)
dev.off()

plotPathwayEnrichments <- function(community){
    
#     subCommunities <- getAllSubCommunities(community)
    subCommunities <- getSubCommunities(community)
    
    if (!is.null(subCommunities) && !NA %in% subCommunities) {

        communitiesOfInterest <- c(community, subCommunities)
        genesOfInterest <- enrichedCommunities[communitiesOfInterest]
        genesOfInterest <- genesOfInterest[!is.na(names(genesOfInterest))]
        
        if (length(genesOfInterest) > 1) {
            res <- compareCluster(genesOfInterest, 
            fun="enrichPathway", organism = "yeast", minGSSize = 5,  pAdjustMethod = "none")

            png(filename=sprintf("community_pathway_enrichment_%s.png", community), 
                width=500 + length(genesOfInterest) * 150)
            print(plot(res))
            dev.off()
        } 
        
    }
    
    print(sprintf("completed %s", community))

}

sapply(communities, plotPathwayEnrichments)

factorClusters <- lapply(enrichedCommunities, function(genes) {
    f <- factor(as.integer(allGenesInDB%in%genes))
    names(f) <- allGenesInDB
    return(f)
})

coms <- factorClusters[getSubCommunities("01")]
m1 <- compareCluster(coms[!is.na(names(coms))], fun="enrichGO", 
                     OrgDb=mapping, minGSSize = 5,  pAdjustMethod = "none")

compareCluster(c(getGenes("01-02-02-02"), getGenes("01-02-02-03")), fun="groupGO", OrgDb=mapping)

# viewPathway(pathName = "Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)", 
#             organism = "yeast", readable = F)

setOntology(ont, loadIC=TRUE)
setEvidenceLevel(evidences="all", organism=org.Sc.sgdORGANISM, gomap=org.Sc.sgdGO)

allGenesORF <- keys(db)


GOenrichmentResults <- sapply(genesInCommunities, function(genesOfInterest) {
    
    conversionTable <- select(db, genesOfInterest, "ORF", "ENTREZID")
    
    GOenrichment(conversionTable$ORF, allGenesORF, cutoff=0.05, method="weight01")
}
)

rownames(GOenrichmentResults) <- c("terms", "p-values", "genes")
colnames(GOenrichmentResults) <- communities

GOenrichmentResults[["terms", "01-05"]]

getSubCommunities("01-05")

GOenrichmentResults[["terms", "01-05-08"]]

getGenes("01-05")

library(GOSemSim)

scGO <- godata(OrgDb = mapping, keytype = ID, ont = ont)

clusterSimFile <- "clusterSimilarity.rda"
if (!file.exists(clusterSimFile)){
    clusterSim <- mclusterSim(clusters=genesInCommunities, semData=scGO)
    save(clusterSim, file=clusterSimFile)
    print("saved")
} else {
    load(clusterSimFile)
    print("loaded")
}

communitiesOfInterest <- c(getSubCommunities("01"), getSubCommunities("01-02"))

length(communitiesOfInterest)

clusterSim <- mclusterSim(clusters=genesInCommunities[communitiesOfInterest], semData=scGO)

clusterSim

library(gridExtra)
grid.table(clusterSim[sort(getSubCommunities("01")),sort(getSubCommunities("01"))])

grid.table(getShortestPath("01"))

grid.table(clusterSim[sort(getSubCommunities("01-02"))[1:6],sort(getSubCommunities("01-02"))[1:6]])

grid.table(getShortestPath("01-02")[1:6, 1:6])


