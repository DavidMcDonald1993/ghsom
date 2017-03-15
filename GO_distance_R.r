
#libraries
library(GO.db)
library(topGO)
library(org.Sc.sgd.db)

#background gene list
setwd('/home/david/Documents/ghsom')
allGenes <- scan("Y2H_union.txt", character())
allGenes <- unique(allGenes) 

##load all community gene lists
setwd("/home/david/Documents/ghsom/union_communities_08")

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

library(ontologyIndex)
go <- get_ontology("/home/david/Documents/ghsom/db/go-basic.obo")

find_representative_term <- function(terms){
    counts <- numeric(length(terms))
    names(counts) <- terms

    for (term in terms) {
        ancestors <- get_term_property(go, "ancestors", term, as_names = FALSE)
        for (ancestor in ancestors) {
            counts[ancestor] <- counts[ancestor] + 1
        }
    }
    return (sort(names(counts), decreasing=TRUE)[1])
}

cutOff <- 0.05

geneLists <- vector("list", numCom) 
GOdataObjects <- vector("list", numCom) 
resultFishers <- vector("list", numCom) 
resultFisher.elims <- vector("list", numCom) 
results <- vector("list", numCom) 
topResults <- vector("list", numCom) 
gos <- vector("list", numCom) 
representativeTerms <- character(length = numCom)

#perform enrichment analyses
for (c in 1:numCom){
    
    #factor of interesting genes
    geneList <- factor(as.integer(allGenes %in% g[[c]]))
    names(geneList) <- allGenes
    geneLists[[c]] <- geneList
    
    #construct topGO object
    GOdata <- new("topGOdata", description=sprintf("topGO object for community %s", c),
                  ontology = "BP", allGenes = geneList,
                  annotationFun = annFUN.org, mapping = "org.Sc.sgd.db", 
                  ID = "ENSEMBL", nodeSize = 10)
    GOdataObjects[[c]] <- GOdata
    
#     #fishers exact test classic
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    resultFishers[[c]] <- resultFisher
    
#     #fishers exact test elimination
    resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
    resultFisher.elims[[c]] <- resultFisher.elim
    
    #tabulate results
    allRes <- GenTable(GOdata, classicFisher = resultFisher,
                  elimFisher = resultFisher.elim,
                  orderBy = "classicFisher", topNodes = 500)
    results[[c]] <- allRes
    
    #go terms <0.01 on both tests
    topResults[[c]] <- subset(allRes, classicFisher < cutOff & elimFisher < cutOff)
    gos[[c]] <- subset(allRes, classicFisher < cutOff & elimFisher < cutOff)$GO.ID
    
    #term that is ancestor of most terms
    representativeTerms[c] <- find_representative_term(gos[[c]])
}

dir.create("/home/david/Documents/ghsom/uetz_go_terms")
setwd("/home/david/Documents/ghsom/uetz_go_terms")
for (c in 1:numCom){
    write.csv(topResults[[c]], sprintf("go_terms_%s", c))
    print(sprintf("saved terms %s", c))
} 

##semantic similarity
library(GOSemSim)
scGO <- godata('org.Sc.sgd.db',  ont="BP", keytype="ENSEMBL")

semSimTable <- mgoSim(representativeTerms, representativeTerms, semData=scGO, measure="Wang", combine=NULL)

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

shortest.path

t

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
        semSims[completed] <- semSimTable[t1, t2]
            
        distances[completed] <- shortest.path[c1, c2]
        
        print(sprintf("Completed: %s", completed))
    }
}
# distances <- distances[distances > 0]
# semSims <- semSims[semSims > 0]

plot(distances, semSims, xlab="Distance on Map", ylab="Semantic Similarity")

cor(distances, semSims)


