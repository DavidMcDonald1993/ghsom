
library(rJava)
source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")
biocLite("RDAVIDWebService")

#libraries
library(GO.db)
library(topGO)
library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(GOSemSim)

# file <- "Uetz_screen"
file <- "yeast_uetz"

p <- 0.2
init <- 1

db <- org.Sc.sgd.db
mapping <- "org.Sc.sgd.db"
ID <- "ENSEMBL"
# db <- org.Hs.eg.db
# mapping <- "org.Hs.eg.db"
# ID <- "ENTREZ"

##load all community gene lists
setwd(sprintf("/home/david/Documents/ghsom/%s_communities_%s_%s", file, p, init))

#background gene list
backgroundFilename <- "all_genes.txt"
allGenes <- scan(backgroundFilename, character())

#communities
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

head(shortest.path)

cutOff <- 0.05

filename <- sprintf("%s-%s-%s.rda", file, p, cutOff)

if (file.exists(filename)){
    
    print(sprintf("loading: %s", filename))
    load(filename)
    print("loaded")
    
} else {
    
    print("creating topGO objects")

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
                      orderBy = "classicFisher")
        results[[c]] <- allRes

        #go terms < cut off 
        gos[[c]] <- score(resultFisher)[score(resultFisher) < cutOff]

        print(sprintf("community %s complete", c))
    }
    
    print(sprintf("Saving data: %s", filename))
    save(geneLists, GOdataObjects, resultFishers, results, gos, file=filename)
    print("saved")
}

##SEMATIC SIMILARITY
#construct gosemsim object
semsimfile <- sprintf("%s-semsimfile.rda", file)
if (file.exists(semsimfile)){
    print(sprintf("loading: %s", semsimfile))
    load(semsimfile)
    print("loaded")
} else {
    print(sprintf("creating %s", semsimfile))
    hsGO <- godata(mapping, ont="BP", keytype=ID)
    save(hsGO, file=semsimfile)
    print(sprintf("saved semsimfile: %s", semsimfile))
}


most_representative_term_ancestor <- function(namedTerms){
    
    counts <- numeric(length(namedTerms))
    names(counts) <- names(namedTerms)

    for (term in names(namedTerms)) {
        ancestors <- as.list(GOBPANCESTOR[term])
        for (ancestor in ancestors[[term]]) {
            if (ancestor %in% names(counts)) {
                counts[ancestor] <- counts[ancestor] + 1
            }
        }

    }
#     return (sort(names(counts), decreasing=TRUE)[1])
    return (names(sort(counts / sum(counts), decreasing=TRUE)[1]))
}

cluster_similarity <- function(c){
    return(mean(mgeneSim(c, semData=hsGO, measure="Wang", verbose=FALSE)))
}

sapply(g, cluster_similarity)

representativeTermsAncestor <- sapply(Filter(length, gos), most_representative_term_ancestor)

representativeTermsAncestor <-representativeTermsAncestor[!is.na(representativeTermsAncestor)] 

representativeTermsAncestor

select(GO.db, keys=representativeTermsAncestor, columns=c("TERM", "DEFINITION"))

simsGOAncestor <- mgoSim(representativeTermsAncestor, representativeTermsAncestor, semData=hsGO, measure="Resnik", combine=NULL)

head(simsGOAncestor)

head(shortest.path)

information_content <- function(term){
    return (goSim(term, term, semData=hsGO, measure="Resnik"))
}

most_representative_term_ic <- function(namedTerms){
    ics <- sapply(names(namedTerms), information_content)
    names(ics) <- names(namedTerms)
    return(names(sort(ics, decreasing=TRUE)[1]))
}

representativeTermsIC <- sapply(Filter(length, gos), most_representative_term_ic)

representativeTermsIC

select(GO.db, keys=representativeTermsIC, columns=c("TERM", "DEFINITION"))

simsGOIC <- mgoSim(representativeTermsIC, representativeTermsIC, semData=hsGO, measure="Resnik", combine=NULL)

head(simsGOIC)

head(shortest.path[lengths(gos) > 0])

simfile <- sprintf("%s-sims.rda", file)
if (file.exists(simfile)){
    print(sprintf("loading: %s", simfile))
    load(simfile)
    print("loaded")
} else {
    sims <- mclusterSim(g, semData=hsGO, measure="Resnik", combine="BMA")
    save(sims, file=simfile)
    print (sprintf("saved sim file: %s", simfile))
}


head(sims)

distances <- numeric(length = (numCom * (numCom - 1)) / 2)
semSims <- numeric(length = (numCom * (numCom - 1)) / 2)

completed <- 0

for (c1 in 1:numCom) {
    
    for (c2 in c1:numCom) {
        
        if (c1 == c2) next   
        
        completed <- completed + 1  
        semSims[completed] <- sims[c1, c2]
            
        distances[completed] <- shortest.path[c1, c2]
        
        print(sprintf("Completed: %s", completed))
    }
}

plot(distances, semSims, xlab="Distance on Map", ylab="Semantic Similarity")

cor(distances, semSims)
