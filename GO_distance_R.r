
#libraries
library(GO.db)
library(topGO)
# library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(GOSemSim)

file <- "yeast_uetz"

ont <- "BP"
p <- 0.8
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

##SEMATIC SIMILARITY
#construct gosemsim object

scGO <- godata(mapping, ont=ont, keytype=ID)
print("DONE")

allGeneNames <- scan(character(), file="../yeast_uetz_communities_0.5_1/all_genes.txt")

g  <- sapply(g, function(i) allGeneNames[as.integer(i)])
allGenes <- allGeneNames[as.integer(allGenes)]

enrichedGOTerms <- function(genes, allGenes, cutoff, correction, ont, mapping, ID, algorithm){
    interestingGenes <- factor(as.integer(allGenes %in% genes))
    names(interestingGenes) <- allGenes
    
    GOdata <- new("topGOdata", description=sprintf("topGO object"),
              ontology = ont, allGenes = interestingGenes,
              annotationFun = annFUN.org, mapping = mapping, 
              ID = ID, nodeSize = 10)
    
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
                      cutoff=0.05, correction=FALSE, ont=ont, mapping=mapping, ID=ID, algorithm="elim")

enrichedGOs[3,]

graphs <- sapply(enrichedGOs[3,], function(i) igraph.from.graphNEL(i$dag))

graphs[[1]]

plot(graphs[[1]])

lengths(enrichedGOs)

enrichedGOs[[1]]

lengths(g)

head(shortest.path)

mgeneSim(g[[1]], semData=scGO, measure="Wang")

mgoSim(names(enrichedGOs[[1]]), names(enrichedGOs[[2]]), semData=scGO, measure="Resnik", combine="BMA")

mgoSim(names(enrichedGOs[[1]]), names(enrichedGOs[[32]]), semData=scGO, measure="Resnik", combine="BMA")

clusterSim(g[[1]], g[[2]], semData=scGO, measure="Wang", combine=NULL)

clusterSim <- mclusterSim(g, semData=scGO, measure="Wang", combine="BMA")

head(clusterSim)

head(shortest.path)

pathways <- read.table("../biochemical_pathways.tab", sep="\t")
cols <- c("pathway_name", "enzyme_name", "E.C._reaction_number", "gene_name", "reference")
colnames(pathways) <- cols

toGene <- function(ORFIdentifiers){
    genes <- character()
    for (identifier in ORFIdentifiers){
        gene <- character()
        try(
            gene <- as.character(org.Sc.sgdGENENAME[identifier])
        )
        genes <- c(genes, gene)
    }
    return(genes)
}

toPath <- function(ORFIdentifiers){
    paths <- character()
    for (identifier in ORFIdentifiers){
        path <- character()
        try(
            path <- as.character(org.Sc.sgdPATH[identifier])
        )
        paths <- c(paths, path)
    }
    return(paths)
}

get_pathways <- function(ORFIdentifiers, pathways) {
    genes  <- toGene(ORFIdentifiers)
    return(subset(pathways, gene_name %in% genes)$pathway_name)
}

get_pathway_genes <- function(ORFIdentifiers, pathways) {
    genes  <- toGene(ORFIdentifiers)
    return(subset(pathways, gene_name %in% genes)$gene_name)
}

pathway_list <- sapply(g, get_pathways, pathways)
pathway_genes <- sapply(g, get_pathway_genes, pathways)

enrichedGOsPathway <- sapply(pathway_genes[lengths(pathway_genes) > 0], enrichedGOTerms, allGenes=allGeneNames, 
                      cutOff=cutOff, correction=correction, ont=ont, mapping=mapping, ID=ID)

range <- 1:length(enrichedGOsPathway)

simsPathway <- sapply(range, function(i) sapply(range, function(j) 
                    mgoSim(names(enrichedGOsPathway[[i]]),
                        names(enrichedGOsPathway[[j]]),
                        semData=scGO, measure="Wang", combine="BMA")))

head(simsPathway)

head(shortest.path)

enrichedGOs[[1]]

geneSimilarities <- sapply(allGenes, function(i) sapply(allGenes, function(j) geneSim(i, j, semData=scGO, combine="BMA")))

geneSimilarities

cutOff <- 0.05

filename <- sprintf("%s-%s-%s-%s.rda", file, p, cutOff, ont)

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
                      ontology = ont, allGenes = geneList,
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
        
        #go terms < cut off  Benjamini-Hochberg multiple hypothesis corrected pval
        gos[[c]] <- score(resultFisher)[which(p.adjust(score(resultFisher), method="BH") <= cutOff)]

        print(sprintf("community %s complete", c))
    }
    
    print(sprintf("Saving data: %s", filename))
    save(geneLists, GOdataObjects, resultFishers, results, gos, file=filename)
    print("saved")
}

print_accession_number <- function(terms, file){
    for (s in strsplit(names(terms), ":")){
        write(s[2], file=file, append=TRUE)
    }
}

###write accession number to file
for (i in 1:length(gos)){
    accessionFile <- sprintf("accession_numbers-%s-%s-%s", cutOff, ont, i)
    print_accession_number(gos[[i]], file=accessionFile)
}

wangAllGeneSim <- mgeneSim(allGenes, semData=scGO, measure="Wang", combine="BMA", verbose=TRUE)

clusters <- hclust(as.dist(-log(wangAllGeneSim)))
clusterCut <- cutree(clusters, numCom)

plot(clusters)

assignedCommunities <- numeric(length(allGenes))
names(assignedCommunities) <- allGenes

for (i in 1:numCom){
    for (geneName in g[[i]]){
        assignedCommunities[geneName] <- i
    }
}

library(NMI)

assignedCommunities <- assignedCommunities[names(assignedCommunities) %in% names(clusterCut)]

assignedCommunitiesDF <- data.frame(assignedCommunities)
assignedCommunitiesDF <- cbind(Row.Names = rownames(assignedCommunitiesDF), assignedCommunitiesDF)

clusterCutDF <- data.frame(clusterCut)
clusterCutDF <- cbind(Row.Names = rownames(clusterCutDF), clusterCutDF)

NMI(assignedCommunitiesDF, clusterCutDF)

most_representative_term_weighted <- function(namedTerms){
    
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
#     return (sort(counts / sum(counts), decreasing=TRUE))
    return (sort(counts / max(counts), decreasing=TRUE))
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
#     return (sort(counts / sum(counts), decreasing=TRUE))
    return (names(sort(counts / sum(counts), decreasing=TRUE)[1]))
}

representativeTermsAncestor <- sapply(Filter(length, gos), most_representative_term_ancestor)

select(GO.db, keys=representativeTermsAncestor, columns=c("TERM", "DEFINITION"))

simsGOAncestor <- mgoSim(representativeTermsAncestor, representativeTermsAncestor, semData=scGO, measure="Wang", combine=NULL)

head(simsGOAncestor)

head(shortest.path)

information_content <- function(term){
    return (goSim(term, term, semData=scGO, measure="Resnik"))
}

most_representative_term_ic <- function(namedTerms){
    ics <- sapply(names(namedTerms), information_content)
    names(ics) <- names(namedTerms)
    return(names(sort(ics, decreasing=TRUE)[1]))
}

representativeTermsIC <- sapply(Filter(length, gos), most_representative_term_ic)

select(GO.db, keys=representativeTermsIC, columns=c("TERM", "DEFINITION"))

simsGOIC <- mgoSim(representativeTermsIC, representativeTermsIC, semData=scGO, measure="Wang", combine=NULL)

head(simsGOIC)

head(shortest.path)

wangClusterSim <- mclusterSim(g, semData=scGO, measure="Wang", combine="BMA")

head(wangClusterSim)

head(shortest.path)

goSims <- matrix(numeric(), nrow=numCom, ncol=numCom)

for (i in 1:numCom){
    for (j in 1:numCom){
        goSims[i, j] = mgoSim(names(gos[[i]]), names(gos[[j]]), measure="Wang", semData=scGO, combine="BMA")
    }
}

head(goSims)

wangGoSims <- sapply(names(enrichedGOs), 
                     function(i) sapply(names(enrichedGOs), 
                                        function(j) mgoSim(i, j, semData=scGO, measure="Wang", combine="BMA")))

wangGoSims

mgeneSim(allGeneNames[as.integer(g[[1]])], semData=scGO, measure="Wang", combine="BMA")

mgoSim(names(enrichedGOs[[1]]), names(enrichedGOs[[2]]), semData=scGO, measure="Wang", combine="BMA")

head(shortest.path)

distances <- numeric(length = (numCom * (numCom - 1)) / 2)
semSims <- numeric(length = (numCom * (numCom - 1)) / 2)

completed <- 0

for (c1 in 1:length(enrichedGOsPathway)) {
    
    for (c2 in c1:length(enrichedGOsPathway)) {
        
        if (c1 == c2) next   
        
        completed <- completed + 1  
        semSims[completed] <- simsPathway[c1, c2]
            
        distances[completed] <- shortest.path[c1, c2]
        
        print(sprintf("Completed: %s", completed))
    }
}

plot(distances, semSims, xlab="Distance on Map", ylab="Shared Paths")

cor(distances, semSims, method="spearman")

library(GOSim)
setOntology(ont, loadIC=FALSE)
setEvidenceLevel(evidences="all",organism=org.Sc.sgdORGANISM, gomap=org.Sc.sgdGO)
e <- GOenrichment(g[[46]], allGenes)

e

goTerms <- e$GOTerms
p.values <- e$p.values

p.values.df <- data.frame(p.values)
p.values.df["go_id"] <- names(p.values)
p.values.df

goTerms <- merge(goTerms, p.values.df, by="go_id")

colnames(goTerms) <- c("GO_ID", "TERM", "DEFINITION", "P_VALUE")
head(goTerms)

library(gridExtra)
grid.table(goTerms[,c("GO_ID", "TERM", "P_VALUE")])

g[[46]]

l <- as.list(org.Sc.sgdGO[["YKL019W"]])

gos <- sapply(l, function(i) i[["GOID"]])

t <- select(GO.db, keys=gos, columns=c("GOID","TERM","ONTOLOGY"))

grid.table(t)


