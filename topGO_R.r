
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("ALL")
biocLite("topGO")
biocLite("hgu95av2.db")
biocLite("Rgraphviz")
biocLite("colMap")

a <- c(1,2,3,4,5,6,7,8,9)
b  <- c(3,2,1,2,2,3,4,6,7,8,9,0,1)
length(intersect(a, b)) / length(union(a, b))

library(clusterProfiler)
x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2",
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1",
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1",
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
eg <- bitr(x, fromType="SYMBOL", toType="GO", OrgDb="org.Hs.eg.db")
head(eg)

#background gene list
setwd('/home/david/Documents/ghsom')
# allGenes <- scan("Y2H_union.txt", character())
allGenes <- scan("Uetz_screen.txt", character())
allGenes <- unique(allGenes)
length(allGenes)

setwd("/home/david/Documents/ghsom/uetz_communities_06")

g <- list()
numCom <- 0
filename <- sprintf("community_%s.txt", numCom)
while (file.exists(filename)) {
    numCom <- numCom + 1
    g[[numCom]] <- scan(filename, character())
    filename <- sprintf("community_%s.txt", numCom)
}
numCom

inter <- function(x, y){
    return(length(intersect(x, y)) / length(union(x, y)))
}

inter(g[[1]], g[[3]])

find_representative_term <- function(terms){
    
    counts <- numeric(length(terms))
    names(counts) <- terms

    for (term in terms) {
        ancestors <- as.list(GOBPANCESTOR[term])
        for (ancestor in ancestors[[term]]) {
            if (ancestor %in% names(counts)) {
                counts[ancestor] <- counts[ancestor] + 1
            }
        }

    }
#     return (sort(names(counts), decreasing=TRUE)[1])
    return (sort(counts / sum(counts), decreasing=TRUE))
}

library(topGO)

filename <- sprintf("uetz-06.rda")

if (file.exists(filename)){
    
    print(sprintf("loading: %s", filename))
    load(filename)
    print("loaded")
    
} else {
    
    cutOff <- 0.05

    geneLists <- vector("list", numCom) 
    GOdataObjects <- vector("list", numCom) 
    resultFishers <- vector("list", numCom) 
    resultFisher.elims <- vector("list", numCom) 
    results <- vector("list", numCom) 
    topResults <- vector("list", numCom) 
    gos <- vector("list", numCom) 
    normalisedRepresentativeTerms <- vector("list", numCom) 
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

        #fishers exact test classic
        resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        resultFishers[[c]] <- resultFisher

        #fishers exact test elimination
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
        normalisedRepresentativeTerms[[c]] <- find_representative_term(gos[[c]])
        representativeTerms[c] <- names(normalisedRepresentativeTerms[[c]][1])
    #     representativeTerms[c] <- gos[[c]][1]

        print(sprintf("community %s complete", c))
    }
    
    print(sprintf("Saving data: %s", filename))
    save(geneLists, GOdataObjects, resultFishers, resultFisher.elims, 
         results, topResults, gos, representativeTerms, file=filename)
}

length(allGenes)

GOdata <- GOdataObjects[[1]]

numGenes(GOdata)

a <- genes(GOdata)
selGenes <- sample(a, 10)
gs <- geneScore(GOdata, whichGenes = selGenes)
gs

length(sigGenes(GOdata))

length(g[[1]])

plot(graph(GOdata))

sel.terms  <- sample(usedGO(GOdata), 10)
num.ann.genes <- countGenesInTerm(GOdata, sel.terms)

num.ann.genes

ann.genes <- genesInTerm(GOdata, sel.terms)

head(ann.genes)

ann.score <- scoresInTerm(GOdata, sel.terms, use.names=TRUE)
head(ann.score)

termStat(GOdata, sel.terms)

test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)

head(sort(score(resultFisher), decreasing=FALSE))

pvalFis <- score(resultFisher)
hist(pvalFis, 50, xlab="p-values")

runTest(GOdata, statistic="fisher", algorithm="classic")

2.3e-05 < 0.05

allRes <- GenTable(GOdata, classic = resultFisher,
                   orderBy = "Significant", topNodes=100)

allRes

names(score(resultFisher)[score(resultFisher) < 0.05])

goID <- allRes[1, "GO.ID"]
print(showGroupDensity(GOdata, goID, ranks = TRUE))

gt <- printGenes(GOdata, whichTerms = goID, chip = "org.Sc.sgd.db", numChar = 40)

?printGenes

scGO <- godata("org.Sc.sgd.db", ont="BP", keytype="ENSEMBL")

weighted_similarity <- function(namedTerms1, namedTerms2) {
    

    s <- mgoSim(names(namedTerms1), 
           names(namedTerms2), semData=scGO, measure="Rel", combine=NULL)

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

shortestPath

cor(t, shortestPath)

mclusterSim(g, semData=scGO, measure="Wang", combine="BMA")

shortestPath <- read.csv("shortest_path.csv", sep=",", header=FALSE)
shortestPath

library(rJava)

source("https://bioconductor.org/biocLite.R")
biocLite("RDAVIDWebService")

library(RDAVIDWebService)
david <- enrichDAVID(gene = gene,
                     idType = "ENSEMBL_GENE_ID",
                     listType = "Gene",
                     annotation = "GOTERM_BP_DIRECT",
                     david.user = "dxm237@cs.bham.ac.uk")

l <- as.character(GOBPANCESTOR["GO:0031023"])
select(GO.db, keys = l, columns = c("TERM", "DEFINITION"))
