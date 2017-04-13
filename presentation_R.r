
#libraries
library(GO.db)
library(topGO)
library(org.Sc.sgd.db)
library(GOSemSim)
library(gridExtra)

file <- "yeast_uetz"

ont <- "BP"
init <- 1

db <- org.Sc.sgd.db
mapping <- "org.Sc.sgd.db"
ID <- "ENSEMBL"

scGO <- godata(mapping, ont=ont, keytype=ID)

ps <- seq(from = 0.1, to = 0.9, by = 0.1)

ps

getMeanSimilarity <- function(p) {
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
    
    clusterSims <- sapply(g, function(i) 
        mean(mgeneSim(genes = i, semData = scGO, measure = "Wang", combine = "BMA", verbose = FALSE)))
    
    return(list(numCom, mean(clusterSims)))
}

results <- sapply(ps, getMeanSimilarity)

results

results[1,]

results[2,]

plot(ps,results[1,],type="l",col="red")
# par(new=TRUE)
plot(ps,results[2,],col="green", add=TRUE)

results[1,]

results[2,]

par(mar = c(5,5,2,5))
plot(ps, results[1,], type="l", col="red3", xlab=expression(e[sg]), xlim=c(0.9,0.1),
     ylab = "Number of Communities Detected")

par(new = TRUE)
plot(ps, results[2,], pch=16, xlim=c(0.9,0.1), axes=FALSE, xlab=NA, ylab=NA, cex=1.2)
axis(side = 4)
mtext(side = 4, line = 3, 'Mean Semantic Similarity')
legend("topleft",
       legend=c("Num Communities", "Semantic Similarity"),
       lty=c(1,0), pch=c(NA, 16), col=c("red3", "black"))


