
library(ReactomePA)
library(org.Sc.sgd.db)
library(clusterProfiler)

reactomeEdgeList <- scan(character(), file="reactome_edgelist.txt")

reactomeEdgeList <- unique(reactomeEdgeList)

head(reactomeEdgeList)

conversionTable <- select(org.Sc.sgd.db, reactomeEdgeList, "ENTREZID", "UNIPROT")

conversionTable  <- conversionTable[!is.na(conversionTable[,"ENTREZID"]),]

head(conversionTable, n=20)

allGenesOfInterest <- conversionTable[,"ENTREZID"]

universe <- keys(org.Sc.sgd.db)
allGenesConversion <- select(org.Sc.sgd.db, keys(org.Sc.sgd.db), columns="ENTREZID", keytype="ORF")

allGenesConversion <- allGenesConversion[!is.na(allGenesConversion[,"ENTREZID"]),]

universe <- allGenesConversion[,"ENTREZID"]
head(universe)

geneList <- factor(as.integer(universe %in% allGenesOfInterest))

head(geneList, n=100)

x <- enrichPathway(gene=conversionTable[,"ENTREZID"], universe = universe,
                   pvalueCutoff=0.05, organism = "yeast")

df <- as.data.frame(x)

head(df)

df[df[,"Description"]=="p53-Independent DNA Damage Response",]

barplot(x, showCategory=10)

dotplot(x, showCategory=15)

enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)

cnetplot(x, categorySize="pvalue", foldChange=geneList)

require(clusterProfiler)
data(gcSample)
res <- compareCluster(gcSample, fun="enrichPathway")

gcSample

png(filename="cluster_pathway_enrichment.png", width = 1500)
plot(res)
dev.off()

tail(geneList, n=100)

data(geneList)
y <- gsePathway(geneList, nPerm=1000,
                minGSSize=120, pvalueCutoff=0.2,
                pAdjustMethod="BH", verbose=FALSE)
result <- as.data.frame(y)
head(result)

viewPathway(pathName = "p53-Independent DNA Damage Response", organism = "yeast", readable = F)
