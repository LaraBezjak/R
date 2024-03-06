library(cluster)
library(ReactomePA)
library(DOSE)

mrnaNorm <- read.table("BRCA_RSEM_genes_normalized.txt", header = F, fill = T, skip = 2)
mrnaIDs <- read.table("BRCA_RSEM_genes_normalized.txt", header = F, fill = T, nrows = 1)
mrnaIDs <- mrnaIDs[, -1][, -1]

samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "-"))[4],0,1212))
sampleType <- as.data.frame(samp)
sampClass <- lapply(samp, function(t) (if (t < 10) return("1") else return("0")))
mrnaClass <- as.data.frame(sampClass)
dim(mrnaNorm)

dim(mrnaClass)

table(unlist(sampClass))

sampClassNum <- lapply(samp, function(t) (if (t < 10) return(1) else return(0)))
mrnaClassNum <- as.data.frame(sampClassNum)

geneNames <- mrnaNorm[1]
dim(geneNames)

mrnaData = t(mrnaNorm[, -1])
rm(samp)
rm(sampClass)
rm(mrnaNorm)
gc()


bssWssFast <- function (X, givenClassArr, numClass=2){
  classVec <- matrix(0, numClass, length(givenClassArr))
  for (k in 1:numClass) {
    temp <- rep(0, length(givenClassArr))
    temp[givenClassArr == (k - 1)] <- 1
    classVec[k, ] <- temp
  }
  classMeanArr <- rep(0, numClass)
  ratio <- rep(0, ncol(X))
  for (j in 1:ncol(X)) {
    overallMean <- sum(X[, j]) / length(X[, j])
    for (k in 1:numClass) {
      classMeanArr[k] <-
      sum(classVec[k, ] * X[, j]) / sum(classVec[k, ])
    }
    classMeanVec <- classMeanArr[givenClassArr + 1]
    bss <- sum((classMeanVec - overallMean)^2)
    wss <- sum((X[, j] - classMeanVec)^2)
    ratio[j] <- bss/wss
  }
  sort(ratio, decreasing = TRUE, index = TRUE)
}

dim(mrnaData)
dim(mrnaClass)
dim(mrnaClassNum)
dim(geneNames)
bss <- bssWssFast(mrnaData, t(mrnaClassNum), 2)
mrnaDataReduced <- mrnaData[,bss$ix[1:100]]
dim(mrnaDataReduced)
trainClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)

#KMeans
set.seed(1)
kmeans.clusters <- kmeans(mrnaDataReduced, 2, nstart = 20)

table(kmeans.clusters$cluster, trainClasses)
clusplot(mrnaDataReduced, kmeans.clusters$cluster, color=TRUE, shade=TRUE,labels=2, lines=0)

table(pam.clusters$clustering, trainClasses)
clusplot(mrnaDataReduced, pam.clusters$clustering, color=TRUE, shade=TRUE,labels=2, lines=0)

#Pathway analysis
genes <- geneNames[bss$ix[1:100],1]
genes <- as.character(genes)
hugoNames <- lapply(genes, function(t) substr(t, 1, regexpr("\\|", t) - 1)) 
entrezNames <- lapply(genes, function(t) substr(t, regexpr("\\|", t) + 1, nchar(t)))
paths = enrichPathway(unlist(entrezNames), pvalueCutoff=1)
head(summary(paths))
enrichMap(paths, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
