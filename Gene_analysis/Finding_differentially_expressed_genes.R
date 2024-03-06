mrnaNorm <- read.table("BRCA_RSEM_genes_normalized.txt", header = F, fill = T, skip = 2)
mrnaIDs <- read.table("BRCA_RSEM_genes_normalized.txt", header = F, fill = T, nrows = 1)

dim(mrnaIDs)
mrnaIDs <- mrnaIDs[, -1][, -1]

dim(mrnaNorm)
dim(mrnaIDs)
mrnaIDs5 = mrnaIDs[,c(1,2,3,4,5)]
mrnaNorm5x5 = mrnaNorm[1:5, 1:5]
head(mrnaIDs5, 2)
head(mrnaNorm, 2)
summary(mrnaNorm5x5)

samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "-"))[4],0,1212))
sampleType <- as.data.frame(samp)
dim(sampleType)

tab <- table(unlist(sampleType))
tab

sampClass <- lapply(samp, function(t) (if (t < 10) return("1") else return("0")))
mrnaClass <- as.data.frame(sampClass)
dim(mrnaClass)
table(unlist(sampClass))

sampClassNum <- lapply(samp, function(t) (if (t < 10) return(1) else return(0)))
mrnaClassNum <- as.data.frame(sampClassNum)

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

bss <- bssWssFast(t(mrnaNorm[, -1]), t(mrnaClassNum), 2)
bss$ix[1:50]

genes <- mrnaNorm[bss$ix[1:50],1]
genes

mrnaSetClass <- rbind(mrnaNorm[bss$ix[1:100],-1], setNames(mrnaClassNum, names(mrnaNorm[, -1])))
dim(mrnaSetClass)

transpose <- as.data.frame(t(mrnaSetClass))
colnames(transpose)[101] <- "class"

normals <- subset(transpose, transpose$class == 0)
dim(normals)

tumors <- subset(transpose, transpose$class == 1)
dim(tumors)

both <- rbind(normals, tumors[1:112,])
dim(both)


color.map <- function(class) { if (class==0) "#FF0000" else "#0000FF" }
groupColors <- unlist(lapply(both$class, color.map))
heatmap(as.matrix(t(both)), scale = "row", col=topo.colors(100), ColSideColors=groupColors)