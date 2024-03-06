library(randomForest)
library(class)

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
rm(sampClass)
rm(mrnaNorm)
gc()

#KNN
trainSet <- mrnaData
testSet <- mrnaData
trainClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)
testClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)
knn.predic <- knn(trainSet, testSet, trainClasses, testClasses,k=1)
cbr.predic = as.vector(knn.predic)
table(cbr.predic, testClasses)
tab <- table(cbr.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)

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

trainSet <- mrnaDataReduced
testSet <- mrnaDataReduced
trainClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)

testClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)

knn.predic <- knn(trainSet, testSet, trainClasses, testClasses,k=1)
table(knn.predic, testClasses)
tab <- table(knn.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab)) 
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)

trainSetClass <- as.data.frame(cbind(trainSet, t(mrnaClassNum[1,]))) 
testSetClass <- as.data.frame(cbind(testSet, t(mrnaClassNum[1,])))
colnames(trainSetClass)[101] <- "class"
trainSetClass$class <- as.factor(trainSetClass$class) 
class(trainSetClass$class) 
rf <- randomForest(class ~., trainSetClass, ntree=100, importance=T) 
colnames(testSetClass)[101] <- "class" 
testSetClass$class <- as.factor(testSetClass$class) 
rf.predic <- predict(rf ,testSetClass)
table(rf.predic, testClasses) 
tab <- table(rf.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)

nbRows <- nrow(mrnaDataReduced)
set.seed(22) 
trainRows <- sample(1:nbRows, .70*nbRows)
trainSet <- mrnaDataReduced[trainRows, ]
testSet <- mrnaDataReduced[-trainRows, ]
dim(trainSet)
dim(testSet)

trainClasses <- unlist(mrnaClassNum[1,trainRows], use.names=FALSE)
testClasses <- unlist(mrnaClassNum[1,-trainRows], use.names=FALSE)
knn.predic <- knn(trainSet, testSet, trainClasses, testClasses,k=1)
knn.predic = as.vector(knn.predic)
table(knn.predic, testClasses)
tab <- table(knn.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)


#RandomForest
trainSetClass <- as.data.frame(cbind(trainSet, t(mrnaClassNum[1,trainRows])))
testSetClass <- as.data.frame(cbind(testSet, t(mrnaClassNum[1,-trainRows])))
colnames(trainSetClass)[101] <- "class"
trainSetClass$class <- as.factor(trainSetClass$class)
class(trainSetClass$class)

rf <- randomForest(class ~., trainSetClass, ntree=100, importance=T)
colnames(testSetClass)[101] <- "class"
testSetClass$class <- as.factor(testSetClass$class)
rf.predic <- predict(rf ,testSetClass)
rf.predic = as.vector(rf.predic)
table(rf.predic, testClasses)
tab <- table(rf.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)