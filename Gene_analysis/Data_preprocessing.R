library(RWeka)
library(mice)
library(Hmisc)
library(VIM)
library(ggplot2)

dataset <- read.table("heart-ch.txt", header = TRUE, sep=",", quote="")
class(dataset)

nrow(dataset) 
ncol(dataset)
dim(dataset)
head(dataset, 0)
summary(dataset)

md.pattern(dataset)
mice_plot <- aggr(dataset, col=c("green","red"),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(dataset), cex.axis=.7,
                  gap=3, ylab=c("Missing data","Pattern"))

dataset[,"chol_imputed"] <- with(dataset, impute(chol, mean))
summary(dataset)

#dataset_n <- Normalize(chol_imputed ~., data = dataset)
#dataset_n <- Normalize(~ chol_imputed , data = dataset)
dataset_n <- Normalize(~. , data = dataset)
summary(dataset_n)

dataset$chol_bin <- as.numeric(cut2(dataset$chol_imputed, g=3))
head(dataset)

dataset[,"chol_bin"] <- as.numeric(cut(dataset[,"chol_imputed"], 5))
head(dataset)

ggplot(dataset,aes(x=age,y=chol, color=num)) + geom_point(size = 4)

ggplot(dataset, aes(chest_pain, fill=factor(num))) + geom_bar()
ggplot(dataset, aes(x=age, fill=factor(num))) + geom_bar()



