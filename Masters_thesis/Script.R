source("https://bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("edgeR")
biocLite("Rsubread")
biocLite("RColorBrewer")
biocLite("org.Hs.eg.db")

library(Rsubread)
library(limma)
library(edgeR)
library(gplots)
library(RColorBrewer)

# Alignment

# QC z uporabo Rsubread
qs1 <- qualityScores(filename="files1/SLO1_1.fq.gz",nreads=100)
png(file="QC_SLO1_1.png")
boxplot(qs1)
dev.off()
qs2 <- qualityScores(filename="files1/SLO3_1.fq.gz",nreads=100)
png(file="QC_SLO3_1.png")
boxplot(qs2)
dev.off()
qs3 <- qualityScores(filename="files1/SLO4_1.fq.gz",nreads=100)
png(file="QC_SLO4_1.png")
boxplot(qs3)
dev.off()
qs4 <- qualityScores(filename="files1/SLO5_1.fq.gz",nreads=100)
png(file="QC_SLO5_1.png")
boxplot(qs4)
dev.off()
qs5 <- qualityScores(filename="files1/SLO7_1.fq.gz",nreads=100)
png(file="QC_SLO7_1.png")
boxplot(qs5)
dev.off()
qs6 <- qualityScores(filename="files1/SLO8_1.fq.gz",nreads=100)
png(file="QC_SLO8_1.png")
boxplot(qs6)
dev.off()

fastq.files1 <- list.files(path = "./files1", pattern = ".fq.gz$", full.names = TRUE)
fastq.files1
fastq.files2 <- list.files(path = "./files2", pattern = ".fq.gz$", full.names = TRUE)
fastq.files2

# sestava indexa referenčnega genoma, če ga nimamo
buildindex(basename="hg19_idx",reference="hg19.fasta")

# sama dejanska poravnava; pazi na vrstni red files1 in files2!
align(index="hg19_idx", readfile1=fastq.files1, readfile2=fastq.files2)

# ponovno zaženi R zaradi kompatibilnosti decimalnih mest

# FEATURE COUNTS
bam.files <- list.files(path = "./files1", pattern = ".BAM$", full.names = TRUE)
bam.files

# delež mapiranih branj po poravnavi

props1 <- propmapped(files=bam.files[1])
props2 <- propmapped(files=bam.files[2]) 
props3 <- propmapped(files=bam.files[3]) 
props4 <- propmapped(files=bam.files[4]) 
props5 <- propmapped(files=bam.files[5]) 
props6 <- propmapped(files=bam.files[6]) 

# štetje fragmentov

fc <- featureCounts(bam.files, annot.inbuilt="hg19", isPairedEnd=TRUE)

names(fc)
fc$stat
dim(fc$counts)
head(fc$counts)
head(fc$annotation)

# XY izključitev zaradi mešanih spolov v skupinah

m<-data.frame(fc$annotation[,c("Chr","GeneID","Length")],fc$counts,stringsAsFactors=FALSE)
m$chroma <- substr(m$Chr, 4, 4)
m$ifelse <- ifelse(m$chroma == "X" | m$chroma == "Y",0,1)
m$true <- (m$ifelse == 1)
selected <- m$true
my.data <- m[selected,]
my.clean.data <- my.data[,c(2:9)] # pazi na število vzorcev (imeti želimo stolpce Chr, GeneID, Length, in svzorce)
View(my.clean.data)

write.table(my.clean.data,file="raw_counts.txt",quote=FALSE,sep="\t", row.names=FALSE)

# PREPROCESIRANJE ZA D_E

seqdata <- read.delim("raw_counts.txt", stringsAsFactors = FALSE)

countdata <- seqdata[,-(1:2)]

rownames(countdata) <- seqdata[,1]
colnames(countdata)

colnames(countdata) <- substr(colnames(countdata), 39, 42)

sampleinfo <- read.delim("targets.txt")
sampleinfo

table(colnames(countdata)==sampleinfo$sample)

# sprememba štetij v CPM in filtriranje na podlagi CPM

myCPM <- cpm(countdata)

png(file="Threshold_correspondence_SLO1.png")
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3), ylab="", xlab="")
abline(v=0.5, col="blue")
abline(h=10, col="red")
dev.off()

png(file="Threshold_correspondence_SLO3.png")
plot(myCPM[,2],countdata[,2],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5, col="blue")
abline(h=10, col="red")
dev.off()

png(file="Threshold_correspondence_SLO4.png")
plot(myCPM[,3],countdata[,3],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5, col="blue")
abline(h=10, col="red")
dev.off()

png(file="Threshold_correspondence_SLO5.png")
plot(myCPM[,4],countdata[,4],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5, col="blue")
abline(h=10, col="red")
dev.off()

png(file="Threshold_correspondence_SLO7.png")
plot(myCPM[,5],countdata[,5],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5, col="blue")
abline(h=10, col="red")
dev.off()

png(file="Threshold_correspondence_SLO8.png")
plot(myCPM[,6],countdata[,6],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5, col="blue")
abline(h=10, col="red")
dev.off()

thresh <- myCPM > 1 # gledamo glede na zgornje grafe pri counts=10
table(rowSums(thresh))
keep <- rowSums(thresh) >= 2
counts.keep <- countdata[keep,]

dgeObj <- DGEList(counts.keep)
names(dgeObj)

png(file="Barplot_of_library_sizes.png")
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2)
dev.off()

logcounts <- cpm(dgeObj,log=TRUE)

png(file="Boxplot_of_logCPMs.png")
boxplot(logcounts, xlab="", ylab="Log2 CPM",las=2)
abline(h=median(logcounts),col="blue")
dev.off()

png(file="PCA_dim1_dim2.png")
levels(sampleinfo$info)
col.cell <- c("red","black")[sampleinfo$info]
data.frame(sampleinfo$info,col.cell)
plotMDS(dgeObj,col=col.cell, ylim=c(-5,5),xlim=c(-5,5))
legend("topleft",fill=c("red","black"),legend=levels(sampleinfo$info))
title("Status")
dev.off()

var_genes <- apply(logcounts, 1, var)
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

png(file="High_var_genes.heatmap.png")
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("blue","red")[sampleinfo$info]
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes\nacross samples",ColSideColors=col.cell,scale="row")
dev.off()

# Normalizacija glede na velikosti knjižnic

dgeObj1 <- calcNormFactors(dgeObj)
dgeObj1

# ponovimo vajo po normalizaciji :)

logcounts.N <- cpm(dgeObj1,log=TRUE)

png(file="Boxplot_of_logCPMs_NORM.png")
boxplot(logcounts.N, xlab="", ylab="Log2 CPM",las=2)
abline(h=median(logcounts.N),col="blue")
dev.off()


png(file="PCA_dim1_dim2_NORM.png")
levels(sampleinfo$info)
col.cell <- c("blue","red")[sampleinfo$info]
data.frame(sampleinfo$info,col.cell)
plotMDS(dgeObj1,col=col.cell, ylim=c(-2.5,2.5),xlim=c(-2.5,2.5))
legend("topleft",fill=c("blue","red"),legend=levels(sampleinfo$info))
title("Status")
dev.off()

var_genes.N <- apply(logcounts.N, 1, var)
head(var_genes.N)
select_var.N <- names(sort(var_genes.N, decreasing=TRUE))[1:500]
head(select_var.N)
highly_variable_lcpm.N <- logcounts.N[select_var.N,]
dim(highly_variable_lcpm.N)
head(highly_variable_lcpm.N)

png(file="High_var_genes.heatmap_NORM.png")
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("blue","red")[sampleinfo$info]
heatmap.2(highly_variable_lcpm.N,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes\nacross samples",ColSideColors=col.cell,scale="row")
dev.off()

png(file="Composition_bias.png")
par(mfrow=c(1,6))
plotMD(logcounts,column = 1)
abline(h=0,col="grey")
plotMD(logcounts,column = 2) 
abline(h=0,col="grey")
plotMD(logcounts,column = 3)
abline(h=0,col="grey")
plotMD(logcounts,column = 4) 
abline(h=0,col="grey")
plotMD(logcounts,column = 5)
abline(h=0,col="grey")
plotMD(logcounts,column = 6) 
abline(h=0,col="grey")
dev.off()

png(file="Composition_bias_NORM.png")
par(mfrow=c(1,6))
plotMD(dgeObj1,column = 1)
abline(h=0,col="grey")
plotMD(dgeObj1,column = 2) 
abline(h=0,col="grey")
plotMD(dgeObj1,column = 3)
abline(h=0,col="grey")
plotMD(dgeObj1,column = 4) 
abline(h=0,col="grey")
plotMD(dgeObj1,column = 5)
abline(h=0,col="grey")
plotMD(dgeObj1,column = 6) 
abline(h=0,col="grey")
dev.off()

# Experimentalni design

group <- sampleinfo$status
group <- factor(group)
design <- model.matrix(~ group)

# VOOM pristop - limma package

y <- voom(dgeObj1,design,plot=TRUE) 


png(file="PCA_dim1_dim2_VOOM.png")
plotMDS(y,xlim=c(-2.5,2.5))
dev.off()

prefit <- lmFit(y,design)
fit <- eBayes(prefit)
topTable(fit, coef=2)

results.limma <- as.data.frame(topTable(fit,n = Inf))
head(results.limma)

write.csv(x=data.frame(results.limma,stringsAsFactors=FALSE),file="raw_differential_expression_limma.csv",quote=FALSE,row.names=TRUE)

# edgeR

dgeObj2 <- estimateCommonDisp(dgeObj1)
dgeObj2 <- estimateGLMTrendedDisp(dgeObj2)
dgeObj2 <- estimateTagwiseDisp(dgeObj2)

png(file="Estimated_dispersions.png")
plotBCV(dgeObj2)
dev.off()

fit.edgeR <- glmFit(dgeObj2, design)
final <- glmLRT(fit.edgeR, coef=2) 
topTags(final)

results.edgeR <- as.data.frame(topTags(final,n = Inf))
head(results.edgeR)

write.csv(x=data.frame(results.edgeR,stringsAsFactors=FALSE),file="raw_differential_expression_edgeR.csv",quote=FALSE,row.names=TRUE)


# Anotacija
library(org.Hs.eg.db)

columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
keys(org.Hs.eg.db, keytype="ENTREZID")[1:10]

my.keys <- c("7453","338758","220213","54625")

my.keys %in% keys(org.Hs.eg.db, keytype="ENTREZID")
all(my.keys %in% keys(org.Hs.eg.db, keytype="ENTREZID"))

annL <- select(org.Hs.eg.db,keys=rownames(results.limma),columns=c("ENTREZID","SYMBOL","GENENAME"))
annE <- select(org.Hs.eg.db,keys=rownames(results.edgeR),columns=c("ENTREZID","SYMBOL","GENENAME"))

table(annL$ENTREZID==rownames(results.limma))
table(annE$ENTREZID==rownames(results.edgeR))

results.annotated.limma <- cbind(annL, results.limma)
results.annotated.edgeR <- cbind(annE, results.edgeR)

write.csv(results.annotated.limma,file="annotated_differential_expression_limma.csv",row.names=FALSE)
write.csv(results.annotated.edgeR,file="annotated_differential_expression_edgeR.csv",row.names=FALSE)

# Shranjevanje R objectov

save(objekt1, objekt2, results.annotated, dgeObj, y, file="ime.Rdata")
load("ime.Rdata")



#compare
l<-subset(results.annotated.limma, results.annotated.limma$P.Value<=0.05)
e<-subset(results.annotated.edgeR, results.annotated.edgeR$PValue<=0.05)
nrow(l)
nrow(e)

same<-c()
for (i in l$ENTREZID) {
  if (is.element(i, e$ENTREZID)) {
    same<-c(same, i)}
}
nrow(results.annotated.limma)
nrow(results.annotated.edgeR)


l.log<-subset(results.annotated.limma, results.annotated.limma$logFC<=-1.5)
l.log1<-subset(results.annotated.limma, results.annotated.limma$logFC>=1.5)
e.log<-subset(results.annotated.limma, results.annotated.edgeR$logFC<=-1.5)
e.log1<-subset(results.annotated.limma, results.annotated.edgeR$logFC>=1.5)
nrow(l.log)+nrow(l.log1)
nrow(e.log)+nrow(e.log1)

same.log<-c()
for (i in l.log$ENTREZID) {
  if (is.element(i, e.log$ENTREZID)) {
    same.log<-c(same.log, i)}
}
for (i in l.log1$ENTREZID) {
  if (is.element(i, e.log1$ENTREZID)) {
    same.log<-c(same.log, i)}
}


summa.fit <- decideTests(fit)
summary(summa.fit)
summa.final<-decideTestsDGE(final)
summary(summa.final)
