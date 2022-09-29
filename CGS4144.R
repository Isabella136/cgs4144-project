library(readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)
library("DESeq2")
library(M3C)

#Getting GEO expression matrix
cts <- read_excel("GSE65683_Required_SRE.xlsx")
for (i in 1:3){
  cts <- cts[,-c(1)]
}
cts <- cts[,-c(2)]

#Preparing density plot
logcts <- data.frame(Element_name = cts$`Element name`)
rangects <- data.frame(Element_name = cts$`Element name`,
                        range = 0)
col = colnames(cts)
for (i in 2:73) {
  logcts[col[i]] <- log(cts[,c(i)],2)
}
for (i in 1:528) {
  rangects[i,2] = max(logcts[,-c(1)][i,]) - min(logcts[,-c(1)][i,])
}
#Density plot
rangects %>% filter(is.finite(range)) %>% ggplot(aes(x=range)) + geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)

#To take care of SREs with no expression
logctsadjusted <- data.frame(Element_name = cts$`Element name`)
rangectsadjusted <- data.frame(Element_name = cts$`Element name`,
                                range = 0)
for (i in 2:73) {
  logctsadjusted[col[i]] <- log(cts[,c(i)]+1,2)
}
for (i in 1:528) {
  rangectsadjusted[i,2] = max(logctsadjusted[,-c(1)][i,]) - min(logctsadjusted[,-c(1)][i,])
}
rangectsadjusted %>% filter(is.finite(range)) %>% ggplot(aes(x=range)) + geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)

#Removing extra columns from sample data matrix
series_matrix <- read.delim("GSE65683_series_matrix2.txt")
series_matrix <- data.frame(t(data.frame(as.data.frame(series_matrix[,-1]))))
series_matrix <- series_matrix[,-c(2:8)]      
series_matrix <- series_matrix[,-c(3:52)]

#Changing group name
for (i in 1:72) {                             
  if (series_matrix$X9[i] == "phenotype group (birth outcome): Group I" || series_matrix$X9[i] == "phenotype group (birth outcome): Group II-i"){
    series_matrix$X9[i] = "Timed Intercourse or Intrauterine Insemination"
  }
  else if (series_matrix$X9[i] == "phenotype group (birth outcome): Group III-i" || series_matrix$X9[i] == "phenotype group (birth outcome): Group III-ii"){
    series_matrix$X9[i] = "Clinic Samples"
  }
  else{
    series_matrix$X9[i] = "Assisted Reproductive Technology"
  }
}

#Setting up for PCA plot
elements = cts[,1]$`Element name`
cts <- data.frame(as.data.frame(cts[,-1]))
j = 1
for (i in 1:528){
  if (elements[i] == "IntronicElmt") {
    elements[i] = paste("IntronicElmt", j, sep="_")
    j = j+1
  }
}
row.names(cts) <- as.character(elements)
cts <- cts[, rownames(series_matrix)]
dds <- DESeqDataSetFromMatrix(countData = cts*1000000,
                              colData = series_matrix,
                              design = ~ X9)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

#PCA plot
pcaData <- plotPCA(vsd, intgroup=c("X9"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pointSize<- 1:72
alpha <- 1:72
for (i in 1:72) {
  if (series_matrix$X9[i] == "Clinic Samples") {
    pointSize[i] = 3
    alpha[i] = 0.2
  }
  else if (series_matrix$X9[i] == "Timed Intercourse or Intrauterine Insemination") {
    pointSize[i] = 5
    alpha[i] = 0.8
  }
  else {
    pointSize[i] = 3
    alpha[i] = 0.8
  }
}
ggplot(pcaData, aes(PC1, PC2, color=X9)) +
  geom_point(size=pointSize, alpha=alpha) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  labs(col='Sample Groups')

#PCA plot w/o clinical samples
i = 1
cts2 = cts
series_matrix2 = series_matrix
while (i < (length(cts2)+1)) {
  if (series_matrix2$X9[i] == "Clinic Samples"){
    cts2 = cts2[,-c(i)]
    series_matrix2 = series_matrix2[-c(i),]
  }
  else {
    i = i+1
  }
}
dds2 <- DESeqDataSetFromMatrix(countData = cts2*1000000,
                              colData = series_matrix2,
                              design = ~ X9)
vsd2 <- varianceStabilizingTransformation(dds2, blind=FALSE)
pcaData2 <- plotPCA(vsd2, intgroup=c("X9"), returnData=TRUE)
percentVar2 <- round(100 * attr(pcaData2, "percentVar"))
pointSize2<- 1:63
alpha2 <- 1:63
for (i in 1:63) {
  if (series_matrix$X9[i] == "Timed Intercourse or Intrauterine Insemination") {
    pointSize2[i] = 5
    alpha2[i] = 0.8
  }
  else {
    pointSize2[i] = 3
    alpha2[i] = 0.8
  }
}
ggplot(pcaData2, aes(PC1, PC2, color=X9)) +
  geom_point(size=pointSize2, alpha=alpha2) +
  xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar2[2],"% variance")) + 
  coord_fixed() + 
  labs(col='Sample Groups')

#t-SNE Plot
labels_tsne <- series_matrix2[,"X9"]
tsne(cts2, labels = as.factor(labels_tsne), perplex = 10, 
     axistextsize = 10, legendtextsize = 15, 
     dotsize = 4, legendtitle = "Sample Groups")

