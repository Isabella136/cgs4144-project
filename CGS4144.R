library(readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)
library("DESeq2")

cts <- read_excel("GSE65683_Required_SRE.xlsx")
for (i in 1:3){
  cts <- cts[,-c(1)]
}
cts <- cts[,-c(2)]
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

#PCA plot
series_matrix <- read.delim("GSE65683_series_matrix2.txt")
series_matrix <- data.frame(t(data.frame(as.data.frame(series_matrix[,-1]))))
series_matrix <- series_matrix[,-c(2:8)]
series_matrix <- series_matrix[,-c(3:52)]

for (i in 1:72) {
  if (series_matrix$X9[i] == "phenotype group (birth outcome): Group I"){
    series_matrix$X9[i] = "control"
  }
  else{
    series_matrix$X9[i] = "experiment"
  }
}

cts2 <- data.frame(as.data.frame(cts[,-1]))
elements = cts[,1]$`Element name`
j = 1
for (i in 1:528){
  if (elements[i] == "IntronicElmt") {
    elements[i] = paste("IntronicElmt", j, sep="_")
    j = j+1
  }
  
}
row.names(cts2) <- as.character(elements)
cts2 <- cts2[, rownames(series_matrix)]

dds <- DESeqDataSetFromMatrix(countData = cts2*1000000,
                              colData = series_matrix,
                              design = ~ X9)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("X9"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pointSize<- 1:72
alpha <- 1:72

for (i in 1:72) {
  if (series_matrix$X9[i] == "control") {
    pointSize[i] = 5
    alpha[i] = 1
  }
  else {
    pointSize[i] = 3
    alpha[i] = 0.5
  }
}
ggplot(pcaData, aes(PC1, PC2, color=X9)) +
  geom_point(size=pointSize, alpha=alpha) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  labs(col='Sample Group')
