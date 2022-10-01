library(readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)
library("DESeq2")
library(M3C)
library(magrittr)
library(org.Hs.eg.db)
library(topGO)
set.seed(12345)

#Getting GEO expression matrix
cts <- read.delim("GSE68086_TEP_data_matrix.txt", row.names=1)

##Preparing density plot
#logcts <- data.frame(cts)
#rangects <- data.frame(Element_name = row.names(cts), range = 0)
#col = colnames(cts)
#for (i in 1:285) {
#  logcts[col[i]] <- log(cts[,c(i)]+1,2)
#}
#for (i in 1:57736) {
#  rangects[i,2] = max(logcts[i,]) - min(logcts[i,])
#}
##Density plot
#rangects %>% filter(is.finite(range)) %>% ggplot(aes(x=range)) + geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)

#Removing extra columns from sample data matrix
series_matrix <- read.delim("GSE68086_series_matrix.txt")
for (i in 247:286) {
  series_matrix[12,i] = series_matrix[10,i]
}
series_matrix <- data.frame(t(data.frame(as.data.frame(series_matrix[,-1]))))
series_matrix <- series_matrix[,-c(2:11)]      
series_matrix <- series_matrix[,-c(3:36)]

#Changing group name
for (i in 1:285) {                             
  if (series_matrix$X12[i] == "cancer type: HC"){
    series_matrix$X12[i] = "Healthy"
  }
  else{
    series_matrix$X12[i] = "Cancer"
  }
}

#PCA plot
colnames(cts) = rownames(series_matrix)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = series_matrix,
                              design = ~ X12)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
#plotPCA(vsd, intgroup=c("X12"))

##t-SNE Plot
#labels_tsne <- series_matrix[,"X12"]
#tsne(cts, labels = as.factor(labels_tsne), perplex = 10, 
#     axistextsize = 10, legendtextsize = 15, 
#     dotsize = 4, legendtitle = "Sample Groups")

#Differential analysis 
deseq_object <- DESeq(dds)
deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(
  deseq_object, 
  coef = 2, 
  res = deseq_results 
)

#Using tidyverse to clean up table (code from example)
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the genes names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("genes") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by threshold
  dplyr::arrange(dplyr::desc(threshold))

mapped_list <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = deseq_df$genes,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "ENTREZID", # The type of gene identifiers you would like to map to
  multiVals = "list"
)
deseq_df2 <- deseq_df
deseq_df2$genes <-mapped_list[deseq_df$genes]

#Volcano Plot
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df2$genes,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

#topGO needs an expression set
pData <- data.frame(group = series_matrix$X12, 
                    row.names = colnames(cts))

cts2 <- as.matrix(cts)
metadata <- data.frame(labelDescription=
                         c("Sample Group"),
                       row.names=c("group"))

phenoData <- new("AnnotatedDataFrame",
                 data=pData, varMetadata=metadata)
exprSet <- ExpressionSet(assayData=as.matrix(cts2),
                         phenoData=phenoData)

#vector of gene names and p-values
temp <- data.frame(p = deseq_df$padj, 
                   zeros = 0,
                   row.names = deseq_df$genes)
temp <- temp[row.names(cts), ]
genelist <- temp$p
names(genelist) <- row.names(temp)
topDiffGenes <- function(allScore) {
  return(allScore < 0.01)
}
sum(topDiffGenes(genelist))
sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = genelist, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ENSEMBL")
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")