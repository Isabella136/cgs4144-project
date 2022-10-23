library(readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)
library("DESeq2")
library(M3C)
library(magrittr)
library(org.Hs.eg.db)
library(topGO)
library(ComplexHeatmap)
library(clusterProfiler)
library(gprofiler2)
library(enrichplot)
library(mclust)
library(cluster)
library(ConsensusClusterPlus)
library(ClusterR)

set.seed(12345)

#Getting GEO expression matrix
cts <- read.delim("GSE68086_TEP_data_matrix.txt", row.names=1)
series_matrix <- read.delim("GSE68086_series_matrix.txt")

density_plot <- function() {
  logcts <- data.frame(cts)
  rangects <- data.frame(Element_name = row.names(cts), range = 0)
  col <- colnames(cts)
  for (i in 1:285)
    logcts[col[i]] <- log(cts[,c(i)]+1,2)
  for (i in 1:57736)
    rangects[i,2] = max(logcts[i,]) - min(logcts[i,])
  rangects %>% filter(is.finite(range)) %>% ggplot(aes(x=range)) + geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)
}

series_matrix_set_up <- function() {
  for (i in 247:286)
    series_matrix[12,i] = series_matrix[10,i]
  series_matrix <- data.frame(t(data.frame(as.data.frame(series_matrix[,-1]))))
  series_matrix <- series_matrix[,-c(2:11)]      
  series_matrix <- series_matrix[,-c(3:36)]
  #Changing group name
  for (i in 1:285) {                             
    if (series_matrix$X12[i] == "cancer type: HC")
      series_matrix$X12[i] = "Healthy"
    else
      series_matrix$X12[i] = "Cancer"
  }
  return(series_matrix)
}

deseqdata <- function() {
  colnames(cts) = rownames(series_matrix)
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = series_matrix,
                                design = ~ X12)
  return(dds)
}

PCA_plot <- function() {
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  plotPCA(vsd, intgroup=c("X12"))
}

tSNE_plot <- function() {
  labels_tsne <- series_matrix[,"X12"]
  tsne(cts, labels = as.factor(labels_tsne), perplex = 10, 
       axistextsize = 10, legendtextsize = 15, 
       dotsize = 4, legendtitle = "Sample Groups")
}

differential_analysis <- function() {
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
    dplyr::arrange(padj)
  return(deseq_df)
}

ensembl_to_hugo <- function(){
  mapped_list <- mapIds(
    org.Hs.eg.db, # Replace with annotation package for your organism
    keys = deseq_df$genes,
    keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
    column = "SYMBOL", # The type of gene identifiers you would like to map to
    multiVals = "list"
  )
  deseq_df2 <- deseq_df
  deseq_df2$genes <-mapped_list[deseq_df$genes]
  return(deseq_df2)
}

volcano_plot <- function() {
  volcano_plot <- EnhancedVolcano::EnhancedVolcano(
    deseq_df,
    lab = deseq_df2$genes,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
  )
}

get_significant_genes <- function() {
  significant_genes <- c()
  for (i in 1:nrow(deseq_df)){
    if(is.na(deseq_df[i, "threshold"]))
      next
    if(deseq_df[i, "threshold"] == TRUE)
      significant_genes <- append(significant_genes, deseq_df[i, "genes"])
  }
  return(significant_genes)
}

gene_list_numbered_vector <- function() {
#vector of gene names and p-values
  genelist <- c()
  j = 1
  for (i in 1:nrow(deseq_df)){
    if (j > length(significant_genes))
      break
    if (deseq_df[i, "genes"] == significant_genes[j]){
      genelist <- append(genelist, deseq_df[i, "padj"])
      j = j+1
    }
  }
  names(genelist) <- significant_genes
  return(genelist)
}

topDiffGenes <- function(allScore) {
  return(allScore < 0.01)
}

topGo <- function() {
  sum(topDiffGenes(genelist))
  sampleGOdata <- new("topGOdata",
                      description = "Simple session", ontology = "BP",
                      allGenes = genelist, geneSel = topDiffGenes,
                      nodeSize = 10,
                      annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ENSEMBL")
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
  resultFisher.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
  allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim,
                     orderBy = "elimFisher", ranksOf = "classicFisher", topNodes = 10)
}

clustProfiler <- function () {
  gene_list <- genelist
  gene_list <- na.omit(gene_list)
  gene_list = sort(gene_list, decreasing=TRUE)
  
  gse <- gseGO(geneList=gene_list, 
               ont ="BP", 
               keyType = "ENSEMBL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "none")
  require(DOSE)
  dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
}


gProfiler2_MF <- function() {
  gProfiler_data_MF <- gost(query = significant_genes, 
                         organism = "hsapiens", ordered_query = FALSE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "g_SCS", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = "GO:MF", as_short_link = FALSE)
  head(gProfiler_data_MF$result, 3)
  gostplot(gProfiler_data_MF, capped = TRUE, interactive = FALSE)
  publish_gosttable(gProfiler_data_MF, highlight_terms = gProfiler_data_MF$result[c(1:10),],
                    use_colors = TRUE, 
                    show_columns = c("source", "term_name", "term_size", "intersection_size"))
}

heatmap <- function() {
  filtered <- deseq_df[deseq_df[2] > 100,]
  filtered <- filtered[abs(filtered[3]) > 2,]
  filtered <- data.frame(filtered)
  filtered <- na.omit(filtered)
  rownames(filtered) <- filtered[,1] #makes the literal names of rows the genes
  filtered <- filtered[,-1] #remove first gene column 
  filtered <- filtered[1:(length(filtered)-1)] #remove true/false
  #create matrix so you can create map
  mat <- counts(dds)[rownames(filtered), ]
  mat <- t(apply(mat, 1, scale))
  coldata <- metadata %>% tibble::column_to_rownames("genes")
  
  Heatmap(mat)
  
  #unused:
  #filtered$genes <- mapIds(org.Hs.eg.db, keys = filtered$genes, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "list")
  #colnames(mat) <- rownames(coldata)
  #map -> Heatmap(mat, cluster_rows = T, cluster_columns = F,
  #               columns_labels = colNames(mat), name = "Heat Bar")
  #Heatmap(mat, cluster_rows=T,cluster_columns=F)
}

gProfiler2_HP <- function() {
  gProfiler_data_HP <- gost(query = significant_genes, 
                         organism = "hsapiens", ordered_query = FALSE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "g_SCS", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = "HP", as_short_link = FALSE)
  head(gProfiler_data_HP$result, 3)
  gostplot(gProfiler_data_HP, capped = TRUE, interactive = FALSE)
  publish_gosttable(gProfiler_data_HP, highlight_terms = gProfiler_data_HP$result[c(1:10),],
                    use_colors = TRUE, 
                    show_columns = c("source", "term_name", "term_size", "intersection_size"))
}

clustering <- function(genesAmt) {
  cluster_data_input <- data.frame(cts[deseq_df$genes[1:genesAmt],], row.names = deseq_df$genes[1:genesAmt])
  gaussian <- Mclust(cluster_data_input)
  pam <- pam(cluster_data_input, 6)
  ccplus <- ConsensusClusterPlus(data.matrix(cluster_data_input))
  kmeans <- kmeans(cluster_data_input, centers=6)
  return(list(gaussian, pam, ccplus, kmeans))
}

series_matrix <- series_matrix_set_up()
dds <- deseqdata()
deseq_df <- differential_analysis()
cluster5000 <- clustering(5000)
cluster10 <- clustering(10)
cluster100 <- clustering(100)
cluster1000 <- clustering(1000)
cluster10000 <- clustering(10000)