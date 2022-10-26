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
library('factoextra')
library(GGally)
library('ggalluvial')
library(viridis)

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
  pam <- pam(dist(t(cluster_data_input)), 7, diss=TRUE)
  ccplus <- ConsensusClusterPlus(data.matrix(cluster_data_input), maxK = 3)
  kmeans <- kmeans(t(cluster_data_input), centers=5)
  return(list(pam, ccplus, kmeans))
}

#done seperately due to taking lots of time by itself
gaussianClusterin <- function(genesAmt) {
  cluster_data_input <- data.frame(cts[deseq_df$genes[1:genesAmt],], row.names = deseq_df$genes[1:genesAmt])
  gaussian <- Mclust(t(cluster_data_input))
}

series_matrix <- series_matrix_set_up()
dds <- deseqdata()
deseq_df <- differential_analysis()
cluster_data_input <- data.frame(cts[deseq_df$genes[1:5000],], row.names = deseq_df$genes[1:5000])

gaussian5000 <- gaussianClusterin(5000)
#plot for gaussian
randProj(t(cluster_data_input), parameters = gaussian5000$parameters, addEllipses= TRUE, main = "Gaussian Cluster", classification = gaussian5000$classification)
gaussian10 <- gaussianClusterin(10)
gaussian100 <- gaussianClusterin(100)
gaussian1000 <- gaussianClusterin(1000)
#gaussian with 10k genes is too expensive
#gaussian10000 <- gaussianClusterin(10000)

cluster5000 <- clustering(5000)
#clusplot for pam
plot(cluster5000[[1]], which.plots = 1, main = "PAM Cluster")
#plot for kmeans clusters
fviz_cluster(cluster5000[[3]], data = t(cluster_data_input), geom = c("point"), axes = c(1,2), main = "Kmeans Cluster, 1v2")
fviz_cluster(cluster5000[[3]], data = t(cluster_data_input), geom = c("point"), axes = c(2,3), main = "Kmeans Cluster, 2v3")
#plot for ccplus
plot(as.dendrogram(cluster5000[[2]][[3]]$consensusTree), ylab = "Height", main = "ConsensusClusterPlus Dendrogram")
rect.hclust(cluster5000[[2]][[3]]$consensusTree, k = 3, border = 2:4)

cluster10 <- clustering(10)
cluster100 <- clustering(100)
cluster1000 <- clustering(1000)
cluster10000 <- clustering(10000)

#pam alluvial
pam_alluvial <- data.frame(pam10 = cluster10[[1]]$clustering, 
                           pam100 = cluster100[[1]]$clustering,
                           pam1000 = cluster1000[[1]]$clustering,
                           pam5000 = cluster5000[[1]]$clustering,
                           pam10000 = cluster10000[[1]]$clustering)
ggplot(pam_alluvial,
    aes(axis1 = pam10 , axis2 = pam100, axis3 = pam1000, axis4 = pam10000)) +
  geom_alluvium(aes(fill = pam5000), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  ggtitle("PAM alluvial") +
  xlab('Gene Count') +
  ylab('Sample Count')

#ccplus alluvial
ccplus_alluvial <- data.frame(cc10 = cluster10[[2]][[3]]$consensusClass, 
                           cc100 = cluster100[[2]][[3]]$consensusClass,
                           cc1000 = cluster1000[[2]][[3]]$consensusClass,
                           cc5000 = cluster5000[[2]][[3]]$consensusClass,
                           cc10000 = cluster10000[[2]][[3]]$consensusClass)
ggplot(ccplus_alluvial,
       aes(axis1 = cc10 , axis2 = cc100, axis3 = cc1000, axis4 = cc10000)) +
  geom_alluvium(aes(fill = cc5000), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  ggtitle("ConsensusClusterPlus alluvial") +
  xlab('Gene Count') +
  ylab('Sample Count')

#kmeans alluvial
kmeans_alluvial <- data.frame(k10 = cluster10[[3]]$cluster, 
                              k100 = cluster100[[3]]$cluster,
                              k1000 = cluster1000[[3]]$cluster,
                              k5000 = cluster5000[[3]]$cluster,
                              k10000 = cluster10000[[3]]$cluster)
ggplot(kmeans_alluvial,
       aes(axis1 = k10 , axis2 = k100, axis3 = k1000, axis4 = k10000)) +
  geom_alluvium(aes(fill = k5000), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  ggtitle("Kmeans alluvial") +
  xlab('Gene Count') +
  ylab('Sample Count')

#gaussian alluvial
gaussian_alluvial <- data.frame(g10 = gaussian10$classification, 
                           g100 = gaussian100$classification,
                           g1000 = gaussian1000$classification,
                           g5000 = gaussian5000$classification)
ggplot(gaussian_alluvial,
       aes(axis1 = g10 , axis2 = g100, axis3 = g1000)) +
  geom_alluvium(aes(fill = g5000), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  ggtitle("Gaussian alluvial") +
  xlab('Gene Count') +
  ylab('Sample Count')



cts_normalized <- cts
for (i in 1:285) {
  cts_normalized[i] <- (cts[i]/sum(cts[i]))*1000000
}
cluster_data_input_normalized <- data.frame(cts_normalized[deseq_df$genes[1:5000],], row.names = deseq_df$genes[1:5000])

#chi-squared test of independence
healthy_cancer <- series_matrix$X12
#pam
cluster_pam_data <- cluster10000[[1]][["clustering"]]
pam_df <- data.frame("Sample" = healthy_cancer, "Number of Clusters" = cluster_pam_data)
table(pam_df)
test_pam <- chisq.test(table(pam_df))

#ccplus
#cluster_ccplus_data <-
#ccplus_df <- data.frame("Sample" = healthy_cancer, "Number of Clusters" = cluster_ccplus_data)
#table(ccplus_df)
#test_ccplus <- chisq.test(table(ccplus_df))

#kmeans
cluster_kmeans_data <- cluster10000[[3]][["cluster"]]
kmeans_df <- data.frame("Sample" = healthy_cancer, "Number of Clusters" = cluster_kmeans_data)
table(kmeans_df)
test_kmeans <- chisq.test(table(kmeans_df))

#gaussian
cluster_gaussian_data <- gaussian5000[["classification"]]
gaussian_df <- data.frame("Sample" = healthy_cancer, "Number of Clusters" = cluster_gaussian_data)
table(gaussian_df)
test_gaussian <- chisq.test(table(gaussian_df))

p_values <- c(test_pam$p.value, #test_ccplus$p.value
              test_kmeans$p.value, test_gaussian$p.value)
p.adjust <- p.adjust(p_values)

statistical_test_results <- data.frame("Cluster Method" = c("PAM", #"ConsensusClusterPlus", 
                                                            "K-Means", "Gaussian"), 
                                       "Chi-Squared" = c(test_pam$statistic, #test_ccplus$statistic,
                                                         test_kmeans$statistic, test_gaussian$statistic),
                                       "P-Values" = p_values, "Adjusted P-Values" = p.adjust)
ggpairs(statistical_test_results)

#Heatmaps Setup
cts_normalized <- cts
for (i in 1:285) {
  cts_normalized[i] <- (cts[i]/sum(cts[i]))*1000000
}

for (i in 1:285) {
  cts_normalized[i] <- cts_normalized[i]+1
}
cts_normalized <- log2(cts_normalized)
cluster_data_input_normalized <- data.frame(cts_normalized[deseq_df$genes[1:5000],], row.names = deseq_df$genes[1:5000])


#Gaussian
gaus_annot <- as.data.frame(gaussian5000$classification)
gaus_annot[,2] <- series_matrix$X12
colnames(gaus_annot) <-c('Cluster Groups', 'Treatment Groups')
colnames(gaus_annot) <-c('Cluster Groups', 'Treatment Groups')
colors <- list('Cluster Groups' =c('1'='red', '2'='blue', '3'='green', '4'='orange', 
                                   '5'='orange', '6'='yellow', '7'='gold', '8'='violet'), 
               'Treatment Groups' =c('Cancer'='pink', 'Healthy'='gray'))
gaus_row_annot <- HeatmapAnnotation(df=gaus_annot, which='row',
                                    col=colors)
m.gaus <- cbind(t(cluster_data_input_normalized), gaussian5000$classification)
Heatmap(m.gaus[,1:5000], name="Scores", left_annotation=gaus_row_annot, show_row_names = F, show_column_names = F, row_title="Genes", column_title = "Samples")

##PAM
pam <- pam(dist(t(cluster_data_input)), 7, diss=TRUE)
pam_annot <- as.data.frame(pam$clustering)
pam_annot[,2] <- series_matrix$X12
colnames(pam_annot) <-c('Cluster Groups', 'Treatment Groups')
colors <- list('Cluster Groups' =c('1'='red', '2'='blue', '3'='green', '4'='orange', 
                                   '5'='orange', '6'='yellow', '7'='gold', '8'='violet'), 
               'Treatment Groups' =c('Cancer'='pink', 'Healthy'='gray'))
pam_row_annot <- HeatmapAnnotation(df=pam_annot, which='row',
                                   col=colors)
m.pam <- cbind(t(cluster_data_input_normalized), pam$clustering)
Heatmap(m.pam[,1:5000], name="Scores", left_annotation=pam_row_annot, show_row_names = F, show_column_names = F, row_title="Genes", column_title = "Samples")

#CCPlus
ccplus <- ConsensusClusterPlus(data.matrix(cluster_data_input), maxK = 3)
ccplus_annot <- as.data.frame(cluster5000[[2]][[3]]$consensusClass)
ccplus_annot[,2] <- series_matrix$X12
colnames(ccplus_annot) <-c('Cluster Groups', 'Treatment Groups')
colors <- list('Cluster Groups' =c('1'='red', '2'='blue', '3'='green', '4'='orange', 
                                   '5'='orange', '6'='yellow', '7'='gold', '8'='violet'),
               'Treatment Groups' =c('Cancer'='pink', 'Healthy'='gray'))
ccplus_row_annot <- HeatmapAnnotation(df=ccplus_annot, which='row',
                                      col=colors)
m.ccplus <- cbind(t(cluster_data_input_normalized),cluster5000[[2]][[3]]$consensusClass)
Heatmap(m.pam[,1:5000], name="Scores", left_annotation=ccplus_row_annot, show_row_names = F, show_column_names = F, row_title="Genes", column_title = "Samples")


#K-means
kmeans_annot <- as.data.frame(cluster10000[[3]][["cluster"]])
kmeans_annot[,2] <- series_matrix$X12
colnames(kmeans_annot) <-c('Cluster Groups', 'Treatment Groups')
colors <- list('Cluster Groups' =c('1'='red', '2'='blue', '3'='green', '4'='orange', 
                                   '5'='orange', '6'='yellow', '7'='gold', '8'='violet'), 
               'Treatment Groups' =c('Cancer'='pink', 'Healthy'='gray'))
kmeans_row_annot <- HeatmapAnnotation(df=kmeans_annot, which='row',
                                      col=colors)
m.kmeans <- cbind(t(cluster_data_input_normalized),cluster10000[[3]][["cluster"]])
Heatmap(m.kmeans[,1:5000], name="Scores", left_annotation=kmeans_row_annot, show_row_names = F, show_column_names = F, row_title="Genes", column_title = "Samples")

