##############################
####    RNA-Seq Analysis  ####
####    Clusterization    ####             
####  Guillermo Paz Lopez ####
####      02/05/2023      ####
##############################

rm(list=ls())

#### Libray load ####
library(edgeR)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(mixOmics)
library(readr)
library(readxl)

#### Files load ####
metadata <- read_excel("metadata_metilacion_rna.xlsx",
                       skip = 1)

counts_matrix <- as.data.frame(readr::read_delim("./rna_data/counts_matrix.txt", 
                                                 delim = "\t", 
                                                 escape_double = FALSE, 
                                                 trim_ws = TRUE))

groups <- factor(c(rep('PD',12),rep('C',14)), levels=c('PD','C'))

# Omitting first 2 rows of metadata file (no RNA-Seq data)
metadata <- metadata[-c(1:2),]

# Gene name as row names
rownames(counts_matrix) <- counts_matrix$gene_name
counts_matrix <- as.matrix(counts_matrix[,-1])

# Changing sample names of counts_matrix
colnames(counts_matrix) <- metadata$Sample

# Omitting rowname (MissingGeneID), exons that has no Gene ID associated
counts <- counts_matrix[-which(rownames(counts_matrix)=='MissingGeneID'),]

#### Counts matrix correction and normalization ####
# Create DGE object with count and gene information
dge_object <- DGEList(counts=counts, 
                      group=groups)

# Identify poorly expresed genes (very few counts in both groups)
omit_genes <- filterByExpr(dge_object)

# Omit poorly expresed genes without keeping sizes
dge_object <- dge_object[omit_genes,,keep.lib.sizes=FALSE]

# Calculate normalization factors (sizes) using the Trimmed Mean of M-values 
# Method proposed by Robinson and Oshlack (2010)
dge_object <- calcNormFactors(dge_object, method = "TMM")

# Normalize using scaling factors and Compute counts per million (CPM) function
counts_norm <- as.data.frame(t(cpm(dge_object,log = T)))

#### Selecting top variable genes ####
genes_mean <- apply(counts_norm, 2, mean)
genes_sd <- apply(counts_norm, 2, sd)

names(genes_mean) <- colnames(counts_norm)
names(genes_sd) <- colnames(counts_norm)

# Selecting genes with sd higher than 10 % of mean
sd_mean_comparison <- function(a=genes_mean, b=genes_sd){
  genes_sel <- c()
  for (i in 1:length(a)){
    if (b[i] > (a[i] * 0.1))
    {
      genes_sel <- c(genes_sel, i)
    }
  }
  return(genes_sel)
}

genes_sel <- sd_mean_comparison()

counts_norm_var <- counts_norm[,genes_sel]

#### K-means clustering analysis ####
# We compute the optimal numver of clusters
kmeans_number <- fviz_nbclust(x = counts_norm_var, 
                              FUNcluster = kmeans,
                              method = 'wss', 
                              k.max = 15, 
                              diss = get_dist(counts_norm_var, method = 'euclidean'), 
                              nstart = 50)

# Clustering observations depending on number of optimal clusters
km_clusters <- kmeans(x = counts_norm_var, 
                      centers = 2, 
                      nstart = 50)

# K-means partitioning graph using previous clusterization
kmeans_cluster <- fviz_cluster(object = km_clusters, 
                               data = counts_norm_var,
                               show.clust.cent = TRUE,
                               ellipse.type = 'euclid', 
                               star.plot = TRUE, 
                               repel = TRUE,
                               geom=c('point', 'text'),
                               show.legend.text = FALSE) +
  labs(title = 'K-means clustering (K = 2)') +
  theme(legend.position = 'right',
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        axis.line = element_line(colour = "black"))

#### Generating and saving graphs ####
ggsave(
  './Graficos/rna_seq/0.kmeans_sd_number.png',
  kmeans_number,
  width = 8,
  height = 3,
  dpi = 1200
)

ggsave(
  './Graficos/rna_seq/0.kmeans_sd_cluster.png',
  kmeans_cluster,
  width = 4,
  height = 4,
  dpi = 1200
)
