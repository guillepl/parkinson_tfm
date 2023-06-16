##############################
#### Methylation Analysis ####
####    Clusterization    ####             
####  Guillermo Paz Lopez ####
####      02/05/2023      ####
##############################

rm(list=ls())

#### Set working directory ####
wd <- getwd()

#### Libray load ####
library(ChAMP)
library(ChAMPdata)
library(factoextra)
library(FactoMineR)

#### Files load ####
# IDAT files
chload <- champ.load(paste0(wd, '/meth_data/data'))

# ChAMP CpG features
data(probe.features)

# Saving metadata as a matrix 
metadata <- chload$pd

# Erasing PD samples 046 and 048 (not in RNA-Seq)
metadata <- metadata[-c(10:11),]

# Creating a vector with group information 
groups <- as.factor(metadata$Sample_Group)

# Saving beta values as a matrix and erasing PD samples 046 and 048
beta_matrix <- chload$beta[,-c(10:11)]

# Changing samples' names
colnames(beta_matrix) <- gsub('MAYO_','',colnames(beta_matrix))
metadata$Sample_Name <- gsub('MAYO_','',colnames(beta_matrix))

# Changing sample names 
colnames(beta_matrix)[1:12] <- paste0(colnames(beta_matrix)[1:12], "_PD")
colnames(beta_matrix)[13:26] <- paste0(colnames(beta_matrix)[13:26], "_C")

#### Normalization ####
beta_norm <- champ.norm(beta=beta_matrix,
                        arraytype="450K",
                        method='BMIQ',
                        cores=2)

#### Selecting top variable CpG ####
beta_norm <- t(beta_norm)
CpG_mean <- apply(beta_norm, 2, mean)
CpG_sd <- apply(beta_norm, 2, sd)

names(CpG_mean) <- colnames(CpG_mean)
names(CpG_sd) <- colnames(CpG_sd)

# Selecting CpG with sd higher than 50 % of mean
sd_mean_comparison <- function(a=CpG_mean, b=CpG_sd, p=0.5){
  genes_sel <- c()
  for (i in 1:length(a)){
    if (b[i] > (a[i] * p))
    {
      genes_sel <- c(genes_sel, i)
    }
  }
  return(genes_sel)
}

CpG_sel <- sd_mean_comparison()
length(CpG_sel)

beta_norm_var <- beta_norm[,CpG_sel]

#### K-means clustering analysis ####
# We compute the optimal numver of clusters
kmeans_number <- fviz_nbclust(x = beta_norm_var, 
                              FUNcluster = kmeans,
                              method = 'wss', 
                              k.max = 15, 
                              diss = get_dist(beta_norm_var, method = 'euclidean'), 
                              nstart = 50)

# Clustering observations depending on number of optimal clusters
km_clusters <- kmeans(x = beta_norm_var, 
                      centers = 4, 
                      nstart = 50)

# K-means partitioning graph using previous clusterization
kmeans_cluster <- fviz_cluster(object = km_clusters, 
                               data = beta_norm_var,
                               show.clust.cent = TRUE,
                               ellipse.type = 'euclid', 
                               star.plot = TRUE, 
                               repel = TRUE,
                               geom=c('point', 'text'),
                               show.legend.text = FALSE) +
  labs(title = 'K-means clustering (K = 4)') +
  theme(legend.position = 'right',
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        axis.line = element_line(colour = "black"))

#### Generating and saving graphs ####
ggsave(
  './kmeans_sd_number.png',
  kmeans_number,
  width = 8,
  height = 3,
  dpi = 1200
)

ggsave(
  './Graficos/metilacion/0.kmeans_sd_cluster.png',
  kmeans_cluster,
  width = 5,
  height = 4,
  dpi = 1200
)
