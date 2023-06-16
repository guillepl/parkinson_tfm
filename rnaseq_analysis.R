##############################
####    RNA-Seq analysis  ####
####  Guillermo Paz Lopez ####
####      02/05/2023      ####
##############################

rm(list=ls())

#### Set working directory ####
setwd("C:/Users/Guillermo Paz/OneDrive/ECAI Bioinformática/Tareas bioinformáticas/Ejercicios Guillermo/TFM")
wd <- getwd()

#### Libray load ####
library(circlize)
library(corrplot)
library(dplyr)
library(edgeR)
library(EnhancedVolcano)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(limma)
library(mixOmics)
library(pheatmap)
library(RColorBrewer)
library(readr)
library(readxl)

#### Color palettes ####
color_duo <- c('royalblue','indianred')
colheatmap <- colheatmap <- colorRampPalette(c("royalblue", "white", "indianred"))(200)

#### Files load ####
metadata <- read_excel("metadata_metilacion_rna.xlsx",
                       skip = 1)

counts_matrix <- as.data.frame(readr::read_delim('rna_data/counts_matrix.txt', 
                                                 delim = "\t", 
                                                 escape_double = FALSE, 
                                                 trim_ws = TRUE))

groups <- as.factor(c(rep('PD',12),rep('Control',14)))
batch_sex <- metadata$Sex
batch_age <- metadata$Age

# Omitting first 2 rows of metadata file (no RNA-Seq data)
metadata <- metadata[-c(1:2),]

# Gene name as row names
rownames(counts_matrix) <- counts_matrix$gene_name
counts_matrix <- as.matrix(counts_matrix[,-1])

# Omitting rowname (MissingGeneID), exons that has no Gene ID associated
counts <- counts_matrix[-which(rownames(counts_matrix)=='MissingGeneID'),]

# Changing sample names of counts_matrix
colnames(counts) <- metadata$Sample

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

# Normalize using scaling factors
counts_norm <- cpm(dge_object,log = T, prior.count=2)
counts_norm <- t(counts_norm)

#### Differential Gene Expression Analysis (DGEA) ####
# DGEA using Wilcoxon
deg_wilcoxon <- function(counts_mx){
  statistic_wilcox <- apply(counts_mx,2,function(x){wilcox.test(x ~ groups)$statistic})
  pvalores_wilcox <- apply(counts_mx,2,function(x){wilcox.test(x ~ groups)$p.value})
  pvalores_wilcox_adj <- p.adjust(pvalores_wilcox, method = "fdr")
  
  comparative_PDvsC <- data.frame('gene' = colnames(counts_mx),
                                  'pvalor' = pvalores_wilcox,
                                  'pvalor_adj' = pvalores_wilcox_adj,
                                  'W_statistic' = statistic_wilcox)
}

deg <- deg_wilcoxon(counts_norm)

# DGEA using edgeR for logFC
d1 <- estimateCommonDisp(dge_object, verbose=T)
d1 <- estimateTagwiseDisp(d1)
et12 <- exactTest(d1, pair=c(1,2))
deg_edge <- et12$table

# Adding logFC from edgeR to Wilcoxon significative genes
deg <- cbind(deg, 'logFC' = deg_edge$logFC[rownames(deg_edge) %in% deg$gene])

# Creating fold-change cutoff
FCcutoff <- 1

##### Volcano-plot of genes ####
volcano_plot <- EnhancedVolcano(toptable = deg,
                                lab = NA,
                                x = 'logFC',
                                y = 'pvalor',
                                ylim = c(0, max(-log10(deg$pvalor))),
                                xlim = c(-2,2),
                                pCutoff = 0.05,
                                FCcutoff = FCcutoff,
                                gridlines.major = FALSE,
                                gridlines.minor = FALSE)

# Selecting DEG from p-value
deg_sig <- deg[deg$pvalor <= 0.05,]
deg_sig_up <- deg_sig[deg_sig$logFC >= 1, ]
deg_sig_down <- deg_sig[deg_sig$logFC <= -1, ]

# Selecting significative normalized counts 
counts_sig <- counts_norm[,deg_sig$gene]
counts_sig_updown <- counts_norm[,c(deg_sig_up$gene, deg_sig_down$gene)]

#### Heatmap of significative DEG ####
# Creating breaks for color scale centered around 0
# Creating breaks for color scale centered around 0
breaksList = seq(min(scale(counts_sig)), 4.1, by = 0.015)
colheatmap2 <- colorRampPalette(colheatmap)(length(breaksList))

heatmap_general <- pheatmap(scale(counts_sig),
                            color=colheatmap2,
                            show_colnames = FALSE,
                            cluster_rows = FALSE,
                            treeheight_col=0, 
                            treeheight_row=30,
                            breaks=breaksList)

breaksList = seq(min(scale(counts_sig_updown)), max(scale(counts_sig_updown)), by = 0.015)
colheatmap2 <- colorRampPalette(colheatmap)(length(breaksList))

heatmap_updown <- pheatmap(scale(counts_sig_updown),
                            color=colheatmap2,
                            show_colnames = FALSE,
                            cluster_rows = FALSE,
                            treeheight_col=0, 
                            treeheight_row=30,
                            breaks=breaksList)

#### Boxplot of top significative DEG by logFC ####
# Funtion to find genes with sample outliers 
find_outlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

# Functions to create dataframe with plot information, and plots
box_genes_data <- function(deg_df, n){
  deg_df <- deg_df[order(abs(deg_df$logFC), decreasing=TRUE),]
  deg_df <- deg_df[1:n,]
  counts_deg <- as.matrix(counts_norm[,deg_df$gene])
  df_graph <- data.frame(samp = rep(rownames(counts_norm), times=nrow(deg_df)),
                         gene = rep(deg_df$gene, each=length(groups)), 
                         value = c(counts_deg),
                         logFC = rep(deg_df$logFC, each=length(groups)), 
                         group = rep(groups, times=nrow(deg_df)))
  df_graph <- df_graph %>%
    group_by(gene, group) %>%
    mutate(outlier = ifelse(find_outlier(value), samp, NA))
  return(df_graph)
}

box_genes_plot <-function(df_graph){
  gene_graph <- ggplot(df_graph, aes(x=factor(gene, levels = unique(df_graph$gene)), 
                                              y=value, 
                                              fill=factor(group))) + 
    geom_boxplot() + 
    stat_boxplot(geom = "errorbar", 
               width=0.2, 
               position=position_dodge(0.75)) +
    scale_fill_manual(values=color_duo) + 
    labs(x = 'Genes', 
       y = 'Nivel de expresion', 
       fill='Grupo') + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  print(gene_graph)
  }

# Creating df and plots of top upregulated and downregulated genes
deg_up_gdf <- box_genes_data(deg_sig_up,6)
deg_up_gplot <- box_genes_plot(deg_up_gdf)

deg_down_gdf <- box_genes_data(deg_sig_down,6)
deg_down_gplot <- box_genes_plot(deg_down_gdf)

#### PCA on significative data ####
pca_fun <- function(mx){
  pca_data <- PCA(mx, scale.unit = TRUE, graph = FALSE)
  pca_plot <- fviz_pca_ind(pca_data, 
                           col.ind = groups,
                           palette = color_duo,
                           addEllipses = TRUE,
                           legend.title = "Grupos",
                           title = paste0('PCA - ',ncol(mx),' genes'),
                           label='none',
                           repel=TRUE,
                           show.legend=FALSE) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

pca_plot_general <- pca_fun(counts_sig)
pca_plot_updown <- pca_fun(counts_sig_updown)


#### Saving data and plots ####
ggsave(
  './Graficos/rna_seq/volcano_plot.png',
  volcano_plot,
  width = 9.5,
  height = 7,
  dpi = 1200
)

ggsave(
  './Graficos/rna_seq/2.heatmap_general.png',
  heatmap_general,
  width = 7.5,
  height = 4,
  dpi = 1200
)

ggsave(
  './Graficos/rna_seq/2.heatmap_updown.png',
  heatmap_updown,
  width = 7.5,
  height = 4,
  dpi = 1200
)

ggsave(
  './Graficos/rna_seq/3.deg_up.png',
  deg_up_gplot,
  width = 6.5,
  height = 4.5,
  dpi = 1200
)

ggsave(
  './Graficos/rna_seq/3.deg_down.png',
  deg_down_gplot,
  width = 6.5,
  height = 4.5,
  dpi = 1200
)

ggsave(
  './Graficos/rna_seq/4.pca_general.png',
  pca_plot_general,
  width = 6,
  height = 4,
  dpi = 1200
)

ggsave(
  './Graficos/rna_seq/4.pca_updown.png',
  pca_plot_updown,
  width = 6,
  height = 4,
  dpi = 1200
)
