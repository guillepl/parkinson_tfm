##############################
####  Enrichment analysis ####
####  Guillermo Paz Lopez ####
####      02/05/2023      ####
##############################

rm(list=ls())

#### Set working directory ####
wd <- getwd()

#### Libraries load ####
library(bedr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(openxlsx)
library(org.Hs.eg.db)
library(pathview)

# Loading data 
genes_rna <- scan('integration_data/genes_enrichment_rna.txt',
                         what='charachter')
genes_rnameth <- scan('integration_data/genes_enrichment_rnameth.txt',
                      what='charachter')

# Transforming ENTREZ ID to numeric and ordering them
gene_conv <- function(genes){
  genes_conversion <- mapIds(org.Hs.eg.db, 
                             genes, 
                             'ENTREZID', 
                             'SYMBOL')
  genes_conversion <- na.omit(genes_conversion)
  genes_conversion <- genes_conversion
  genes_conversion <- as.numeric(genes_conversion)
  genes_conversion <- sort(genes_conversion,decreasing = T)
  return(genes_conversion)
}

genes_rna <- gene_conv(genes_rna)
genes_rnameth <- gene_conv(genes_rnameth)

#### GOEA of General Genes from (RNA) ####
go_fun <- function(genes, feature){
  go_res <- enrichGO(gene          =  genes,
                     OrgDb         = 'org.Hs.eg.db',
                     keyType       = 'ENTREZID',
                     ont           = feature, 
                     pAdjustMethod = 'fdr',
                     pvalueCutoff  = 0.05)
  return(go_res)
}

go_mf_rna <- go_fun(genes_rna, 'MF')
go_bp_rna <- go_fun(genes_rna, 'BP')
go_cc_rna <- go_fun(genes_rna, 'CC')

go_mf_rna_sig <- go_mf_rna@result[go_mf_rna@result$p.adjust<=0.05,]
go_bp_rna_sig <- go_bp_rna@result[go_bp_rna@result$p.adjust<=0.05,]
go_cc_rna_sig <- go_cc_rna@result[go_cc_rna@result$p.adjust<=0.05,]

##### Molecular function ####
mf_rgo<- setReadable(go_mf_rna, OrgDb = 'org.Hs.eg.db', 'ENTREZID')
mf_simgo <- pairwise_termsim(mf_rgo)
  
# Dot plot of significant GO terms
mf_dotplot <- dotplot(mf_rgo, showCategory=20)
  
# Significant GO terms Directed Acyclic Graph (DAG)
mf_dagplot <- goplot(mf_rgo)
  
# Gene-concept Network
mf_cnetplot <- cnetplot(mf_simgo,
                        categorySize="pvalue",
                        cex_label_category = 0.6, 
                        cex_label_gene=0.6)

##### Biological processes ####
bp_rgo<- setReadable(go_bp_rna, OrgDb = 'org.Hs.eg.db', 'ENTREZID')
bp_simgo <- pairwise_termsim(bp_rgo)

# Dot plot of significant GO terms
bp_dotplot <- dotplot(bp_rgo, showCategory=20)

# Significant GO terms Directed Acyclic Graph (DAG)
bp_dagplot <- goplot(bp_rgo)

# Treeplot
bp_treeplot <- treeplot(bp_simgo)

# Gene-concept Network
bp_cnetplot <- cnetplot(bp_simgo, 
                        showCategory = c('proteasome-mediated ubiquitin-dependent protein catabolic process',
                                         'protein monoubiquitination',
                                         'internal protein amino acid acetylation',
                                         'peptidyl-lysine acetylation',
                                         'peptidyl-lysine modification',
                                         'regulation of protein stability'),
                        categorySize="pvalue",
                        cex_label_category = 0.6, 
                        cex_label_gene=0.6)

##### Cellular components ####
cc_rgo<- setReadable(go_cc_rna, OrgDb = 'org.Hs.eg.db', 'ENTREZID')
cc_simgo <- pairwise_termsim(cc_rgo)

# Dot plot of significant GO terms
cc_dotplot <- dotplot(cc_rgo, showCategory=20)

# Significant GO terms Directed Acyclic Graph (DAG)
cc_dagplot <- goplot(cc_rgo)

# Treeplot
cc_treeplot <- treeplot(cc_simgo)

# Gene-concept Network
cc_cnetplot <- cnetplot(cc_simgo, 
                        showCategory = c('inner mitochondrial membrane protein complex',
                                         'mitochondrial inner membrane',
                                         'mitochondrial protein-containing complex',
                                         'mitochondrial matrix',
                                         'ATPase complex'),
                        categorySize="pvalue",
                        cex_label_category = 0.6, 
                        cex_label_gene=0.6)

##### Saving plots #### 
# Molecular functions
ggsave(
    './Graficos/enriquecimiento/genes_rna/mf_dotplot.png',
    mf_dotplot,
    width = 6,
    height = 4,
    dpi = 1200
)
  
ggsave(
    './Graficos/enriquecimiento/genes_rna/mf_dagplot.png',
    mf_dagplot,
    width = 6,
    height = 4,
    dpi = 1200
)
  
ggsave(
    './Graficos/enriquecimiento/genes_rna/mf_cnetplot.png',
    mf_cnetplot,
    width = 9,
    height = 7,
    dpi = 1200
)

# Biological processes
ggsave(
  './Graficos/enriquecimiento/genes_rna/bp_dotplot.png',
  bp_dotplot,
  width = 6,
  height = 12,
  dpi = 1200
)

ggsave(
  './Graficos/enriquecimiento/genes_rna/bp_dagplot.png',
  bp_dagplot,
  width = 12,
  height = 8,
  dpi = 1200
)

ggsave(
  './Graficos/enriquecimiento/genes_rna/bp_treeplot.png',
  bp_treeplot,
  width = 16,
  height = 8,
  dpi = 1200
)

ggsave(
  './Graficos/enriquecimiento/genes_rna/bp_cnetplot.png',
  bp_cnetplot,
  width = 12,
  height = 10,
  dpi = 1200
)

# Cellular components
ggsave(
  './Graficos/enriquecimiento/genes_rna/cc_dotplot.png',
  cc_dotplot,
  width = 6,
  height = 9,
  dpi = 1200
)

ggsave(
  './Graficos/enriquecimiento/genes_rna/cc_dagplot.png',
  cc_dagplot,
  width = 12,
  height = 8,
  dpi = 1200
)

ggsave(
  './Graficos/enriquecimiento/genes_rna/cc_treeplot.png',
  cc_treeplot,
  width = 13,
  height = 8,
  dpi = 1200
)

ggsave(
  './Graficos/enriquecimiento/genes_rna/cc_cnetplot.png',
  cc_cnetplot,
  width = 11,
  height = 9,
  dpi = 1200
)
