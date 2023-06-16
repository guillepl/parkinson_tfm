##############################
#### Methylation analysis ####
####  Guillermo Paz Lopez ####
####      05/04/2023      ####
##############################

rm(list=ls())

#### Set working directory ####
wd <- getwd()

#### Libray load ####
library(car)
library(ChAMP)
library(ChAMPdata)
library(EnhancedVolcano)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(limma)
library(minfi)
library(missMethyl)
library(pheatmap)
library(RColorBrewer)

#### Color palettes ####
color_duo <- c('royalblue','indianred')
colheatmap <- colheatmap <- colorRampPalette(c("royalblue", "white", "indianred"))(200)

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

#### Correction ####
beta_norm <- champ.norm(beta=beta_matrix,
                        arraytype="450K",
                        method='BMIQ',
                        cores=2)

#### Quality Control step #### 
##### Density plot of beta values - QCstep ####
png("Graficos/metilacion/1.beta_distribution.png", width = 150, height =100, units='mm', res = 300)
minfi::densityPlot(beta_matrix,
                   sampGroups = groups, 
                   xlim = c(-0.2,1.2), 
                   legend = FALSE,
                   main=expression(paste('Distribución de los valores ', beta,' en las muestras')),
                   xlab = 'Valores beta',
                   cex = 1.5)
legend('topright', 
       legend=c('PD','Control'), 
       title='Grupos',
       fill = c('#1B9E77', '#D95F02'), 
       cex = 0.7)
dev.off()

#### Batch effect and Singular Value Decomposition (SVD) plot ####
metadata_sel <- metadata[,c('Sample_Name','Sample_Group','Age','Sex')]

# We use beta_norm matrix as a dataframe to ensure champ.SVD works
batch_info <- champ.SVD(beta=as.data.frame(beta_norm), 
                        pd=metadata_sel)

# Plotting batch effect with drawheatmap ChAMP customized function
drawheatmap <- function(svdPV.m){
  myPalette <- c("darkred","red","orange","pink","white");
  breaks.v <- c(-10000,-10,-5,-2,log10(0.05),0);
  image(x=1:nrow(svdPV.m), y=1:ncol(svdPV.m), z=log10(svdPV.m), col=myPalette, breaks=breaks.v, xlab="", ylab="", axes=FALSE, main= "Singular Value Decomposition Analysis (SVD)");
  axis(1,at=1:nrow(svdPV.m),labels=paste("PC-",1:nrow(svdPV.m),sep=""),las=2);
  suppressWarnings(axis(2,at=1:ncol(svdPV.m),labels=colnames(svdPV.m),las=2));
  legend(x=-(topPCA/2.5),
         y=2.5,
         cex=0.7,
         legend=c(expression("p < 1x"~10^{-10}),
                  expression("p < 1x"~10^{-5}),"p < 0.01", "p < 0.05", "p > 0.05"), 
         fill=myPalette,
         par('usr')[2], 
         par('usr')[4], 
         xpd=NA);
}
topPCA <- 6
png("Graficos/metilacion/2.batch_effect.png", width = 200, height =100, units='mm', res = 300)
par(mar=c(5,15,2,1));
batch_heatmap <- drawheatmap(batch_info)
dev.off()

champ.SVD(beta=as.data.frame(beta_norm),
          pd=metadata_sel)

#### Normality and Homocedasticity tests ####
norm_test <- function(beta_norm){
  normality_pvalue_PD <- p.adjust(apply(beta_norm, 1, function(cg) shapiro.test(cg[1:14])$p.value), 'fdr')
  normality_pvalue_Control <- p.adjust(apply(beta_norm, 1, function(cg) shapiro.test(cg[15:28])$p.value), 'fdr')
  normality_df <- data.frame('CpG' = rownames(beta_norm),
                             'PD' = normality_pvalue_PD,
                             'Control' = normality_pvalue_Control,
                             'Normality' = normality_pvalue_PD > 0.05 & normality_pvalue_Control > 0.05)
  return(normality_df)
}

var_test <- function(beta_norm){
  variance_pvalue <- p.adjust(apply(beta_norm, 1, function(cg) leveneTest(cg ~ groups)$`Pr(>F)`[[1]]), 'fdr')
  variance_df <- data.frame('CpG' = rownames(beta_norm),
                            'levene_pvalue' = variance_pvalue,
                            'Homocedasticity' = variance_pvalue > 0.05)
  return(variance_df)
}

# We save the results in a dataframe
normality_df <- norm_test(beta_norm)
variance_df <- var_test(beta_norm)

# Proportion of CpG sites that are normal in both groups, PD and Control
prop.table(table(normality_df$Normality))
prop.table(table(variance_df$Homocedasticity))
# --> More than 50 % of CpG sites can be assumed as normal --> Parametric tests will be performed

#### Differentially Methylated Probes/Positions (DMPs) #### 
# Differential Methylation Analysis using Limma
design <- model.matrix(~ groups - 1)
colnames(design) <- c("Control", "PD")

fit <- lmFit(beta_norm, design)
contrast.matrix <- makeContrasts("PD-Control", levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
dmp <- topTable(fit2, number=nrow(beta_norm))

# Creating fold-change/delta(beta) cutoff
FCcutoff <- 0.1

##### Volcano-plot of CpGs ####
volcano_plot <- EnhancedVolcano(toptable = dmp,
                                lab = NA,
                                x = 'logFC',
                                y = 'P.Value',
                                ylim = c(0, max(-log10(dmp$P.Value))),
                                xlim = c(-0.3,0.3),
                                pCutoff = 0.05,
                                FCcutoff = FCcutoff,
                                legendLabels = c('NS', 
                                                 expression(paste(Delta,'(',beta,')')),
                                                 'p-value', 
                                                 expression(paste('p-value and ',Delta,'(',beta,')'))),
                                xlab = expression(paste(Delta,'(',beta,')')),
                                gridlines.major = FALSE,
                                gridlines.minor = FALSE)

# Selecting significative CpGs depending on pvalue
dmp_sig <- dmp[dmp$adj.P.Val <= 0.05,]

beta_norm_sig <- beta_norm[rownames(dmp_sig),]

# Adding information about each CpG
dmp_sig <- cbind(dmp_sig, probe.features[rownames(dmp_sig),])

# Creation of two new datasets: Hypermethylated and Hypomethylated CpG sites
dmp_sig_hyper <- dmp_sig[(dmp_sig$logFC > FCcutoff),]
dmp_sig_hypo <- dmp_sig[(dmp_sig$logFC < -FCcutoff),]

beta_norm_hyperhyposig <- beta_norm[c(rownames(dmp_sig_hyper),
                                      rownames(dmp_sig_hypo)),]

#### Checking distribution of CpG sites ####
cpg_distribution <- function(dmp_df, relation){
  if(relation == 'cgi' | relation == 'feature'){
    CpGposition <- data.frame(table(dmp_df[,relation]))
    CpGposition <- cbind(CpGposition, 
                         'Percentage'=round((as.numeric(CpGposition$Freq)/sum(CpGposition$Freq))*100,digits = 1))
    return(CpGposition)
  }else{print('Wrong relation, must be feature or cgi')}
}

# Distribution depending on neighborhood relation
CpGposition <- cpg_distribution(dmp_sig, 'cgi')
CpGpositionhyper <- cpg_distribution(dmp_sig_hyper, 'cgi')
CpGpositionhypo <- cpg_distribution(dmp_sig_hypo, 'cgi')

# Distribution depending on location within gene relation
CpGlocation <- cpg_distribution(dmp_sig, 'feature')
CpGlocationhyper <- cpg_distribution(dmp_sig_hyper, 'feature')
CpGlocationhypo <- cpg_distribution(dmp_sig_hypo, 'feature')

#####  Pie charts of CpG island features by neighborhood ####
# Representation of general pie chart, hypermethylated pie chart and hypomethylated pie chart
cpg_distribution_graph <- function(CpGposition_df){
  graph <- ggplot(CpGposition_df, aes(x='', y=Freq, fill=Var1)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    geom_text(aes(x = 1.7, label = paste0(Percentage, " %\n",'(',Freq,')')), position = position_stack(vjust=0.5)) +
    labs(x = NULL, y = NULL, fill = 'Region') +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_manual(values=c("#a3e797","#78aed9","#e3967d","#e3e797","#97e7ca","gray","palevioletred1"))
  print(graph)
}
cpgposition_general_pie <- cpg_distribution_graph(CpGposition)
cpgposition_hyper_pie <- cpg_distribution_graph(CpGpositionhyper)
cpgposition_hypo_pie <- cpg_distribution_graph(CpGpositionhypo)
cpglocation_general_pie <- cpg_distribution_graph(CpGlocation)
cpglocation_hypo_pie <- cpg_distribution_graph(CpGlocationhypo)

cpglocation_hyper_pie <-graph <- ggplot(CpGlocationhyper, aes(x='', y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(x = 1.7, label = paste0(Percentage, " %\n",'(',Freq,')')), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = 'Region') +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values=c("#78aed9","#e3967d","#e3e797","#97e7ca","gray","palevioletred1"))

#### Transposing beta values matrixes ####
beta_norm_sig_t <- t(beta_norm_sig)
beta_norm_hyperhyposig_t <- t(beta_norm_hyperhyposig)

#### Heatmap of significative CpGs ####
# Creating breaks for color scale centered around 0
breaksList = seq(min(scale(beta_norm_sig_t)), max(scale(beta_norm_sig_t)), by = 0.015)
colheatmap2 <- colorRampPalette(colheatmap)(length(breaksList))

# Generating heatmaps
heatmap_general <- pheatmap(scale(beta_norm_sig_t), 
                            color=colheatmap2,
                            show_colnames = FALSE,
                            cluster_rows = FALSE,
                            treeheight_col=0, 
                            treeheight_row=30,
                            breaks=breaksList)

breaksList = seq(min(scale(beta_norm_hyperhyposig_t)), max(scale(beta_norm_hyperhyposig_t)), by = 0.015)
colheatmap2 <- colorRampPalette(colheatmap)(length(breaksList))

heatmap_hyperhypo <- pheatmap(scale(beta_norm_hyperhyposig_t), 
                              color=colheatmap2,
                              show_colnames = FALSE,
                              cluster_rows = FALSE, 
                              treeheight_col=0, 
                              treeheight_row=30,
                              breaks=breaksList)

#### PCA on significative data ####
pca_fun <- function(beta_mx){
  pca_data <- PCA(beta_mx, scale.unit = TRUE, graph = FALSE)
  pca_plot <- fviz_pca_ind(pca_data, 
                           col.ind = groups,
                           palette = color_duo,
                           addEllipses = TRUE,
                           legend.title = "Grupos",
                           title = paste0('PCA - ',ncol(beta_mx),' CpG'),
                           label='none',
                           repel=TRUE,
                           show.legend=FALSE) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

pca_plot_general <- pca_fun(beta_norm_sig_t)
pca_plot_hyperhypo <- pca_fun(beta_norm_hyperhyposig_t)

#### Annotation: grouping CpG sites by Gene name ####
cpg_gene_count <- function(dmp_list){
  
  gene_list <- unique(na.omit(dmp_list$gene))
  gene_list <- gene_list[! gene_list %in% '']
  
  df <- data.frame('genes' = gene_list,
                   'dmp.num' = 0,
                   'dmp.list' = NA)
  
  for(i in 1:length(gene_list)){
    
    print(paste0(i,'/',length(gene_list)))#
    
    gene <- gene_list[i] 
    probes <- dmp_list[dmp_list$gene == gene,]
    dmp_gene <- rownames(dmp_list)[rownames(dmp_list) %in% rownames(probes)]
    if(length(dmp_gene) != 0){
      df[i,2] <- length(dmp_gene)
      df[i,3] <- paste0(dmp_gene, collapse=',')
    }
  }
  return(df)
}

gene_hypo <- cpg_gene_count(dmp_sig_hypo)
gene_hyper <- cpg_gene_count(dmp_sig_hyper)
gene_general <- cpg_gene_count(dmp_sig)


#### Differentially Methylated Regions (DMRs) *** #### 
set.seed(1234)
dmr <- champ.DMR(beta = beta_norm,
                 pheno = groups,
                 method = "Bumphunter")

#### Saving data and plots ####
save(beta_norm_sig, dmp_sig, dmp_sig_hyper, dmp_sig_hypo,
     gene_general, gene_hyper, gene_hypo, metadata, 
     probe.features, file='methylation_data.RData')

ggsave(
  './Graficos/metilacion/3.volcano_plot.png',
  volcano_plot,
  width = 6,
  height = 7,
  dpi = 1200
)

ggsave(
  './Graficos/metilacion/cpg_position_general_pie.png',
  cpgposition_general_pie,
  width = 6,
  height = 7,
  dpi = 1200
)

ggsave(
  './Graficos/metilacion/cpg_position_hyper_pie.png',
  cpgposition_hyper_pie,
  width = 6,
  height = 7,
  dpi = 1200
)

ggsave(
  './Graficos/metilacion/cpg_position_hypo_pie.png',
  cpgposition_hypo_pie,
  width = 6,
  height = 7,
  dpi = 1200
)

ggsave(
  './Graficos/metilacion/cpg_location_general_pie.png',
  cpglocation_general_pie,
  width = 6,
  height = 7,
  dpi = 1200
)

ggsave(
  './Graficos/metilacion/cpg_location_hyper_pie.png',
  cpglocation_hyper_pie,
  width = 6,
  height = 7,
  dpi = 1200
)

ggsave(
  './Graficos/metilacion/cpg_location_hypo_pie.png',
  cpglocation_hypo_pie,
  width = 6,
  height = 7,
  dpi = 1200
)

ggsave(
  './Graficos/metilacion/heatmap_general.png',
  heatmap_general,
  width = 7.5,
  height = 4,
  dpi = 1200
)

ggsave(
  './Graficos/metilacion/heatmap_hyperhypo.png',
  heatmap_hyperhypo,
  width = 6,
  height = 4,
  dpi = 1200
)

ggsave(
  './Graficos/metilacion/pca_general.png',
  pca_plot_general,
  width = 6,
  height = 4,
  dpi = 1200
)

ggsave(
  './Graficos/metilacion/pca_hyperhypo.png',
  pca_plot_hyperhypo,
  width = 6,
  height = 4,
  dpi = 1200
)