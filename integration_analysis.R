##############################
#### Integration analysis ####
####  Guillermo Paz Lopez ####
####      02/05/2023      ####
##############################

rm(list=ls())

#### Set working directory ####
wd <- getwd()

#### Libray load ####
library(biomaRt)
library(MatrixEQTL)
library(mixOmics)

#### Data load #### 
load('methylation_data.RData')
load('rna_data.RData')

#### Preparing covariates data ####
covariates <- rbind(metadata$Age, metadata$Sex)
covariates[2,] <- ifelse(covariates[2,] == 'F', 0, 1)
covariates <- apply(covariates, 2, as.numeric)
colnames(covariates) <- metadata$Sample
rownames(covariates) <- c('Age','Sex') 

#### Preparing data format for MatrixEQTL analysis ####
# Preparing matrixes colnames and order
counts_sig <- t(counts_sig)
counts_sig <- counts_sig[,colnames(beta_norm_sig)]

##### Preparing DMP name and location format ####
dmp_location <- data.frame('snpid' = rownames(dmp_sig),
                           'chr' = paste0('chr',dmp_sig$CHR),
                           'pos' = dmp_sig$MAPINFO)

dmp_value <- beta_norm_sig

##### Preparing DEG name and location format ####
human38 = useMart(biomart = 'ensembl',
                  dataset = 'hsapiens_gene_ensembl')

deg_location <- getBM(attributes=c('hgnc_symbol', 
                               'chromosome_name', 
                               'start_position', 
                               'end_position'),
                  filters='hgnc_symbol',
                  values=list(deg_sig$gene),
                  mart=human38)

deg_location <- deg_location[!grepl('_', deg_location$chromosome_name) ,]

colnames(deg_location) <- c('geneid','chr','left','right')

deg_location$chr <- paste0('chr',deg_location$chr)

deg_value <- counts_sig

##### Using matrixEQTL ####
base.dir = find.package('MatrixEQTL')

# Reading CpG information
cpgs = SlicedData$new()
cpgs$fileDelimiter = "\t"      
cpgs$fileOmitCharacters = "NA"
cpgs$fileSkipRows = 1          
cpgs$fileSkipColumns = 1       
cpgs$fileSliceSize = 2000      
cpgs$LoadFile('integration_data/dmp_value.txt')

# Reading gene information
gene = SlicedData$new()
gene$fileDelimiter = "\t"      
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1          
gene$fileSkipColumns = 1       
gene$fileSliceSize = 2000      
gene$LoadFile('integration_data/deg_value.txt')

# Reading covariates information
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      
cvrt$fileOmitCharacters = "NA"; 
cvrt$fileSkipRows = 1;          
cvrt$fileSkipColumns = 1;       
if(length('integration_data/covariates.txt')>0) {
  cvrt$LoadFile('integration_data/covariates.txt');
}

# Running the analysis
cpgspos = read.table('integration_data/dmp_location.txt', header = TRUE, stringsAsFactors = FALSE);
genepos = read.table('integration_data/deg_location.txt', header = TRUE, stringsAsFactors = FALSE);

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR 

output_file_name_cis = tempfile()
output_file_name_tra = tempfile()

pvOutputThreshold_cis = 0.05
pvOutputThreshold_tra = 0.05

errorCovariance = numeric()

cisDist = 2000

me = Matrix_eQTL_main(
  snps = cpgs, 
  gene = gene, 
  cvrt = SlicedData$new(),
  output_file_name      = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis  = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = cpgspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

# Saving results as dataframes
dmp_cis <- me$cis$eqtls
dmp_trans <- me$trans$eqtls

dmp_cis_sig <- dmp_cis[dmp_cis$FDR <= 0.05,]
dmp_trans_sig <- dmp_trans[dmp_trans$FDR <= 0.05,]

dmp_trans_sig_gene <- table(dmp_trans_sig$gene)[order(table(dmp_trans_sig$gene), decreasing=TRUE)]

#### Finding DEGs with DMPs #### 
length(deg_sig$gene[deg_sig$gene %in% gene_general$genes]) # 369
length(deg_sig_up$gene[deg_sig_up$gene %in% gene_hyper$genes])
length(deg_sig_down$gene[deg_sig_down$gene %in% gene_hyper$genes])
length(deg_sig_up$gene[deg_sig_up$gene %in% gene_hypo$genes]) # 1
length(deg_sig_down$gene[deg_sig_down$gene %in% gene_hypo$genes])

genes_enrichment <- deg_sig$gene[deg_sig$gene %in% gene_general$genes]

#### Integration using multiblock sPLS (DIABLO) ####
# Creating data
counts_sig <- counts_sig[colnames(beta_norm_sig),]
omics_list <- list('meth'=t(beta_norm_sig), 'rna'=counts_sig)

groups <- as.factor(metadata$Sample_Group)

# Calculating correlation between the two omics
pls_cor <- pls(omics_list$rna, omics_list$meth, ncomp = 1)
cor_value <- round(cor(pls_cor$variates$X, pls_cor$variates$Y),
                   digits=2)

##### Testing for best variable number ####
# Creating a design matrix
design <- matrix(cor_value, 
                ncol = length(omics_list), 
                nrow = length(omics_list), 
                dimnames = list(names(omics_list), names(omics_list)))
diag(design) <- 0

# Training a model to select min variables that reduce error rate
set.seed(123) 
test.keepX <- list(meth = c(seq(100, 200, 25)),
                   rna = c(seq(50, 100, 15)))

blockpls_testing <- tune.block.splsda(omics_list, 
                                      groups, 
                                      ncomp = 2, 
                                      design = design,
                                      test.keepX = test.keepX,
                                      validation = 'Mfold', 
                                      folds = 10, 
                                      nrepeat = 5, 
                                      BPPARAM = BiocParallel::SnowParam(workers = 2),
                                      dist = 'centroids.dist')

list.keepX <- blockpls_testing$choice.keepX

# Creating sPLS-DA (DIABLO) model with variables selected
blockpls_result <- block.splsda(omics_list, 
                                groups, 
                                ncomp = 2,
                                keepX = list.keepX,
                                design = design)

#### Plotting and saving componentes and variables' correlation ####
col_duo <- c('royalblue','indianred')
names(col_duo) <- c('PD','Control')
colheatmap <- colheatmap <- colorRampPalette(c("royalblue", "white", "indianred"))(200)

# First component correlation
png("Graficos/integracion/diablo_cor_comp1.png", width = 150, height = 100, units='mm', res = 300)
plotDiablo(blockpls_result, ncomp = 1, col.per.group = col_duo)
dev.off()

# Second component correlation
png("Graficos/integracion/diablo_cor_comp2.png", width = 150, height = 100, units='mm', res = 300)
plotDiablo(blockpls_result, ncomp = 2, col.per.group = col_duo)
dev.off()

# Components and variables' correlation
png("Graficos/integracion/diablo_var_circ.png", width = 150, height = 150, units='mm', res = 300)
plotVar(blockpls_result, var.names = FALSE, 
        style = 'graphics', legend = TRUE, 
        pch = c(16, 17), cex = c(1,1), 
        col = col_duo)
dev.off()

# Loadings of top contributive variables ***
png("Graficos/integracion/var_contrib.png", width = 200, height = 100, units='mm', res = 300)
plotLoadings(blockpls_result, 
             comp = 1, 
             ndisplay = 20,
             contrib = 'max', 
             method = 'median',
             legend.color = col_duo)
dev.off()

# Heatmap - Clustered Image Map (CIM)
png("Graficos/integracion/diablo_hm_cluster.png", width = 400, height = 150, units='mm', res = 300)
cimDiablo(blockpls_result, 
          color = colheatmap,
          color.blocks = c('deeppink4','gold4'),
          color.Y = c('cyan4','darkorange3'),
          comp = 1, 
          margin=c(8,20), 
          legend.position = "right")
dev.off()

#### Saving data ####
write.table(selectVar(blockpls_result, block = 'rna', comp = 1), 
            file='genes_blocksplsda.txt',
            row.names=FALSE,
            quote=FALSE)

ggsave(
  './Graficos/integracion/diablo_variable_selection.png',
  plot(blockpls_testing, col = col_duo),
  width = 6,
  height = 4,
  dpi = 1200
)
