# Parkinson_TFM
This is the code I used to perform the different analyses to study differential expression and methylation between control and PD samples.

In order of use:
- genome_index.sh: Bash. It's used to create the genome index using STAR aligner.
- align_script.sh: Bash. It's used to align the paired-end FASTQ files to the reference genome.
- rnaseq_kcluster: R. It's used as a control step to study differences between samples and check for non identified groups using RNA-Seq data.
- methylation_kcluster: R. Same as previoues but using methylation data (450k array).
- rnaseq_analysis: R. Steps used to analyze expression data and detect differentially expressed genes.
- methylation_analysis: R. Steps used to analyze methylation data an detect differentially methylated positions/probes.
- enrichment_analysis: R. Steps used to perform gene ontology enrichment analysis using previous selected genes.
- integration_analysis: R. Steps used to perform eQTM analysis and integration analysis using DIABLO tool from mixOmics R package.
