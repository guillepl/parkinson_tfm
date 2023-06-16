#!/usr/bin/env bash
# The name to show in queue lists for this job:
#SBATCH -J genome_index

# Number of desired cpus:
#SBATCH --cpus-per-task=8

# Amount of RAM needed for this job:
#SBATCH --mem=48gb

# The time the job will be running:
#SBATCH --time=100:00:00

# To use GPUs you have to request them:
##SBATCH --gres=gpu:1

# Set output and error files
#SBATCH --error=genome_index.%J.err
#SBATCH --output=genome_index.%J.out

# To load some software (you can show the list with 'module avail'):
module load star/2.7.9a

# the program to execute with its parameters:
hostname

STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir hg38_index \
--genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile Homo_sapiens.GRCh38.109.gtf \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentGene gene_name \
--sjdbGTFtagExonParentTranscript transcript_name \
--sjdbOverhang 50


