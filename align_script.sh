#!/usr/bin/env bash
# The name to show in queue lists for this job:
#SBATCH -J align_TFM

# Number of desired cpus:
#SBATCH --cpus-per-task=8

# Amount of RAM needed for this job:
#SBATCH --mem=48gb

# The time the job will be running:
#SBATCH --time=100:00:00

# To use GPUs you have to request them:
##SBATCH --gres=gpu:1

# Set output and error files
#SBATCH --error=align_TFM.%J.err
#SBATCH --output=align_TFM.%J.out

# To load some software (you can show the list with 'module avail'):
module load star/2.7.9a

# the program to execute with its parameters:
hostname

FILES=()
n=0
for file in data_trim/*.fastq.gz; do sample=${file:10:11};
    FILES[n]+="$sample"
    n=$((n+1))
done
echo ${FILES[*]}

FILES_u=($(for file in "${FILES[@]}"; do echo "${file}"; done | sort -u))
echo ${FILES_u[*]}
 
for file in "${FILES_u[@]}"; do
    echo "$file"

  STAR --runMode alignReads \
  --outFilterType BySJout \
  --outFilterMultimapNmax 20 \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \
  --runThreadN 8 \
  --genomeDir hg38_index/ \
  --readFilesIn data_trim/${file}_1_trim.fastq.gz data_trim/${file}_2_trim.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --alignEndsType Local \
  --quantMode GeneCounts \
  --outFileNamePrefix alignment/${file}

done