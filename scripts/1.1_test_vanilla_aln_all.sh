#!/bin/bash
#SBATCH -A ansari.prj
#SBATCH -J aln
#SBATCH -p short
#SBATCH -o logs/test_aln.out
#SBATCH -e logs/test_aln.err

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
#conda activate lcwgs_17
conda activate aln17

#!/usr/bin/env bash
set -euo pipefail

home_dir="/well/ansari/users/osr869/"
samples=(
  "28EI-003-125"
  "28EI-003-045"
  "28EI-003-118"
  "28EI-003-084"
  "28EI-003-079"
  "28EI-003-069"
)

reference="data/references/GRCh38.fa"
mkdir -p data/bams

for sample in "${samples[@]}"; do
  fastq1="data/fastq_cleaned/${sample}_1.fastq.gz"
  fastq2="data/fastq_cleaned/${sample}_2.fastq.gz"

  bam="data/bams/${sample}.bam"
  bai="data/bams/${sample}.bam.bai"
  tmp1="data/bams/${sample}.tmp1.bam"
  metric="data/bams/${sample}.metrics.txt"

  picard="java -Xmx40G -Xms20G -jar ${home_dir}conda/skylake/envs/vani_aln17/share/picard-slim-2.27.4-0/picard.jar"

  bwa mem -t 6 "$reference" "$fastq1" "$fastq2" | samtools view -b -o "$tmp1"

  samtools sort -@6 -m 1G -o "$bam" "$tmp1"
  samtools index "$bam"

  picard AddOrReplaceReadGroups \
    -VERBOSITY ERROR \
    -I "$bam" \
    -O "$tmp1" \
    -RGLB OGC \
    -RGPL Illumina \
    -RGPU unknown \
    -RGSM "$sample"

  $picard FixMateInformation -I "$tmp1"

  samtools sort -@6 -m 1G -o "$bam" "$tmp1"

  $picard MarkDuplicates \
    -I "$bam" \
    -O "$tmp1" \
    -M "$metric" \
    --REMOVE_DUPLICATES

  samtools sort -@6 -m 1G -o "$bam" "$tmp1"
  samtools index "$bam"

  rm -f "$tmp1"
done