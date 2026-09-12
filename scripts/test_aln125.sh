#!/bin/bash
#SBATCH -A ansari.prj
#SBATCH -J 125
#SBATCH -p short
#SBATCH --mem=90G
#SBATCH -o logs/aln_125.out
#SBATCH -e logs/aln_125.err

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
#conda activate lcwgs_17
conda activate aln17

#snakemake -s test_vanilla_aln.smk --use-conda -c 2

home_dir="/well/ansari/users/osr869/"
sample="28EI-003-125"

fastq1="data/fastq_cleaned/${sample}_1.fastq.gz"
fastq2="data/fastq_cleaned/${sample}_2.fastq.gz"
reference="data/references/GRCh38.fa"

bam="data/bams/${sample}.bam"
bai="data/bams/${sample}.bam.bai"
tmp1="data/bams/${sample}.tmp1.bam"
metric="data/bams/${sample}.metrics.txt"

picard="java -Xmx80G -Xms20G -jar ${home_dir}conda/skylake/envs/aln17/share/picard-slim-2.27.4-0/picard.jar"

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