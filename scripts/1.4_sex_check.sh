#!/usr/bin/env bash
#SBATCH --job-name=sexcheck
#SBATCH --array=1-30%7
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/sexcheck_%A_%a.out
#SBATCH --error=logs/sexcheck_%A_%a.err

module load BCFtools
module load PLINK

samples_tsv="/well/ansari/users/osr869/host/lcwgs_edits/results/X_check/test14.tsv"
ref_fasta="/well/ansari/users/osr869/host/lcwgs_edits/data/references/GRCh38.fa"
outdir="/well/ansari/users/osr869/host/lcwgs_edits/results/X_check"

chrX_nonPAR="chrX:2781480-155701382"

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samples_tsv")
sample_id=$(printf '%s' "$line" | cut -f1)
bam="/well/ansari/users/osr869/host/lcwgs_edits/data/bams/${sample_id}.bam"

bcf="${outdir}/${sample_id}_chrX_nonPAR.bcf"
plink_prefix="${outdir}/${sample_id}_sexcheck"

bcftools mpileup -q 20 -Q 20 -r "$chrX_nonPAR" -f "$ref_fasta" "$bam" |
bcftools call -mv -f GQ -Ob -o "$bcf"

bcftools index -f "$bcf"

plink --bcf "$bcf" --check-sex --allow-extra-chr --out "$plink_prefix"