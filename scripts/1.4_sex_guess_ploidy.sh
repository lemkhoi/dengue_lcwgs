#!/usr/bin/env bash
#SBATCH --job-name=sexcheck
#SBATCH --array=1-30%7
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/X_ploidy/sexcheck_%A_%a.out
#SBATCH --error=logs/X_ploidy/sexcheck_%A_%a.err

module load BCFtools
module load PLINK

samples_tsv="/well/ansari/users/osr869/host/lcwgs_edits/results/X_ploidy/test14.tsv"
ref_fasta="/well/ansari/users/osr869/host/lcwgs_edits/data/references/GRCh38.fa"
outdir="/well/ansari/users/osr869/host/lcwgs_edits/results/X_ploidy"

chrX_nonPAR="chrX:2781480-155701382"

mapfile -t lines < "$samples_tsv"

for ((i=0; i<${#lines[@]}; i++)); do
    line="${lines[$i]}"
    sample_id=$(printf '%s' "$line" | cut -f1)
    bam="/well/ansari/users/osr869/host/lcwgs_edits/data/bams/${sample_id}.bam"

    bcf="${outdir}/${sample_id}_chrX_nonPAR.bcf"
    ploidy_out="${outdir}/${sample_id}_guess_ploidy.txt"

    bcftools mpileup -q 20 -Q 20 -r "$chrX_nonPAR" -f "$ref_fasta" "$bam" |
    bcftools call -m -Ob -o "$bcf"

    bcftools index -f "$bcf"

    bcftools +guess-ploidy -v "$bcf" > "$ploidy_out"
done

# indir="/well/ansari/users/osr869/host/lcwgs_edits/results/X_ploidy"
# out_sum="${outdir}/guess_ploidy_combined.tsv"

# printf "Sample\tPredicted_sex\tlogP_Haploid_perSite\tlogP_Diploid_perSite\tnSites\tScore\n" > "$out"

# for f in "${outdir}"/*_guess_ploidy.txt; do
#     [[ -e "$f" ]] || continue

#     awk -v OFS="\t" '
#         $1 == "SEX" && $2 != "Sample" {
#             print $2, $3, $4, $5, $6, $7
#         }
#     ' "$f" >> "$out_sum"
# done

# echo "$out_sum"