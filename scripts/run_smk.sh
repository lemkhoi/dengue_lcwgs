#!/bin/bash
#SBATCH -A ansari.prj
#SBATCH -J tm_filter
#SBATCH -p short
#SBATCH --mem=80G
#SBATCH -o logs/snakemake_plot.out
#SBATCH -e logs/snakemake_plot.err

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate lcwgs

#cd /gpfs3/well/ansari/users/osr869/host/lcwgs_edits/scripts

#snakemake --snakefile vis_imputation_plots_bcf.smk -j 10
snakemake --snakefile prepare_plink.smk -j 7

