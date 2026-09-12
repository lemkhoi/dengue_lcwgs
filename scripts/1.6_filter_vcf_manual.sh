#!/bin/bash

# module load Miniforge3/24.1.2-0
# conda activate lcwgs
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
ldd /well/ansari/users/osr869/host/lcwgs_edits/software/QCTool/qctool/build/release/apps/qctool_v2.2.2 | grep stdc++

cd /well/ansari/users/osr869/host/lcwgs_edits/

for chr in {1..22}; do

    samples="results/wip_vcfs/GAsP_b38/merged/filtered/quilt.chr${chr}.samples"
    vcf="results/wip_vcfs/GAsP_b38/merged/maf/quilt.chr${chr}.vcf.gz"
    bgen="results/wip_vcfs/GAsP_b38/merged/filtered/quilt.chr${chr}.bgen"
    snpstat="results/wip_vcfs/GAsP_b38/merged/filtered/quilt.chr${chr}.snpstat.tsv"

    echo "ID" > "$samples"
    echo "0" >> "$samples"
    bcftools query -l "$vcf" >> "$samples"

    /well/ansari/users/osr869/host/lcwgs_edits/software/QCTool/qctool/build/release/apps/qctool_v2.2.2 \
        -g "$vcf" \
        -s "$samples" \
        -vcf-genotype-field GP \
        -og "$bgen" \
        -bgen-bits 8 \
        -bgen-compression zstd

    /well/ansari/users/osr869/host/lcwgs_edits/software/QCTool/qctool/build/release/apps/qctool_v2.2.2 \
        -g "$bgen" \
        -s "$samples" \
        -snp-stats \
        -snp-stats-columns info \
        -osnp "$snpstat"

done