#!/bin/bash

#remove samples with fraction uncovered > 60%
#remove INFO score < 0.9 (variants still stick around after zp's filter)
module load BCFtools

data_dir="/well/ansari/users/osr869/host/lcwgs_edits/results/wip_vcfs/GAsP_b38/merged/gw_omni_filtered"
filtered_dir="/well/ansari/users/osr869/host/lcwgs_edits/topmed/covered40_samples"
sample_list="/well/ansari/users/osr869/host/lcwgs_edits/topmed/covered40_samples.txt"

set -euo pipefail

mkdir -p "${filtered_dir}"

for i in {1..22}; do
    echo "Processing chr${i}"

    bcftools view \
        -S "${sample_list}" \
        -Ou \
        "${data_dir}/lc.chr${i}.vcf.gz" | \
    bcftools filter \
        -i 'INFO/INFO_SCORE > 0.9' \
        -Oz \
        -o "${filtered_dir}/lc.chr${i}.vcf.gz"

    bcftools index -t "${filtered_dir}/lc.chr${i}.vcf.gz"
done