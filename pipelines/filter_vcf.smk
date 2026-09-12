include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"
home_dir = config['home_dir']

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append(f"{home_dir}software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

chromosome = [i for i in range(1,23)]
batches = config["seq_batches"]

rule generate_topmed_manifest_maf:
    input:
        af = "data/gnomAD_MAFs/eas/gnomAD_MAF_eas_chr{chr}.txt"
    output:
        manifest = "results/pre_gwas/variant_manifest_by_maf/chr{chr}.tsv"
    resources: mem = '40G'
    threads: 4
    params:
        maf = config['maf_filter']
    shell: """
        awk '$5 > {params.maf} {{OFS="\t"; print $1, $2}}' {input.af} > {output.manifest}
    """

rule merge_batch_vcf_and_filter_by_gnomAD_maf:
    input:
        lc_vcf = expand("results/imputation/vcfs/GAsP_b38/{batch}/quilt.chr{chr}.vcf.gz", batch = batches, allow_missing = True),
        manifest = "results/pre_gwas/variant_manifest_by_maf/chr{chr}.tsv"
    output:
        merged = "results/wip_vcfs/GAsP_b38/merged/maf/quilt.chr{chr}.vcf.gz",
        tmp_vcf = "results/wip_vcfs/GAsP_b38/merged/raw_merged/chr{chr}.vcf.gz"
    resources: mem = '40G'
    threads: 4
    shell: """
        mkdir -p results/wip_vcfs/GAsP_b38/merged/maf/

        mkdir -p results/wip_vcfs/GAsP_b38/merged/raw_merged/

        bcftools merge -Oz -o {output.tmp_vcf} {input.lc_vcf}
        tabix -f {output.tmp_vcf}

        bcftools view -R {input.manifest} {output.tmp_vcf} | \
        bcftools annotate -x ID -I +'%CHROM\_%POS\_%REF\_%ALT' -Oz -o {output.merged}
        tabix -f {output.merged}
        rm {output.tmp_vcf}.tbi
    """

rule generate_merged_info:
    input:
        merged = "results/wip_vcfs/GAsP_b38/merged/maf/quilt.chr{chr}.vcf.gz"
    output:
        sample = temp("results/wip_vcfs/GAsP_b38/merged/filtered/quilt.chr{chr}.samples"), 
        bgen = temp("results/wip_vcfs/GAsP_b38/merged/filtered/quilt.chr{chr}.bgen"),
        snp_stat = temp("results/wip_vcfs/GAsP_b38/merged/filtered/quilt.chr{chr}.snpstat.tsv"),
    resources: mem = '40G'
    threads: 4
    params:
        qctool = tools['qctool']
    shell: """
        mkdir -p results/wip_vcfs/GAsP_b38/merged/filtered/

        echo "ID" > {output.sample}
        echo "0" >> {output.sample}
        bcftools query -l {input.merged} >> {output.sample}

        {params.qctool} \
        -g {input.merged} \
        -s {output.sample} \
        -vcf-genotype-field GP \
        -og {output.bgen} -bgen-bits 8 -bgen-compression zstd

        {params.qctool} \
        -g {output.bgen} \
        -s {output.sample} \
        -snp-stats -snp-stats-columns info \
        -osnp {output.snp_stat}
    """

rule filter_merged_info:
    input:
        merged = "results/wip_vcfs/GAsP_b38/merged/maf/quilt.chr{chr}.vcf.gz",
        snp_stat = "results/wip_vcfs/GAsP_b38/merged/filtered/quilt.chr{chr}.snpstat.tsv"
    output:
        manifest = temp("results/wip_vcfs/GAsP_b38/merged/filtered/quilt.chr{chr}.gw.manifest.tsv"),
        filtered = "results/wip_vcfs/GAsP_b38/merged/gw_filtered/lc.chr{chr}.vcf.gz",
    resources: mem = '40G'
    threads: 4
    params:
        info = config['gw_info_filter']
    shell: """
        mkdir -p results/wip_vcfs/GAsP_b38/merged/gw_filtered/

        grep -v '^#' {input.snp_stat} | tail -n +2 | \
        awk '$8 > {params.info} {{OFS="\t"; print $3, $4}}' > {output.manifest}

        bcftools view -R {output.manifest} -Oz -o {output.filtered} {input.merged}
        tabix -f {output.filtered}
    """

rule filter_lc_sites:
    input:
        filtered = "results/wip_vcfs/GAsP_b38/merged/gw_filtered/lc.chr{chr}.vcf.gz",
        sites = "data/chip/omni5m/omni5m_sites.tsv"
    output:
        vcf = "results/wip_vcfs/GAsP_b38/merged/gw_omni_filtered/lc.chr{chr}.vcf.gz",
        sites = temp("results/wip_vcfs/GAsP_b38/merged/gw_omni_filtered/site.chr{chr}.vcf.gz"),
    resources: mem = '40G'
    threads: 4
    shell: """
        mkdir -p results/wip_vcfs/GAsP_b38/merged/gw_omni_filtered/

        awk -F'\t' '$1 ~ /^[0-9]+$/ && $1 >= 1 && $1 <= 22 {{ $1="chr"$1; print }}' \
        OFS='\t' {input.sites} > {output.sites}

        bcftools view -R {output.sites} {input.filtered} | \
        bcftools sort -Oz -o {output.vcf}
        tabix -f {output.vcf}
    """

rule filter_chr6_for_hla_imputation:
    input:
        merged = "results/wip_vcfs/GAsP_b38/merged/maf/quilt.chr6.vcf.gz",
        snp_stat = "results/wip_vcfs/GAsP_b38/merged/filtered/quilt.chr6.snpstat.tsv",
        b38_scaffold = "results/hla/reference/multiEth_sites.b38.vcf.gz"
    output:
        manifest = temp("results/wip_vcfs/GAsP_b38/merged/filtered/quilt.chr6.hla.manifest.tsv"),
        scaffold = temp("results/wip_vcfs/GAsP_b38/merged/filtered/multiEth_sites.b38.txt"),
        vcf = "results/wip_vcfs/GAsP_b38/merged/hla_filtered/lc.chr6.vcf.gz",
        tmp_vcf = temp("results/wip_vcfs/GAsP_b38/merged/hla_filtered/tmp.chr6.vcf.gz")
    localrule: True
    params:
        info = config['hla_info_filter']
    shell: """
        mkdir -p results/wip_vcfs/GAsP_b38/merged/hla_filtered/

        grep -v '^#' {input.snp_stat} | tail -n +2 | \
        awk '$8 > {params.info} {{OFS="\t"; print $3, $4}}' > {output.manifest}

        bcftools query -f '%CHROM\t%POS\n' {input.b38_scaffold} > {output.scaffold}

        bcftools view -R {output.scaffold} -m2 -M2 -v snps -Oz -o {output.tmp_vcf} {input.merged}
        tabix -f {output.tmp_vcf}

        bcftools view -R {output.manifest} -Oz -o {output.vcf} {output.tmp_vcf}
        rm {output.tmp_vcf}.tbi
        tabix -f {output.vcf}
    """

