configfile: "pipelines/config.json"
include: "software.smk"
include: "auxiliary.smk"
home_dir = config['home_dir']

import os
import json
import pandas as pd
import numpy as np
import sys
sys.path.append(f"{home_dir}software/lcwgsus/")
import lcwgsus
chromosome = [i for i in range(1,23)]

RECOMB_POP = config["RECOMB_POP"]
WINDOWSIZE = config["WINDOWSIZE"]
BUFFER = config["BUFFER"]
ref_panel = config['ref_panel']

samples_lc = read_tsv_as_lst(config['samples_lc'])

def get_bams_from_batch(wildcards):
    file = f"data/sample_tsvs/batches/samples_{wildcards.batch}.tsv"
    return [f"data/bams/{s}.bam" for s in read_tsv_as_lst(file)]

rule prepare_bamlist:
    input:
        bams = get_bams_from_batch
    output:
        bamlist = "results/imputation/bamlists/{batch}/bamlist.txt"
    localrule: True
    shell: """
        mkdir -p results/imputation/bamlists/{wildcards.batch}/

        ls {input.bams} > {output.bamlist}
    """

rule convert_recomb:
    input:
        m = f"results/imputation/{RECOMB_POP}/{RECOMB_POP}-{{chr}}-final.txt.gz",
        chain = "results/imputation/hg19ToHg38.over.chain.gz",
        liftover = "results/imputation/liftOver"
    output:
        f"results/imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    params:
        script = "scripts/quilt_accessories/make_b38_recomb_map.R",
        outdir = f"{home_dir}hcv/results/imputation/"
    threads: 1
    wildcard_constraints: chr='\d{1,2}'
    localrule: True
    shell: """
        R -f {params.script} --args {params.outdir} {RECOMB_POP} {wildcards.chr}
    """

rule convert_ref:
    input:
        vcf = f"data/ref_panel/{ref_panel}/{ref_panel}.chr{{chr}}.vcf.gz",
        tbi = f"data/ref_panel/{ref_panel}/{ref_panel}.chr{{chr}}.vcf.gz.tbi"
    output:
        tmp_vcf = temp(f"results/imputation/refs/{ref_panel}/tmp.{ref_panel}.chr{{chr}}.vcf.gz"),
        hap = f"results/imputation/refs/{ref_panel}/{ref_panel}.chr{{chr}}.hap.gz",
        legend = f"results/imputation/refs/{ref_panel}/{ref_panel}.chr{{chr}}.legend.gz",
        samples = f"results/imputation/refs/{ref_panel}/{ref_panel}.chr{{chr}}.samples"
    wildcard_constraints:
        chr='\d{1,2}'
    threads: 4
    resources: mem = '30G'
    params: 
        outdir = f"results/imputation/refs/{ref_panel}/",
        prefix = f"results/imputation/refs/{ref_panel}/{ref_panel}"
    shell: """
        mkdir -p {params.outdir}

        bcftools norm -m+ {input.vcf} | bcftools view -m2 -M2 -v snps | bcftools sort -Oz -o {output.tmp_vcf}
        tabix {output.tmp_vcf}

        bcftools convert -h \
        {params.prefix}.chr{wildcards.chr} {output.tmp_vcf}
    """

rule determine_chunks:
    input:
        legend = expand(f"results/imputation/refs/{ref_panel}/{ref_panel}.chr{{chr}}.legend.gz", chr = chromosome)
    output:
        json = f"results/imputation/refs/{ref_panel}/regions.json"
    params:
        analysis_dir = f"results/imputation/refs/{ref_panel}/",
        code = "scripts/quilt_accessories/determine_chunks.R",
        json = "data/imputation_accessories/5Mb_chunks.json"
    shell: """
        mkdir -p {params.analysis_dir}

        if [ -f {params.json} ]; then
            cp {params.json} {output.json}
        else
            Rscript {params.code} {params.analysis_dir} {WINDOWSIZE} {BUFFER} {ref_panel}
        fi
    """