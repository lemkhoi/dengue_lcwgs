include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("/well/ansari/users/gjx698/software/lcwgsus/")
import lcwgsus
from lcwgsus.variables import *

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]
batches = config["seq_batches"]

coverotron_all_dict = {}
coverotron_all_ary = []
region_file = "data/imputation_accessories/5Mb_chunks_for_coverage.json"
for b in batches:
    coverotron_output_prefix = f'results/coverage/coverotron/{b}/chr'
    bam_prefix = f'results/coverage/tmp/{b}/chr'
    samples_b = read_tsv_as_lst(f'data/sample_tsvs/batches/samples_{b}.tsv')
    chunk_bams, cov_files_all = get_bam_concat_lst(region_file, samples_b, bam_prefix, coverotron_output_prefix)
    coverotron_all_ary = coverotron_all_ary + cov_files_all
    coverotron_all_dict[b] = chunk_bams

rule keep_primary_reads:
    input:
        bam = 'data/bams/{id}.bam'
    output:
        bam = temp('results/coverage/tmp/{batch}/chr{chr}.{regionStart}.{regionEnd}/{id}.bam')
    resources:
        mem = '30G'
    localrule: True
    shell: """
        mkdir -p results/coverage/tmp/{wildcards.batch}/chr{wildcards.chr}.{wildcards.regionStart}.{wildcards.regionEnd}/

        samtools view -h -f 2 -F 2304 {input.bam} chr{wildcards.chr}:{wildcards.regionStart}-{wildcards.regionEnd} | \
        samtools sort - -o {output.bam}
        samtools index {output.bam}
    """
    
def get_bams_from_batch(wildcards):
    return coverotron_all_dict[wildcards.batch][f'{wildcards.chr}.{wildcards.regionStart}.{wildcards.regionEnd}']

rule calculate_bin_coverage_per_chunk:
    input:
        bam = get_bams_from_batch
    output:
        gyp_coverage = 'results/coverage/coverotron/{batch}/chr{chr}.{regionStart}.{regionEnd}.tsv'
    resources:
        mem = '40G'
    threads: 4
    params: coverotron = tools['coverotron']
    shell: """
        mkdir -p results/coverage/coverotron/{wildcards.batch}/

        {params.coverotron} -bin 1000 \
        -range chr{wildcards.chr}:{wildcards.regionStart}-{wildcards.regionEnd} \
        -reads {input.bam} > \
        {output.gyp_coverage}
    """




rule keep_primary_reads_hla:
    input:
        bam = 'data/bams/{id}.bam'
    output:
        bam = temp('results/coverage/tmp/{batch}/chr6.25000000.34000000/{id}.bam')
    resources:
        mem = '30G'
    localrule: True
    shell: """
        mkdir -p results/coverage/tmp/{wildcards.batch}/chr6.25000000.34000000/

        samtools view -h -f 2 -F 2304 {input.bam} chr6:25000000-34000000 | \
        samtools sort - -o {output.bam}
        samtools index {output.bam}
    """

def get_bams_from_batch1(wildcards):
    samples = read_tsv_as_lst(f'data/sample_tsvs/batches/samples_{wildcards.batch}.tsv')
    bamlist = []
    for s in samples:
        bamlist.append(f'results/coverage/tmp/{wildcards.batch}/chr6.25000000.34000000/{s}.bam')
    return bamlist

rule calculate_bin_coverage_hla:
    input:
        bam = get_bams_from_batch1
    output:
        gyp_coverage = 'results/coverage/coverotron/hla/hla_{batch}.tsv'
    resources:
        mem = '40G'
    threads: 4
    params: coverotron = tools['coverotron']
    shell: """
        mkdir -p results/coverage/coverotron/hla/

        {params.coverotron} -bin 1000 \
        -range chr6:25000000-34000000 \
        -reads {input.bam} > \
        {output.gyp_coverage}
    """
