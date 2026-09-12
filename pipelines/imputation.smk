configfile: "pipelines/config.json"
include: "software.smk"
include: "auxiliary.smk"
home_dir = config['home_dir']

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append(f"{home_dir}software/lcwgsus/")
import lcwgsus

chromosome = [i for i in range(1,23)]

RECOMB_POP = config["RECOMB_POP"]
NGEN = config["NGEN"]
WINDOWSIZE = config["WINDOWSIZE"]
BUFFER = config["BUFFER"]
batches = config["seq_batches"]
ref_panel = config['ref_panel']

all_RData = {}
all_vcf_lst = {}
all_vcf_dict = {}

for p in batches:
    region = f"results/imputation/refs/{ref_panel}/regions.json"
    ref_prefix = f"results/imputation/refs/{ref_panel}/RData/ref_package.chr"
    vcf_prefix = f"results/imputation/vcfs/{ref_panel}/{p}/regions/quilt.chr"
    all_RData[p], all_vcf_lst[p], all_vcf_dict[p] = get_vcf_concat_lst(region, ref_prefix, vcf_prefix)

rule prepare_ref:
    input:
        json = f"results/imputation/refs/{ref_panel}/regions.json",
        hap = f"results/imputation/refs/{ref_panel}/{ref_panel}.chr{{chr}}.hap.gz",
        legend = f"results/imputation/refs/{ref_panel}/{ref_panel}.chr{{chr}}.legend.gz",
        recomb =  f"results/imputation/{RECOMB_POP}/{RECOMB_POP}-chr{{chr}}-final.b38.txt.gz"
    output:
        RData = f"results/imputation/refs/{ref_panel}/RData/ref_package.chr{{chr}}.{{regionStart}}.{{regionEnd}}.RData"
    resources:
        mem = '30G'
    params:
        threads = 4
    shell: """
        mkdir -p results/imputation/refs/{ref_panel}/RData/other/
        R -e 'library("data.table"); library("QUILT"); QUILT_prepare_reference( \
        outputdir="results/imputation/refs/{ref_panel}/RData/other/", \
        chr="chr{wildcards.chr}", \
        nGen={NGEN}, \
        reference_haplotype_file="{input.hap}" ,\
        reference_legend_file="{input.legend}", \
        genetic_map_file="{input.recomb}", \
        regionStart={wildcards.regionStart}, \
        regionEnd={wildcards.regionEnd}, \
        buffer=0, \
        output_file="{output.RData}")'
    """

rule quilt_imputation:
    input:
        bamlist = "results/imputation/bamlists/{batch}/bamlist.txt",
        RData = rules.prepare_ref.output.RData
    output:
        vcf = temp(f"results/imputation/vcfs/{ref_panel}/{{batch}}/regions/quilt.chr{{chr}}.{{regionStart}}.{{regionEnd}}.vcf.gz")
    resources:
        mem = '30G'
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell: """
        ## set a seed here, randomly, so can try to reproduce if it fails
        SEED=7626
        mkdir -p "results/imputation/vcfs/{ref_panel}/{wildcards.batch}/regions/"
        R -e 'library("data.table"); library("QUILT"); QUILT( \
        outputdir="results/imputation/refs/{ref_panel}/RData/other/", \
        chr="chr{wildcards.chr}", \
        regionStart={wildcards.regionStart}, \
        regionEnd={wildcards.regionEnd}, \
        buffer=0, \
        bamlist="{input.bamlist}", \
        prepared_reference_filename="{input.RData}", \
        output_filename="{output.vcf}", \
        seed='${{SEED}}')'
    """

def get_input_vcfs_as_list(wildcards):
    return(all_vcf_dict[wildcards.batch][str(wildcards.chr)])

def get_input_vcfs_as_string(wildcards):
    return(" ".join(map(str, all_vcf_dict[wildcards.batch][str(wildcards.chr)])))

rule concat_quilt_vcf:
    input:
        vcfs = get_input_vcfs_as_list
    output:
        vcf = f"results/imputation/vcfs/{ref_panel}/{{batch}}/quilt.chr{{chr}}.vcf.gz"
    resources:
        mem = '30G'
    params:
        threads = 1,
        input_string = get_input_vcfs_as_string
    wildcard_constraints:
        chr='\w{1,2}',
        regionStart='\d{1,9}',
        regionEnd='\d{1,9}'
    shell: """
        if [ -e {output.vcf}.temp1.vcf.gz ]; then
            rm {output.vcf}.temp1.vcf.gz
        fi
        if [ -e {output.vcf}.temp2.vcf.gz ]; then
            rm {output.vcf}.temp2.vcf.gz
        fi

        bcftools concat --ligate-force -Oz -o {output.vcf}.temp1.vcf.gz {params.input_string}

        gunzip -c {output.vcf}.temp1.vcf.gz | grep '#' > {output.vcf}.temp2.vcf
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT:GP:DS\t[%GT:%GP:%DS\t]\n' \
        {output.vcf}.temp1.vcf.gz  >> {output.vcf}.temp2.vcf

        bcftools sort -Oz -o {output.vcf} {output.vcf}.temp2.vcf
        
        tabix {output.vcf}
        rm {output.vcf}.temp*
    """
