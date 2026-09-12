configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"
home_dir = config['home_dir']

import io
import os
import re
import json
import pandas as pd
import numpy as np
import math
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append(f"{home_dir}software/lcwgsus/")
sys.path.append(f"{home_dir}software/QUILT_sus/QUILT/Python/")
import lcwgsus

from lcwgsus.variables import *
from hla_phase_functions import *
from hla_align_functions import *

hla_ref_panel_indir = "results/hla/imputation/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
IPD_IMGT_version = config['IPD_IMGT_version']
ref_panel = config["ref_panel"]
RECOMB_POP = config["RECOMB_POP"]

rule prepare_hla_reference_panel_method:
    input:
        hla_types_panel = f"{hla_ref_panel_indir}20181129_HLA_types_full_1000_Genomes_Project_panel.txt",
        ipd_igmt = f"{hla_ref_panel_indir}IPD-IMGT-HLA_v{IPD_IMGT_version}.zip",
        fasta = "data/references/GRCh38.fa",
        genetic_map = f"data/imputation_accessories/maps/{RECOMB_POP}-chr6-final.b38.txt.gz",
        hap = f"{hla_ref_panel_indir}oneKG.hap.gz",
        legend = f"{hla_ref_panel_indir}oneKG.legend.gz",
        sample = f"{hla_ref_panel_indir}oneKG.samples"
    output:
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{gene}fullallelesfilledin.RData", gene = hla_genes)
    resources:
        mem = '50G'
    threads: 5
    conda: "sus2"
    params:
        quilt_hla_prep = tools['quilt_hla_prep'],
        refseq = f"{home_dir}software/QUILT_sus/hla_ancillary_files/refseq.hg38.chr6.26000000.34000000.txt.gz",
        region_exclude_file = f"{home_dir}software/QUILT_sus/hla_ancillary_files/hlagenes.txt",
        hla_ref_panel_outdir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/"
    shell: """
        {params.quilt_hla_prep} \
        --outputdir={params.hla_ref_panel_outdir} \
        --nGen=100 \
        --hla_types_panel={input.hla_types_panel} \
        --ipd_igmt_alignments_zip_file={input.ipd_igmt} \
        --ref_fasta={input.fasta} \
        --refseq_table_file={params.refseq} \
        --full_regionStart=25587319 \
        --full_regionEnd=33629686 \
        --buffer=500000 \
        --region_exclude_file={params.region_exclude_file} \
        --genetic_map_file={input.genetic_map} \
        --reference_haplotype_file={input.hap} \
        --reference_legend_file={input.legend} \
        --reference_sample_file={input.sample} \
        --hla_regions_to_prepare="c('A','B','C','DQB1','DRB1')" \
        --nCores=5
    """

rule hla_alignment_matrix:
    input:
        bam = "data/bams/{id}.bam",
        db_files = expand(f"{home_dir}recyclable_files/hla_reference_files/v{IPD_IMGT_version}_oneKG_only/{{gene}}.ssv", gene = HLA_GENES_ALL_EXPANDED)
    output:
        matrix = f"results/hla/imputation/WFA_alignments/v{IPD_IMGT_version}/{{id}}/{{gene}}/AS_matrix.ssv"
    resources:
        mem = '120G'
    conda: "sus2"
    threads: 16
    params:
        hla_gene_information_file = f"{home_dir}software/QUILT_sus/hla_ancillary_files/hla_gene_information_expanded.tsv",
        db_dir = f"{home_dir}recyclable_files/hla_reference_files/v{IPD_IMGT_version}_oneKG_only/",
        strict = True,
        align_all = False
    script:
        f"{home_dir}software/QUILT_sus/QUILT/Python/hla_calculate_alignment_matrix.py"

rule hla_imputation_method:
    input:
        bam = "data/bams/{id}.bam",
        ref_panel = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{gene}fullallelesfilledin.RData",
        prepared_db = f'{home_dir}recyclable_files/hla_reference_files/v{IPD_IMGT_version}_oneKG_only/{{gene}}.ssv',
        matrix = f"results/hla/imputation/WFA_alignments/v{IPD_IMGT_version}/{{id}}/{{gene}}/AS_matrix.ssv"
    output:
        bamlist = temp("results/hla/imputation/bamlists_fv/{id}.{gene}.txt"),
        imputed = "results/hla/imputation/QUILT_HLA_result_method/{id}/{gene}/quilt.hla.output.combined.all.txt"
    resources:
        mem = '20G'
    threads: 2
    params:
        quilt_hla = tools['quilt_hla'],
        fa_dict = "data/references/GRCh38.dict",
        ref_dir = "results/hla/imputation/ref_panel/QUILT_prepared_reference_method/"
    conda: "sus2"
    shell: """
        mkdir -p results/hla/imputation/QUILT_HLA_result_method/{wildcards.id}/{wildcards.gene}/
        ls {input.bam} > {output.bamlist}

        {params.quilt_hla} \
        --outputdir="results/hla/imputation/QUILT_HLA_result_method/{wildcards.id}/{wildcards.gene}/" \
        --bamlist={output.bamlist} \
        --region={wildcards.gene} \
        --prepared_hla_reference_dir={params.ref_dir} \
        --quilt_hla_haplotype_panelfile={params.ref_dir}/quilt.hrc.hla.{wildcards.gene}.haplotypes.RData \
        --dict_file={params.fa_dict}
    """