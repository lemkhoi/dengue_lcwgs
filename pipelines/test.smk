configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
import pyreadr
home_dir = config['home_dir']
sys.path.append(f"{home_dir}software/lcwgsus/")
sys.path.append(f'{home_dir}software/QUILT_test/QUILT/Python/')
import lcwgsus
from lcwgsus.variables import *
from hla_phase_functions import *
from hla_align_functions import *

hla_ref_panel_indir = "results/hla/imputation/ref_panel/auxiliary_files/"
hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
samples_fv = read_tsv_as_lst('data/sample_tsvs/fv_idt_names.tsv')
read_lengths = config["sr_read_lengths"]

rule all:
    input:
        uncoverage_rate_all = "results/sr_coverage/uncov_all.tsv",
        uncoverage_rate = expand("results/sr_coverage/uncov/samples/{rl}_uncov.tsv", rl = read_lengths)


rule compute_bedgraph_nozero:
    input:
        fq1_R1 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250060_repeat/A250074/Samples/DefaultProject/{id}/{id}_R1.fastq.gz",
        fq1_R2 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250060_repeat/A250074/Samples/DefaultProject/{id}/{id}_R2.fastq.gz",

        fq2_R1 = ""
    output:
        bedgraph = temp("results/coverage/bedgraphs/{rl}_bedgraph_nozero.bed")
    resources: mem = '50G'
    shell: """
        zcat sample_R1_L001.fastq.gz sample_R1_L002.fastq.gz | gzip > sample_R1_merged.fastq.gz

    """
