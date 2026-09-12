configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import json
import pandas as pd
import numpy as np
import sys
import os
home_dir = config['home_dir']


samples_A300000 = read_tsv_as_lst('data/sample_tsvs/batches/samples_A300000.tsv')
rule merge1_all:
    input:
        merge1_R1 = expand("results/merge_boost/merge1/{id}/{id}_R1.fastq.gz", id = samples_A300000),
        merge1_R2 = expand("results/merge_boost/merge1/{id}/{id}_R2.fastq.gz", id = samples_A300000)

rule merge1:
    input:
        fq1_R1 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250090/A250090/Samples/DefaultProject/{id}/{id}_R1.fastq.gz",
        fq1_R2 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250090/A250090/Samples/DefaultProject/{id}/{id}_R2.fastq.gz",

        fq2_R1 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250091/A250091/Samples/DefaultProject/{id}/{id}_R1.fastq.gz",
        fq2_R2 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250091/A250091/Samples/DefaultProject/{id}/{id}_R2.fastq.gz",
    output:
        R1 = "results/merge_boost/merge1/{id}/{id}_R1.fastq.gz",
        R2 = "results/merge_boost/merge1/{id}/{id}_R2.fastq.gz"
    resources: mem = '50G'
    shell: """
        mkdir -p results/merge_boost/merge1/{wildcards.id}/

        zcat {input.fq1_R1} {input.fq2_R1} | gzip > {output.R1}
        zcat {input.fq1_R2} {input.fq2_R2} | gzip > {output.R2}
    """



samples_A100000 = read_tsv_as_lst('data/sample_tsvs/batches/samples_A100000.tsv')
samples_A210000 = read_tsv_as_lst('data/sample_tsvs/batches/samples_A200000.tsv')[:4]
rule merge2_all:
    input:
        merge21_R1 = expand("results/merge_boost/merge21/{id}/{id}_R1.fastq.gz", id = samples_A100000),
        merge21_R2 = expand("results/merge_boost/merge21/{id}/{id}_R2.fastq.gz", id = samples_A100000),

        merge22_R1 = expand("results/merge_boost/merge22/{id}/{id}_R1.fastq.gz", id = samples_A210000),
        merge22_R2 = expand("results/merge_boost/merge22/{id}/{id}_R2.fastq.gz", id = samples_A210000)

rule merge21:
    input:
        fq1_R1 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250049/Samples/DefaultProject/{id}/{id}_R1.fastq.gz",
        fq1_R2 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250049/Samples/DefaultProject/{id}/{id}_R2.fastq.gz",

        fq2_R1 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250065/A250065/Samples/DefaultProject/{id}/{id}_R1.fastq.gz",
        fq2_R2 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250065/A250065/Samples/DefaultProject/{id}/{id}_R2.fastq.gz",
    output:
        R1 = "results/merge_boost/merge21/{id}/{id}_R1.fastq.gz",
        R2 = "results/merge_boost/merge21/{id}/{id}_R2.fastq.gz"
    resources: mem = '50G'
    shell: """
        mkdir -p results/merge_boost/merge21/{wildcards.id}/

        zcat {input.fq1_R1} {input.fq2_R1} | gzip > {output.R1}
        zcat {input.fq1_R2} {input.fq2_R2} | gzip > {output.R2}
    """

rule merge22:
    input:
        fq1_R1 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250049/Samples/DefaultProject/{id}/{id}_R1.fastq.gz",
        fq1_R2 = "/well/ansari/projects/human/LCWGS/Vietnam_main_phase/raw_data/dengue/A250049/Samples/DefaultProject/{id}/{id}_R2.fastq.gz",

        fq2_R1 = "results/merge_boost/merge1/{id}/{id}_R1.fastq.gz",
        fq2_R2 = "results/merge_boost/merge1/{id}/{id}_R2.fastq.gz"
    output:
        R1 = "results/merge_boost/merge22/{id}/{id}_R1.fastq.gz",
        R2 = "results/merge_boost/merge22/{id}/{id}_R2.fastq.gz"
    resources: mem = '50G'
    shell: """
        mkdir -p results/merge_boost/merge22/{wildcards.id}/

        zcat {input.fq1_R1} {input.fq2_R1} | gzip > {output.R1}
        zcat {input.fq1_R2} {input.fq2_R2} | gzip > {output.R2}
    """