configfile: "pipelines/config.json"
include: "auxiliary.smk"
home_dir = config['home_dir']

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append(f"{home_dir}software/lcwgsus/")
import lcwgsus

rule trimmomatic:
    input:
        fastq1 = "data/fastq/{id}_1.fastq.gz",
        fastq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        fwd_pair = temp("data/fastq_cleaned/{id}_1.fastq.gz"),
        rev_pair = temp("data/fastq_cleaned/{id}_2.fastq.gz"),
        fwd_unpair = temp("data/fastq_cleaned/{id}_unpaired_1.fastq.gz"),
        rev_unpair = temp("data/fastq_cleaned/{id}_unpaired_2.fastq.gz")
    params:
        adapters = config['adapter']
    threads: 4
    resources:
        mem = '30G'
    shell: """
        trimmomatic PE -phred33 \
        {input.fastq1} {input.fastq2} \
        {output.fwd_pair} {output.fwd_unpair} \
        {output.rev_pair} {output.rev_unpair} \
        ILLUMINACLIP:{params.adapters}:2:30:10:2:true MINLEN:50
    """
