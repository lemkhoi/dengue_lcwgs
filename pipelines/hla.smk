configfile: "pipelines/config.json"
include: "auxiliary.smk"
include: "software.smk"

import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append(f"{home_dir}software/lcwgsus/")
import lcwgsus

samples_lc = read_tsv_as_lst(config['samples_lc'])

rule hla_la_calling:
    input:
        bam = "data/bams/{id}.bam"
    output:
        called = "results/hla/call/{id}/hla/R1_bestguess.txt"
    resources: mem = '40G'
    threads: 4
    params:
        hla_la = tools['hla-la'],
        working_dir = f'{home_dir}hcv/results/hla/call/'
    shell: """
        mkdir -p results/hla/call/{wildcards.id}/
        module load Java/17

        HLA-LA.pl \
        --BAM {input.bam} \
        --graph PRG_MHC_GRCh38_withIMGT \
        --workingDir {params.working_dir} \
        --sampleID {wildcards.id} \
        --maxThreads {threads}

        mv results/hla/call/{wildcards.id}/reads_per_level.txt \
        results/hla/call/{wildcards.id}/hla/

        rm results/hla/call/{wildcards.id}/*bam*
        rm results/hla/call/{wildcards.id}/*fastq
    """