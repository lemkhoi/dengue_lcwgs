include: "auxiliary.smk"
include: "software.smk"
configfile: "pipelines/config.json"
home_dir = config['home_dir']

import pandas as pd
import numpy as np
import sys
sys.path.append(f"{home_dir}software/lcwgsus/")
import lcwgsus

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]

rule compute_bedgraph_nozero:
    input:
        bam = "data/bams/{id}.bam"
    output:
        bedgraph = temp("results/coverage/bedgraphs/{id}_bedgraph_nozero.bed")
    resources: mem = '50G'
    shell: """
        bedtools genomecov -ibam {input.bam} -bg | \
        awk '$1 ~ /^chr(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22)$/' \
        > {output.bedgraph}
    """

rule calculate_avg_coverage:
    input:
        bedgraph = "results/coverage/bedgraphs/{id}_bedgraph_nozero.bed"
    output:
        access_coverage = temp("results/coverage/tmp/{id}_access_coverage.txt"),
        avg_coverage = "results/coverage/depth/samples/{id}_depth.tsv"
    params:
        access_bed = config['access_bed'],
        access_bed_length = 2526390487
    shell: """
        mkdir -p results/coverage/tmp/
        mkdir -p results/coverage/depth/samples/

        bedtools intersect -a {input.bedgraph} -b {params.access_bed} -wb | cut -f1-4 > {output.access_coverage}
        sum_product=$(awk '{{ sum += ($3-$2)*$4 }} END {{ print sum }}' {output.access_coverage})
        result=$(echo "scale=4; ($sum_product/{params.access_bed_length})" | bc)
        echo "{wildcards.id}\t$result" > {output.avg_coverage}
    """

rule aggregate_avg_coverage:
    input:
        files = expand("results/coverage/depth/samples/{id}_depth.tsv", id = samples_lc)
    output:
        avg_coverage = "results/coverage/depth_all.tsv"
    localrule: True
    shell: """
        cat {input.files} >> {output.avg_coverage}
    """

rule calculate_hla_coverage:
    input:
        bedgraph = "results/coverage/bedgraphs/{id}_bedgraph_nozero.bed"
    output:
        access_coverage = temp("results/coverage/hla_tmp/{id}_access_coverage.txt"),
        avg_coverage = "results/coverage/hla_depth/samples/{id}_depth.tsv"
    params:
        hla_bed = config['hla_bed'],
        hla_bed_length = 25500
    shell: """
        mkdir -p results/coverage/hla_tmp/
        mkdir -p results/coverage/hla_depth/samples/

        bedtools intersect -a {input.bedgraph} -b {params.hla_bed} -wb | cut -f1-4 > {output.access_coverage}
        sum_product=$(awk '{{ sum += ($3-$2)*$4 }} END {{ print sum }}' {output.access_coverage})
        result=$(echo "scale=4; ($sum_product/{params.hla_bed_length})" | bc)
        echo "{wildcards.id}\t$result" > {output.avg_coverage}
    """

rule aggregate_hla_coverage:
    input:
        files = expand("results/coverage/hla_depth/samples/{id}_depth.tsv", id = samples_lc)
    output:
        avg_coverage = "results/coverage/hla_depth_all.tsv"
    localrule: True
    shell: """
        cat {input.files} >> {output.avg_coverage}
    """

rule calculate_uncoverage_rate:
    input:
        bam = "data/bams/{id}.bam",
        bedgraph = "results/coverage/bedgraphs/{id}_bedgraph_nozero.bed"
    output:
        uncoverage_rate = "results/coverage/uncov/samples/{id}_uncov.tsv"
    params:
        access_bed = config['access_bed']
    resources:
        mem = '60G'
    threads: 4
    shell: """
        mkdir -p results/coverage/uncov/samples/

        result=$(bedtools coverage -a {params.access_bed} -b {input.bedgraph} -hist | grep all | head -n 1 | cut -f5)
        echo "{wildcards.id}\t$result" > {output.uncoverage_rate}
    """

rule aggregate_uncoverage_rate:
    input:
        files = expand("results/coverage/uncov/samples/{id}_uncov.tsv", id = samples_lc)
    output:
        uncoverage_rate = "results/coverage/uncov_all.tsv"
    localrule: True
    shell: """
        cat {input.files} >> {output.uncoverage_rate}
    """

rule compute_bedgraph_sex_chr:
    input:
        bam = "data/bams/{id}.bam"
    output:
        bedgraph = temp("results/coverage/bedgraphs/sex/{id}_bedgraph_nozero.bed")
    resources: mem = '50G'
    shell: """
        mkdir -p results/coverage/bedgraphs/sex/

        bedtools genomecov -ibam {input.bam} -bg | \
        awk '$1 ~ /^chr(X|Y)$/' \
        > {output.bedgraph}
    """

rule calculate_sexchr_coverage:
    input:
        bedgraph = "results/coverage/bedgraphs/sex/{id}_bedgraph_nozero.bed"
    output:
        cx = temp("results/coverage/sex/{id}.chrXY.tsv"),
        cov = "results/coverage/sex/{id}.sexchr.tsv",
    params:
        lenX = 156040895,
        lenY = 57227415
    shell: """
        mkdir -p results/coverage/sex/

        grep 'chrX' {input.bedgraph} > {output.cx}
        sum_product=$(awk '{{ sum += ($3-$2)*$4 }} END {{ print sum }}' {output.cx})
        resultX=$(echo "scale=4; ($sum_product/{params.lenX})" | bc)

        grep 'chrY' {input.bedgraph} > {output.cx}
        sum_product=$(awk '{{ sum += ($3-$2)*$4 }} END {{ print sum }}' {output.cx})
        resultY=$(echo "scale=4; ($sum_product/{params.lenY})" | bc)
        echo "{wildcards.id}\t$resultX\t$resultY" > {output.cov}
    """

rule aggregate_sex_cov:
    input:
        files = expand("results/coverage/sex/{id}.sexchr.tsv", id = samples_lc)
    output:
        sex_cov = "results/coverage/sex_all.tsv"
    localrule: True
    shell: """
        cat {input.files} >> {output.sex_cov}
    """
