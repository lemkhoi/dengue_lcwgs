# include: "prepare.smk"

include: "preprocess.smk"
# include: "reference.smk"
include: "alignment.smk"

include: "dup_rate.smk"
include: "coverage.smk"
include: "qc.smk"

# include: "imputation_prep.smk"
include: "imputation.smk"
include: "filter_vcf.smk"

include: "hla.smk"
include: "hla_imputation_method.smk"
include: "regions.smk"

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

samples_lc = read_tsv_as_lst(config['samples_lc'])
chromosome = [i for i in range(1,23)]

RECOMB_POP = config["RECOMB_POP"]
ref_panel = config['ref_panel']
batches = config["seq_batches"]

hla_genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
IPD_IMGT_version = config['IPD_IMGT_version']

region_file = "data/imputation_accessories/5Mb_chunks.json"
mGen_vcf_prefix = "data/ref_panel/GAsP_b38/regions/chr"
mGen_chunk_RData, mGen_chunk_vcf_lst, mGen_chunk_vcf_dict = get_vcf_concat_lst(region_file, '', mGen_vcf_prefix)

all_RData = {}
all_vcf_lst = {}
all_vcf_dict = {}
for p in batches:
    region = f"results/imputation/refs/{ref_panel}/regions.json"
    ref_prefix = f"results/imputation/refs/{ref_panel}/RData/ref_package.chr"
    vcf_prefix = f"results/imputation/vcfs/{ref_panel}/{p}/regions/quilt.chr"
    all_RData[p], all_vcf_lst[p], all_vcf_dict[p] = get_vcf_concat_lst(region, ref_prefix, vcf_prefix)

# rule prepare_all:
#     input:
#         lifted = "results/hla/reference/multiEth_sites.b38.vcf.gz",

#         csv = "data/chip/omni5m/omni5m.csv",
#         pos = "data/chip/omni5m/omni5m_sites.tsv",

#         db = expand(f'{home_dir}recyclable_files/hla_reference_files/v{IPD_IMGT_version}_aligners/{{gene}}.ssv', gene = HLA_GENES_ALL_EXPANDED),
#         db_filtered = expand(f'{home_dir}recyclable_files/hla_reference_files/v{IPD_IMGT_version}_oneKG_only/{{gene}}.ssv', gene = HLA_GENES_ALL_EXPANDED),

rule alignment_and_qc_all:
    input:        
        bams = expand("data/bams/{id}.bam", id = samples_lc),
        bamlists = expand("results/imputation/bamlists/{batch}/bamlist.txt", batch = batches),

        duprate = expand("results/dup_rate/samples/dup_rate_{id}.tsv", id = samples_lc),
        dup_rate_all = "results/dup_rate/dup_rate_all.tsv",
        depth = expand("results/coverage/depth/samples/{id}_depth.tsv", id = samples_lc),
        depth_all = "results/coverage/depth_all.tsv",
        hla_depth = expand("results/coverage/hla_depth/samples/{id}_depth.tsv", id = samples_lc),
        hla_depth_all = "results/coverage/hla_depth_all.tsv",
        uncov = expand("results/coverage/uncov/samples/{id}_uncov.tsv", id = samples_lc),
        uncov_all = "results/coverage/uncov_all.tsv",
        sex = expand("results/coverage/sex/{id}.sexchr.tsv", id = samples_lc),
        sex_all = "results/coverage/sex_all.tsv",

        metrics = "results/metrics.tsv",
        hist_cov = "graphs/qc/hist_coverage.png",
        hist_hla_cov = "graphs/qc/hist_hla_coverage.png",
        qc_metrics = "graphs/qc/qc_metrics.png",
        sex_cov = "graphs/qc/sex_coverage.png",

        # hist_cov_per_batch = expand("graphs/qc/{batch}/hist_coverage.png", batch = batches),
        # hist_hla_cov_per_batch = expand("graphs/qc/{batch}/hist_hla_coverage.png", batch = batches),
        # qc_metrics_per_batch = expand("graphs/qc/{batch}/qc_metrics.png", batch = batches),
        # sex_cov_per_batch = expand("graphs/qc/{batch}/sex_coverage.png", batch = batches),

rule qc_all:
    input:
        bamlists = expand("results/imputation/bamlists/{batch}/bamlist.txt", batch = batches),

        duprate = expand("results/dup_rate/samples/dup_rate_{id}.tsv", id = samples_lc),
        dup_rate_all = "results/dup_rate/dup_rate_all.tsv",
        depth = expand("results/coverage/depth/samples/{id}_depth.tsv", id = samples_lc),
        depth_all = "results/coverage/depth_all.tsv",
        hla_depth = expand("results/coverage/hla_depth/samples/{id}_depth.tsv", id = samples_lc),
        hla_depth_all = "results/coverage/hla_depth_all.tsv",        
        uncov = expand("results/coverage/uncov/samples/{id}_uncov.tsv", id = samples_lc),
        uncov_all = "results/coverage/uncov_all.tsv",
        sex = expand("results/coverage/sex/{id}.sexchr.tsv", id = samples_lc),
        sex_all = "results/coverage/sex_all.tsv",

        metrics = "results/metrics.tsv",
        hist_cov = "graphs/qc/hist_coverage.png",
        hist_hla_cov = "graphs/qc/hist_hla_coverage.png",
        qc_metrics = "graphs/qc/qc_metrics.png",
        sex_cov = "graphs/qc/sex_coverage.png",
        
        # hist_cov_per_batch = expand("graphs/qc/{batch}/hist_coverage.png", batch = batches),
        # hist_hla_cov_per_batch = expand("graphs/qc/{batch}/hist_hla_coverage.png", batch = batches),
        # qc_metrics_per_batch = expand("graphs/qc/{batch}/qc_metrics.png", batch = batches),
        # sex_cov_per_batch = expand("graphs/qc/{batch}/sex_coverage.png", batch = batches),

rule imp_all:
    input:
        RData = [convert_dict_to_lst(all_RData)],
        vcfs = expand(f"results/imputation/vcfs/{ref_panel}/{{batch}}/quilt.chr{{chr}}.vcf.gz", chr = chromosome, batch = batches),

        confident_lc_vcf = expand("results/wip_vcfs/GAsP_b38/merged/gw_filtered/lc.chr{chr}.vcf.gz", chr = chromosome),
        two_stage_lc_vcf = expand("results/wip_vcfs/GAsP_b38/merged/gw_omni_filtered/lc.chr{chr}.vcf.gz", chr = chromosome),
        two_stage_hla_vcf = "results/wip_vcfs/GAsP_b38/merged/hla_filtered/lc.chr6.vcf.gz"

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

rule hla_all:
    input:
        hla_called = expand("results/hla/call/{id}/hla/R1_bestguess.txt", id = samples_lc),

        AS_matrices = expand(f"results/hla/imputation/WFA_alignments/v{IPD_IMGT_version}/{{id}}/{{gene}}/AS_matrix.ssv", gene = hla_genes, id = samples_lc),
        ref_panel = expand("results/hla/imputation/ref_panel/QUILT_prepared_reference_method/HLA{gene}fullallelesfilledin.RData", gene = hla_genes),
        hla_imputed = expand("results/hla/imputation/QUILT_HLA_result_method/{id}/{gene}/quilt.hla.output.combined.all.txt", gene = hla_genes, id = samples_lc),

rule sv_all:
    input:
        coverage_files = coverotron_all_ary

rule hla_coverage_all:
    input:
        gyp_coverage = expand('results/coverage/coverotron/hla/hla_{batch}.tsv', batch = batches)
