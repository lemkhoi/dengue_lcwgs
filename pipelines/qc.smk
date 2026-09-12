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

rule merge_metrics:
    input:
        dup_rate_all = "results/dup_rate/dup_rate_all.tsv",
        depth_all = "results/coverage/depth_all.tsv",
        uncov_all = "results/coverage/uncov_all.tsv",
        sex_all = "results/coverage/sex_all.tsv",
        hla_all = "results/coverage/hla_depth_all.tsv"
    output:
        metrics = "results/metrics.tsv"
    localrule: True
    run:
        def merge_metrics_outer(cov_all, uncov_all, dup_all, sex_all, hla_all, save=False, outdir=None, save_name=None):
            coverage = pd.read_csv(cov_all, sep = '\t', header = None)
            coverage.columns = ['sample', 'coverage']

            hla_coverage = pd.read_csv(hla_all, sep = '\t', header = None)
            hla_coverage.columns = ['sample', 'hla_coverage']

            uncoverage = pd.read_csv(uncov_all, sep = '\t', header = None)
            uncoverage.columns = ['sample', 'uncoverage']
            
            dup_rate = pd.read_csv(dup_all, sep = '\t', header = None)
            dup_rate.columns = ['sample', 'dup_rate']

            sex = pd.read_csv(sex_all, sep = '\t', header = None)
            sex.columns = ['sample', 'chrX', 'chrY']
            
            metrics = pd.merge(coverage, hla_coverage, on = 'sample', how='outer')
            metrics = pd.merge(metrics, dup_rate, on = 'sample', how='outer')
            metrics = pd.merge(metrics, uncoverage, on = 'sample', how='outer')
            metrics['skew'] = 0
            
            def calculate_skew(r):
                uncov = r['uncoverage']
                cov = r['coverage']
                r['skew'] = uncov/poisson.pmf(0, cov)
                return r
            
            metrics = metrics.apply(calculate_skew, axis = 1)
            metrics = pd.merge(metrics, sex, on = 'sample', how='outer')
            
            save_tsv(metrics, save, outdir, save_name)
            return metrics
        merge_metrics_outer(input.depth_all, input.uncov_all, input.dup_rate_all,
              input.sex_all, input.hla_all, True, 'results/', 'metrics.tsv')

rule plot_hist_coverage:
    input:
        metrics = "results/metrics.tsv"
    output:
        hist_cov = "graphs/qc/hist_coverage.png"
    localrule: True
    run:
        lcwgsus.plot_hist_coverage(input.metrics, None, True, 'graphs/qc/', 'hist_coverage.png')

rule plot_hist_hla_coverage:
    input:
        metrics = "results/metrics.tsv"
    output:
        hist_cov = "graphs/qc/hist_hla_coverage.png"
    localrule: True
    run:
        lcwgsus.plot_hist_hla_coverage(input.metrics, None, True, 'graphs/qc/', 'hist_hla_coverage.png')

rule plot_qc_metrics:
    input:
        metrics = "results/metrics.tsv"
    output:
        qc_metrics = "graphs/qc/qc_metrics.png"
    localrule: True
    run:
        lcwgsus.plot_qc_metrics(input.metrics, None, True, 'graphs/qc/', 'qc_metrics.png')

rule plot_sex_coverage:
    input:
        metrics = "results/metrics.tsv"
    output:
        sex_cov = "graphs/qc/sex_coverage.png"
    localrule: True
    run:
        lcwgsus.plot_sex_coverage(input.metrics, None, None, True, 'graphs/qc/', 'sex_coverage.png')

rule plot_hist_coverage_per_batch:
    input:
        metrics = "results/metrics.tsv"
    output:
        hist_cov = "graphs/qc/{batch}/hist_coverage.png"
    localrule: True
    params:
        sample_file = 'data/sample_tsvs/batches/samples_{batch}.tsv'
    run:
        samples = lcwgsus.read_tsv_as_lst(params.sample_file)
        lcwgsus.plot_hist_coverage(input.metrics, samples, True, f'graphs/qc/{wildcards.batch}/', 'hist_coverage.png')

rule plot_hist_hla_coverage_per_batch:
    input:
        metrics = "results/metrics.tsv"
    output:
        hist_cov = "graphs/qc/{batch}/hist_hla_coverage.png"
    localrule: True
    params:
        sample_file = 'data/sample_tsvs/batches/samples_{batch}.tsv'
    run:
        samples = lcwgsus.read_tsv_as_lst(params.sample_file)
        lcwgsus.plot_hist_hla_coverage(input.metrics, samples, True, f'graphs/qc/{wildcards.batch}/', 'hist_hla_coverage.png')

rule plot_qc_metrics_per_batch:
    input:
        metrics = "results/metrics.tsv"
    output:
        qc_metrics = "graphs/qc/{batch}/qc_metrics.png"
    localrule: True
    params:
        sample_file = 'data/sample_tsvs/batches/samples_{batch}.tsv'
    run:
        samples = lcwgsus.read_tsv_as_lst(params.sample_file)
        lcwgsus.plot_qc_metrics(input.metrics, samples, True, f'graphs/qc/{wildcards.batch}/', 'qc_metrics.png')

rule plot_sex_coverage_per_batch:
    input:
        metrics = "results/metrics.tsv"
    output:
        sex_cov = "graphs/qc/{batch}/sex_coverage.png"
    localrule: True
    params:
        sample_file = 'data/sample_tsvs/batches/samples_{batch}.tsv'
    run:
        samples = lcwgsus.read_tsv_as_lst(params.sample_file)
        lcwgsus.plot_sex_coverage(input.metrics, None, samples, True, f'graphs/qc/{wildcards.batch}/', 'sex_coverage.png')

