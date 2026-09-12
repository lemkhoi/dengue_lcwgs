home_dir = "/well/ansari/users/osr869/host/lcwgs_edits/"
import pandas as pd
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.append(f"{home_dir}software/lcwgsus/")
import lcwgsus
chrs = [i for i in range(1,2)]

rule all:
    input:
        f"{home_dir}graphs/info_vs_maf_sum.png",
        f"{home_dir}graphs/sum_INFO.tsv",
        f"{home_dir}graphs/quit_raw_combined.tsv"

rule generate_vcf_df:
    input:
        vcfs = home_dir + "results/wip_vcfs/GAsP_b38/merged/raw_merged/chr{chr}.vcf.gz"
    output:
        vcf_dfs = home_dir + "results/wip_vcfs/GAsP_b38/merged/raw_merged_df/df.chr{chr}.txt"
    resources:
        mem = '40G'
    run:
        df = lcwgsus.read_vcf(input.vcfs)
        df = lcwgsus.extract_info(df, info_cols = ['EAF', 'INFO_SCORE', 'HWE', 'ERC', 'EAC', 'PAF'], attribute = 'INFO', drop_attribute = True)
        df[['chr', 'pos', 'ref', 'alt', 'EAF', 'INFO_SCORE', 'HWE', 'ERC', 'EAC', 'PAF']].to_csv(output.vcf_dfs, sep="\t", index=None)

rule merge_vcf_af:
    input:
        vcf_dfs = home_dir + "results/wip_vcfs/GAsP_b38/merged/raw_merged_df/df.chr{chr}.txt",
        afs = home_dir + "data/gnomAD_MAFs/eas/gnomAD_MAF_eas_chr{chr}.txt"
    output:
        vcf_af = home_dir + "results/wip_vcfs/GAsP_b38/merged/raw_merged_maf/vcf.af.chr{chr}.txt"
    run:
        vcf = pd.read_csv(input.vcf_dfs, sep="\t")
        afs = lcwgsus.read_af(input.afs)

        df = pd.merge(vcf,
                afs,
                on=['chr', 'pos', 'ref', 'alt'],
                how="left").dropna()
        df.to_csv(output.vcf_af, sep="\t", index=None)

rule combine_vcf_af:
    input:
        vcf_af = expand(home_dir + "results/wip_vcfs/GAsP_b38/merged/raw_merged_maf/vcf.af.chr{chr}.txt", chr=chrs)
    output:
        combined_file = f"{home_dir}graphs/quit_raw_combined.tsv"
    run:
        df_list = [pd.read_csv(f, sep="\t") for f in input.vcf_af]
        df = pd.concat(df_list, ignore_index=True)
        df.to_csv(output.combined_file, index=False, sep='\t')

rule generate_info_plot_af:
    input:
        vcf_af = expand(home_dir + "results/wip_vcfs/GAsP_b38/merged/raw_merged_maf/vcf.af.chr{chr}.txt", chr=chrs)
    output:
        plot_file = f"{home_dir}graphs/info_vs_maf_sum.png",
        sum_file = f"{home_dir}graphs/sum_INFO.tsv"
    params:
        MAF_ary = lcwgsus.MAF_ARY
    run:
        # vcf
        df_list = [pd.read_csv(f, sep="\t") for f in input.vcf_af]
        df = pd.concat(df_list, ignore_index=True)
        #df = df[df.INFO_SCORE > 0.9]

        #figure
        df['classes'] = np.digitize(df['MAF'], params.MAF_ary)

        # Calculate boxplot statistics
        summary = (
            df.groupby('classes')['INFO_SCORE']
            .agg(
                q1=lambda x: x.quantile(0.25),
                median='median',
                q3=lambda x: x.quantile(0.75),
                minimum='min',
                maximum='max'
            )
            .reset_index()
        )

        # Draw boxplots manually
        fig, ax = plt.subplots(figsize=(12, 8))
        for _, row in summary.iterrows():
            x = row['classes']
            # whiskers
            ax.vlines(
                x,
                row['minimum'],
                row['maximum'],
                color='black'
            )
            # box
            ax.add_patch(
                plt.Rectangle(
                    (x - 0.3, row['q1']),
                    0.6,
                    row['q3'] - row['q1'],
                    fill=False,
                    edgecolor='black',
                    linewidth=1.5
                )
            )
            # median
            ax.hlines(
                row['median'],
                x - 0.3,
                x + 0.3,
                color='black',
                linewidth=2
            )
        plt.xlabel('Allele Frequencies (%)')
        plt.ylabel('INFO SCORE')
        plt.title('INFO Score vs Allele Frequencies')
        ax.set_xticklabels(params.MAF_ary[np.sort(df['classes'].unique()) - 1])
        plt.savefig(output.plot_file)
        summary.to_csv(output.sum_file, sep='\t', index=False)
